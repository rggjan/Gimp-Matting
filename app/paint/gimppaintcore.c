/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <string.h>

#include <gegl.h>

#include "libgimpbase/gimpbase.h"
#include "libgimpmath/gimpmath.h"

#include "paint-types.h"

#include "base/pixel-region.h"
#include "base/temp-buf.h"
#include "base/tile-manager.h"
#include "base/tile.h"

#include "paint-funcs/paint-funcs.h"

#include "core/gimp.h"
#include "core/gimp-utils.h"
#include "core/gimpdrawable.h"
#include "core/gimpimage.h"
#include "core/gimpimage-undo.h"
#include "core/gimppickable.h"
#include "core/gimpprojection.h"

#include "gimppaintcore.h"
#include "gimppaintcoreundo.h"
#include "gimppaintoptions.h"

#include "gimpairbrush.h"

#include "gimp-intl.h"

#define STROKE_BUFFER_INIT_SIZE      2000

enum
{
  PROP_0,
  PROP_UNDO_DESC
};


/*  local function prototypes  */

static void      gimp_paint_core_finalize            (GObject          *object);
static void      gimp_paint_core_set_property        (GObject          *object,
                                                      guint             property_id,
                                                      const GValue     *value,
                                                      GParamSpec       *pspec);
static void      gimp_paint_core_get_property        (GObject          *object,
                                                      guint             property_id,
                                                      GValue           *value,
                                                      GParamSpec       *pspec);

static gboolean  gimp_paint_core_real_start          (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *paint_options,
                                                      const GimpCoords *coords,
                                                      GError          **error);
static gboolean  gimp_paint_core_real_pre_paint      (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *options,
                                                      GimpPaintState    paint_state,
                                                      guint32           time);
static void      gimp_paint_core_real_paint          (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *options,
                                                      const GimpCoords *coords,
                                                      GimpPaintState    paint_state,
                                                      guint32           time);
static void      gimp_paint_core_real_post_paint     (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *options,
                                                      GimpPaintState    paint_state,
                                                      guint32           time);
static void      gimp_paint_core_real_interpolate    (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *options,
                                                      guint32           time);
static TempBuf * gimp_paint_core_real_get_paint_area (GimpPaintCore    *core,
                                                      GimpDrawable     *drawable,
                                                      GimpPaintOptions *options,
                                                      const GimpCoords *coords);
static GimpUndo* gimp_paint_core_real_push_undo      (GimpPaintCore    *core,
                                                      GimpImage        *image,
                                                      const gchar      *undo_desc);

static void      paint_mask_to_canvas_tiles          (GimpPaintCore    *core,
                                                      PixelRegion      *paint_maskPR,
                                                      gdouble           paint_opacity);
static void      paint_mask_to_canvas_buf            (GimpPaintCore    *core,
                                                      PixelRegion      *paint_maskPR,
                                                      gdouble           paint_opacity);
static void      canvas_tiles_to_canvas_buf          (GimpPaintCore    *core);


G_DEFINE_TYPE (GimpPaintCore, gimp_paint_core, GIMP_TYPE_OBJECT)

#define parent_class gimp_paint_core_parent_class

static gint global_core_ID = 1;


static void
gimp_paint_core_class_init (GimpPaintCoreClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->finalize     = gimp_paint_core_finalize;
  object_class->set_property = gimp_paint_core_set_property;
  object_class->get_property = gimp_paint_core_get_property;

  klass->start               = gimp_paint_core_real_start;
  klass->pre_paint           = gimp_paint_core_real_pre_paint;
  klass->paint               = gimp_paint_core_real_paint;
  klass->post_paint          = gimp_paint_core_real_post_paint;
  klass->interpolate         = gimp_paint_core_real_interpolate;
  klass->get_paint_area      = gimp_paint_core_real_get_paint_area;
  klass->push_undo           = gimp_paint_core_real_push_undo;

  g_object_class_install_property (object_class, PROP_UNDO_DESC,
                                   g_param_spec_string ("undo-desc", NULL, NULL,
                                                        _("Paint"),
                                                        GIMP_PARAM_READWRITE |
                                                        G_PARAM_CONSTRUCT_ONLY));
}

static void
gimp_paint_core_init (GimpPaintCore *core)
{
  core->ID               = global_core_ID++;

  core->undo_desc        = NULL;

  core->distance         = 0.0;
  core->pixel_dist       = 0.0;
  core->x1               = 0;
  core->y1               = 0;
  core->x2               = 0;
  core->y2               = 0;

  core->use_saved_proj   = FALSE;

  core->undo_tiles       = NULL;
  core->saved_proj_tiles = NULL;
  core->canvas_tiles     = NULL;

  core->orig_buf         = NULL;
  core->orig_proj_buf    = NULL;
  core->canvas_buf       = NULL;
}

static void
gimp_paint_core_finalize (GObject *object)
{
  GimpPaintCore *core = GIMP_PAINT_CORE (object);

  gimp_paint_core_cleanup (core);

  g_free (core->undo_desc);
  core->undo_desc = NULL;

  G_OBJECT_CLASS (parent_class)->finalize (object);
}

static void
gimp_paint_core_set_property (GObject      *object,
                              guint         property_id,
                              const GValue *value,
                              GParamSpec   *pspec)
{
  GimpPaintCore *core = GIMP_PAINT_CORE (object);

  switch (property_id)
    {
    case PROP_UNDO_DESC:
      g_free (core->undo_desc);
      core->undo_desc = g_value_dup_string (value);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_paint_core_get_property (GObject    *object,
                              guint       property_id,
                              GValue     *value,
                              GParamSpec *pspec)
{
  GimpPaintCore *core = GIMP_PAINT_CORE (object);

  switch (property_id)
    {
    case PROP_UNDO_DESC:
      g_value_set_string (value, core->undo_desc);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static gboolean
gimp_paint_core_real_start (GimpPaintCore    *core,
                            GimpDrawable     *drawable,
                            GimpPaintOptions *paint_options,
                            const GimpCoords *coords,
                            GError          **error)
{
  return TRUE;
}

static gboolean
gimp_paint_core_real_pre_paint (GimpPaintCore    *core,
                                GimpDrawable     *drawable,
                                GimpPaintOptions *paint_options,
                                GimpPaintState    paint_state,
                                guint32           time)
{
  return TRUE;
}

static void
gimp_paint_core_real_paint (GimpPaintCore    *core,
                            GimpDrawable     *drawable,
                            GimpPaintOptions *paint_options,
                            const GimpCoords *coords,
                            GimpPaintState    paint_state,
                            guint32           time)
{
}

static void
gimp_paint_core_real_post_paint (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 GimpPaintState    paint_state,
                                 guint32           time)
{
}

static void
gimp_paint_core_real_interpolate (GimpPaintCore    *core,
                                  GimpDrawable     *drawable,
                                  GimpPaintOptions *paint_options,
                                  guint32           time)
{
  gimp_paint_core_paint (core, drawable, paint_options,
                         GIMP_PAINT_STATE_MOTION, time);

  core->last_coords = core->cur_coords;
}

static TempBuf *
gimp_paint_core_real_get_paint_area (GimpPaintCore    *core,
                                     GimpDrawable     *drawable,
                                     GimpPaintOptions *paint_options,
                                     const GimpCoords *coords)
{
  return NULL;
}

static GimpUndo *
gimp_paint_core_real_push_undo (GimpPaintCore *core,
                                GimpImage     *image,
                                const gchar   *undo_desc)
{
  return gimp_image_undo_push (image, GIMP_TYPE_PAINT_CORE_UNDO,
                               GIMP_UNDO_PAINT, undo_desc,
                               0,
                               "paint-core", core,
                               NULL);
}


/*  public functions  */

void
gimp_paint_core_paint (GimpPaintCore    *core,
                       GimpDrawable     *drawable,
                       GimpPaintOptions *paint_options,
                       GimpPaintState    paint_state,
                       guint32           time)
{
  GimpPaintCoreClass *core_class;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_DRAWABLE (drawable));
  g_return_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)));
  g_return_if_fail (GIMP_IS_PAINT_OPTIONS (paint_options));

  core_class = GIMP_PAINT_CORE_GET_CLASS (core);

  if (core_class->pre_paint (core, drawable,
                             paint_options,
                             paint_state, time))
    {

      if (paint_state == GIMP_PAINT_STATE_MOTION)
        {
          /* Save coordinates for gimp_paint_core_interpolate() */
          core->last_paint.x = core->cur_coords.x;
          core->last_paint.y = core->cur_coords.y;
        }

      core_class->paint (core, drawable,
                         paint_options,
                         &core->cur_coords,
                         paint_state, time);

      core_class->post_paint (core, drawable,
                              paint_options,
                              paint_state, time);
    }
}

gboolean
gimp_paint_core_start (GimpPaintCore     *core,
                       GimpDrawable      *drawable,
                       GimpPaintOptions  *paint_options,
                       const GimpCoords  *coords,
                       GError           **error)
{
  GimpItem *item;

  g_return_val_if_fail (GIMP_IS_PAINT_CORE (core), FALSE);
  g_return_val_if_fail (GIMP_IS_DRAWABLE (drawable), FALSE);
  g_return_val_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)), FALSE);
  g_return_val_if_fail (GIMP_IS_PAINT_OPTIONS (paint_options), FALSE);
  g_return_val_if_fail (coords != NULL, FALSE);
  g_return_val_if_fail (error == NULL || *error == NULL, FALSE);

  item = GIMP_ITEM (drawable);

  if (core->stroke_buffer)
    {
      g_array_free (core->stroke_buffer, TRUE);
      core->stroke_buffer = NULL;
    }

  core->stroke_buffer = g_array_sized_new (TRUE, TRUE,
                                           sizeof (GimpCoords),
                                           STROKE_BUFFER_INIT_SIZE);

  core->cur_coords = *coords;

  if (! GIMP_PAINT_CORE_GET_CLASS (core)->start (core, drawable,
                                                 paint_options,
                                                 coords, error))
    {
      return FALSE;
    }

  /*  Allocate the undo structure  */
  if (core->undo_tiles)
    tile_manager_unref (core->undo_tiles);

  core->undo_tiles = tile_manager_new (gimp_item_get_width  (item),
                                       gimp_item_get_height (item),
                                       gimp_drawable_bytes (drawable));

  /*  Allocate the saved proj structure  */
  if (core->saved_proj_tiles)
    tile_manager_unref (core->saved_proj_tiles);

  core->saved_proj_tiles = NULL;

  if (core->use_saved_proj)
    {
      GimpImage    *image    = gimp_item_get_image (item);
      GimpPickable *pickable = GIMP_PICKABLE (gimp_image_get_projection (image));
      TileManager  *tiles    = gimp_pickable_get_tiles (pickable);

      core->saved_proj_tiles = tile_manager_new (tile_manager_width (tiles),
                                                 tile_manager_height (tiles),
                                                 tile_manager_bpp (tiles));
    }

  /*  Allocate the canvas blocks structure  */
  if (core->canvas_tiles)
    tile_manager_unref (core->canvas_tiles);

  core->canvas_tiles = tile_manager_new (gimp_item_get_width  (item),
                                         gimp_item_get_height (item),
                                         1);

  /*  Get the initial undo extents  */

  core->x1 = core->x2 = core->cur_coords.x;
  core->y1 = core->y2 = core->cur_coords.y;

  core->last_paint.x = -1e6;
  core->last_paint.y = -1e6;

  /*  Freeze the drawable preview so that it isn't constantly updated.  */
  gimp_viewable_preview_freeze (GIMP_VIEWABLE (drawable));

  return TRUE;
}

void
gimp_paint_core_finish (GimpPaintCore *core,
                        GimpDrawable  *drawable,
                        gboolean       push_undo)
{
  GimpImage *image;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_DRAWABLE (drawable));
  g_return_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)));

  if (core->stroke_buffer)
    {
      g_array_free (core->stroke_buffer, TRUE);
      core->stroke_buffer = NULL;
    }

  image = gimp_item_get_image (GIMP_ITEM (drawable));

  /*  Determine if any part of the image has been altered--
   *  if nothing has, then just return...
   */
  if ((core->x2 == core->x1) || (core->y2 == core->y1))
    {
      gimp_viewable_preview_thaw (GIMP_VIEWABLE (drawable));
      return;
    }

  if (push_undo)
    {
      gimp_image_undo_group_start (image, GIMP_UNDO_GROUP_PAINT,
                                   core->undo_desc);

      GIMP_PAINT_CORE_GET_CLASS (core)->push_undo (core, image, NULL);

      gimp_drawable_push_undo (drawable, NULL,
                               core->x1, core->y1,
                               core->x2 - core->x1, core->y2 - core->y1,
                               core->undo_tiles,
                               TRUE);

      gimp_image_undo_group_end (image);
    }

  tile_manager_unref (core->undo_tiles);
  core->undo_tiles = NULL;

  if (core->saved_proj_tiles)
    {
      tile_manager_unref (core->saved_proj_tiles);
      core->saved_proj_tiles = NULL;
    }

  gimp_viewable_preview_thaw (GIMP_VIEWABLE (drawable));
}

static void
gimp_paint_core_copy_valid_tiles (TileManager *src_tiles,
                                  TileManager *dest_tiles,
                                  gint         x,
                                  gint         y,
                                  gint         w,
                                  gint         h)
{
  Tile *src_tile;
  gint  i, j;

  for (i = y; i < (y + h); i += (TILE_HEIGHT - (i % TILE_HEIGHT)))
    {
      for (j = x; j < (x + w); j += (TILE_WIDTH - (j % TILE_WIDTH)))
        {
          src_tile = tile_manager_get_tile (src_tiles,
                                            j, i, FALSE, FALSE);

          if (tile_is_valid (src_tile))
            {
              src_tile = tile_manager_get_tile (src_tiles,
                                                j, i, TRUE, FALSE);

              tile_manager_map_tile (dest_tiles, j, i, src_tile);

              tile_release (src_tile, FALSE);
            }
        }
    }
}

void
gimp_paint_core_cancel (GimpPaintCore *core,
                        GimpDrawable  *drawable)
{
  gint x, y;
  gint width, height;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_DRAWABLE (drawable));
  g_return_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)));

  /*  Determine if any part of the image has been altered--
   *  if nothing has, then just return...
   */
  if ((core->x2 == core->x1) || (core->y2 == core->y1))
    return;

  if (gimp_rectangle_intersect (core->x1, core->y1,
                                core->x2 - core->x1,
                                core->y2 - core->y1,
                                0, 0,
                                gimp_item_get_width  (GIMP_ITEM (drawable)),
                                gimp_item_get_height (GIMP_ITEM (drawable)),
                                &x, &y, &width, &height))
    {
      gimp_paint_core_copy_valid_tiles (core->undo_tiles,
                                        gimp_drawable_get_tiles (drawable),
                                        x, y, width, height);
    }

  tile_manager_unref (core->undo_tiles);
  core->undo_tiles = NULL;

  if (core->saved_proj_tiles)
    {
      tile_manager_unref (core->saved_proj_tiles);
      core->saved_proj_tiles = NULL;
    }

  gimp_drawable_update (drawable, x, y, width, height);

  gimp_viewable_preview_thaw (GIMP_VIEWABLE (drawable));
}

void
gimp_paint_core_cleanup (GimpPaintCore *core)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));

  if (core->undo_tiles)
    {
      tile_manager_unref (core->undo_tiles);
      core->undo_tiles = NULL;
    }

  if (core->saved_proj_tiles)
    {
      tile_manager_unref (core->saved_proj_tiles);
      core->saved_proj_tiles = NULL;
    }

  if (core->canvas_tiles)
    {
      tile_manager_unref (core->canvas_tiles);
      core->canvas_tiles = NULL;
    }

  if (core->orig_buf)
    {
      temp_buf_free (core->orig_buf);
      core->orig_buf = NULL;
    }

  if (core->orig_proj_buf)
    {
      temp_buf_free (core->orig_proj_buf);
      core->orig_proj_buf = NULL;
    }

  if (core->canvas_buf)
    {
      temp_buf_free (core->canvas_buf);
      core->canvas_buf = NULL;
    }
}

void
gimp_paint_core_interpolate (GimpPaintCore    *core,
                             GimpDrawable     *drawable,
                             GimpPaintOptions *paint_options,
                             const GimpCoords *coords,
                             guint32           time)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_DRAWABLE (drawable));
  g_return_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)));
  g_return_if_fail (GIMP_IS_PAINT_OPTIONS (paint_options));
  g_return_if_fail (coords != NULL);

  core->cur_coords = *coords;

  GIMP_PAINT_CORE_GET_CLASS (core)->interpolate (core, drawable,
                                                 paint_options, time);
}

void
gimp_paint_core_set_current_coords (GimpPaintCore    *core,
                                    const GimpCoords *coords)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (coords != NULL);

  core->cur_coords = *coords;
}

void
gimp_paint_core_get_current_coords (GimpPaintCore    *core,
                                    GimpCoords       *coords)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (coords != NULL);

  *coords = core->cur_coords;

}

void
gimp_paint_core_set_last_coords (GimpPaintCore    *core,
                                 const GimpCoords *coords)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (coords != NULL);

  core->last_coords = *coords;
}

void
gimp_paint_core_get_last_coords (GimpPaintCore *core,
                                 GimpCoords    *coords)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (coords != NULL);

  *coords = core->last_coords;
}

/**
 * gimp_paint_core_round_line:
 * @core:                 the #GimpPaintCore
 * @options:              the #GimpPaintOptions to use
 * @constrain_15_degrees: the modifier state
 *
 * Adjusts core->last_coords and core_cur_coords in preparation to
 * drawing a straight line. If @center_pixels is TRUE the endpoints
 * get pushed to the center of the pixels. This avoids artefacts
 * for e.g. the hard mode. The rounding of the slope to 15 degree
 * steps if ctrl is pressed happens, as does rounding the start and
 * end coordinates (which may be fractional in high zoom modes) to
 * the center of pixels.
 **/
void
gimp_paint_core_round_line (GimpPaintCore    *core,
                            GimpPaintOptions *paint_options,
                            gboolean          constrain_15_degrees)
{
  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_PAINT_OPTIONS (paint_options));

  if (gimp_paint_options_get_brush_mode (paint_options) == GIMP_BRUSH_HARD)
    {
      core->last_coords.x = floor (core->last_coords.x) + 0.5;
      core->last_coords.y = floor (core->last_coords.y) + 0.5;
      core->cur_coords.x  = floor (core->cur_coords.x ) + 0.5;
      core->cur_coords.y  = floor (core->cur_coords.y ) + 0.5;
    }

  if (constrain_15_degrees)
    gimp_constrain_line (core->last_coords.x, core->last_coords.y,
                         &core->cur_coords.x, &core->cur_coords.y,
                         GIMP_CONSTRAIN_LINE_15_DEGREES);
}


/*  protected functions  */

TempBuf *
gimp_paint_core_get_paint_area (GimpPaintCore    *core,
                                GimpDrawable     *drawable,
                                GimpPaintOptions *paint_options,
                                const GimpCoords *coords)
{
  g_return_val_if_fail (GIMP_IS_PAINT_CORE (core), NULL);
  g_return_val_if_fail (GIMP_IS_DRAWABLE (drawable), NULL);
  g_return_val_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)), NULL);
  g_return_val_if_fail (GIMP_IS_PAINT_OPTIONS (paint_options), NULL);
  g_return_val_if_fail (coords != NULL, NULL);

  return GIMP_PAINT_CORE_GET_CLASS (core)->get_paint_area (core, drawable,
                                                           paint_options,
                                                           coords);
}

TempBuf *
gimp_paint_core_get_orig_image (GimpPaintCore *core,
                                GimpDrawable  *drawable,
                                gint           x,
                                gint           y,
                                gint           width,
                                gint           height)
{
  PixelRegion   srcPR;
  PixelRegion   destPR;
  Tile         *undo_tile;
  gboolean      release_tile;
  gint          h;
  gint          pixelwidth;
  gint          drawable_width;
  gint          drawable_height;
  const guchar *s;
  guchar       *d;
  gpointer      pr;

  g_return_val_if_fail (GIMP_IS_PAINT_CORE (core), NULL);
  g_return_val_if_fail (GIMP_IS_DRAWABLE (drawable), NULL);
  g_return_val_if_fail (core->undo_tiles != NULL, NULL);

  core->orig_buf = temp_buf_resize (core->orig_buf,
                                    gimp_drawable_bytes (drawable),
                                    x, y, width, height);

  drawable_width  = gimp_item_get_width  (GIMP_ITEM (drawable));
  drawable_height = gimp_item_get_height (GIMP_ITEM (drawable));

  gimp_rectangle_intersect (x, y,
                            width, height,
                            0, 0,
                            drawable_width, drawable_height,
                            &x, &y,
                            &width, &height);

  /*  configure the pixel regions  */
  pixel_region_init (&srcPR, gimp_drawable_get_tiles (drawable),
                     x, y, width, height,
                     FALSE);

  pixel_region_init_temp_buf (&destPR, core->orig_buf,
                              x - core->orig_buf->x,
                              y - core->orig_buf->y,
                              width, height);

  for (pr = pixel_regions_register (2, &srcPR, &destPR);
       pr != NULL;
       pr = pixel_regions_process (pr))
    {
      /*  If the undo tile corresponding to this location is valid, use it  */
      undo_tile = tile_manager_get_tile (core->undo_tiles,
                                         srcPR.x, srcPR.y,
                                         FALSE, FALSE);

      if (tile_is_valid (undo_tile))
        {
          release_tile = TRUE;

          undo_tile = tile_manager_get_tile (core->undo_tiles,
                                             srcPR.x, srcPR.y,
                                             TRUE, FALSE);
          s = tile_data_pointer (undo_tile, srcPR.x, srcPR.y);
        }
      else
        {
          release_tile = FALSE;

          s = srcPR.data;
        }

      d = destPR.data;

      pixelwidth = srcPR.w * srcPR.bytes;

      h = srcPR.h;
      while (h --)
        {
          memcpy (d, s, pixelwidth);

          s += srcPR.rowstride;
          d += destPR.rowstride;
        }

      if (release_tile)
        tile_release (undo_tile, FALSE);
    }

  return core->orig_buf;
}

TempBuf *
gimp_paint_core_get_orig_proj (GimpPaintCore *core,
                               GimpPickable  *pickable,
                               gint           x,
                               gint           y,
                               gint           width,
                               gint           height)
{
  TileManager  *src_tiles;
  PixelRegion   srcPR;
  PixelRegion   destPR;
  Tile         *saved_tile;
  gboolean      release_tile;
  gint          h;
  gint          pixelwidth;
  gint          pickable_width;
  gint          pickable_height;
  const guchar *s;
  guchar       *d;
  gpointer      pr;

  g_return_val_if_fail (GIMP_IS_PAINT_CORE (core), NULL);
  g_return_val_if_fail (GIMP_IS_PICKABLE (pickable), NULL);
  g_return_val_if_fail (core->saved_proj_tiles != NULL, NULL);

  core->orig_proj_buf = temp_buf_resize (core->orig_proj_buf,
                                         gimp_pickable_get_bytes (pickable),
                                         x, y, width, height);

  src_tiles = gimp_pickable_get_tiles (pickable);

  pickable_width  = tile_manager_width  (src_tiles);
  pickable_height = tile_manager_height (src_tiles);

  gimp_rectangle_intersect (x, y,
                            width, height,
                            0, 0,
                            pickable_width, pickable_height,
                            &x, &y,
                            &width, &height);

  /*  configure the pixel regions  */
  pixel_region_init (&srcPR, src_tiles,
                     x, y, width, height,
                     FALSE);

  pixel_region_init_temp_buf (&destPR, core->orig_proj_buf,
                              x - core->orig_proj_buf->x,
                              y - core->orig_proj_buf->y,
                              width, height);

  for (pr = pixel_regions_register (2, &srcPR, &destPR);
       pr != NULL;
       pr = pixel_regions_process (pr))
    {
      /*  If the saved tile corresponding to this location is valid, use it  */
      saved_tile = tile_manager_get_tile (core->saved_proj_tiles,
                                          srcPR.x, srcPR.y,
                                          FALSE, FALSE);

      if (tile_is_valid (saved_tile))
        {
          release_tile = TRUE;

          saved_tile = tile_manager_get_tile (core->saved_proj_tiles,
                                              srcPR.x, srcPR.y,
                                              TRUE, FALSE);
          s = tile_data_pointer (saved_tile, srcPR.x, srcPR.y);
        }
      else
        {
          release_tile = FALSE;

          s = srcPR.data;
        }

      d = destPR.data;

      pixelwidth = srcPR.w * srcPR.bytes;

      h = srcPR.h;
      while (h --)
        {
          memcpy (d, s, pixelwidth);

          s += srcPR.rowstride;
          d += destPR.rowstride;
        }

      if (release_tile)
        tile_release (saved_tile, FALSE);
    }

  return core->orig_proj_buf;
}

void
gimp_paint_core_paste (GimpPaintCore            *core,
                       PixelRegion              *paint_maskPR,
                       GimpDrawable             *drawable,
                       gdouble                   paint_opacity,
                       gdouble                   image_opacity,
                       GimpLayerModeEffects      paint_mode,
                       GimpPaintApplicationMode  mode)
{
  TileManager *alt = NULL;
  PixelRegion  srcPR;

  /*  set undo blocks  */
  gimp_paint_core_validate_undo_tiles (core, drawable,
                                       core->canvas_buf->x,
                                       core->canvas_buf->y,
                                       core->canvas_buf->width,
                                       core->canvas_buf->height);

  if (core->use_saved_proj)
    {
      GimpImage      *image      = gimp_item_get_image (GIMP_ITEM (drawable));
      GimpProjection *projection = gimp_image_get_projection (image);
      gint            off_x;
      gint            off_y;
      gint            x, y;
      gint            w, h;

      gimp_item_get_offset (GIMP_ITEM (drawable), &off_x, &off_y);

      if (gimp_rectangle_intersect (core->canvas_buf->x + off_x,
                                    core->canvas_buf->y + off_y,
                                    core->canvas_buf->width,
                                    core->canvas_buf->height,
                                    0, 0,
                                    tile_manager_width (core->saved_proj_tiles),
                                    tile_manager_height (core->saved_proj_tiles),
                                    &x, &y, &w, &h))
        {
          gimp_paint_core_validate_saved_proj_tiles (core,
                                                     GIMP_PICKABLE (projection),
                                                     x, y, w, h);
        }
    }

  /*  If the mode is CONSTANT:
   *   combine the canvas buf, the paint mask to the canvas tiles
   */
  if (mode == GIMP_PAINT_CONSTANT)
    {
      /* Some tools (ink) paint the mask to paint_core->canvas_tiles
       * directly. Don't need to copy it in this case.
       */
      if (paint_maskPR->tiles != core->canvas_tiles)
        {
          /*  initialize any invalid canvas tiles  */
          gimp_paint_core_validate_canvas_tiles (core,
                                                 core->canvas_buf->x,
                                                 core->canvas_buf->y,
                                                 core->canvas_buf->width,
                                                 core->canvas_buf->height);

          paint_mask_to_canvas_tiles (core, paint_maskPR, paint_opacity);
        }

      canvas_tiles_to_canvas_buf (core);
      alt = core->undo_tiles;
    }
  /*  Otherwise:
   *   combine the canvas buf and the paint mask to the canvas buf
   */
  else
    {
      paint_mask_to_canvas_buf (core, paint_maskPR, paint_opacity);
    }

  /*  intialize canvas buf source pixel regions  */
  pixel_region_init_temp_buf (&srcPR, core->canvas_buf,
                              0, 0,
                              core->canvas_buf->width,
                              core->canvas_buf->height);

  /*  apply the paint area to the image  */
  gimp_drawable_apply_region (drawable, &srcPR,
                              FALSE, NULL,
                              image_opacity, paint_mode,
                              alt,  /*  specify an alternative src1  */
                              NULL,
                              core->canvas_buf->x,
                              core->canvas_buf->y);

  /*  Update the undo extents  */
  core->x1 = MIN (core->x1, core->canvas_buf->x);
  core->y1 = MIN (core->y1, core->canvas_buf->y);
  core->x2 = MAX (core->x2, core->canvas_buf->x + core->canvas_buf->width);
  core->y2 = MAX (core->y2, core->canvas_buf->y + core->canvas_buf->height);

  /*  Update the drawable  */
  gimp_drawable_update (drawable,
                        core->canvas_buf->x,
                        core->canvas_buf->y,
                        core->canvas_buf->width,
                        core->canvas_buf->height);
}

/* This works similarly to gimp_paint_core_paste. However, instead of
 * combining the canvas to the paint core drawable using one of the
 * combination modes, it uses a "replace" mode (i.e. transparent
 * pixels in the canvas erase the paint core drawable).

 * When not drawing on alpha-enabled images, it just paints using
 * NORMAL mode.
 */
void
gimp_paint_core_replace (GimpPaintCore            *core,
                         PixelRegion              *paint_maskPR,
                         GimpDrawable             *drawable,
                         gdouble                   paint_opacity,
                         gdouble                   image_opacity,
                         GimpPaintApplicationMode  mode)
{
  PixelRegion  srcPR;

  if (! gimp_drawable_has_alpha (drawable))
    {
      gimp_paint_core_paste (core, paint_maskPR, drawable,
                             paint_opacity,
                             image_opacity, GIMP_NORMAL_MODE,
                             mode);
      return;
    }

  /*  set undo blocks  */
  gimp_paint_core_validate_undo_tiles (core, drawable,
                                       core->canvas_buf->x,
                                       core->canvas_buf->y,
                                       core->canvas_buf->width,
                                       core->canvas_buf->height);

  if (mode == GIMP_PAINT_CONSTANT)
    {
      /* Some tools (ink) paint the mask to paint_core->canvas_tiles
       * directly. Don't need to copy it in this case.
       */
      if (paint_maskPR->tiles != core->canvas_tiles)
        {
          /*  initialize any invalid canvas tiles  */
          gimp_paint_core_validate_canvas_tiles (core,
                                                 core->canvas_buf->x,
                                                 core->canvas_buf->y,
                                                 core->canvas_buf->width,
                                                 core->canvas_buf->height);

          /* combine the paint mask and the canvas tiles */
          paint_mask_to_canvas_tiles (core, paint_maskPR, paint_opacity);

          /* initialize the maskPR from the canvas tiles */
          pixel_region_init (paint_maskPR, core->canvas_tiles,
                             core->canvas_buf->x,
                             core->canvas_buf->y,
                             core->canvas_buf->width,
                             core->canvas_buf->height,
                             FALSE);
        }
    }
  else
    {
      /* The mask is just the paint_maskPR */
    }

  /*  intialize canvas buf source pixel regions  */
  pixel_region_init_temp_buf (&srcPR, core->canvas_buf,
                              0, 0,
                              core->canvas_buf->width,
                              core->canvas_buf->height);

  /*  apply the paint area to the image  */
  gimp_drawable_replace_region (drawable, &srcPR,
                                FALSE, NULL,
                                image_opacity,
                                paint_maskPR,
                                core->canvas_buf->x,
                                core->canvas_buf->y);

  /*  Update the undo extents  */
  core->x1 = MIN (core->x1, core->canvas_buf->x);
  core->y1 = MIN (core->y1, core->canvas_buf->y);
  core->x2 = MAX (core->x2, core->canvas_buf->x + core->canvas_buf->width) ;
  core->y2 = MAX (core->y2, core->canvas_buf->y + core->canvas_buf->height) ;

  /*  Update the drawable  */
  gimp_drawable_update (drawable,
                        core->canvas_buf->x,
                        core->canvas_buf->y,
                        core->canvas_buf->width,
                        core->canvas_buf->height);
}

/**
 * Smooth and store coords in the stroke buffer
 */

void
gimp_paint_core_smooth_coords (GimpPaintCore    *core,
                               GimpPaintOptions *paint_options,
                               GimpCoords       *coords)
{
  GimpSmoothingOptions *smoothing_options = paint_options->smoothing_options;
  GArray               *history           = core->stroke_buffer;

  if (core->stroke_buffer == NULL)
    return; /* Paint core has not initalized yet */

  if (smoothing_options->use_smoothing &&
      smoothing_options->smoothing_quality > 0)
    {
      gint       i;
      guint      length;
      gint       min_index;
      gdouble    gaussian_weight  = 0.0;
      gdouble    gaussian_weight2 = SQR (smoothing_options->smoothing_factor);
      gdouble    velocity_sum     = 0.0;
      gdouble    scale_sum        = 0.0;

      g_array_append_val (history, *coords);

      if (history->len < 2)
        return; /* Just dont bother, nothing to do */

      coords->x = coords->y = 0.0;

      length = MIN (smoothing_options->smoothing_quality, history->len);

      min_index = history->len - length;

      if (gaussian_weight2 != 0.0)
        gaussian_weight = 1 / (sqrt (2 * G_PI) * smoothing_options->smoothing_factor);

      for (i = history->len - 1; i >= min_index; i--)
        {
          gdouble     rate        = 0.0;
          GimpCoords *next_coords = &g_array_index (history,
                                                    GimpCoords, i);

          if (gaussian_weight2 != 0.0)
            {
              /* We use gaussian function with velocity as a window function */
              velocity_sum += next_coords->velocity * 100;
              rate = gaussian_weight * exp (-velocity_sum*velocity_sum / (2 * gaussian_weight2));
            }

          scale_sum += rate;
          coords->x += rate * next_coords->x;
          coords->y += rate * next_coords->y;
        }

      if (scale_sum != 0.0)
        {
          coords->x /= scale_sum;
          coords->y /= scale_sum;
        }

    }

}


static void
canvas_tiles_to_canvas_buf (GimpPaintCore *core)
{
  PixelRegion srcPR;
  PixelRegion maskPR;

  /*  combine the canvas tiles and the canvas buf  */
  pixel_region_init_temp_buf (&srcPR, core->canvas_buf,
                              0, 0,
                              core->canvas_buf->width,
                              core->canvas_buf->height);

  pixel_region_init (&maskPR, core->canvas_tiles,
                     core->canvas_buf->x,
                     core->canvas_buf->y,
                     core->canvas_buf->width,
                     core->canvas_buf->height,
                     FALSE);

  /*  apply the canvas tiles to the canvas buf  */
  apply_mask_to_region (&srcPR, &maskPR, OPAQUE_OPACITY);
}

static void
paint_mask_to_canvas_tiles (GimpPaintCore *core,
                            PixelRegion   *paint_maskPR,
                            gdouble        paint_opacity)
{
  PixelRegion srcPR;

  /*   combine the paint mask and the canvas tiles  */
  pixel_region_init (&srcPR, core->canvas_tiles,
                     core->canvas_buf->x,
                     core->canvas_buf->y,
                     core->canvas_buf->width,
                     core->canvas_buf->height,
                     TRUE);

  /*  combine the mask to the canvas tiles  */
  combine_mask_and_region (&srcPR, paint_maskPR,
                           paint_opacity * 255.999, GIMP_IS_AIRBRUSH (core));
}

static void
paint_mask_to_canvas_buf (GimpPaintCore *core,
                          PixelRegion   *paint_maskPR,
                          gdouble        paint_opacity)
{
  PixelRegion srcPR;

  /*  combine the canvas buf and the paint mask to the canvas buf  */
  pixel_region_init_temp_buf (&srcPR, core->canvas_buf,
                              0, 0,
                              core->canvas_buf->width,
                              core->canvas_buf->height);

  /*  apply the mask  */
  apply_mask_to_region (&srcPR, paint_maskPR, paint_opacity * 255.999);
}

void
gimp_paint_core_validate_undo_tiles (GimpPaintCore *core,
                                     GimpDrawable  *drawable,
                                     gint           x,
                                     gint           y,
                                     gint           w,
                                     gint           h)
{
  gint i, j;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_DRAWABLE (drawable));
  g_return_if_fail (core->undo_tiles != NULL);

  for (i = y; i < (y + h); i += (TILE_HEIGHT - (i % TILE_HEIGHT)))
    {
      for (j = x; j < (x + w); j += (TILE_WIDTH - (j % TILE_WIDTH)))
        {
          Tile *dest_tile = tile_manager_get_tile (core->undo_tiles,
                                                   j, i, FALSE, FALSE);

          if (! tile_is_valid (dest_tile))
            {
              Tile *src_tile =
                tile_manager_get_tile (gimp_drawable_get_tiles (drawable),
                                       j, i, TRUE, FALSE);
              tile_manager_map_tile (core->undo_tiles, j, i, src_tile);
              tile_release (src_tile, FALSE);
            }
        }
    }
}

void
gimp_paint_core_validate_saved_proj_tiles (GimpPaintCore *core,
                                           GimpPickable  *pickable,
                                           gint           x,
                                           gint           y,
                                           gint           w,
                                           gint           h)
{
  gint i, j;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (GIMP_IS_PICKABLE (pickable));
  g_return_if_fail (core->saved_proj_tiles != NULL);

  for (i = y; i < (y + h); i += (TILE_HEIGHT - (i % TILE_HEIGHT)))
    {
      for (j = x; j < (x + w); j += (TILE_WIDTH - (j % TILE_WIDTH)))
        {
          Tile *dest_tile = tile_manager_get_tile (core->saved_proj_tiles,
                                                   j, i, FALSE, FALSE);

          if (! tile_is_valid (dest_tile))
            {
              Tile *src_tile =
                tile_manager_get_tile (gimp_pickable_get_tiles (pickable),
                                       j, i, TRUE, FALSE);

              tile_manager_map_tile (core->saved_proj_tiles, j, i, src_tile);
              tile_release (src_tile, FALSE);
            }
        }
    }
}

void
gimp_paint_core_validate_canvas_tiles (GimpPaintCore *core,
                                       gint           x,
                                       gint           y,
                                       gint           w,
                                       gint           h)
{
  gint i, j;

  g_return_if_fail (GIMP_IS_PAINT_CORE (core));
  g_return_if_fail (core->canvas_tiles != NULL);

  for (i = y; i < (y + h); i += (TILE_HEIGHT - (i % TILE_HEIGHT)))
    {
      for (j = x; j < (x + w); j += (TILE_WIDTH - (j % TILE_WIDTH)))
        {
          Tile *tile = tile_manager_get_tile (core->canvas_tiles, j, i,
                                              FALSE, FALSE);

          if (! tile_is_valid (tile))
            {
              tile = tile_manager_get_tile (core->canvas_tiles, j, i,
                                            TRUE, TRUE);
              memset (tile_data_pointer (tile, 0, 0), 0, tile_size (tile));
              tile_release (tile, TRUE);
            }
        }
    }
}

