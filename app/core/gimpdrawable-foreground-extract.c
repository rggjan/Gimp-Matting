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

#include <gegl.h>

#include "libgimpbase/gimpbase.h"

#include "core-types.h"

#include "base/pixel-region.h"
#include "base/siox.h"
#include "base/tile-manager.h"


#include "gimpchannel.h"
#include "gimpdrawable.h"
#include "gimpdrawable-foreground-extract.h"
#include "gimpimage.h"
#include "gimpprogress.h"

#include "gimp-intl.h"

/*  public functions  */

void
gimp_drawable_foreground_extract (GimpDrawable              *drawable,
                                  GimpForegroundExtractMode  mode,
                                  GimpDrawable              *mask,
                                  GimpProgress              *progress)
{ 
  g_return_if_fail (GIMP_IS_DRAWABLE (mask));
  g_return_if_fail (mode == GIMP_FOREGROUND_EXTRACT_SIOX);

  /*state =
    gimp_drawable_foreground_extract_siox_init (drawable,
                                                0, 0,
                                                gimp_item_get_width  (GIMP_ITEM (mask)),
                                                gimp_item_get_height (GIMP_ITEM (mask)));

  if (state)
    {
      gimp_drawable_foreground_extract_siox (mask, state,
                                             SIOX_REFINEMENT_RECALCULATE,
                                             SIOX_DEFAULT_SMOOTHNESS,
                                             sensitivity,
                                             FALSE,
                                             progress);

      gimp_drawable_foreground_extract_siox_done (state);
    }
  */
  

}

SioxState *
gimp_drawable_foreground_extract_siox_init (GimpDrawable *drawable,
                                            gint          x,
                                            gint          y,
                                            gint          width,
                                            gint          height)
{
  const guchar *colormap = NULL;
  gboolean      intersect;
  gint          offset_x;
  gint          offset_y;
  //SioxState    *state;
  

  g_return_val_if_fail (GIMP_IS_DRAWABLE (drawable), NULL);
  g_return_val_if_fail (gimp_item_is_attached (GIMP_ITEM (drawable)), NULL);

  //tile_manager_get
  //tile_manager_get_tile()
  //pixel_region_get_col()
  //gimp_tile_cache_size()
  //gimp_pixel_rgn_get_pixel()
  //gimp_pixel_rgn_get_pixel
/*
  
    PixelRegion region;
  gpointer    pr;
  gint        row, col;
  
    pixel_region_init (&region, gimp_drawable_get_tiles(drawable), 0, 0, 300, 300, TRUE);

  g_printerr ("fgextract step #2 -> %d clusters\n", region.bytes);
  g_print("This is a test\n");
  printf("abc");
  for (pr = pixel_regions_register (1, &region);
       pr != NULL; pr = pixel_regions_process (pr))
    {
      guchar *data = region.data;

      for (row = 0; row < region.h; row++)
        {
          guchar *d = data;

          // everything that fits the mask is in the image
          for (col = 0; col < region.w; col++, d+=3)
            {
              d[0]=0;
              d[1]=255;
              d[2]=0;
            }

          data += region.rowstride;
        }
    }
  
  
  gimp_drawable_update (drawable, 0, 0, 300, 300);
*/
  
  if (gimp_drawable_is_indexed (drawable))
    colormap = gimp_drawable_get_colormap (drawable);

  gimp_item_get_offset (GIMP_ITEM (drawable), &offset_x, &offset_y);

  intersect = gimp_rectangle_intersect (offset_x, offset_y,
                                        gimp_item_get_width  (GIMP_ITEM (drawable)),
                                        gimp_item_get_height (GIMP_ITEM (drawable)),
                                        x, y, width, height,
                                        &x, &y, &width, &height);


  /* FIXME:
   * Clear the mask outside the rectangle that we are working on?
   */

  if (! intersect)
    return NULL;

  return siox_init (gimp_drawable_get_tiles (drawable), colormap,
                    offset_x, offset_y,
                    x, y, width, height);
}

void
gimp_drawable_foreground_extract_siox (GimpDrawable       *mask,
                                       GimpLayer          *result_layer,
                                       SioxState          *state,
                                       SioxRefinementType  refinement,
                                       gint                smoothness,
                                       const gdouble       sensitivity[3],
                                       gboolean            multiblob,
                                       GimpProgress       *progress)
{
  gint x1, y1;
  gint x2, y2;

  g_return_if_fail (GIMP_IS_DRAWABLE (mask));
  g_return_if_fail (gimp_drawable_bytes (mask) == 1);

  g_return_if_fail (state != NULL);

  g_return_if_fail (progress == NULL || GIMP_IS_PROGRESS (progress));

  if (progress)
    gimp_progress_start (progress, _("Foreground Extraction"), FALSE);

  if (GIMP_IS_CHANNEL (mask))
    {
      gimp_channel_bounds (GIMP_CHANNEL (mask), &x1, &y1, &x2, &y2);
    }
  else
    {
      x1 = 0;
      y1 = 0;
      x2 = gimp_item_get_width  (GIMP_ITEM (mask));
      y2 = gimp_item_get_height (GIMP_ITEM (mask));
    }


  siox_foreground_extract (state, refinement,
                           gimp_drawable_get_tiles (mask), x1, y1, x2, y2,
                           smoothness, sensitivity, multiblob,
                           (SioxProgressFunc) gimp_progress_set_value,
                           progress, gimp_drawable_get_tiles(GIMP_DRAWABLE(result_layer)));


  if (progress)
    gimp_progress_end (progress);

/* start */


/*
  pixel_region_init (&region, gimp_drawable_get_tiles(mask), 0, 0, 300, 300, TRUE);

  for (pr = pixel_regions_register (1, &region);
       pr != NULL; pr = pixel_regions_process (pr))
    {
      guchar *data = region.data;

      for (row = 0; row < region.h; row++)
        {
          guchar *d = data;

          // everything that fits the mask is in the image
          for (col = 0; col < region.w; col++, d++)
            {
              *d = col*2;
            }

          data += region.rowstride;
        }
    }

  
  g_print("This is a test 2\n");
  printf("abc");*/
/* end */

  gimp_drawable_update (GIMP_DRAWABLE(result_layer), x1, y1, x2, y2);
  gimp_drawable_update (mask, x1, y1, x2, y2);
}

void
gimp_drawable_foreground_extract_siox_done (SioxState *state)
{
  g_return_if_fail (state != NULL);

  siox_done (state);
}
