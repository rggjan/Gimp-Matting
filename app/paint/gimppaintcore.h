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

#ifndef __GIMP_PAINT_CORE_H__
#define __GIMP_PAINT_CORE_H__


#include "core/gimpobject.h"


#define GIMP_TYPE_PAINT_CORE            (gimp_paint_core_get_type ())
#define GIMP_PAINT_CORE(obj)            (G_TYPE_CHECK_INSTANCE_CAST ((obj), GIMP_TYPE_PAINT_CORE, GimpPaintCore))
#define GIMP_PAINT_CORE_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST ((klass), GIMP_TYPE_PAINT_CORE, GimpPaintCoreClass))
#define GIMP_IS_PAINT_CORE(obj)         (G_TYPE_CHECK_INSTANCE_TYPE ((obj), GIMP_TYPE_PAINT_CORE))
#define GIMP_IS_PAINT_CORE_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE ((klass), GIMP_TYPE_PAINT_CORE))
#define GIMP_PAINT_CORE_GET_CLASS(obj)  (G_TYPE_INSTANCE_GET_CLASS ((obj), GIMP_TYPE_PAINT_CORE, GimpPaintCoreClass))


typedef struct _GimpPaintCoreClass GimpPaintCoreClass;

struct _GimpPaintCore
{
  GimpObject   parent_instance;

  gint         ID;               /*  unique instance ID                  */

  gchar       *undo_desc;        /*  undo description                    */

  GimpCoords   start_coords;     /*  starting coords (for undo only)     */

  GimpCoords   cur_coords;       /*  current coords                      */
  GimpCoords   last_coords;      /*  last coords                         */

  GimpVector2  last_paint;       /*  last point that was painted         */

  gdouble      distance;         /*  distance traveled by brush          */
  gdouble      pixel_dist;       /*  distance in pixels                  */

  gint         x1, y1;           /*  undo extents in image coords        */
  gint         x2, y2;           /*  undo extents in image coords        */

  gboolean     use_saved_proj;   /*  keep the unmodified proj around     */

  TileManager *undo_tiles;       /*  tiles which have been modified      */
  TileManager *saved_proj_tiles; /*  proj tiles which have been modified */
  TileManager *canvas_tiles;     /*  the buffer to paint the mask to     */

  TempBuf     *orig_buf;         /*  the unmodified drawable pixels      */
  TempBuf     *orig_proj_buf;    /*  the unmodified projection pixels    */
  TempBuf     *canvas_buf;       /*  the buffer to paint pixels to       */

  GArray      *stroke_buffer;
};

struct _GimpPaintCoreClass
{
  GimpObjectClass  parent_class;

  /*  virtual functions  */
  gboolean   (* start)          (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 const GimpCoords *coords,
                                 GError          **error);

  gboolean   (* pre_paint)      (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 GimpPaintState    paint_state,
                                 guint32           time);
  void       (* paint)          (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 const GimpCoords *coords,
                                 GimpPaintState    paint_state,
                                 guint32           time);
  void       (* post_paint)     (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 GimpPaintState    paint_state,
                                 guint32           time);

  void       (* interpolate)    (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 guint32           time);

  TempBuf  * (* get_paint_area) (GimpPaintCore    *core,
                                 GimpDrawable     *drawable,
                                 GimpPaintOptions *paint_options,
                                 const GimpCoords *coords);

  GimpUndo * (* push_undo)      (GimpPaintCore    *core,
                                 GimpImage        *image,
                                 const gchar      *undo_desc);
};


GType     gimp_paint_core_get_type                  (void) G_GNUC_CONST;

void      gimp_paint_core_paint                     (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     GimpPaintOptions *paint_options,
                                                     GimpPaintState    state,
                                                     guint32           time);

gboolean  gimp_paint_core_start                     (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     GimpPaintOptions *paint_options,
                                                     const GimpCoords *coords,
                                                     GError          **error);
void      gimp_paint_core_finish                    (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     gboolean          push_undo);
void      gimp_paint_core_cancel                    (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable);
void      gimp_paint_core_cleanup                   (GimpPaintCore    *core);

void      gimp_paint_core_interpolate               (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     GimpPaintOptions *paint_options,
                                                     const GimpCoords *coords,
                                                     guint32           time);

void      gimp_paint_core_set_current_coords        (GimpPaintCore    *core,
                                                     const GimpCoords *coords);
void      gimp_paint_core_get_current_coords        (GimpPaintCore    *core,
                                                     GimpCoords       *coords);

void      gimp_paint_core_set_last_coords           (GimpPaintCore    *core,
                                                     const GimpCoords *coords);
void      gimp_paint_core_get_last_coords           (GimpPaintCore    *core,
                                                     GimpCoords       *coords);

void      gimp_paint_core_round_line                (GimpPaintCore    *core,
                                                     GimpPaintOptions *options,
                                                     gboolean          constrain_15_degrees);


/*  protected functions  */

TempBuf * gimp_paint_core_get_paint_area            (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     GimpPaintOptions *options,
                                                     const GimpCoords *coords);
TempBuf * gimp_paint_core_get_orig_image            (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     gint              x,
                                                     gint              y,
                                                     gint              width,
                                                     gint              height);
TempBuf * gimp_paint_core_get_orig_proj             (GimpPaintCore    *core,
                                                     GimpPickable     *pickable,
                                                     gint              x,
                                                     gint              y,
                                                     gint              width,
                                                     gint              height);

void      gimp_paint_core_paste             (GimpPaintCore            *core,
                                             PixelRegion              *paint_maskPR,
                                             GimpDrawable             *drawable,
                                             gdouble                   paint_opacity,
                                             gdouble                   image_opacity,
                                             GimpLayerModeEffects      paint_mode,
                                             GimpPaintApplicationMode  mode);
void      gimp_paint_core_replace           (GimpPaintCore            *core,
                                             PixelRegion              *paint_maskPR,
                                             GimpDrawable             *drawable,
                                             gdouble                   paint_opacity,
                                             gdouble                   image_opacity,
                                             GimpPaintApplicationMode  mode);

void      gimp_paint_core_validate_undo_tiles       (GimpPaintCore    *core,
                                                     GimpDrawable     *drawable,
                                                     gint              x,
                                                     gint              y,
                                                     gint              w,
                                                     gint              h);
void      gimp_paint_core_validate_saved_proj_tiles (GimpPaintCore    *core,
                                                     GimpPickable     *pickable,
                                                     gint              x,
                                                     gint              y,
                                                     gint              w,
                                                     gint              h);
void      gimp_paint_core_validate_canvas_tiles     (GimpPaintCore    *core,
                                                     gint              x,
                                                     gint              y,
                                                     gint              w,
                                                     gint              h);

void      gimp_paint_core_smooth_coords             (GimpPaintCore    *core,
                                                     GimpPaintOptions *paint_options,
                                                     GimpCoords       *coords);


#endif  /*  __GIMP_PAINT_CORE_H__  */
