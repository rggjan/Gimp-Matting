/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * gimpbrush-scale.h
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

#ifndef __GIMP_BRUSH_SCALE_H__
#define __GIMP_BRUSH_SCALE_H__


/*  virtual functions of GimpBrush, don't call directly  */

void      gimp_brush_real_transform_size   (GimpBrush *brush,
                                            gdouble    scale,
                                            gdouble    angle,
                                            gint      *scaled_width,
                                            gint      *scaled_height);
TempBuf * gimp_brush_real_transform_mask   (GimpBrush *brush,
                                            gdouble    scale,
                                            gdouble    angle);
TempBuf * gimp_brush_real_transform_pixmap (GimpBrush *brush,
                                            gdouble    scale,
                                            gdouble    angle);


#endif  /*  __GIMP_BRUSH_SCALE_H__  */