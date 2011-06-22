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

#ifndef  __GIMP_DRAWABLE_FOREGROUND_EXTRACT_H__
#define  __GIMP_DRAWABLE_FOREGROUND_EXTRACT_H__


/*  general API (as seen from the PDB)  */

void       gimp_drawable_foreground_extract (GimpDrawable              *drawable,
                                             GimpForegroundExtractMode  mode,
                                             GimpDrawable              *mask,
                                             GimpProgress              *progress);

/*  SIOX specific API  */

MattingState * gimp_drawable_foreground_extract_siox_init   (GimpDrawable *drawable,
                                                          gint          x,
                                                          gint          y,
                                                          gint          width,
                                                          gint          height);
void        gimp_drawable_foreground_extract_siox  (GimpDrawable       *mask,
						    GimpLayer          *result_layer,
						    GimpLayer          *working_layer,
                                                    MattingState          *state,
                                                    SioxRefinementType  refinemane,
                                                    gfloat start_percentage,
                                                    const gdouble       sensitivity[3],
                                                    gboolean            multiblob,
                                                    GimpProgress       *progress);
void        gimp_drawable_foreground_extract_siox_done (MattingState      *state);


#endif  /*  __GIMP_DRAWABLE_FOREGROUND_EXTRACT_H__  */
