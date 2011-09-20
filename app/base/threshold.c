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

#include <glib-object.h>

#include "base-types.h"

#include "pixel-region.h"
#include "threshold.h"


void
threshold (Threshold   *tr,
           PixelRegion *srcPR,
           PixelRegion *destPR)
{
  const guchar *src, *s;
  guchar       *dest, *d;
  gboolean      has_alpha;
  gint          alpha;
  gint          w, h, b;
  gint          value;

  h         = srcPR->h;
  src       = srcPR->data;
  dest      = destPR->data;
  has_alpha = pixel_region_has_alpha (srcPR);
  alpha     = has_alpha ? srcPR->bytes - 1 : srcPR->bytes;

  while (h--)
    {
      w = srcPR->w;
      s = src;
      d = dest;

      while (w--)
        {
          if (tr->color)
            {
              value = MAX (s[RED], s[GREEN]);
              value = MAX (value, s[BLUE]);

              value = (value >= tr->low_threshold &&
                       value <= tr->high_threshold ) ? 255 : 0;
            }
          else
            {
              value = (s[GRAY] >= tr->low_threshold &&
                       s[GRAY] <= tr->high_threshold) ? 255 : 0;
            }

          for (b = 0; b < alpha; b++)
            d[b] = value;

          if (has_alpha)
            d[alpha] = s[alpha];

          s += srcPR->bytes;
          d += destPR->bytes;
        }

      src  += srcPR->rowstride;
      dest += destPR->rowstride;
    }
}
