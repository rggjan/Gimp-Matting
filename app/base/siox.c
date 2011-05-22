/*
 * SIOX: Simple Interactive Object Extraction
 *
 * For algorithm documentation refer to:
 * G. Friedland, K. Jantz, L. Knipping, R. Rojas:
 * "Image Segmentation by Uniform Color Clustering
 *  -- Approach and Benchmark Results",
 * Technical Report B-05-07, Department of Computer Science,
 * Freie Universitaet Berlin, June 2005.
 * http://www.inf.fu-berlin.de/inst/pubs/tr-b-05-07.pdf
 *
 * See http://www.siox.org/ for more information.
 *
 * Algorithm idea by Gerald Friedland.
 * This implementation is Copyright (C) 2005
 * by Gerald Friedland <fland@inf.fu-berlin.de>
 * and Kristian Jantz <jantz@inf.fu-berlin.de>
 * and Tobias Lenz <tlenz@inf.fu-berlin.de>.
 *
 * Adapted for GIMP by Sven Neumann <sven@gimp.org>
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

#include "libgimpbase/gimpbase.h"
#include "libgimpmath/gimpmath.h"

#include "base-types.h"
#include "core/core-types.h"

#include "paint-funcs/paint-funcs.h"

#include "cpercep.h"
#include "pixel-region.h"
#include "tile.h"
#include "tile-manager.h"
#include "siox.h"


/* Thresholds in the mask:
 *   pixels < SIOX_LOW  are known background
 *   pixels > SIOX_HIGH are known foreground
 */
#define SIOX_LOW  1
#define SIOX_HIGH 254

/* When clustering:
 *   use LAB for color images (3 dims),
 *   use L only for grayscale images (1 dim)
 */
#define SIOX_COLOR_DIMS 3
#define SIOX_GRAY_DIMS  1

/* For findmaxblob:
 * Find all blobs with area not smaller than sizefactor of biggest blob
 * CHECKME: Should the user decide this with a slider?
 */
#define MULTIBLOB_DEFAULT_SIZEFACTOR 4
#define MULTIBLOB_ONE_BLOB_ONLY      0

#define BIG_CACHE_CHANNELS 4
#define BIG_CACHE_W 64*3*BIG_CACHE_CHANNELS
#define BIG_CACHE_H 64*3
#define BIG_CACHE_SIZE BIG_CACHE_W*BIG_CACHE_H

#define GET_PIXEL(big_cache, x, y, color) big_cache[(64+y)*BIG_CACHE_W\
                                           +(64+x)*BIG_CACHE_CHANNELS+color]

/* #define SIOX_DEBUG  */

typedef struct
{
  gfloat l;
  gfloat a;
  gfloat b;
  gint   cardinality;
} lab;


/* A struct that holds SIOX current state */
struct _SioxState
{
  TileManager  *pixels;
  const guchar *colormap;
  gint          bpp;
  gint          offset_x;
  gint          offset_y;
  gint          x;
  gint          y;
  gint          width;
  gint          height;
  GHashTable   *cache;
  lab          *bgsig;
  lab          *fgsig;
  gint          bgsiglen;
  gint          fgsiglen;
  gint          xsbpp;
};

/* A struct that holds the classification result */
typedef struct
{
  gfloat bgdist;
  gfloat fgdist;
} classresult;


static void
siox_cache_entry_free (gpointer entry)
{
  g_slice_free (classresult, entry);
}

/* Progressbar update callback */
static inline void
siox_progress_update (SioxProgressFunc  progress_callback,
                      gpointer          progress_data,
                      gdouble           value)
{
  if (progress_data)
    progress_callback (progress_data, value);
}

/*  assumes that lab starts with an array of floats (l,a,b)  */
#define CURRENT_VALUE(points, i, dim) (((const gfloat *) (points + i))[dim])

/* squared euclidean distance */
static inline float
euklid (const lab *p,
        const lab *q)
{
  return (SQR (p->l - q->l) + SQR (p->a - q->a) + SQR (p->b - q->b));
}

/* a struct that contains information about a blob */
struct blob
{
  gint     seedx;
  gint     seedy;
  gint     size;
  gboolean mustkeep;
};

/* Creates a key for the hashtable from a given pixel color value */
static inline gint
create_key (const guchar *src,
            gint          bpp,
            const guchar *colormap)
{
  switch (bpp)
    {
    case 3:                     /* RGB  */
    case 4:                     /* RGBA */
      return (src[RED] << 16 | src[GREEN] << 8 | src[BLUE]);
    case 2:
    case 1:
      if (colormap)             /* INDEXED(A) */
        {
          gint i = *src * 3;

          return (colormap[i + RED]   << 16 |
                  colormap[i + GREEN] << 8  |
                  colormap[i + BLUE]);
        }
      else                      /* GRAY(A) */
        {
          return *src;
        }
    default:
      return 0;
    }
}

/**
 * siox_init:
 * @pixels:   the tiles to extract the foreground from
 * @colormap: colormap in case @pixels are indexed, %NULL otherwise
 * @offset_x: horizontal offset of @pixels with respect to the @mask
 * @offset_y: vertical offset of @pixels with respect to the @mask
 * @x:        horizontal offset into the mask
 * @y:        vertical offset into the mask
 * @width:    width of working area on mask
 * @height:   height of working area on mask
 *
 * Initializes the SIOX segmentator.
 * Creates and returns a SioxState struct that has to be passed to all
 * function calls of this module as it maintaines the state.
 *
'* Returns: a new siox state structure.
 */
SioxState *
siox_init (TileManager  *pixels,
           const guchar *colormap,
           gint          offset_x,
           gint          offset_y,
           gint          x,
           gint          y,
           gint          width,
           gint          height)
{
  SioxState *state;

  g_return_val_if_fail (pixels != NULL, NULL);
  g_return_val_if_fail (x >= 0, NULL);
  g_return_val_if_fail (y >= 0, NULL);

  state = g_slice_new (SioxState);

  state->pixels   = pixels;
  state->colormap = colormap;
  state->offset_x = offset_x;
  state->offset_y = offset_y;
  state->x        = x;
  state->y        = y;
  state->width    = width;
  state->height   = height;
  state->bgsig    = NULL;
  state->fgsig    = NULL;
  state->bgsiglen = 0;
  state->fgsiglen = 0;
  state->bpp      = tile_manager_bpp (pixels);

  state->cache = g_hash_table_new_full (g_direct_hash,
                                        NULL, NULL,
                                        (GDestroyNotify) siox_cache_entry_free);

  cpercep_init ();

#ifdef SIOX_DEBUG
  g_printerr ("siox.c: siox_init (bpp=%d, "
              "x=%d, y=%d, width=%d, height=%d, offset_x=%d, offset_y=%d)\n",
              state->bpp, x, y, width, height, offset_x, offset_y);
#endif

  return state;
}

void
load_big_cache (TileManager *source, guchar *big_cache, gint tx, gint ty)
{
  gint    xdiff;
  gint    ydiff;
  gint    x, y;
  gint    bx, by;
  gint    width_tile, height_tile;
  
  Tile   *src_tile;
  guchar *pointer;

  g_return_if_fail (tile_manager_bpp(source) == 4);

  for (xdiff = -1; xdiff <= 1; xdiff++)
    {
      for (ydiff = -1; ydiff <= 1; ydiff++)
        {
          // No idea why we have tx + ydiff here, but otherwise it doesn't work!
          src_tile = tile_manager_get_at (source, tx + ydiff, ty + xdiff, TRUE, FALSE);
          width_tile = 0;
          height_tile = 0;

          if (src_tile)
            {
              pointer = tile_data_pointer (src_tile, 0, 0);
              width_tile = tile_ewidth (src_tile);
              height_tile = tile_eheight (src_tile);

              for (x = 0; x < width_tile; x++) // TODO what if src tile not that big?
                {
                  bx = (xdiff + 1)*64 + x;
                  for (y = 0; y < height_tile; y++)
                    {
                      by = (ydiff + 1)*64 + y;

                      big_cache[by * BIG_CACHE_W + bx * 4] = *pointer;
                      big_cache[by * BIG_CACHE_W + bx * 4 + 1] = *(pointer + 1);
                      big_cache[by * BIG_CACHE_W + bx * 4 + 2] = *(pointer + 2);
                      big_cache[by * BIG_CACHE_W + bx * 4 + 3] = *(pointer + 3);

                      pointer += 4;
                    }
                }
              
              tile_release (src_tile, FALSE);
            }

          for (x = width_tile; x < 64; x++)
            {
              bx = (xdiff + 1)*64 + x;
              for (y = height_tile; y < 64; y++)
                {
                  by = (ydiff + 1)*64 + y;

                  big_cache[by * BIG_CACHE_W + bx * 4] = 0;
                  big_cache[by * BIG_CACHE_W + bx * 4 + 1] = 0;
                  big_cache[by * BIG_CACHE_W + bx * 4 + 2] = 0;
                  big_cache[by * BIG_CACHE_W + bx * 4 + 3] = 128;
                }
            }
        }
    }
}

void
initialize_new_layer (TileManager* source_layer,
                      TileManager* destination_layer,
                      TileManager* mask_layer)
{
  PixelRegion src, dest, mask;
  PixelRegionIterator *pr;
  gint row, col;
  int width, height;
  
  width = tile_manager_width (source_layer);
  height = tile_manager_height(source_layer);

  pixel_region_init (&src, source_layer, 0, 0, width, height, FALSE);
  pixel_region_init (&dest, destination_layer, 0, 0, width, height, TRUE);
  pixel_region_init (&mask, mask_layer, 0, 0, width, height, FALSE);

  g_return_if_fail (src.bytes == 3 && dest.bytes == 4 && mask.bytes == 1); // TODO check if indexed etc...

  for (pr = pixel_regions_register (3, &src, &dest, &mask);
          pr != NULL;
          pr = pixel_regions_process (pr))
    {
      const guchar *mask_data = mask.data;
      const guchar *src_data = src.data;
      guchar *dest_data = dest.data;

      for (row = 0; row < src.h; row++)
        {
          const guchar *m = mask_data;
          const guchar *s = src_data;
          guchar *d = dest_data;

          for (col = 0; col < src.w; col++, s += src.bytes, d += dest.bytes, ++m) // TODO check if 4 is ok
            {
              d[0] = s[0];
              d[1] = s[1];
              d[2] = s[2];
              if (m[0] == MATTING_USER_FOREGROUND)
                d[3] = 255;
              else if (m[0] == MATTING_USER_BACKGROUND)
                d[3] = 0;
              else
                d[3] = 128;
            }

          src_data += src.rowstride;
          dest_data += dest.rowstride;
          mask_data += mask.rowstride;
        }
    }
}

/**
 * siox_foreground_extract:
 * @state:       current state struct as constructed by siox_init
 * @refinement:  #SioxRefinementType
 * @mask:        a mask indicating sure foreground (255), sure background (0)
 *               and undecided regions ([1..254]).
 * @x1:          region of interest
 * @y1:          region of interest
 * @x2:          region of interest
 * @y2:          region of interest
 * @sensitivity: a double array with three entries specifing the accuracy,
 *               a good value is: { 0.64, 1.28, 2.56 }
 * @smoothness:  boundary smoothness (a good value is 3)
 * @multiblob:   allow multiple blobs (true) or only one (false)
 * @progress_callback: a progress callback
 * @progress_data: data passed to @progress_callback
 *
 * Writes the resulting segmentation into @mask. The region of
 * interest as specified using @x1, @y1, @x2 and @y2 defines the
 * bounding box of the background and undecided areas. No changes to
 * the mask are done outside this rectangle.
 */
void
siox_foreground_extract (SioxState          *state,
                         SioxRefinementType  refinement,
                         TileManager        *mask,
                         gint                x1,
                         gint                y1,
                         gint                x2,
                         gint                y2,
                         gint                smoothness,
                         const gdouble       sensitivity[3],
                         gboolean            multiblob,
                         SioxProgressFunc    progress_callback,
                         gpointer            progress_data,
                         TileManager        *result_layer,
                         TileManager        *working_layer)
{
  gint         width, height;
  guchar      *big_cache;
  gint         tiles_x, tiles_y;
  gint         i;

  Tile        *tile;

  gint         tx, ty, x, y;
  guchar      *pointer;

  gint         radius, n;
  guchar       value;
  
  g_return_if_fail (state != NULL);
  g_return_if_fail (mask != NULL && tile_manager_bpp (mask) == 1);
  g_return_if_fail (x1 >= 0);
  g_return_if_fail (x2 > x1 && x2 <= tile_manager_width (mask));
  g_return_if_fail (y1 >= 0);
  g_return_if_fail (y2 > y1 && y2 <= tile_manager_height (mask));
  g_return_if_fail (smoothness >= 0);
  g_return_if_fail (progress_data == NULL || progress_callback != NULL);

  width  = state->width;
  height = state->height;

  g_return_if_fail (TILE_WIDTH == 64 && TILE_HEIGHT == 64);

  g_return_if_fail (tile_manager_bpp (state->pixels) == 3);
  g_return_if_fail (tile_manager_bpp (result_layer) == 4);

  big_cache = g_malloc (BIG_CACHE_SIZE);

  tiles_x = tile_manager_tiles_per_col (state->pixels);
  tiles_y = tile_manager_tiles_per_row (state->pixels);
  //tiles_x = state->pixels->ntile_cols;
  //tiles_y = state->pixels->ntile_rows;

  initialize_new_layer(state->pixels, working_layer, mask);

  for (tx = 0; tx < tiles_x-1; tx++)
    {
      for (ty = 0; ty < tiles_y-1; ty++)
        {
          load_big_cache(working_layer, big_cache, tx, ty);

          tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);
          pointer = tile_data_pointer (tile, 0, 0);

          for (x=0; x<64; x++)
            {
              for (y=0; y<64; y++, pointer += 4)
                {
                  if (GET_PIXEL (big_cache, x, y, 3) != 128)
                    {
                      pointer[0] = GET_PIXEL (big_cache, x, y, 0);
                      pointer[1] = GET_PIXEL (big_cache, x, y, 1);
                      pointer[2] = GET_PIXEL (big_cache, x, y, 2);
                      pointer[3] = GET_PIXEL (big_cache, x, y, 3);

                      continue;
                    }

                  for (radius = 0; radius <= SEARCH_RADIUS; radius++)
                    {
                      for (n = -radius; n < radius; n++)
                        {
                          value = GET_PIXEL (big_cache, x + n, y + radius, 3);
                          if (value != 128)
                            {
                              pointer[0] = 255;
                              pointer[1] = 0;
                              pointer[2] = 0;
                              pointer[3] = 255;
                            }

                          value = GET_PIXEL (big_cache, x + n, y - radius, 3);
                          if (value != 128)
                            {
                              pointer[0] = 255;
                              pointer[1] = 0;
                              pointer[2] = 0;
                              pointer[3] = 255;
                            }

                          value = GET_PIXEL (big_cache, x + radius, y + n, 3);
                          if (value != 128)
                            {
                              pointer[0] = 255;
                              pointer[1] = 0;
                              pointer[2] = 0;
                              pointer[3] = 255;
                            }

                          value = GET_PIXEL (big_cache, x - radius, y + n, 3);
                          if (value != 128)
                            {
                              pointer[0] = 255;
                              pointer[1] = 0;
                              pointer[2] = 0;
                              pointer[3] = 255;
                            }
                        }
                    }
                }
            }

            tile_release (tile, TRUE);
        }
    }

  g_free(big_cache);
}

/*
 // Could set to FALSE, TRUE, but then we get a warning...
          dst_tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);
          g_return_if_fail (dst_tile);
          pointer = tile_data_pointer (dst_tile, 0, 0);

          for (x = 64; x < 64*2; x++)
            {
              for (y = 64; y < 64*2; y++)
                {
                  *pointer = big_cache[y * BIG_CACHE_W + x*4];
                  *(pointer+1) = big_cache[y * BIG_CACHE_W + x*4 + 1];
                  *(pointer+2) = big_cache[y * BIG_CACHE_W + x*4 + 2];
                  *(pointer+3) = 255;

                  pointer += 4;
                }
            }

          tile_release (dst_tile, TRUE);

 */


/**
 * siox_done:
 * @state: The state of this tool.
 *
 * Frees the memory assciated with the state.
 */
void
siox_done (SioxState *state)
{
  g_return_if_fail (state != NULL);

  g_free (state->fgsig);
  g_free (state->bgsig);
  g_hash_table_destroy (state->cache);

  g_slice_free (SioxState, state);

#ifdef SIOX_DEBUG
  g_printerr ("siox.c: siox_done()\n");
#endif
}
