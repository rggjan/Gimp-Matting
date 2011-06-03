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

#define IMAGE_DEBUG_PPM

#ifdef IMAGE_DEBUG_PPM
#include "stdio.h"
#include "stdlib.h"
#endif

#define DEBUG_EXTENSION

#define SEARCH_RADIUS 10
#define MATTING_SQUARED_COLOR_DISTANCE 25

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

#define CACHE_CHANNELS 4

#define BIG_CACHE_W 64*3*CACHE_CHANNELS
#define BIG_CACHE_H 64*3
#define BIG_CACHE_SIZE BIG_CACHE_W*BIG_CACHE_H

#define GET_PIXEL(big_cache, x, y, color) big_cache[(64+y)*BIG_CACHE_W\
                                           +(64+x)*CACHE_CHANNELS+color]

#define BIGGER_CACHE_W 64*7*CACHE_CHANNELS
#define BIGGER_CACHE_H 64*7
#define BIGGER_CACHE_SIZE BIGGER_CACHE_W*BIGGER_CACHE_H

#define GET_PIXEL_BIGGER(bigger_cache, x, y, color) bigger_cache[(192+y)*BIGGER_CACHE_W\
                                           +(192+x)*CACHE_CHANNELS+color]

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

#ifdef IMAGE_DEBUG_PPM
static void
debug_image (const char* filename, int dimx, int dimy, guchar *data, int padding, int num_colors)
{
  int i, j;
  FILE *fp = fopen (filename, "wb"); /* b - binary mode */

  if (num_colors == 3)
    {
      fprintf (fp, "P6\n%d %d\n255\n", dimx, dimy);
    }
  else if (num_colors == 1)
    {
      fprintf (fp, "P5\n%d %d\n255\n", dimx, dimy);
    }
  else
    {
      printf ("Problem!\n");
      exit (1);
    }

  {
    guchar* current = data;

    for (j = 0; j < dimy; ++j)
      {
        for (i = 0; i < dimx; ++i)
          {
            if (num_colors == 3)
              {
                static unsigned char color[3];
                color[0] = current[0];
                color[1] = current[1];
                color[2] = current[2];
                (void) fwrite (color, 1, 3, fp);
              }
            else
              {
                static unsigned char color;
                color = current[0];
                (void) fwrite (&color, 1, 1, fp);
              }

            current += padding;
          }
      }
  }
  fclose (fp);
}
#endif

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

static void
load_big_cache (TileManager *source, guchar *big_cache, gint tx, gint ty,
                gint radius)
{
  gint    xdiff;
  gint    ydiff;
  gint    x, y;
  gint    bx, by;
  gint    width_tile, height_tile;
  gint   cache_width = (radius*2+1)*64*CACHE_CHANNELS;
  
  Tile   *src_tile;
  guchar *pointer;

  g_return_if_fail (tile_manager_bpp(source) == 4);

  for (ydiff = -radius; ydiff <= radius; ydiff++)
    {
      for (xdiff = -radius; xdiff <= radius; xdiff++)
        {
          src_tile = tile_manager_get_at (source, tx + xdiff, ty + ydiff, TRUE, FALSE);
          width_tile = 0;
          height_tile = 0;

          if (src_tile)
            {
              pointer = tile_data_pointer (src_tile, 0, 0);
              width_tile = tile_ewidth (src_tile);
              height_tile = tile_eheight (src_tile);

              for (y = 0; y < height_tile; y++)
                {
                  by = (ydiff + radius)*64 + y;
                  for (x = 0; x < width_tile; x++)
                    {
                      bx = (xdiff + radius)*64 + x;

                      big_cache[by * cache_width + bx * 4] = *pointer;
                      big_cache[by * cache_width + bx * 4 + 1] = *(pointer + 1);
                      big_cache[by * cache_width + bx * 4 + 2] = *(pointer + 2);
                      big_cache[by * cache_width + bx * 4 + 3] = *(pointer + 3);

                      pointer += 4;
                    }
                }
              
              tile_release (src_tile, FALSE);
            }

          for (y = 0; y < 64; y++)
            {
              by = (ydiff + radius)*64 + y;
              for (x = width_tile; x < 64; x++)
                {
                  bx = (xdiff + radius)*64 + x;

                  big_cache[by * cache_width + bx * 4] = 255;
                  big_cache[by * cache_width + bx * 4 + 1] = 0;
                  big_cache[by * cache_width + bx * 4 + 2] = 0;
                  big_cache[by * cache_width + bx * 4 + 3] = 128;
                }
            }
          
           for (y = height_tile; y < 64; y++)
            {
              by = (ydiff + radius)*64 + y;
              for (x = 0; x < width_tile; x++)
                {
                  bx = (xdiff + radius)*64 + x;

                  big_cache[by * cache_width + bx * 4] = 255;
                  big_cache[by * cache_width + bx * 4 + 1] = 0;
                  big_cache[by * cache_width + bx * 4 + 2] = 0;
                  big_cache[by * cache_width + bx * 4 + 3] = 128;
                }
            }
          
        }
    }
}

// retrieves the x/y coordinates form a key

static inline void
get_pos_from_key (gint64 *key,
                  gint *x,
                  gint *y)
{
  gint64 temp = *key;
  *x = (temp>>32 & 0xffffffff);
  *y = temp & 0xffffffff;
}

// evaluate the energy function for found color fg/bg for pixel situated at x/y

static float
objective_function (guchar *fg,
                    guchar *bg,
                    gint x,
                    gint y,
                    guchar *bigger_cache)
{
  //gint en = 3;
  //gint ea = 2;
  gint ef = 1;
  //gint eb = 4;
  gint xi, yi;
  float np, ap, dpb, dpf, *pointer;
  float newAlpha, pfp, r, g, b;

  dpb = pow (bg[3], ef);
  dpf = pow (fg[3], ef);

  for (xi = -1; xi < 2; xi++)
    {
      for (yi = -1; yi < 2; yi++)
        {
          if (xi != 0 || yi != 0)
            {
              r = GET_PIXEL_BIGGER (bigger_cache, x + xi, y + yi, 0);
              g = GET_PIXEL_BIGGER (bigger_cache, x + xi, y + yi, 1);
              b = GET_PIXEL_BIGGER (bigger_cache, x + xi, y + yi, 2);
              newAlpha = (fg[0] - bg[0]) * r + (fg[1] - bg[1]) * g + (fg[2] - bg[2]) * b;
              newAlpha = newAlpha / sqrt ((fg[0] - bg[0])*(fg[0] - bg[0])+(fg[1] - bg[1])*(fg[1] - bg[1])+(fg[2] - bg[2])*(fg[2] - bg[2]));
              // TODO: check if it should be 1 - newAlpha
              newAlpha = (newAlpha > 1 ? 1 : (newAlpha < 0 ? 0 : newAlpha));
              np += (r - newAlpha * fg[0] + (1 - newAlpha) * bg[0])*(r - newAlpha * fg[0] + (1 - newAlpha) * bg[0])
                      + (g - newAlpha * fg[1] + (1 - newAlpha) * bg[1])*(g - newAlpha * fg[0] + (1 - newAlpha) * bg[1])
                      + (b - newAlpha * fg[2] + (1 - newAlpha) * bg[2])*(b - newAlpha * fg[0] + (1 - newAlpha) * bg[2]);
            }
        }
    }

  // TODO: this calculation will happen twice (also when storing the final value)
  r = GET_PIXEL_BIGGER (bigger_cache, x, y, 0);
  g = GET_PIXEL_BIGGER (bigger_cache, x, y, 1);
  b = GET_PIXEL_BIGGER (bigger_cache, x, y, 2);
  newAlpha = (fg[0] - bg[0]) * r + (fg[1] - bg[1]) * g + (fg[2] - bg[2]) * b;
  newAlpha = newAlpha / sqrt ((fg[0] - bg[0])*(fg[0] - bg[0])+(fg[1] - bg[1])*(fg[1] - bg[1])+(fg[2] - bg[2])*(fg[2] - bg[2]));
  // TODO: check if it should be 1 - newAlpha
  newAlpha = (newAlpha > 1 ? 1 : (newAlpha < 0 ? 0 : newAlpha));

  // TODO uncomment things below... commented them because of compiler errors!
  
  pointer = &fg[4];
  pfp = *pointer;
  pointer = &bg[4];
  pfp = *pointer / (pfp + *pointer);
  // TODO: should newAlpha be 1 or 255 here?
  ap = pfp + (1 - 2 * pfp) * newAlpha;

  return np + ap + dpb + dpf;
}

// searches for known regions for a given unknown pixel in a hash table

static void
search_neighborhood (gpointer key,
                     guchar *value,
                     gpointer *args)
{
  gint pos_x, pos_y;
  gint tx, ty;
  gint x, y, distance, direction, toggle;
  guchar values[8 * 4 * 2]; //(rgb + distance + float magn. gradient) * 4 directions * fground/bground
  guchar prevval[3 * 4];
  guchar a;
  double angle;

  // Structure: [fg1, bg1, fg2, bg2, ...]
  // 1 = right, 2 = up, ...
  gboolean found[8];

  gint permutation[4];
  Tile *tile;
  guchar *pointer;
  float *pointertemp;

  for (direction = 0; direction < 8; direction++)
    {
      found[direction] = FALSE;
    }

  // initialize to original value
  for (distance = 0; distance < 4; distance++)
    {
      prevval[distance] = value[0];
      prevval[distance + 1] = value[1];
      prevval[distance + 2] = value[2];
    }
  permutation[0] = 1;
  permutation[1] = 1;
  permutation[2] = -1;
  permutation[3] = -1;


  //g_printf("key: %d\n", key);

  // TODO do this with shift, (little endian problems...)
  get_pos_from_key (key, &pos_x, &pos_y); //TODO: x and y swapped here!
  
  //g_printf("key: x = %i, y = %i\n", pos_x, pos_y);

  // TODO add pi somehow
  angle = (pos_x % 3) + (pos_y % 3) * 3;

  tx = (pos_x - (pos_x % 64)) / 64;
  ty = (pos_y - (pos_y % 64)) / 64;

  pos_x = pos_x - 64 * tx;
  pos_y = pos_y - 64 * ty;

  // TODO: add list of sorted unknown regions to traverse so that same tiles are not loaded several times
  /*if (args[3] != tx && args[4] != ty)
    {
      load_big_cache ((TileManager*) args[2], (guchar*) args[0], tx, ty, 3);
      g_printf ("Cache loaded! for tiles %i %i\n", tx, ty);
      args[3] = tx;
      args[4] = ty;
    }*/ // TODO uncomment this!

  for (distance = 6; distance < 3 * 64; distance += 6)
    {
      x = floor (cos (angle) * distance);
      y = floor (sin (angle) * distance);
      for (direction = 0; direction < 4; direction++)
        {
          if (found[direction * 2] && found[direction * 2 + 1])
            {
              break;
            }
          else
            {
              guchar r, g, b;
              
              gint xtmp = permutation[direction] * x + pos_x;
              gint ytmp = permutation[(direction + 1) % 4] * y + pos_y;

              r = GET_PIXEL_BIGGER (((guchar*) args[0]), xtmp, ytmp, 0);
              g = GET_PIXEL_BIGGER (((guchar*) args[0]), xtmp, ytmp, 1);
              b = GET_PIXEL_BIGGER (((guchar*) args[0]), xtmp, ytmp, 2);
              a = GET_PIXEL_BIGGER (((guchar*) args[0]), xtmp, ytmp, 3);

              // check if it's foreground or background, depending on a
              // toggle = 0 equals foreground
              for (toggle = 0; toggle < 2; toggle++)
                {
                  // add gradient distance (as float)
                  if (!found[direction + toggle])
                    {
                      // TODO: Check if values is initialized to zero
                      pointertemp = &values[direction * 16 + (toggle * 8) + 4];
                      *pointertemp += sqrt ((prevval[direction] - r)*(prevval[direction] - r)
                                            + (prevval[direction + 1] - g)*(prevval[direction + 1] - g)
                                            + (prevval[direction + 2] - b)*(prevval[direction + 2] - b));
                    }

                  if (a == (toggle == 0 ? 255 : 0) && !found[direction*2 + toggle])
                    {

                      if (xtmp > 100 || ytmp > 100)
                        {
                          //g_printf ("found for pixel pos: %i %i, pixel location %i %i\n", pos_x, pos_y, xtmp, ytmp);

                          //g_printf ("rgba = %i, %i, %i, %i\n", r, g, b, a);
                        }
                      
                      values[direction * 16 + (toggle * 8)] = r;
                      values[direction * 16 + 1 + (toggle * 8)] = g;
                      values[direction * 16 + 2 + (toggle * 8)] = b;
                      values[direction * 16 + 3 + (toggle * 8)] = distance;
                      found[direction + toggle] = TRUE;
                    }
                }
            }
        }
    }

  {
    float min = -1;
    gint minindexf = -1;
    gint minindexb = -1;
    float temp;
    float newAlpha;

    // calculate energy function for every fg/bg pair
    for (distance = 0; distance < 4; distance++)
      {
        if (found[distance * 2])
          {
            for (direction = 0; direction < 4; direction++)
              {
                if (found[direction * 2 + 1])
                  {
                    temp = objective_function (&values[distance * 16],
                                               &values[direction * 16 + 8],
                                               pos_x,
                                               pos_y,
                                               args[0]);
                    if (temp < min || min < 0)
                      {
                        min = temp;
                        minindexf = distance;
                        minindexb = direction;
                      }
                  }

              }
          }
      }
    tile = tile_manager_get_at (args[1], tx, ty, TRUE, TRUE);
    pointer = tile_data_pointer (tile, pos_x, pos_y);

    newAlpha = (values[minindexf] - values[minindexb]) * value[0] + (values[minindexf + 1] - values[minindexb + 1]) * value[1] + (values[minindexf + 2] - values[minindexb + 2]) * value[2];
    newAlpha = newAlpha / sqrt ((values[minindexf] - values[minindexb])*(values[minindexf] - values[minindexb])+(values[minindexf + 1] - values[minindexb + 1])*(values[minindexf + 1] - values[minindexb + 1])+(values[minindexf + 2] - values[minindexb + 2])*(values[minindexf + 2] - values[minindexb + 2]));
    // TODO: check if it should be 1 - newAlpha
    newAlpha = (newAlpha > 1 ? 1 : (newAlpha < 0 ? 0 : newAlpha));

    // test: do combination of best match
    // should return a very similar image as before
    pointer[0] = floor (values[minindexf] * newAlpha) + floor (values[minindexb]*(1 - newAlpha));
    pointer[1] = floor (values[minindexf + 1] * newAlpha) + floor (values[minindexb + 1]*(1 - newAlpha));
    pointer[2] = floor (values[minindexf + 2] * newAlpha) + floor (values[minindexb + 2]*(1 - newAlpha));
    pointer[3] = 255;
    tile_release (tile, TRUE);
    //printf("values: %i %i %i | %i %i %i | %i %i %i | %i %i %i\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
  }
}


static void
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

static void
update_mask (TileManager* result_layer,
             TileManager* mask_layer)
{
  PixelRegion result, mask;
  PixelRegionIterator *pr;
  gint row, col;
  int width, height;

  width = tile_manager_width (result_layer);
  height = tile_manager_height(result_layer);

  pixel_region_init (&result, result_layer, 0, 0, width, height, FALSE);
  pixel_region_init (&mask, mask_layer, 0, 0, width, height, TRUE);

  g_return_if_fail (result.bytes == 4 && mask.bytes == 1); // TODO check if indexed etc...

  for (pr = pixel_regions_register (2, &result, &mask);
          pr != NULL;
          pr = pixel_regions_process (pr))
    {
      const guchar *result_data = result.data;
      guchar *mask_data = mask.data;

      for (row = 0; row < result.h; row++)
        {
          const guchar *r = result_data;
          guchar *m = mask_data;

          for (col = 0; col < result.w; col++, r += result.bytes, ++m) // TODO check if 4 is ok
            {
              guchar mask;
              mask = m[0];
              if (mask != MATTING_USER_FOREGROUND &&
                  mask != MATTING_USER_BACKGROUND)
                {
                  guchar result = r[3];
                  if (result == 255)
                    m[0] = MATTING_ALGO_FOREGROUND;
                  else if (result == 0)
                    m[0] = MATTING_ALGO_BACKGROUND;
                  else
                    m[0] = MATTING_ALGO_UNDEFINED;
                }
            }

          result_data += result.rowstride;
          mask_data += mask.rowstride;
        }
    }
}

static inline gboolean
check_closeness (guchar color[3], guchar *big_cache, gint x, gint y, guchar* result)
{
  guchar value;
  gint   color_distance_sum;
  gint   color_distance;
  gint   i;

  value = GET_PIXEL (big_cache, x, y, 3);

  if (value != 128)
    {
      color_distance_sum = 0;
      for (i = 0; i < 3; i++)
        {
          color_distance = color[i] -
                  GET_PIXEL (big_cache, x, y, i);
          color_distance_sum +=
                  color_distance*color_distance;
        }

      if (color_distance_sum < MATTING_SQUARED_COLOR_DISTANCE)
        {
          *result = value;
          return TRUE;
        }
    }

  return FALSE;
}

static inline void
search_for_neighbours (guchar* big_cache, gint x, gint y, guchar* result)
{
  guchar color[3];
  gint   i;
  gint   radius;
  gint   n;

  guchar alpha;

  alpha = GET_PIXEL (big_cache, x, y, 3);
  if (alpha != 128)
    {
      *result = alpha;
      return;
    }

  for (i = 0; i < 3; i++)
    {
      color[i] = GET_PIXEL (big_cache, x, y, i);
    }

  for (radius = 0; radius <= SEARCH_RADIUS; radius++)
    {
      for (n = -radius; n < radius; n++)
        {
          if (check_closeness (color, big_cache, x+radius, y+n, result))
            return;

          if (check_closeness (color, big_cache, x-radius, y+n, result))
            return;

          if (check_closeness (color, big_cache, x+n, y+radius, result))
            return;

          if (check_closeness (color, big_cache, x+n, y-radius, result))
            return;
        }
    }

  *result = 128;
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
  //gint         width, height;
  guchar      *big_cache;
  guchar      *bigger_cache;
  gint         tiles_x, tiles_y;

  Tile        *tile;

  gint         tx, ty, x, y;
  guchar      *pointer;
  gpointer     foreach_args[5];
  gint         loaded_tile_x, loaded_tile_y;
  gint         i;

  gboolean     unknown;
  
  static GHashTable *unknown_hash = NULL;
  guchar      *unknown_pixel;
  //guchar      *resulttest;
  
  //resulttest = g_malloc (8);

  unknown_hash = g_hash_table_new(g_int64_hash, g_int64_equal);
  
  g_return_if_fail (state != NULL);
  g_return_if_fail (mask != NULL && tile_manager_bpp (mask) == 1);
  g_return_if_fail (x1 >= 0);
  g_return_if_fail (x2 > x1 && x2 <= tile_manager_width (mask));
  g_return_if_fail (y1 >= 0);
  g_return_if_fail (y2 > y1 && y2 <= tile_manager_height (mask));
  g_return_if_fail (smoothness >= 0);
  g_return_if_fail (progress_data == NULL || progress_callback != NULL);

  //width  = state->width;
  //height = state->height;

  g_return_if_fail (TILE_WIDTH == 64 && TILE_HEIGHT == 64);

  g_return_if_fail (tile_manager_bpp (state->pixels) == 3);
  g_return_if_fail (tile_manager_bpp (result_layer) == 4);

  big_cache = g_malloc (BIG_CACHE_SIZE);

  tiles_x = tile_manager_tiles_per_col (state->pixels);
  tiles_y = tile_manager_tiles_per_row (state->pixels);
  //tiles_x = state->pixels->ntile_cols;
  //tiles_y = state->pixels->ntile_rows;

  initialize_new_layer (state->pixels, working_layer, mask);

  i = 0;
  {
    TileManager* tmp;
    i++;

    for (tx = 0; tx < tiles_x - 1; tx++)
      {
        for (ty = 0; ty < tiles_y - 1; ty++)
          {
            static char buffer[100];

            load_big_cache (working_layer, big_cache, tx, ty, 1);

            snprintf (buffer, 100, "big_cache_tx_%i_ty_%i.ppm", tx, ty);
            debug_image (buffer, 64 * 3, 64 * 3, big_cache, 4, 3);

            snprintf (buffer, 100, "big_cache_tx_%i_ty_%i_alpha.ppm", tx, ty);
            debug_image (buffer, 64 * 3, 64 * 3, big_cache + 3, 4, 1);

            tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);
            pointer = tile_data_pointer (tile, 0, 0);

            for (x = 0; x < 64; x++)
              {
                for (y = 0; y < 64; y++, pointer += 4)
                  {
                    pointer[0] = GET_PIXEL (big_cache, x, y, 0);
                    pointer[1] = GET_PIXEL (big_cache, x, y, 1);
                    pointer[2] = GET_PIXEL (big_cache, x, y, 2);

                    search_for_neighbours (big_cache, x, y, pointer + 3);

                    unknown = (GET_PIXEL (big_cache, x, y, 3) == 128);
                    
                    if (unknown)
                      {
                        gint64 *addr;

                        unknown_pixel = g_malloc (8 * sizeof (guchar));
                        unknown_pixel[0] = pointer[0];
                        unknown_pixel[1] = pointer[1];
                        unknown_pixel[2] = pointer[2];
                        unknown_pixel[3] = pointer[3];
                        unknown_pixel[4] = 0;
                        unknown_pixel[5] = 0;
                        unknown_pixel[6] = 0;
                        unknown_pixel[7] = 0;
                        // TODO: free this memory
                        addr = g_malloc (sizeof (gint64));
                        *addr = (((gint64) x) << 32) + y;

                        //g_printf("Inserting xy %i, %i\n", x, y);

                        g_hash_table_insert (unknown_hash, addr, unknown_pixel);
                      }
                  }
              }

            tile_release (tile, TRUE);
          }
      }
    tmp = working_layer;
    working_layer = result_layer;
    result_layer = tmp;
  }

#ifdef DEBUG_EXTENSION
  update_mask (result_layer, mask);
  return;
#endif

  bigger_cache = g_malloc (BIGGER_CACHE_SIZE);
  foreach_args[0] = bigger_cache;
  foreach_args[1] = working_layer;
  foreach_args[2] = result_layer;
  loaded_tile_x = -1;
  loaded_tile_y = -1;
  foreach_args[3] = &loaded_tile_x;
  foreach_args[4] = &loaded_tile_y;

  // TODO replace foreach
  g_hash_table_foreach (unknown_hash, (GHFunc) search_neighborhood, foreach_args);

  // TODO do this only once
  g_free (big_cache);
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
