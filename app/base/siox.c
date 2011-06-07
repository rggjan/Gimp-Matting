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
#include <glib-2.0/glib/gprintf.h>

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

//#define IMAGE_DEBUG_PPM

//#define DEBUG_PHASE1
//#define DEBUG_PHASE2

#define DEBUG_PREDEFINED_MASK
#define DEBUG_PREDEFINED_MASK_WRITE FALSE

#ifdef IMAGE_DEBUG_PPM
#include "stdio.h"
#include "stdlib.h"
#endif

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
#define SIOX_GRAY_DIMS  2

/* For findmaxblob:
 * Find all blobs with area not smaller than sizefactor of biggest blob
 * CHECKME: Should the user decide this with a slider?
 */
#define MULTIBLOB_DEFAULT_SIZEFACTOR 4
#define MULTIBLOB_ONE_BLOB_ONLY      0


typedef union
{
  guint64 value;
  struct
  {
    guint32 x;
    guint32 y;
  } coords;
} HashAddress;

struct HashEntry_
{
  guchar color[3];

  guchar foreground[3];
  guchar background[3];

  guchar foreground_refined[3];
  guchar background_refined[3];

  guchar alpha;
  guchar alpha_refined;

  gfloat sigma_f_squared;
  gfloat sigma_b_squared;

  gboolean pair_found;
  HashAddress this;
  HashAddress next;
};

typedef struct HashEntry_ HashEntry;

// TODO should be 300*6/64, actually...
#define BIG_CACHE_CHANNELS 4
#define BIG_CACHE_RADIUS 3
#define BIG_CACHE_SIZE ((BIG_CACHE_RADIUS*2+1)*64)
#define GET_PIXEL(big_cache, x, y, color) (big_cache[BIG_CACHE_RADIUS*64+y][BIG_CACHE_RADIUS*64+x][color])
typedef guchar BigCache[BIG_CACHE_SIZE][BIG_CACHE_SIZE][BIG_CACHE_CHANNELS];

#define HASH_CACHE_EXTRA 14
#define HASH_CACHE_SIZE (64+HASH_CACHE_EXTRA*2)
#define GET_ENTRY(hash_cache, x, y) (hash_cache[HASH_CACHE_EXTRA+y][HASH_CACHE_EXTRA+x])
typedef HashEntry* HashCache[HASH_CACHE_SIZE][HASH_CACHE_SIZE];

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

static void debug_cache (const char* filename, BigCache cache, int radius)
{
  static char filename_alpha[100];
  static char filename_rgb[100];
  const int size = (radius * 2 + 1) * 64;

  FILE *f_color, *f_alpha;

  snprintf (filename_alpha, 100, "%s_alpha.ppm", filename);
  snprintf (filename_rgb, 100, "%s.ppm", filename);

  f_color = fopen (filename_rgb, "wb"); /* b - binary mode */
  f_alpha = fopen (filename_alpha, "wb"); /* b - binary mode */

  fprintf (f_color, "P6\n%d %d\n255\n", size, size);
  fprintf (f_alpha, "P5\n%d %d\n255\n", size, size);
  {
    int x, y;
    for (y = -radius * 64; y < (radius + 1) * 64; y++)
      {
        for (x = -radius * 64; x < (radius + 1) * 64; x++)
          {
            fwrite (&GET_PIXEL(cache, x, y, 0), 1, 3, f_color);
            fwrite (&GET_PIXEL(cache, x, y, 3), 1, 1, f_alpha);
          }
      }
  }
  fclose (f_color);
  fclose (f_alpha);
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
load_big_cache (TileManager *source, BigCache big_cache, gint tx, gint ty,
                gint radius)
{
  gint    xdiff;
  gint    ydiff;
  gint    x, y;
  gint    bx, by;
  gint    width_tile, height_tile;

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
                  by = ydiff * 64 + y;
                  for (x = 0; x < width_tile; x++)
                    {
                      bx = xdiff * 64 + x;

                      GET_PIXEL(big_cache, bx, by, 0) = pointer[0];
                      GET_PIXEL(big_cache, bx, by, 1) = pointer[1];
                      GET_PIXEL(big_cache, bx, by, 2) = pointer[2];
                      GET_PIXEL(big_cache, bx, by, 3) = pointer[3];

                      pointer += 4;
                    }
                }

              tile_release (src_tile, FALSE);
            }

          for (y = 0; y < 64; y++)
            {
              by = ydiff * 64 + y;
              for (x = width_tile; x < 64; x++)
                {
                  bx = xdiff * 64 + x;

                  GET_PIXEL(big_cache, bx, by, 0) = 0;
                  GET_PIXEL(big_cache, bx, by, 1) = 0;
                  GET_PIXEL(big_cache, bx, by, 2) = 0;
                  GET_PIXEL(big_cache, bx, by, 3) = 128;
                }
            }

          for (y = height_tile; y < 64; y++)
            {
              by = ydiff * 64 + y;
              for (x = 0; x < width_tile; x++)
                {
                  bx = xdiff * 64 + x;

                  GET_PIXEL(big_cache, bx, by, 0) = 0;
                  GET_PIXEL(big_cache, bx, by, 1) = 0;
                  GET_PIXEL(big_cache, bx, by, 2) = 0;
                  GET_PIXEL(big_cache, bx, by, 3) = 128;
                }
            }

        }
    }
}

// TODO only use rgb in certain places
typedef struct
{
  guchar color[3];
  gboolean found;
  gint distance;
  gfloat gradient;
  gint x;
  gint y;
} SearchStructure;

// Project the Point P onto the line from A-B.
// alpha_pointer receives the calculated alpha value between 0 and 1
// where (0=A, 1=B)
// returns the squared distance of the best projected point to P
static gfloat projection (guchar A[3], guchar B[3], guchar P[3], float* alpha_pointer)
{
  gfloat ABx = B[0] - A[0];
  gfloat ABy = B[1] - A[1];
  gfloat ABz = B[2] - A[2];

  gfloat APx = P[0] - A[0];
  gfloat APy = P[1] - A[1];
  gfloat APz = P[2] - A[2];

  gfloat dot_AB = ABx * ABx + ABy * ABy + ABz * ABz;
  gfloat alpha = (ABx * APx + ABy * APy + ABz * APz) / dot_AB;

  gfloat PPx, PPy, PPz;

  alpha = alpha > 1 ? 1 : (alpha < 0 ? 0 : alpha);

  PPx = P[0] - (A[0] + alpha * ABx);
  PPy = P[1] - (A[1] + alpha * ABy);
  PPz = P[2] - (A[2] + alpha * ABz);

  if (alpha_pointer)
    *alpha_pointer = alpha;

  return PPx * PPx + PPy * PPy + PPz * PPz;
}

// evaluate the energy function for found color fg/bg for pixel situated at x/y
static float
objective_function (SearchStructure *fg,
                    SearchStructure *bg,
                    gint x,
                    gint y,
                    BigCache big_cache,
                    float* best_alpha)
{
  gint xi, yi;
  float ap, pfp;

  float Np = 0;

  float finalAlpha;

  for (yi = -1; yi < 2; yi++)
    {
      for (xi = -1; xi < 2; xi++)
        {
          // TODO check borders
          guchar P[3];

          P[0] = GET_PIXEL (big_cache, x + xi, y + yi, 0);
          P[1] = GET_PIXEL (big_cache, x + xi, y + yi, 1);
          P[2] = GET_PIXEL (big_cache, x + xi, y + yi, 2);
          if (xi != 0 || yi != 0)
            {
              Np += projection (fg->color, bg->color, P, NULL);
            }
          else
            {
              Np += projection (fg->color, bg->color, P, &finalAlpha);
            }
        }
    }

  // TODO calculate pfp, this is a GLOBAL value!!!
  // pfp = bg->gradient / (fg->gradient + bg->gradient);
  pfp = 0.5;
  ap = pfp + (1 - 2 * pfp) * (1 - finalAlpha);

  *best_alpha = finalAlpha;
  return pow(Np, 3) * pow(bg->distance, 4) * pow(fg->distance, 1) * pow(ap, 2);
}

// searches for known regions for a given unknown pixel in a hash table

static inline gfloat dist_squared(guchar x1, guchar y1, guchar z1, guchar x2, guchar y2, guchar z2)
{
  gint diff_x = (gint)x1 - (gint)x2;
  gint diff_y = (gint)y1 - (gint)y2;
  gint diff_z = (gint)z1 - (gint)z2;

  return (gfloat) (diff_x * diff_x + diff_y * diff_y + diff_z * diff_z);
}

static gfloat calculate_variance (guchar P[3], gint x, gint y, BigCache big_cache)
{
  gint yi;
  gint xi;
  guchar my_alpha = GET_PIXEL (big_cache, x, y, 3);

  gfloat sum = 0;
  gint num_found = 0;

  for (yi = -2; yi <= 2; yi++)
    {
      for (xi = -2; xi <= 2; xi++)
        {
          // TODO: got a SEGFAULT once!!!
          // TODO check borders
          if (GET_PIXEL (big_cache, x, y, 3) == my_alpha)
            {
              sum += dist_squared(GET_PIXEL (big_cache, x + xi, y + yi, 0),
                                  GET_PIXEL (big_cache, x + xi, y + yi, 1),
                                  GET_PIXEL (big_cache, x + xi, y + yi, 2),
                                  P[0],
                                  P[1],
                                  P[2]);

              num_found++;
            }
        }
    }

  return sum / num_found;
}

typedef struct
{
  gint x;
  gint y;
  gfloat diff;
} TopColor;

static void load_hash_cache (GHashTable* unknown_hash, HashCache hash_cache,
                             gint tx, gint ty)
{

  gint x, y;
  gint xoff = tx * 64;
  gint yoff = ty * 64;

  for (y = -HASH_CACHE_EXTRA; y < HASH_CACHE_EXTRA + 64; y++)
    {
      for (x = -HASH_CACHE_EXTRA; x < HASH_CACHE_EXTRA + 64; x++)
        {
          HashAddress address;
          address.coords.x = x + xoff;
          address.coords.y = y + yoff;

          GET_ENTRY(hash_cache, x, y) = g_hash_table_lookup (unknown_hash, &address);
        }
    }
}

static void inline
compare_neighborhood (HashEntry* entry, gint *current_tx, gint* current_ty,
                      HashCache hash_cache, GHashTable* unknown_hash,
                      BigCache big_cache)
{
  gint pos_x, pos_y;
  gint tx, ty;

  gint xdiff, ydiff;
  gint num;
  gint matches = 0;

  HashEntry *current;
  TopColor top3[3];

  // Load coordinates from entry
  pos_x = entry->this.coords.x;
  pos_y = entry->this.coords.y;

  tx = pos_x / 64;
  ty = pos_y / 64;

  pos_x = pos_x - 64 * tx;
  pos_y = pos_y - 64 * ty;

  if (*current_tx != tx || *current_ty != ty)
    {
      load_hash_cache (unknown_hash, hash_cache, tx, ty);

      *current_tx = tx;
      *current_ty = ty;
    }

  for (num = 0; num < 3; num++)
    {
      top3[num].diff = INFINITY;
    }

  for (ydiff = -HASH_CACHE_EXTRA; ydiff <= HASH_CACHE_EXTRA; ydiff++)
    {
      for (xdiff = -HASH_CACHE_EXTRA; xdiff <= HASH_CACHE_EXTRA; xdiff++)
        {
          current = GET_ENTRY(hash_cache, pos_x + xdiff, pos_y + ydiff);

          if (current && current->pair_found)
            {
              gfloat current_alpha;
              // TODO: check if this projection is right! do we have to use other color values maybe?
              gfloat temp = projection (current->foreground,
                                        current->background,
                                        entry->color,
                                        &current_alpha);

              // check if color is better than least best of colors, add the color and sort the list
              if (temp < top3[2].diff)
                {
                  gint i = 2;
                  for (num = 1; num >= 0; num--)
                    {
                      if (temp < top3[num].diff)
                        {
                          i = num;
                        }
                    }
                  for (num = 2; num > i; num--)
                    {
                      top3[num] = top3[num-1];
                    }
                  top3[i].diff = temp;
                  top3[i].x = pos_x + xdiff;
                  top3[i].y = pos_y + ydiff;
                  matches++;
                }
            }
        }
    }

  if (matches >= 1)
    {
      gfloat colordiff, lower, upper;

      guchar new_fg[3] = {0, 0, 0};
      guchar new_bg[3] = {0, 0, 0};
      gfloat new_sigma_f_squared = 0;
      gfloat new_sigma_b_squared = 0;

      gint index;

      // TODO check if we have to treat all colors in unit cube

      matches = matches > 3 ? 3 : matches;

      for (num = 0; num < matches; num++)
        {
          current = GET_ENTRY(hash_cache, top3[num].x, top3[num].y);
          for (index = 0; index < 3; index++)
            {
              // TODO: check if we should use ints here!
              new_fg[index] += floor(current->foreground[index] / matches);
              new_bg[index] += floor(current->background[index] / matches);
              new_sigma_f_squared += current->sigma_f_squared / matches;
              new_sigma_b_squared += current->sigma_b_squared / matches;
            }
        }

      colordiff = dist_squared(entry->color[0],
                               entry->color[1],
                               entry->color[2],
                               new_fg[0],
                               new_fg[1],
                               new_fg[2]);

      if (colordiff > new_sigma_f_squared)
        {
          for (index = 0; index < 3; index++)
            {
              entry->foreground_refined[index] = new_fg[index];
            }
        }
      else
        {
          for (index = 0; index < 3; index++)
            {
              entry->foreground_refined[index] = entry->color[index];
            }
        }

      colordiff = dist_squared(entry->color[0],
                               entry->color[1],
                               entry->color[2],
                               new_bg[0],
                               new_bg[1],
                               new_bg[2]);

      if (colordiff > new_sigma_b_squared)
        {
          for (index = 0; index < 3; index++)
            {
              entry->background_refined[index] = new_bg[index];
            }
        }
      else
        {
          for (index = 0; index < 3; index++)
            {
              entry->background_refined[index] = entry->color[index];
            }
        }

      gfloat current_alpha;
      // TODO: check if this projection is right! do we have to use other color values maybe?
      projection (entry->foreground_refined,
                  entry->background_refined,
                  entry->color,
                  &current_alpha);

      entry->alpha_refined = (1 - current_alpha) * 255;

    }
  else
    {
      entry->foreground_refined[0] = 255;
      entry->foreground_refined[1] = 0;
      entry->foreground_refined[2] = 0;
      entry->alpha_refined = 255;
    }

  /*
    if (min == INFINITY)
      {
        entry->foreground_refined[0] = 255;
        entry->foreground_refined[1] = 0;
        entry->foreground_refined[2] = 0;
        entry->alpha_refined = 255;
      }*/
}

static void inline
search_neighborhood (HashEntry* entry, gint *current_tx, gint *current_ty,
                     BigCache big_cache, TileManager* layer)
{
  gint pos_x, pos_y;
  gint tx, ty;
  gint xdiff, ydiff;

  gint direction, toggle, distance;

  SearchStructure found[2][4];

  guchar prevval[3 * 4];
  double angle;
  gint permutation[4];

  /*
  Tile *tile;
  guchar *pointer;
  float *pointertemp;
  */
  for (toggle = 0; toggle < 2; toggle++)
    {
      for (direction = 0; direction < 4; direction++)
        {
          found[toggle][direction].found = FALSE;
        }
    }


  // initialize to original value
  for (distance = 0; distance < 4; distance++)
    {
      prevval[distance] = entry->color[0];
      prevval[distance + 1] = entry->color[1];
      prevval[distance + 2] = entry->color[2];
    }

  // Load coordinates from entry
  pos_x = entry->this.coords.x;
  pos_y = entry->this.coords.y;

  tx = pos_x / 64;
  ty = pos_y / 64;

  pos_x = pos_x - 64 * tx;
  pos_y = pos_y - 64 * ty;

  if (*current_tx != tx || *current_ty != ty)
    {
      load_big_cache (layer, big_cache, tx, ty, 3);
      //g_printf ("Cache loaded! for tiles %i %i\n", tx, ty);
      *current_tx = tx;
      *current_ty = ty;

#ifdef IMAGE_DEBUG_PPM
      {
        static gchar buffer[100];
        snprintf (buffer, 100, "bigger_cache_tx_%i_ty_%i", tx, ty);
        debug_cache (buffer, big_cache, 3);
      }
#endif
    }


  // in a 9x9 window, we want to have values in a 90° window
  // TODO check if this really works!
  angle = ((pos_x % 3) + (pos_y % 3) * 3) * 2.*G_PI / 9. / 4.;

  for (distance = 6; distance < 3 * 64; distance += 6)
    {
      // TODO precalculate these!
      xdiff = cos (angle) * distance;
      ydiff = sin (angle) * distance;

      permutation[0] = xdiff;
      permutation[1] = ydiff;
      permutation[2] = -xdiff;
      permutation[3] = -ydiff;

      for (direction = 0; direction < 4; direction++)
        {
          if (found[0][direction].found && found[1][direction].found)
            {
              continue;
            }
          else
            {
              guchar r, g, b, a;

              gint xtmp = permutation[direction] + pos_x;
              gint ytmp = permutation[(direction + 1) % 4] + pos_y;

              r = GET_PIXEL (big_cache, xtmp, ytmp, 0);
              g = GET_PIXEL (big_cache, xtmp, ytmp, 1);
              b = GET_PIXEL (big_cache, xtmp, ytmp, 2);
              a = GET_PIXEL (big_cache, xtmp, ytmp, 3);

              // check if it's foreground or background, depending on alpha
              // toggle = 0 equals foreground
              for (toggle = 0; toggle < 2; toggle++)
                {

                  // add gradient distance
                  if (!found[toggle][direction].found)
                    {
                      // TODO: Check if value is initialized to zero
                      // TODO: Is thir correct? I think we must use sobel...
                      found[toggle][direction].gradient += sqrt ((prevval[direction] - r) * (prevval[direction] - r)
                                                           + (prevval[direction + 1] - g) * (prevval[direction + 1] - g)
                                                           + (prevval[direction + 2] - b) * (prevval[direction + 2] - b));
                    }
                  if (a == (toggle == 0 ? 255 : 0) &&
                      !found[toggle][direction].found)
                    {
                      found[toggle][direction].color[0] = r;
                      found[toggle][direction].color[1] = g;
                      found[toggle][direction].color[2] = b;
                      found[toggle][direction].distance = distance;
                      found[toggle][direction].found = TRUE;
                      found[toggle][direction].x = xtmp;
                      found[toggle][direction].y = ytmp;
                    }
                }
            }
        }
    }

  {
    gfloat min = -1;
    gint minindexf = -1;
    gint minindexb = -1;

    gint foreground_direction, background_direction;

    gfloat best_alpha = -1;

    // calculate energy function for every fg/bg pair
    for (foreground_direction = 0; foreground_direction < 4; foreground_direction++)
      {
        if (found[0][foreground_direction].found)
          {
            for (background_direction = 0; background_direction < 4; background_direction++)
              {
                if (found[1][background_direction].found)
                  {
                    float current_alpha;

                    gfloat temp = objective_function (&(found[0][foreground_direction]),
                                                      &(found[1][background_direction]),
                                                      pos_x,
                                                      pos_y,
                                                      big_cache,
                                                      &current_alpha);
                    if (temp < min || min < 0)
                      {
                        best_alpha = current_alpha;
                        min = temp;
                        minindexf = foreground_direction;
                        minindexb = background_direction;
                      }
                  }

              }
          }
      }

    if (minindexf != -1 && minindexb != -1 && best_alpha != -1)
      {
        SearchStructure best_foreground = found[0][minindexf];
        SearchStructure best_background = found[1][minindexb];

        // test: do combination of best match
        // should return a very similar image as before
        entry->foreground[0] = best_foreground.color[0];
        entry->foreground[1] = best_foreground.color[1];
        entry->foreground[2] = best_foreground.color[2];

        entry->background[0] = best_background.color[0];
        entry->background[1] = best_background.color[1];
        entry->background[2] = best_background.color[2];

        entry->alpha = (1 - best_alpha) * 255;
        entry->pair_found = TRUE;

        entry->sigma_b_squared = calculate_variance (best_background.color, best_background.x, best_background.y, big_cache);
        entry->sigma_f_squared = calculate_variance (best_foreground.color, best_foreground.x, best_foreground.y, big_cache);
      }

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

#ifdef DEBUG_PREDEFINED_MASK
static void
read_write_mask (TileManager* mask_layer, gboolean write)
{
  PixelRegion mask;
  PixelRegionIterator *pr;
  gint row, col;
  int width, height;

  FILE *fp;
  char* filename = "mask_out.jj";

  if (write)
    fp = fopen (filename, "wb"); /* b - binary mode */
  else
    fp = fopen (filename, "rb"); /* b - binary mode */

  width = tile_manager_width (mask_layer);
  height = tile_manager_height(mask_layer);

  pixel_region_init (&mask, mask_layer, 0, 0, width, height, write);

  g_return_if_fail (mask.bytes == 1); // TODO check if indexed etc...

  for (pr = pixel_regions_register (1, &mask);
       pr != NULL;
       pr = pixel_regions_process (pr))
    {
      guchar *mask_data = mask.data;

      for (row = 0; row < mask.h; row++)
        {
          guchar *m = mask_data;

          for (col = 0; col < mask.w; col++, ++m)
            {
              if (write)
                fwrite (m, 1, 1, fp);
              else
                fread (m, 1, 1, fp);
            }

          mask_data += mask.rowstride;
        }
    }

  fclose(fp);
}
#endif



#ifdef DEBUG_EXTENSION
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
#endif

static inline gboolean
check_closeness (guchar color[3], BigCache big_cache, gint x, gint y, guchar* result)
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
            color_distance * color_distance;
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
search_for_neighbours (BigCache big_cache, gint x, gint y, guchar* result)
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
          if (check_closeness (color, big_cache, x + radius, y + n, result))
            return;

          if (check_closeness (color, big_cache, x - radius, y + n, result))
            return;

          if (check_closeness (color, big_cache, x + n, y + radius, result))
            return;

          if (check_closeness (color, big_cache, x + n, y - radius, result))
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
  BigCache     big_cache;
  gint         tiles_x, tiles_y;
  Tile        *tile;

  gint         tx, ty, x, y;
  guchar      *pointer;
  //gpointer     foreach_args[5];
  //gint         loaded_tile_x, loaded_tile_y;

  gboolean     unknown;

  HashEntry   *first_entry = NULL;
  HashEntry   *previous_entry = NULL;

  static       GHashTable *unknown_hash = NULL;

  unknown_hash = g_hash_table_new(g_int64_hash, g_int64_equal);

  g_return_if_fail (state != NULL);
  g_return_if_fail (mask != NULL && tile_manager_bpp (mask) == 1);
  g_return_if_fail (x1 >= 0);
  g_return_if_fail (x2 > x1 && x2 <= tile_manager_width (mask));
  g_return_if_fail (y1 >= 0);
  g_return_if_fail (y2 > y1 && y2 <= tile_manager_height (mask));
  g_return_if_fail (smoothness >= 0);
  g_return_if_fail (progress_data == NULL || progress_callback != NULL);
  g_return_if_fail (TILE_WIDTH == 64 && TILE_HEIGHT == 64);
  g_return_if_fail (tile_manager_bpp (state->pixels) == 3);
  g_return_if_fail (tile_manager_bpp (result_layer) == 4);

  tiles_x = tile_manager_tiles_per_col (state->pixels);
  tiles_y = tile_manager_tiles_per_row (state->pixels);

#ifdef DEBUG_PREDEFINED_MASK
  read_write_mask (mask, DEBUG_PREDEFINED_MASK_WRITE);
#endif

  initialize_new_layer (state->pixels, working_layer, mask);

  for (ty = 0; ty < tiles_y; ty++)
    {
      for (tx = 0; tx < tiles_x; tx++)
        {
          guint   height_tile;
          guint   width_tile;

          load_big_cache (working_layer, big_cache, tx, ty, 1);

#ifdef IMAGE_DEBUG_PPM
          {
            static char buffer[100];

            snprintf (buffer, 100, "big_cache_tx_%i_ty_%i", tx, ty);
            debug_cache (buffer, big_cache, 1);
          }
#endif

          tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);
          pointer = tile_data_pointer (tile, 0, 0);

          width_tile = tile_ewidth (tile);
          height_tile = tile_eheight (tile);

          for (y = 0; y < height_tile; y++)
            {
              for (x = 0; x < width_tile; x++, pointer += 4)
                {
                  pointer[0] = GET_PIXEL (big_cache, x, y, 0);
                  pointer[1] = GET_PIXEL (big_cache, x, y, 1);
                  pointer[2] = GET_PIXEL (big_cache, x, y, 2);

                  search_for_neighbours (big_cache, x, y, pointer + 3);

                  unknown = (GET_PIXEL (big_cache, x, y, 3) == 128);

                  if (unknown)
                    {
                      HashEntry *entry = g_slice_new (HashEntry);
                      // TODO: free this memory

                      entry->foreground[0] = 255;
                      entry->foreground[1] = 0;
                      entry->foreground[2] = 0;
                      entry->alpha = 255;
                      entry->pair_found = FALSE;
                      entry->this.coords.x = tx * 64 + x;
                      entry->this.coords.y = ty * 64 + y;

                      // TODO maybe look this up in image, instead of
                      // saving it redundantly in cache
                      entry->color[0] = pointer[0];
                      entry->color[1] = pointer[1];
                      entry->color[2] = pointer[2];

                      if (previous_entry != NULL)
                        previous_entry->next = entry->this;
                      previous_entry = entry;

                      if(first_entry == NULL)
                        first_entry = entry;

                      g_hash_table_insert (unknown_hash, &(entry->this), entry);
                    }
                }
            }

          tile_release (tile, TRUE);
        }
    }

  g_return_if_fail(previous_entry != NULL);

  // End the list with a null pointer
  previous_entry->next.value = 0;

#ifdef DEBUG_PHASE1
  update_mask (result_layer, mask);
  return;
#endif

  // Phase 2, loop over values in hash and fill them in
  {
    HashEntry *current = first_entry;
    gint current_tx = -1;
    gint current_ty = -1;

    while (current != NULL && current->next.value != 0)
      {
        search_neighborhood (current, &current_tx, &current_ty,
                             big_cache, result_layer);

        current = g_hash_table_lookup (unknown_hash, &(current->next));
      }
  }

#ifndef DEBUG_PHASE2
  // Phase 3, get better values from neighbours
  {
    HashEntry *current = first_entry;
    HashCache hash_cache;
    gint current_tx = -1;
    gint current_ty = -1;

    while (current != NULL && current->next.value != 0)
      {
        compare_neighborhood (current, &current_tx, &current_ty,
                              hash_cache, unknown_hash, big_cache);

        current = g_hash_table_lookup (unknown_hash, &(current->next));
      }
  }
#endif

  // Last phase, fill values from hash back into result layer
  for (ty = 0; ty < tiles_y; ty++)
    {
      for (tx = 0; tx < tiles_x; tx++)
        {
          guint height_tile;
          guint width_tile;

          tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);
          pointer = tile_data_pointer (tile, 0, 0);

          width_tile = tile_ewidth (tile);
          height_tile = tile_eheight (tile);

          for (y = 0; y < height_tile; y++)
            {
              for (x = 0; x < width_tile; x++, pointer += 4)
                {
                  HashEntry *current;
                  HashAddress address;
                  address.coords.x = tx * 64 + x;
                  address.coords.y = ty * 64 + y;

                  current = g_hash_table_lookup (unknown_hash, &address);

                  if (current != NULL)
                    {
#ifndef DEBUG_PHASE2
                      pointer[0] = current->foreground_refined[0];
                      pointer[1] = current->foreground_refined[1];
                      pointer[2] = current->foreground_refined[2];
                      pointer[3] = current->alpha_refined;
#else
                      pointer[0] = current->foreground[0];
                      pointer[1] = current->foreground[1];
                      pointer[2] = current->foreground[2];
                      pointer[3] = current->alpha;
#endif
                      //pointer[3] = 255;
                    }
                }
            }
        }
    }
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
