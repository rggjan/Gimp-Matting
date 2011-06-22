/*
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

#include <glib-object.h>
#include <glib-2.0/glib/gprintf.h>

#include "libgimpmath/gimpmath.h"

#include "core/core-types.h"

#include "pixel-region.h"
#include "tile.h"

#include "tile-manager.h"
#include "siox.h"

//#define IMAGE_DEBUG_PPM

// Stop after certain phase
#define DEBUG_PHASE 3

// 0 = foreground, 1 = background, 2 = alpha
// #define DEBUG_SHOW_SPECIAL 0

// TRUE for writing
// FALSE for reading
// undefined for normal mode
#define DEBUG_PREDEFINED_MASK_WRITE FALSE

#ifdef IMAGE_DEBUG_PPM
#include "stdio.h"
#include "stdlib.h"
#endif

// TODO move slider before painting doesn't work!

#define LAMBDA (10)
#define SEARCH_RADIUS (10)
#define MATTING_SQUARED_COLOR_DISTANCE (25)

const gfloat gauss[13][13] = {{1.7102e-06, 8.0986e-06, 2.8906e-05, 7.7761e-05, 0.00015767, 0.00024095, 0.00027754, 0.00024095, 0.00015767, 7.7761e-05, 2.8906e-05, 8.0986e-06, 1.7102e-06},
  {8.0986e-06, 3.8351e-05, 0.00013688, 0.00036824, 0.00074664, 0.001141, 0.0013143, 0.001141, 0.00074664, 0.00036824, 0.00013688, 3.8351e-05, 8.0986e-06},
  {2.8906e-05, 0.00013688, 0.00048856, 0.0013143, 0.0026649, 0.0040726, 0.0046911, 0.0040726, 0.0026649, 0.0013143, 0.00048856, 0.00013688, 2.8906e-05},
  {7.7761e-05, 0.00036824, 0.0013143, 0.0035357, 0.007169, 0.010956, 0.01262, 0.010956, 0.007169, 0.0035357, 0.0013143, 0.00036824, 7.7761e-05},
  {0.00015767, 0.00074664, 0.0026649, 0.007169, 0.014536, 0.022215, 0.025588, 0.022215, 0.014536, 0.007169, 0.0026649, 0.00074664, 0.00015767},
  {0.00024095, 0.001141, 0.0040726, 0.010956, 0.022215, 0.033949, 0.039104, 0.033949, 0.022215, 0.010956, 0.0040726, 0.001141, 0.00024095},
  {0.00027754, 0.0013143, 0.0046911, 0.01262, 0.025588, 0.039104, 0.045042, 0.039104, 0.025588, 0.01262, 0.0046911, 0.0013143, 0.00027754},
  {0.00024095, 0.001141, 0.0040726, 0.010956, 0.022215, 0.033949, 0.039104, 0.033949, 0.022215, 0.010956, 0.0040726, 0.001141, 0.00024095},
  {0.00015767, 0.00074664, 0.0026649, 0.007169, 0.014536, 0.022215, 0.025588, 0.022215, 0.014536, 0.007169, 0.0026649, 0.00074664, 0.00015767},
  {7.7761e-05, 0.00036824, 0.0013143, 0.0035357, 0.007169, 0.010956, 0.01262, 0.010956, 0.007169, 0.0035357, 0.0013143, 0.00036824, 7.7761e-05},
  {2.8906e-05, 0.00013688, 0.00048856, 0.0013143, 0.0026649, 0.0040726, 0.0046911, 0.0040726, 0.0026649, 0.0013143, 0.00048856, 0.00013688, 2.8906e-05},
  {8.0986e-06, 3.8351e-05, 0.00013688, 0.00036824, 0.00074664, 0.001141, 0.0013143, 0.001141, 0.00074664, 0.00036824, 0.00013688, 3.8351e-05, 8.0986e-06},
  {1.7102e-06, 8.0986e-06, 2.8906e-05, 7.7761e-05, 0.00015767, 0.00024095, 0.00027754, 0.00024095, 0.00015767, 7.7761e-05, 2.8906e-05, 8.0986e-06, 1.7102e-06}
};

#define GAUSS(x, y) (gauss[x+6][y+6])

typedef union
{
  guint32 value;
  struct
  {
    guint16 x; // TODO assert x not bigger than this!
    guint16 y;
  } coords;
} HashAddress;

struct HashEntry_
{
  guchar color[3]; // TODO remove

  guchar foreground[3];
  guchar background[3];

  guchar foreground_refined[3];
  guchar background_refined[3];

  guchar alpha;
  guchar alpha_refined;

  gfloat sigma_f_squared;
  gfloat sigma_b_squared;

  gfloat confidence;

  gboolean valid;

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

/* A struct that holds the MATTING current state */
struct _MattingState
{
  TileManager *pixels;
  TileManager *result_layer;
  TileManager *mask;

  gboolean enough_pixels;

  BigCache big_cache;
  HashCache hash_cache;

  gint x1, y1, x2, y2;
  gint tx, ty;

  gint width, height;
  //HashCache    hash_cache;
};

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

MattingState *
siox_init (TileManager  *pixels,
           const guchar *colormap,
           gint          offset_x,
           gint          offset_y,
           gint          x,
           gint          y,
           gint          width,
           gint          height)
{
  MattingState *state;

  g_return_val_if_fail (pixels != NULL, NULL);
  state = g_slice_new (MattingState);

  state->pixels   = pixels;
  state->enough_pixels = FALSE;

  return state;
}

static void
load_big_cache (TileManager *source, TileManager *mask, BigCache big_cache,
                gint tx, gint ty, gint radius)
{
  gint    xdiff;
  gint    ydiff;
  gint    x, y;
  gint    bx, by;
  gint    width_tile, height_tile;

  Tile   *src_tile = NULL;
  Tile   *mask_tile = NULL;
  guchar *src_pointer = NULL;
  guchar *mask_pointer = NULL;

  guchar mask_bpp = tile_manager_bpp(mask);

  g_return_if_fail ((tile_manager_bpp(source) == 4 && mask == NULL) ||
                    (tile_manager_bpp(source) == 3));

  for (ydiff = -radius; ydiff <= radius; ydiff++)
    {
      for (xdiff = -radius; xdiff <= radius; xdiff++)
        {
          src_tile = tile_manager_get_at (source, tx + xdiff, ty + ydiff, TRUE, FALSE);
          if (mask)
            mask_tile = tile_manager_get_at (mask, tx + xdiff, ty + ydiff, TRUE, FALSE);

          width_tile = 0;
          height_tile = 0;

          if (src_tile)
            {
              src_pointer = tile_data_pointer (src_tile, 0, 0);
              if (mask)
                mask_pointer = tile_data_pointer (mask_tile, 0, 0);

              width_tile = tile_ewidth (src_tile);
              height_tile = tile_eheight (src_tile);

              for (y = 0; y < height_tile; y++)
                {
                  by = ydiff * 64 + y;
                  for (x = 0; x < width_tile; x++)
                    {
                      bx = xdiff * 64 + x;

                      GET_PIXEL(big_cache, bx, by, 0) = src_pointer[0];
                      GET_PIXEL(big_cache, bx, by, 1) = src_pointer[1];
                      GET_PIXEL(big_cache, bx, by, 2) = src_pointer[2];
                      if (mask)
                        GET_PIXEL(big_cache, bx, by, 3) = mask_pointer[mask_bpp-1];
                      else
                        GET_PIXEL(big_cache, bx, by, 3) = src_pointer[3];

                      if (mask)
                        {
                          src_pointer += 3;
                          mask_pointer += mask_bpp;
                        }
                      else
                        {
                          src_pointer += 4;
                        }
                    }
                }

              tile_release (src_tile, FALSE);
              if (mask)
                tile_release (mask_tile, FALSE);
            }

          for (y = 0; y < 64; y++)
            {
              by = ydiff * 64 + y;
              for (x = width_tile; x < 64; x++)
                {
                  bx = xdiff * 64 + x;

                  //GET_PIXEL(big_cache, bx, by, 0) = 0;
                  //GET_PIXEL(big_cache, bx, by, 1) = 0;
                  //GET_PIXEL(big_cache, bx, by, 2) = 0;
                  GET_PIXEL(big_cache, bx, by, 3) = 128;
                }
            }

          for (y = height_tile; y < 64; y++)
            {
              by = ydiff * 64 + y;
              for (x = 0; x < width_tile; x++)
                {
                  bx = xdiff * 64 + x;

                  //GET_PIXEL(big_cache, bx, by, 0) = 0;
                  //GET_PIXEL(big_cache, bx, by, 1) = 0;
                  //GET_PIXEL(big_cache, bx, by, 2) = 0;
                  GET_PIXEL(big_cache, bx, by, 3) = 128;
                }
            }
        }
    }
}

typedef struct
{
  guchar color[3];
  gboolean found;
  gint distance;
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

  // Normalize, so that it is in the unit cube of colors
  return (PPx * PPx + PPy * PPy + PPz * PPz) / (255. * 255.);
}

// evaluate the energy function for found color fg/bg for pixel situated at x/y
static inline float
objective_function (SearchStructure *fg,
                    SearchStructure *bg,
                    gint x, gint xmin, gint xmax,
                    gint y, gint ymin, gint ymax,
                    float* best_alpha,
                    float pfp,
                    MattingState *state)
{
  gint xi, yi;
  float ap;

  float Np = 0;

  float finalAlpha;

  for (yi = ymin; yi <= ymax; yi++)
    {
      for (xi = xmin; xi <= ymax; xi++)
        {
          // TODO check borders
          guchar *P = &(GET_PIXEL (state->big_cache, x + xi, y + yi, 0));

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

static gfloat calculate_variance (guchar P[3], gint x, gint y, MattingState *state)
{
  gint yi;
  gint xi;
  guchar my_alpha = GET_PIXEL (state->big_cache, x, y, 3);

  gfloat sum = 0;
  gint num_found = 0;

  for (yi = -2; yi <= 2; yi++)
    {
      for (xi = -2; xi <= 2; xi++)
        {
          if (GET_PIXEL (state->big_cache, x + xi, y + yi, 3) == my_alpha &&
              state->tx * 64 + x + xi < state->width  && state->tx * 64 + x >= 0 &&
              state->ty * 64 + y + yi < state->height && state->ty * 64 + y >= 0)
            {
              sum += dist_squared(GET_PIXEL (state->big_cache, x + xi, y + yi, 0),
                                  GET_PIXEL (state->big_cache, x + xi, y + yi, 1),
                                  GET_PIXEL (state->big_cache, x + xi, y + yi, 2),
                                  P[0],
                                  P[1],
                                  P[2]);

              num_found++;
            }
        }
    }

  if (num_found == 0)
    {
      g_printf("num_found = 0... this should not happen!\n");
      return 0;
    }
  else
    {
      return sum / num_found;
    }
}

typedef struct
{
  gint x;
  gint y;
  gfloat diff;
} TopColor;

static void load_hash_cache (GHashTable* unknown_hash, HashCache hash_cache,
                             gint tx, gint ty, gint bordersize)
{

  gint x, y;
  gint xoff = tx * 64;
  gint yoff = ty * 64;

  for (y = -bordersize; y < bordersize + 64; y++)
    {
      for (x = -bordersize; x < bordersize + 64; x++)
        {
          HashAddress address;
          address.coords.x = x + xoff;
          address.coords.y = y + yoff;

          GET_ENTRY(hash_cache, x, y) = g_hash_table_lookup (unknown_hash, &address);
        }
    }
}

static void inline
compare_neighborhood (HashEntry* entry, GHashTable* unknown_hash,
                      MattingState *state)
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

  if (state->tx != tx || state->ty != ty)
    {
      load_hash_cache (unknown_hash, state->hash_cache, tx, ty, HASH_CACHE_EXTRA);

      state->tx = tx;
      state->ty = ty;
    }

  for (num = 0; num < 3; num++)
    {
      top3[num].diff = INFINITY;
    }

  // TODO: change this to radius! (not square)
  for (ydiff = -HASH_CACHE_EXTRA; ydiff <= HASH_CACHE_EXTRA; ydiff++)
    {
      for (xdiff = -HASH_CACHE_EXTRA; xdiff <= HASH_CACHE_EXTRA; xdiff++)
        {
          current = GET_ENTRY(state->hash_cache, pos_x + xdiff, pos_y + ydiff);

          if (current && current->valid)
            {
              gfloat current_alpha;

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
      gfloat colordiff, current_alpha;

      guchar new_fg[3] = {0, 0, 0};
      guchar new_bg[3] = {0, 0, 0};
      gfloat new_sigma_f_squared = 0;
      gfloat new_sigma_b_squared = 0;

      gint index;

      matches = matches > 3 ? 3 : matches;

      for (num = 0; num < matches; num++)
        {
          current = GET_ENTRY(state->hash_cache, top3[num].x, top3[num].y);
          for (index = 0; index < 3; index++)
            {
              // TODO: check if we should use ints here!
              new_fg[index] += floor(current->foreground[index] / matches);
              new_bg[index] += floor(current->background[index] / matches);
            }

          new_sigma_f_squared += current->sigma_f_squared / matches;
          new_sigma_b_squared += current->sigma_b_squared / matches;
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

      {
        gfloat mp;
        gint i;
        gboolean same = TRUE;

        mp = projection (new_fg,
                         new_bg,
                         entry->color,
                         &current_alpha);

        entry->alpha_refined = (1 - current_alpha) * 255;

        for (i = 0; i < 3; i++)
          {
            if (entry->foreground_refined[i] != entry->background_refined[i])
              {
                same = FALSE;
                break;
              }
          }

        if (same)
          {
            entry->confidence = 1e-8;
          }
        else
          {
            mp = projection (entry->foreground_refined,
                             entry->background_refined,
                             entry->color,
                             NULL);

            entry->confidence = exp(-LAMBDA * sqrt(mp));
          }
      }
      entry->valid = TRUE;
    }
}

static void inline
calculate_final_colors (gfloat gauss, gint x, gint y, gfloat alpha_orig,
                        gdouble fq[3], gdouble fd[3], gdouble bq[3], gdouble bd[3],
                        gdouble *meandiff_q, gdouble *meandiff_d,
                        gdouble *low_freq_alpha_q, gdouble *low_freq_alpha_d,
                        HashCache hash_cache, BigCache big_cache, gboolean is_middle)
{
  guchar fg[3], bg[3];
  gint i;
  gfloat alpha, confidence;
  gdouble weight;
  HashEntry* current;
// TODO final_confidence should not get NAN!!!
  current = GET_ENTRY(hash_cache, x, y);
  if(current == NULL)
    {
      alpha = GET_PIXEL (big_cache, x, y, 3);
      alpha = (alpha == 255 ? 1 : 0);

      for (i = 0; i < 3; i++)
        {
          if (alpha == 1)
            {
              fg[i] = GET_PIXEL (big_cache, x, y, i);
            }
          else
            {
              bg[i] = GET_PIXEL (big_cache, x, y, i);
            }
        }
      confidence = 1;
    }
  else
    {
      gfloat diff_weight;
      for (i = 0; i < 3; i++)
        {
          fg[i] = current->foreground_refined[i];
          bg[i] = current->background_refined[i];
        }
      alpha = current->alpha_refined / 255.;
      confidence = current->confidence;

      diff_weight = alpha * (1 - alpha) * confidence;
      *meandiff_d += diff_weight;
      *meandiff_q += diff_weight * sqrt(dist_squared(
                                          current->foreground_refined[0],
                                          current->foreground_refined[1],
                                          current->foreground_refined[2],
                                          current->background_refined[0],
                                          current->background_refined[1],
                                          current->background_refined[2]));
    }

  // TODO: Special Case for original pixel
  if (is_middle)
    {
      weight = gauss * confidence;
    }
  else
    {
      gfloat diff = alpha_orig - alpha;
      if (diff < 0)
        diff = -diff;

      weight = gauss * confidence * diff;
    }

  for (i = 0; i < 3; i++)
    {
      gdouble low_freq_weight = confidence * gauss + (current == NULL);

      fq[i] += weight * alpha * fg[i];
      fd[i] += weight * alpha;
      bq[i] += weight * (1 - alpha) * bg[i];
      bd[i] += weight * (1 - alpha);

      *low_freq_alpha_q += low_freq_weight * alpha;
      *low_freq_alpha_d += low_freq_weight;
    }
}

static void inline
local_smoothing (HashEntry* entry, gint *current_tx, gint* current_ty,
                 HashCache hash_cache, GHashTable* unknown_hash,
                 TileManager *working_layer, MattingState *state)
{
  gint pos_x, pos_y;
  gint tx, ty;
  gint i;

  gint xdiff, ydiff;
  gdouble fq[3] = {0, 0, 0};
  gdouble fd[3] = {0, 0, 0};
  gdouble bq[3] = {0, 0, 0};
  gdouble bd[3] = {0, 0, 0};

  gdouble meandiff_q = 0;
  gdouble meandiff_d = 0;

  gdouble low_freq_alpha_q = 0;
  gdouble low_freq_alpha_d = 0;

  // Load coordinates from entry
  pos_x = entry->this.coords.x;
  pos_y = entry->this.coords.y;

  tx = pos_x / 64;
  ty = pos_y / 64;

  pos_x = pos_x - 64 * tx;
  pos_y = pos_y - 64 * ty;

  if (*current_tx != tx || *current_ty != ty)
    {
      load_hash_cache (unknown_hash, hash_cache, tx, ty, 6);
//      load_big_cache (working_layer, big_cache, tx, ty, 1);

      *current_tx = tx;
      *current_ty = ty;
    }


  for (ydiff = -6; ydiff <= 6; ydiff++)
    {
      for (xdiff = -6; xdiff <= 6; xdiff++)
        {
          calculate_final_colors (GAUSS(xdiff, ydiff),
                                  pos_x + xdiff, pos_y + ydiff,
                                  entry->alpha_refined / 255.,
                                  fq, fd, bq, bd,
                                  &meandiff_q, &meandiff_d,
                                  &low_freq_alpha_q, &low_freq_alpha_d,
                                  hash_cache, state->big_cache,
                                  xdiff == 0 && ydiff == 0);
        }
    }

  for (i = 0; i < 3; i++)
    {
      entry->foreground[i] = fq[i] / fd[i];
      entry->background[i] = bq[i] / bd[i];
    }

  {
    gfloat current_alpha, mp;

    gdouble low_freq_alpha = low_freq_alpha_q / low_freq_alpha_d;

    gdouble final_confidence = sqrt(dist_squared(entry->foreground[0],
                                    entry->foreground[1],
                                    entry->foreground[2],
                                    entry->background[0],
                                    entry->background[1],
                                    entry->background[2]));

    if (meandiff_d != 0)
      {
        gdouble meandiff = meandiff_q / meandiff_d;
        final_confidence /= meandiff;
      }
    else
      {
        final_confidence = 1;
      }

    if (final_confidence > 1)
      final_confidence = 1;

    mp = projection (entry->foreground,
                     entry->background,
                     entry->color,
                     &current_alpha);
    final_confidence *= exp(-LAMBDA * sqrt(mp));

    //g_printf("low_freq_alpha: %f\n", low_freq_alpha);

    // TODO here can be NANS!!!
    if (!(final_confidence <= 1 || final_confidence >= 0))
      g_printf("Problem!: final_confidence: %f\n", final_confidence);

    entry->alpha = (final_confidence * (1 - current_alpha) + (1 - final_confidence) * low_freq_alpha) * 255;
    //entry->alpha = (1-current_alpha)*255;
  }
}

static void inline
search_neighborhood (HashEntry* entry, MattingState *state)
{
  gint pos_x, pos_y;
  gint orig_pos_x, orig_pos_y;
  gint tx, ty;
  gint xdiff, ydiff;
  gfloat min_gradients[2] = {INFINITY, INFINITY};
  gfloat gradients[4] = {0, 0, 0, 0};

  gint direction, toggle, distance;

  SearchStructure found[2][4];

  guchar prevval[4][3];
  double angle;
  gint permutation[4];

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
      prevval[distance][0] = entry->color[0];
      prevval[distance][1] = entry->color[1];
      prevval[distance][2] = entry->color[2];
    }

  // Load coordinates from entry
  orig_pos_x = entry->this.coords.x;
  orig_pos_y = entry->this.coords.y;

  tx = orig_pos_x / 64;
  ty = orig_pos_y / 64;

  pos_x = orig_pos_x - 64 * tx;
  pos_y = orig_pos_y - 64 * ty;

  if (state->tx != tx || state->ty != ty)
    {
      load_big_cache (state->pixels, state->result_layer, state->big_cache, tx, ty, 3);

#ifdef IMAGE_DEBUG_PPM
      {
        static char buffer[100];

        snprintf (buffer, 100, "pixels_result_%i_%i", tx, ty);
        debug_cache (buffer, state->big_cache, 3);
      }
#endif

      state->tx = tx;
      state->ty = ty;
    }

  // in a 9x9 window, we want to have values in a 90° window
  angle = ((orig_pos_x % 3) + (orig_pos_y % 3) * 3) * 2.*G_PI / 9. / 4.;

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

              r = GET_PIXEL (state->big_cache, xtmp, ytmp, 0);
              g = GET_PIXEL (state->big_cache, xtmp, ytmp, 1);
              b = GET_PIXEL (state->big_cache, xtmp, ytmp, 2);
              a = GET_PIXEL (state->big_cache, xtmp, ytmp, 3);

              gradients[direction] += sqrt(dist_squared(prevval[direction][0],
                                           prevval[direction][1],
                                           prevval[direction][2],
                                           r, g, b));

              prevval[direction][0] = r;
              prevval[direction][1] = g;
              prevval[direction][2] = b;

              // check if it's foreground or background, depending on alpha
              // toggle = 0 equals foreground
              for (toggle = 0; toggle < 2; toggle++)
                {
                  if (a == (toggle == 0 ? 255 : 0) &&
                      !found[toggle][direction].found)
                    {
                      if (gradients[direction] < min_gradients[toggle])
                        min_gradients[toggle] = gradients[direction];

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
    gfloat min = INFINITY;
    gint minindexf = -1;
    gint minindexb = -1;

    gint foreground_direction, background_direction;

    gfloat best_alpha = -1;
    gfloat pfp = 0.5;
    if (min_gradients[1] > 0 || min_gradients[0] > 0)
      pfp = min_gradients[1] / (min_gradients[0] + min_gradients[1]);

    // calculate energy function for every fg/bg pair
    for (foreground_direction = 0; foreground_direction < 4; foreground_direction++)
      {
        if (found[0][foreground_direction].found)
          {
            for (background_direction = 0; background_direction < 4; background_direction++)
              {
                if (found[1][background_direction].found)
                  {
                    float current_alpha, cost;
                    gint xmin = -1, ymin = -1, xmax = 1, ymax = 1;

                    if (orig_pos_x <= 0)
                      xmin++;
                    else if (orig_pos_x >= state->width - 1)
                      xmax--;
                    if (orig_pos_y <= 0)
                      ymin++;
                    else if (orig_pos_y >= state->height - 1)
                      ymax++;

                    cost = objective_function (&(found[0][foreground_direction]),
                                               &(found[1][background_direction]),
                                               pos_x, xmin, xmax,
                                               pos_y, ymin, ymax,
                                               &current_alpha,
                                               pfp,
                                               state);
                    if (cost < min)
                      {
                        best_alpha = current_alpha;
                        min = cost;
                        minindexf = foreground_direction;
                        minindexb = background_direction;
                      }
                  }

              }
          }
      }

    if (min != INFINITY)
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
        entry->valid = TRUE;

        entry->sigma_b_squared = calculate_variance (best_background.color, best_background.x, best_background.y, state);
        entry->sigma_f_squared = calculate_variance (best_foreground.color, best_foreground.x, best_foreground.y, state);
      }

    //printf("values: %i %i %i | %i %i %i | %i %i %i | %i %i %i\n", values[0], values[1], values[2], values[3], values[4], values[5], values[6], values[7], values[8], values[9], values[10], values[11]);
  }
}

static void
initialize_new_layer (TileManager* destination_layer,
                      TileManager* mask_layer,
                      MattingState *state)
{
  PixelRegion src, dest, mask;
  PixelRegionIterator *pr;
  gint row, col;
  gint width, height;

  width = state->x2 - state->x1;
  height = state->y2 - state->y1;

  pixel_region_init (&src, state->pixels, state->x1, state->y1, width, height, FALSE);
  pixel_region_init (&dest, destination_layer, state->x1, state->y1, width, height, TRUE);
  pixel_region_init (&mask, mask_layer, state->x1, state->y1, width, height, FALSE);

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
                {
                  d[3] = 255;
                }
              else if (m[0] == MATTING_USER_BACKGROUND)
                {
                  d[3] = 0;
                }
              else
                {
                  d[3] = 128;
                }
            }

          src_data += src.rowstride;
          dest_data += dest.rowstride;
          mask_data += mask.rowstride;
        }
    }
}

static gfloat
mask_percent_unknown (TileManager* mask_layer, MattingState* state)
{
  PixelRegion mask;
  PixelRegionIterator *pr;
  gint row, col;
  gint pixels_unknown = 0;
  guint width = state->x2 - state->x1;
  guint height = state->y2 - state->y1;
  guint pixels_total = width * height;

  pixel_region_init (&mask, mask_layer, state->x1, state->y1, width, height, TRUE);

  if (mask.bytes != 1)
    return 1;

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
              if (m[0] != MATTING_USER_FOREGROUND && m[0] != MATTING_USER_BACKGROUND)
                {
                  pixels_unknown++;
                }
            }

          mask_data += mask.rowstride;
        }
    }
  return (gfloat) pixels_unknown / pixels_total;
}

#ifdef DEBUG_PREDEFINED_MASK_WRITE
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

static void
update_mask (TileManager* result_layer,
             TileManager* mask_layer,
             MattingState *state)
{
  PixelRegion result, mask;
  PixelRegionIterator *pr;
  gint row, col;
  gint width, height;

  width = state->x2 - state->x1;
  height = state->y2 - state->y1; // TODO give this as argument

  pixel_region_init (&result, result_layer, state->x1, state->y1, width, height, FALSE);
  pixel_region_init (&mask, mask_layer, state->x1, state->y1, width, height, TRUE);

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
                  if (result > 128)
                    m[0] = MATTING_ALGO_FOREGROUND;
                  else
                    m[0] = MATTING_ALGO_BACKGROUND;
                }
            }

          result_data += result.rowstride;
          mask_data += mask.rowstride;
        }
    }
}

static inline gboolean
check_closeness (guchar color[3], BigCache big_cache, gint x, gint y, guchar* result)
{
  guchar value;
  gint   color_distance_sum;
  gint   color_distance;
  gint   i;

  value = GET_PIXEL (big_cache, x, y, 3);

  if (value == MATTING_USER_FOREGROUND || value == MATTING_USER_BACKGROUND)
    {
      color_distance_sum = 0;
      for (i = 0; i < 3; i++)
        {
          color_distance = (gint)color[i] - GET_PIXEL (big_cache, x, y, i);
          color_distance_sum += color_distance * color_distance;
        }

      if (color_distance_sum < MATTING_SQUARED_COLOR_DISTANCE)
        {
          if (value == MATTING_USER_FOREGROUND)
            *result = 255;
          else
            *result = 0;
          return TRUE;
        }
    }

  return FALSE;
}

static inline guchar
search_for_neighbours (gint x, gint y, guchar* color, MattingState *state)
{
  gint   radius, n;
  guchar alpha;

  alpha = GET_PIXEL (state->big_cache, x, y, 3);

  if (alpha == MATTING_USER_BACKGROUND)
    {
      return 0;
    }

  if (alpha == MATTING_USER_FOREGROUND)
    {
      return 255;
    }

  for (radius = 0; radius <= SEARCH_RADIUS; radius++)
    {
      for (n = -radius; n < radius; n++)
        {
          if (check_closeness (color, state->big_cache, x + radius, y + n, &alpha))
            return alpha;

          if (check_closeness (color, state->big_cache, x - radius, y + n, &alpha))
            return alpha;

          if (check_closeness (color, state->big_cache, x + n, y + radius, &alpha))
            return alpha;

          if (check_closeness (color, state->big_cache, x + n, y - radius, &alpha))
            return alpha;
        }
    }
  return 128;
}

void
siox_foreground_extract (MattingState       *state,
                         SioxRefinementType  refinement,
                         TileManager        * mask,
                         gint                x1,
                         gint                y1,
                         gint                x2,
                         gint                y2,
                         gfloat              start_percentage,
                         const gdouble       sensitivity[3],
                         gboolean            multiblob,
                         SioxProgressFunc    progress_callback,
                         gpointer            progress_data,
                         TileManager        * result_layer,
                         TileManager        * working_layer)
{
  Tile        *tile;

  gint         tx, ty, x, y;
  guchar      *pointer;

  HashEntry   *first_entry = NULL;
  HashEntry   *previous_entry = NULL;

  static       GHashTable *unknown_hash = NULL;

  unknown_hash = g_hash_table_new(g_int_hash, g_int_equal); // TODO assert int = int32

  state->width = tile_manager_width (mask);
  state->height = tile_manager_height (mask);

  g_return_if_fail (state != NULL);
  g_return_if_fail (mask != NULL && tile_manager_bpp (mask) == 1);
  g_return_if_fail (x1 >= 0);
  g_return_if_fail (x2 > x1 && x2 <= state->width);
  g_return_if_fail (y1 >= 0);
  g_return_if_fail (y2 > y1 && y2 <= state->height);
  g_return_if_fail (start_percentage >= 0 && start_percentage <= 1);
  g_return_if_fail (progress_data == NULL || progress_callback != NULL);

  g_return_if_fail (tile_manager_bpp (state->pixels) == 3);
  g_return_if_fail (tile_manager_bpp (result_layer) == 4);

#ifdef DEBUG_PREDEFINED_MASK_WRITE
  read_write_mask (mask, DEBUG_PREDEFINED_MASK_WRITE);
#endif

  state->x1 = x1;
  state->y1 = y1;
  state->x2 = x2;
  state->y2 = y2;

  state->result_layer = result_layer;
  state->mask = mask;

  if (!state->enough_pixels)
    {
      gfloat unknown_percent = mask_percent_unknown (mask, state);

      g_printf("Unknown pixels: %f (%f)\n", unknown_percent, start_percentage);
      if (unknown_percent > (1 - start_percentage))
        {
          return;
        }
      state->enough_pixels = TRUE;
    }

  //initialize_new_layer (result_layer, mask, state);

  for (ty = y1 / 64; ty <= y2 / 64; ty++)
    {
      for (tx = x1 / 64; tx <= x2 / 64; tx++)
        {
          guint   height_tile;
          guint   width_tile;

          load_big_cache (state->pixels, mask, state->big_cache, tx, ty, 1);

#ifdef IMAGE_DEBUG_PPM
          {
            static char buffer[100];

            snprintf (buffer, 100, "pixels_mask_%i_%i", tx, ty);
            debug_cache (buffer, state->big_cache, 1);
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
                  guchar alpha;

                  pointer[0] = GET_PIXEL (state->big_cache, x, y, 0);
                  pointer[1] = GET_PIXEL (state->big_cache, x, y, 1);
                  pointer[2] = GET_PIXEL (state->big_cache, x, y, 2);
                  alpha = search_for_neighbours (x, y, pointer, state);
                  pointer[3] = alpha;

                  if (alpha == 128)
                    {
                      HashEntry *entry = g_slice_new (HashEntry);
                      // TODO: free this memory

                      // TODO: check if this can really be uncomented
                      //entry->foreground[0] = 255;
                      //entry->foreground[1] = 0;
                      //entry->foreground[2] = 0;
                      //entry->alpha = 255;
                      entry->valid = FALSE;
                      entry->this.coords.x = tx * 64 + x;
                      entry->this.coords.y = ty * 64 + y;

                      // TODO maybe look this up in image, instead of
                      // saving it redundantly in cache
                      entry->color[0] = pointer[0];
                      entry->color[1] = pointer[1];
                      entry->color[2] = pointer[2];
                      // TODO remove color from entry

                      // Linked list creation
                      if(first_entry == NULL)
                        first_entry = entry;
                      if (previous_entry != NULL)
                        previous_entry->next = entry->this;
                      previous_entry = entry;

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

  if (DEBUG_PHASE > 1)
    {
      // Phase 2, loop over values in hash and fill them in
      {
        HashEntry *current = first_entry;
        state->tx = -1;
        state->ty = -1;

        while (current != NULL && current->next.value != 0)
          {
            search_neighborhood (current, state);

            current = g_hash_table_lookup (unknown_hash, &(current->next));
          }
      }

      if (DEBUG_PHASE > 2)
        {
          // Phase 3, get better values from neighbours
          {
            HashEntry *current = first_entry;
            state->tx = -1;
            state->ty = -1;

            while (current != NULL && current->next.value != 0)
              {
                compare_neighborhood (current, unknown_hash, state);

                current = g_hash_table_lookup (unknown_hash, &(current->next));
              }
          }

          if (DEBUG_PHASE > 3)
            // Phase 4, get final color values
            {
              HashEntry *current = first_entry;
              HashCache hash_cache;
              gint current_tx = -1;
              gint current_ty = -1;

              while (current != NULL && current->next.value != 0)
                {
                  local_smoothing (current, &current_tx, &current_ty,
                                   hash_cache, unknown_hash, working_layer, state);

                  current = g_hash_table_lookup (unknown_hash, &(current->next));
                }
            }
        }
    }

  // Last phase, fill values from hash back into result layer
  {
    HashEntry *current = first_entry;
    HashEntry *tmp;
    Tile *tile = NULL;
    gint current_tx = -1;
    gint current_ty = -1;
    gint tx, ty, pos_x, pos_y;
    gint i;

    while (current != NULL && current->next.value != 0)
      {
#ifndef DEBUG_SHOW_SPECIAL
        if (current->valid)
#endif
          {
            tx = current->this.coords.x / 64;
            ty = current->this.coords.y / 64;

            pos_x = current->this.coords.x - 64 * tx;
            pos_y = current->this.coords.y - 64 * ty;

            if (current_tx != tx || current_ty != ty)
              {
                if (tile != NULL)
                  tile_release(tile, TRUE);

                tile = tile_manager_get_at (result_layer, tx, ty, TRUE, TRUE);

                current_tx = tx;
                current_ty = ty;
              }

            pointer = tile_data_pointer (tile, pos_x, pos_y);

#ifdef DEBUG_SHOW_SPECIAL
            if (current->valid)
              {
                pointer[3] = 255;

                for (i = 0; i < 3; i++)
                  {
                    switch (DEBUG_SHOW_SPECIAL)
                      {
                      case 0:
                        pointer[i] = current->foreground[i];
                        break;
                      case 1:
                        pointer[i] = current->background[i];
                        break;
                      case 2:
                        pointer[i] = current->alpha;
                        break;
                      }
                  }
              }
            else
              {
                pointer[0] = 255;
                pointer[1] = 0;
                pointer[2] = 0;
                pointer[3] = 255;
              }
#else
            for (i = 0; i < 3; i++)
              {
                if (DEBUG_PHASE == 3)
                  pointer[i] = current->foreground_refined[i];
                else
                  pointer[i] = current->foreground[i];
              }
            if (DEBUG_PHASE == 3)
              pointer[3] = current->alpha_refined;
            else
              pointer[3] = current->alpha;
#endif
          }

        tmp = g_hash_table_lookup (unknown_hash, &(current->next));
        g_slice_free(HashEntry, current);
        current = tmp;
      }
  }

  /*
    for (ty = y1 / 64; ty <= y2 / 64; ty++)
      {
        for (tx = x1 / 64; tx <= x2 / 64; tx++)
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
  #ifdef DEBUG_PHASE3
                        pointer[0] = current->foreground_refined[0];
                        pointer[1] = current->foreground_refined[1];
                        pointer[2] = current->foreground_refined[2];
  #ifdef DEBUG_SHOW_SPECIAL
                        if (DEBUG_SHOW_SPECIAL == 1)
                          {
                            pointer[0] = current->foreground_refined[0];
                            pointer[1] = current->foreground_refined[1];
                            pointer[2] = current->foreground_refined[2];
                          }
                        else if (DEBUG_SHOW_SPECIAL == 2)
                          {
                            pointer[0] = current->background_refined[0];
                            pointer[1] = current->background_refined[1];
                            pointer[2] = current->background_refined[2];
                          }
                        else
                          {
                            pointer[0] = current->alpha_refined;
                            pointer[1] = current->alpha_refined;
                            pointer[2] = current->alpha_refined;
                          }

                        pointer[3] = 255;
  #else
    pointer[3] = current->alpha_refined;
  #endif // DEBUG_SHOW_SPECIAL
  #else
    pointer[0] = current->foreground[0];
    pointer[1] = current->foreground[1];
    pointer[2] = current->foreground[2];
  #ifdef DEBUG_SHOW_SPECIAL
    if (DEBUG_SHOW_SPECIAL == 1)
      {
        pointer[0] = current->foreground[0];
        pointer[1] = current->foreground[1];
        pointer[2] = current->foreground[2];
      }
    else if (DEBUG_SHOW_SPECIAL == 2)
      {
        pointer[0] = current->background[0];
        pointer[1] = current->background[1];
        pointer[2] = current->background[2];
      }
    else
      {
        pointer[0] = current->alpha;
        pointer[1] = current->alpha;
        pointer[2] = current->alpha;
      }
    pointer[3] = 255;
  #else
    pointer[3] = current->alpha;
  #endif // DEBUG_SHOW_SPECIAL
  #endif
                      }
                    else
                      {
  #ifdef DEBUG_SHOW_SPECIAL
                        pointer[3] = 0;
  #endif
                      }
                  }
              }
          }
      }*/

  update_mask (result_layer, mask, state);
}

void
siox_done (MattingState * state)
{
  g_return_if_fail (state != NULL);

  g_slice_free (MattingState, state);
}
