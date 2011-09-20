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

#include <glib-object.h>

#include "base-types.h"

#include "tile.h"
#include "tile-cache.h"
#include "tile-manager.h"
#include "tile-manager-private.h"
#include "tile-rowhints.h"
#include "tile-swap.h"
#include "tile-private.h"


static void  tile_manager_allocate_tiles (TileManager *tm);

#ifdef TILE_PROFILING
extern gint tile_exist_peak;
extern gint tile_exist_count;
#endif

#ifdef GIMP_UNSTABLE
GList *tile_managers = NULL;
#endif


GType
gimp_tile_manager_get_type (void)
{
  static GType type = 0;

  if (! type)
    type = g_boxed_type_register_static ("TileManager",
                                         (GBoxedCopyFunc) tile_manager_ref,
                                         (GBoxedFreeFunc) tile_manager_unref);

  return type;
}

#ifdef GIMP_UNSTABLE
void
tile_manager_exit (void)
{
  if (tile_managers)
    {
      g_warning ("%d tile managers leaked", g_list_length (tile_managers));

      while (tile_managers)
        {
          g_printerr ("unref tile manager %p (%d x %d)\n",
                      tile_managers->data,
                      tile_manager_width (tile_managers->data),
                      tile_manager_height (tile_managers->data));

          tile_manager_unref (tile_managers->data);
        }
    }
}
#endif

static inline gint
tile_manager_get_tile_num (TileManager *tm,
                           gint         xpixel,
                           gint         ypixel)
{
  if ((xpixel < 0) || (xpixel >= tm->width) ||
      (ypixel < 0) || (ypixel >= tm->height))
    return -1;

  return (ypixel / TILE_HEIGHT) * tm->ntile_cols + (xpixel / TILE_WIDTH);
}

TileManager *
tile_manager_new (gint width,
                  gint height,
                  gint bpp)
{
  TileManager *tm;

  g_return_val_if_fail (width > 0 && height > 0, NULL);
  g_return_val_if_fail (bpp > 0 && bpp <= 4, NULL);

  tm = g_slice_new0 (TileManager);

  tm->ref_count   = 1;
  tm->width       = width;
  tm->height      = height;
  tm->bpp         = bpp;
  tm->ntile_rows  = (height + TILE_HEIGHT - 1) / TILE_HEIGHT;
  tm->ntile_cols  = (width  + TILE_WIDTH  - 1) / TILE_WIDTH;
  tm->cached_num  = -1;

#ifdef GIMP_UNSTABLE
  tile_managers = g_list_prepend (tile_managers, tm);
#endif

  return tm;
}

TileManager *
tile_manager_ref (TileManager *tm)
{
  g_return_val_if_fail (tm != NULL, NULL);

  tm->ref_count++;

  return tm;
}

void
tile_manager_unref (TileManager *tm)
{
  g_return_if_fail (tm != NULL);

  tm->ref_count--;

  if (tm->ref_count < 1)
    {
#ifdef GIMP_UNSTABLE
      tile_managers = g_list_remove (tile_managers, tm);
#endif

      if (tm->cached_tile)
        tile_release (tm->cached_tile, FALSE);

      if (tm->tiles)
        {
          gint ntiles = tm->ntile_rows * tm->ntile_cols;
          gint i;

          for (i = 0; i < ntiles; i++)
            tile_detach (tm->tiles[i], tm, i);

          g_free (tm->tiles);
        }

      g_slice_free (TileManager, tm);
    }
}

TileManager *
tile_manager_duplicate (TileManager *tm)
{
  TileManager *copy;
  gint         n_tiles;
  gint         i;

  g_return_val_if_fail (tm != NULL, NULL);

  copy = tile_manager_new (tm->width, tm->height, tm->bpp);

  tile_manager_allocate_tiles (copy);

  n_tiles = tm->ntile_rows * tm->ntile_cols;

  for (i = 0; i < n_tiles; i++)
    {
      Tile *tile;

      tile = tile_manager_get (tm, i, TRUE, FALSE);
      tile_manager_map (copy, i, tile);
      tile_release (tile, FALSE);
    }

  return copy;
}

void
tile_manager_set_validate_proc (TileManager      *tm,
                                TileValidateProc  proc,
                                gpointer          user_data)
{
  g_return_if_fail (tm != NULL);

  tm->validate_proc = proc;
  tm->user_data     = user_data;
}

Tile *
tile_manager_get_tile (TileManager *tm,
                       gint         xpixel,
                       gint         ypixel,
                       gboolean     wantread,
                       gboolean     wantwrite)
{
  g_return_val_if_fail (tm != NULL, NULL);

  return tile_manager_get (tm,
                           tile_manager_get_tile_num (tm, xpixel, ypixel),
                           wantread, wantwrite);
}

Tile *
tile_manager_get (TileManager *tm,
                  gint         tile_num,
                  gboolean     wantread,
                  gboolean     wantwrite)
{
  Tile *tile;
  gint  ntiles;

  g_return_val_if_fail (tm != NULL, NULL);

  ntiles = tm->ntile_rows * tm->ntile_cols;

  if ((tile_num < 0) || (tile_num >= ntiles))
    return NULL;

  if (! tm->tiles)
    tile_manager_allocate_tiles (tm);

  tile = tm->tiles[tile_num];

  if (G_UNLIKELY (wantwrite && ! wantread))
    g_warning ("WRITE-ONLY TILE... UNTESTED!");

#ifdef DEBUG_TILE_MANAGER
  if (G_UNLIKELY (tile->share_count && tile->write_count))
    g_printerr (">> MEEPITY %d,%d <<\n", tile->share_count, tile->write_count);
#endif

  if (wantread)
    {
      if (wantwrite)
        {
          if (tile_num == tm->cached_num)
            {
              tile_release (tm->cached_tile, FALSE);

              tm->cached_tile = NULL;
              tm->cached_num  = -1;
            }

          if (tile->share_count > 1)
            {
              /* Copy-on-write required */
              Tile *new = tile_new (tile->bpp);

              new->ewidth  = tile->ewidth;
              new->eheight = tile->eheight;
              new->valid   = tile->valid;

              new->size    = new->ewidth * new->eheight * new->bpp;
              new->data    = g_new (guchar, new->size);

#ifdef TILE_PROFILING
              tile_exist_count++;
              if (tile_exist_count > tile_exist_peak)
                tile_exist_peak = tile_exist_count;
#endif

              if (tile->rowhint)
                {
                  tile_allocate_rowhints (new);

                  memcpy (new->rowhint, tile->rowhint,
                          new->eheight * sizeof (TileRowHint));
                }

              if (tile->data)
                {
                  memcpy (new->data, tile->data, new->size);
                }
              else
                {
                  tile_lock (tile);
                  memcpy (new->data, tile->data, new->size);
                  tile_release (tile, FALSE);
                }

              tile_detach (tile, tm, tile_num);
              tile_attach (new, tm, tile_num);

              tile = new;
              tm->tiles[tile_num] = tile;
            }

	  /* must lock before marking dirty */
	  tile_lock (tile);
          tile->write_count++;
          tile->dirty = TRUE;
        }
      else
        {
#ifdef DEBUG_TILE_MANAGER
          if (G_UNLIKELY (tile->write_count))
            g_printerr ("STINK! r/o on r/w tile (%d)\n", tile->write_count);
#endif
          tile_lock (tile);
        }
    }

  return tile;
}

Tile *
tile_manager_get_at (TileManager *tm,
                     gint         tile_col,
                     gint         tile_row,
                     gboolean     wantread,
                     gboolean     wantwrite)
{
  g_return_val_if_fail (tm != NULL, NULL);

  if (tile_col < 0 || tile_col >= tm->ntile_cols ||
      tile_row < 0 || tile_row >= tm->ntile_rows)
    return NULL;

  return tile_manager_get (tm,
                           tile_row * tm->ntile_cols + tile_col,
                           wantread, wantwrite);
}

void
tile_manager_validate_tile (TileManager *tm,
                            Tile        *tile)
{
  g_return_if_fail (tm != NULL);
  g_return_if_fail (tile != NULL);

  tile->valid = TRUE;

  if (tm->validate_proc)
    {
      (* tm->validate_proc) (tm, tile, tm->user_data);
    }
  else
    {
      /*  Set the contents of the tile to empty  */
      memset (tile->data, 0, tile_size (tile));
    }

#ifdef DEBUG_TILE_MANAGER
  g_printerr ("%c", tm->user_data ? 'V' : 'v');
#endif
}

static void
tile_manager_allocate_tiles (TileManager *tm)
{
  Tile       **tiles;
  const gint   nrows       = tm->ntile_rows;
  const gint   ncols       = tm->ntile_cols;
  const gint   right_tile  = tm->width  - ((ncols - 1) * TILE_WIDTH);
  const gint   bottom_tile = tm->height - ((nrows - 1) * TILE_HEIGHT);
  gint         i, j, k;

  g_assert (tm->tiles == NULL);

  tiles = g_new (Tile *, nrows * ncols);

  for (i = 0, k = 0; i < nrows; i++)
    {
      for (j = 0; j < ncols; j++, k++)
        {
          Tile *new = tile_new (tm->bpp);

          tile_attach (new, tm, k);

          if (j == (ncols - 1))
            new->ewidth = right_tile;

          if (i == (nrows - 1))
            new->eheight = bottom_tile;

          new->size = new->ewidth * new->eheight * new->bpp;

          tiles[k] = new;
        }
    }

  tm->tiles = tiles;
}

static void
tile_manager_invalidate_tile (TileManager  *tm,
                              gint          tile_num)
{
  Tile *tile = tm->tiles[tile_num];

  if (! tile->valid)
    return;

  if (tile_num == tm->cached_num)
    {
      tile_release (tm->cached_tile, FALSE);

      tm->cached_tile = NULL;
      tm->cached_num  = -1;
    }

  if (tile->cached)
    tile_cache_flush (tile);

  if (G_UNLIKELY (tile->share_count > 1))
    {
      /* This tile is shared.  Replace it with a new invalid tile. */
      Tile *new = tile_new (tile->bpp);

      new->ewidth  = tile->ewidth;
      new->eheight = tile->eheight;
      new->size    = tile->size;

      tile_detach (tile, tm, tile_num);
      tile_attach (new, tm, tile_num);

      tile = new;
      tm->tiles[tile_num] = tile;
    }

  tile->valid = FALSE;

  if (tile->data)
    {
      g_free (tile->data);
      tile->data = NULL;
#ifdef TILE_PROFILING
      tile_exist_count--;
#endif
    }

  if (tile->swap_offset != -1)
    {
      /* If the tile is on disk, then delete its
       *  presence there.
       */
      tile_swap_delete (tile);
    }
}

static void
tile_manager_invalidate_pixel (TileManager  *tm,
                               gint          xpixel,
                               gint          ypixel)
{
  gint num = tile_manager_get_tile_num (tm, xpixel, ypixel);

  if (num < 0)
    return;

  tile_manager_invalidate_tile (tm, num);
}

void
tile_manager_map_tile (TileManager *tm,
                       gint         xpixel,
                       gint         ypixel,
                       Tile        *srctile)
{
  g_return_if_fail (tm != NULL);
  g_return_if_fail (srctile != NULL);

  tile_manager_map (tm,
                    tile_manager_get_tile_num (tm, xpixel, ypixel),
                    srctile);
}

void
tile_manager_map (TileManager *tm,
                  gint         tile_num,
                  Tile        *srctile)
{
  Tile *tile;

  g_return_if_fail (tm != NULL);
  g_return_if_fail (srctile != NULL);
  g_return_if_fail (tile_num >= 0);
  g_return_if_fail (tile_num < tm->ntile_rows * tm->ntile_cols);

  if (G_UNLIKELY (! tm->tiles))
    {
      g_warning ("%s: empty tile level - initializing", G_STRLOC);

      tile_manager_allocate_tiles (tm);
    }

  tile = tm->tiles[tile_num];

#ifdef DEBUG_TILE_MANAGER
  g_printerr (")");
#endif

  if (G_UNLIKELY (! srctile->valid))
    g_warning("%s: srctile not validated yet!  please report", G_STRLOC);

  if (G_UNLIKELY (tile->ewidth  != srctile->ewidth  ||
                  tile->eheight != srctile->eheight ||
                  tile->bpp     != srctile->bpp))
    {
      g_warning ("%s: nonconformant map (%p -> %p)",
                 G_STRLOC, srctile, tile);
    }

  tile_detach (tile, tm, tile_num);

#ifdef DEBUG_TILE_MANAGER
  g_printerr (">");
#endif

#ifdef DEBUG_TILE_MANAGER
  g_printerr (" [src:%p tm:%p tn:%d] ", srctile, tm, tile_num);
#endif

  tile_attach (srctile, tm, tile_num);

  tm->tiles[tile_num] = srctile;

#ifdef DEBUG_TILE_MANAGER
  g_printerr ("}\n");
#endif
}

void
tile_manager_invalidate_area (TileManager *tm,
                              gint         x,
                              gint         y,
                              gint         w,
                              gint         h)
{
  gint  i;
  gint  j;

  /*  if no tiles have been allocated, there's no need to invalidate any  */
  if (! tm->tiles)
    return;

  for (i = y; i < (y + h); i += (TILE_HEIGHT - (i % TILE_HEIGHT)))
    for (j = x; j < (x + w); j += (TILE_WIDTH - (j % TILE_WIDTH)))
      {
        tile_manager_invalidate_pixel (tm, j, i);
      }
}

gint
tile_manager_width (const TileManager *tm)
{
  g_return_val_if_fail (tm != NULL, 0);

  return tm->width;
}

gint
tile_manager_height (const TileManager *tm)
{
  g_return_val_if_fail (tm != NULL, 0);

  return tm->height;
}

gint
tile_manager_bpp (const TileManager *tm)
{
  g_return_val_if_fail (tm != NULL, 0);

  return tm->bpp;
}

gint64
tile_manager_get_memsize (const TileManager *tm,
                          gboolean           sparse)
{
  /*  the tile manager itself  */
  gint64 memsize = sizeof (TileManager);

  if (! tm)
    return 0;

  /*  the array of tiles  */
  memsize += (gint64) tm->ntile_rows * tm->ntile_cols * (sizeof (Tile) +
                                                         sizeof (gpointer));

  /*  the memory allocated for the tiles  */
  if (sparse)
    {
      if (tm->tiles)
        {
          Tile   **tiles = tm->tiles;
          gint64   size  = TILE_WIDTH * TILE_HEIGHT * tm->bpp;
          gint     i, j;

          for (i = 0; i < tm->ntile_rows; i++)
            for (j = 0; j < tm->ntile_cols; j++, tiles++)
              {
                if (tile_is_valid (*tiles))
                  memsize += size;
              }
        }
    }
  else
    {
      memsize += (gint64) tm->width * tm->height * tm->bpp;
    }

  return memsize;
}

static inline gint
tile_manager_locate_tile (TileManager *tm,
                          Tile        *tile)
{
  TileLink *tl;

  for (tl = tile->tlink; tl; tl = tl->next)
    {
      if (tl->tm == tm)
        break;
    }

  if (G_UNLIKELY (tl == NULL))
    {
      g_warning ("%s: tile not attached to manager", G_STRLOC);
      return 0;
    }

  return tl->tile_num;
}

void
tile_manager_get_tile_col_row (TileManager *tm,
                               Tile        *tile,
                               gint        *tile_col,
                               gint        *tile_row)
{
  gint tile_num;

  g_return_if_fail (tm != NULL);
  g_return_if_fail (tile != NULL);
  g_return_if_fail (tile_col != NULL && tile_row != NULL);

  tile_num = tile_manager_locate_tile (tm, tile);

  *tile_col = tile_num % tm->ntile_cols;
  *tile_row = tile_num / tm->ntile_cols;
}

void
tile_manager_get_tile_coordinates (TileManager *tm,
                                   Tile        *tile,
                                   gint        *x,
                                   gint        *y)
{
  gint tile_col;
  gint tile_row;

  g_return_if_fail (tm != NULL);
  g_return_if_fail (tile != NULL);
  g_return_if_fail (x != NULL && y != NULL);

  tile_manager_get_tile_col_row (tm, tile, &tile_col, &tile_row);

  *x = TILE_WIDTH  * tile_col;
  *y = TILE_HEIGHT * tile_row;
}

void
tile_manager_map_over_tile (TileManager *tm,
                            Tile        *tile,
                            Tile        *srctile)
{
  TileLink *tl;

  g_return_if_fail (tm != NULL);
  g_return_if_fail (tile != NULL);
  g_return_if_fail (srctile != NULL);

  for (tl = tile->tlink; tl; tl = tl->next)
    {
      if (tl->tm == tm)
        break;
    }

  if (G_UNLIKELY (tl == NULL))
    {
      g_warning ("%s: tile not attached to manager", G_STRLOC);
      return;
    }

  tile_manager_map (tm, tl->tile_num, srctile);
}

void
tile_manager_read_pixel_data (TileManager *tm,
                              gint         x1,
                              gint         y1,
                              gint         x2,
                              gint         y2,
                              guchar      *buffer,
                              guint        stride)
{
  guint x, y;

  for (y = y1; y <= y2; y += TILE_HEIGHT - (y % TILE_HEIGHT))
    for (x = x1; x <= x2; x += TILE_WIDTH - (x % TILE_WIDTH))
      {
        Tile         *tile = tile_manager_get_tile (tm, x, y, TRUE, FALSE);
        const guchar *s    = TILE_DATA_POINTER (tile, x, y);
        guchar       *d    = buffer + stride * (y - y1) + tm->bpp * (x - x1);
        guint         rows, cols;
        guint         srcstride;

        rows = tile->eheight - y % TILE_HEIGHT;
        if (rows > (y2 - y + 1))
          rows = y2 - y + 1;

        cols = tile->ewidth - x % TILE_WIDTH;
        if (cols > (x2 - x + 1))
          cols = x2 - x + 1;

        srcstride = tile->ewidth * tile->bpp;

        while (rows--)
          {
            memcpy (d, s, cols * tm->bpp);

            s += srcstride;
            d += stride;
          }

        tile_release (tile, FALSE);
      }
}

void
tile_manager_write_pixel_data (TileManager  *tm,
                               gint          x1,
                               gint          y1,
                               gint          x2,
                               gint          y2,
                               const guchar *buffer,
                               guint         stride)
{
  guint x, y;

  for (y = y1; y <= y2; y += TILE_HEIGHT - (y % TILE_HEIGHT))
    for (x = x1; x <= x2; x += TILE_WIDTH - (x % TILE_WIDTH))
      {
        Tile         *tile = tile_manager_get_tile (tm, x, y, TRUE, TRUE);
        const guchar *s    = buffer + stride * (y - y1) + tm->bpp * (x - x1);
        guchar       *d    = TILE_DATA_POINTER (tile, x, y);
        guint         rows, cols;
        guint         dststride;

        rows = tile->eheight - y % TILE_HEIGHT;
        if (rows > (y2 - y + 1))
          rows = y2 - y + 1;

        cols = tile->ewidth - x % TILE_WIDTH;
        if (cols > (x2 - x + 1))
          cols = x2 - x + 1;

        dststride = tile->ewidth * tile->bpp;

        while (rows--)
          {
            memcpy (d, s, cols * tm->bpp);

            s += stride;
            d += dststride;
          }

        tile_release (tile, TRUE);
      }
}

void
tile_manager_read_pixel_data_1 (TileManager *tm,
                                gint         x,
                                gint         y,
                                guchar      *buffer)
{
  const gint num = tile_manager_get_tile_num (tm, x, y);

  if (num < 0)
    return;

  if (num != tm->cached_num)    /* must fetch a new tile */
    {
      Tile *tile;

      if (tm->cached_tile)
        tile_release (tm->cached_tile, FALSE);

      tm->cached_num  = -1;
      tm->cached_tile = NULL;

      /*  use a temporary variable instead of assigning to
       *  tm->cached_tile directly to make sure tm->cached_num
       *  and tm->cached_tile are always in a consistent state.
       *  (the requested tile might be invalid and needs to be
       *  validated, which would call tile_manager_get() recursively,
       *  which in turn would clear the cached tile) See bug #472770.
       */
      tile = tile_manager_get (tm, num, TRUE, FALSE);

      tm->cached_num  = num;
      tm->cached_tile = tile;
    }

  {
    const guchar *src = TILE_DATA_POINTER (tm->cached_tile, x, y);

    switch (tm->bpp)
      {
      case 4:
        *buffer++ = *src++;
      case 3:
        *buffer++ = *src++;
      case 2:
        *buffer++ = *src++;
      case 1:
        *buffer++ = *src++;
      }
  }
}

void
tile_manager_write_pixel_data_1 (TileManager  *tm,
                    gint          x,
                    gint          y,
                    const guchar *buffer)
{
  Tile   *tile = tile_manager_get_tile (tm, x, y, TRUE, TRUE);
  guchar *dest = TILE_DATA_POINTER (tile, x, y);

  switch (tm->bpp)
    {
    case 4:
      *dest++ = *buffer++;
    case 3:
      *dest++ = *buffer++;
    case 2:
      *dest++ = *buffer++;
    case 1:
      *dest++ = *buffer++;
    }

  tile_release (tile, TRUE);
}
