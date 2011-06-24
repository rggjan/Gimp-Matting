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

#include <gtk/gtk.h>

#include "widgets-types.h"

#include "gimpdeviceinfo.h"
#include "gimpdeviceinfo-coords.h"


static const GimpCoords default_coords = GIMP_COORDS_DEFAULT_VALUES;


/*  public functions  */

gboolean
gimp_device_info_get_event_coords (GimpDeviceInfo *info,
                                   GdkWindow      *window,
                                   const GdkEvent *event,
                                   GimpCoords     *coords)
{
  gdouble x;

  if (gdk_event_get_axis (event, GDK_AXIS_X, &x))
    {
      *coords = default_coords;

      coords->x = x;
      gdk_event_get_axis (event, GDK_AXIS_Y, &coords->y);

      if (gdk_event_get_axis (event, GDK_AXIS_PRESSURE, &coords->pressure))
        {
          coords->pressure = gimp_device_info_map_axis (info,
                                                        GDK_AXIS_PRESSURE,
                                                        coords->pressure);
        }

      if (gdk_event_get_axis (event, GDK_AXIS_XTILT, &coords->xtilt))
        {
          coords->xtilt = gimp_device_info_map_axis (info,
                                                     GDK_AXIS_XTILT,
                                                     coords->xtilt);
        }

      if (gdk_event_get_axis (event, GDK_AXIS_YTILT, &coords->ytilt))
        {
          coords->ytilt = gimp_device_info_map_axis (info,
                                                     GDK_AXIS_YTILT,
                                                     coords->ytilt);
        }

      if (gdk_event_get_axis (event, GDK_AXIS_WHEEL, &coords->wheel))
        {
          coords->wheel = gimp_device_info_map_axis (info,
                                                     GDK_AXIS_WHEEL,
                                                     coords->wheel);
        }

      return TRUE;
    }

  gimp_device_info_get_device_coords (info, window, coords);

  return FALSE;
}

void
gimp_device_info_get_device_coords (GimpDeviceInfo *info,
                                    GdkWindow      *window,
                                    GimpCoords     *coords)
{
  gdouble axes[GDK_AXIS_LAST];

  *coords = default_coords;

  gdk_device_get_state (info->device, window, axes, NULL);

  gdk_device_get_axis (info->device, axes, GDK_AXIS_X, &coords->x);
  gdk_device_get_axis (info->device, axes, GDK_AXIS_Y, &coords->y);

  if (gdk_device_get_axis (info->device,
                           axes, GDK_AXIS_PRESSURE, &coords->pressure))
    {
      coords->pressure = gimp_device_info_map_axis (info,
                                                    GDK_AXIS_PRESSURE,
                                                    coords->pressure);
    }

  if (gdk_device_get_axis (info->device,
                           axes, GDK_AXIS_XTILT, &coords->xtilt))
    {
      coords->xtilt = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_XTILT,
                                                 coords->xtilt);
    }

  if (gdk_device_get_axis (info->device,
                           axes, GDK_AXIS_YTILT, &coords->ytilt))
    {
      coords->ytilt = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_YTILT,
                                                 coords->ytilt);
    }

  if (gdk_device_get_axis (info->device,
                           axes, GDK_AXIS_WHEEL, &coords->wheel))
    {
      coords->wheel = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_WHEEL,
                                                 coords->wheel);
    }
}

void
gimp_device_info_get_time_coords (GimpDeviceInfo *info,
                                  GdkTimeCoord   *event,
                                  GimpCoords     *coords)
{
  *coords = default_coords;

  gdk_device_get_axis (info->device, event->axes, GDK_AXIS_X, &coords->x);
  gdk_device_get_axis (info->device, event->axes, GDK_AXIS_Y, &coords->y);

  /*  CLAMP() the return value of each *_get_axis() call to be safe
   *  against buggy XInput drivers.
   */

  if (gdk_device_get_axis (info->device,
                           event->axes, GDK_AXIS_PRESSURE, &coords->pressure))
    {
      coords->pressure = gimp_device_info_map_axis (info,
                                                    GDK_AXIS_PRESSURE,
                                                    coords->pressure);
    }

  if (gdk_device_get_axis (info->device,
                           event->axes, GDK_AXIS_XTILT, &coords->xtilt))
    {
      coords->xtilt = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_XTILT,
                                                 coords->xtilt);
    }

  if (gdk_device_get_axis (info->device,
                           event->axes, GDK_AXIS_YTILT, &coords->ytilt))
    {
      coords->ytilt = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_YTILT,
                                                 coords->ytilt);
    }

  if (gdk_device_get_axis (info->device,
                           event->axes, GDK_AXIS_WHEEL, &coords->wheel))
    {
      coords->wheel = gimp_device_info_map_axis (info,
                                                 GDK_AXIS_WHEEL,
                                                 coords->wheel);
    }
}

gboolean
gimp_device_info_get_event_state (GimpDeviceInfo  *info,
                                  GdkWindow       *window,
                                  const GdkEvent  *event,
                                  GdkModifierType *state)
{
  if (gdk_event_get_state (event, state))
    return TRUE;

  gimp_device_info_get_device_state (info, window, state);

  return FALSE;
}

void
gimp_device_info_get_device_state (GimpDeviceInfo  *info,
                                   GdkWindow       *window,
                                   GdkModifierType *state)
{
  gdk_device_get_state (info->device, window, NULL, state);
}
