/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995 Spencer Kimball and Peter Mattis
 *
 * gimpdisplayshell-grab.c
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

#include "display-types.h"

#include "widgets/gimpdeviceinfo.h"
#include "widgets/gimpdevices.h"
#include "widgets/gimpdevicemanager.h"

#include "gimpdisplay.h"
#include "gimpdisplayshell.h"
#include "gimpdisplayshell-grab.h"


gboolean
gimp_display_shell_pointer_grab (GimpDisplayShell *shell,
                                 const GdkEvent   *event,
                                 GdkEventMask      event_mask)
{
  GdkGrabStatus status;

  g_return_val_if_fail (GIMP_IS_DISPLAY_SHELL (shell), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);
  g_return_val_if_fail (shell->pointer_grabbed == FALSE, FALSE);

  status = gdk_pointer_grab (gtk_widget_get_window (shell->canvas),
                             FALSE, event_mask, NULL, NULL,
                             gdk_event_get_time (event));

  if (status == GDK_GRAB_SUCCESS)
    {
      shell->pointer_grabbed = TRUE;

      return TRUE;
    }
  else if (status == GDK_GRAB_ALREADY_GRABBED)
    {
      GimpDeviceManager *manager;
      GdkDisplay        *gdk_display;

      manager = gimp_devices_get_manager (shell->display->gimp);
      gdk_display = gtk_widget_get_display (GTK_WIDGET (shell));

      /*  EEK: trying to grab an extended device always returns
       *  ALREADY_GRABBED, so simply assume the grab succeeded anyway
       */
      if (gimp_device_manager_get_current_device (manager)->device !=
          gdk_display_get_core_pointer (gdk_display))
        {
          shell->pointer_grabbed = TRUE;

          return TRUE;
        }
    }

  g_printerr ("%s: gdk_pointer_grab failed with status %d\n",
              G_STRFUNC, status);

  return FALSE;
}

void
gimp_display_shell_pointer_ungrab (GimpDisplayShell *shell,
                                   const GdkEvent   *event)
{
  g_return_if_fail (GIMP_IS_DISPLAY_SHELL (shell));
  g_return_if_fail (event != NULL);
  g_return_if_fail (shell->pointer_grabbed == TRUE);

  gdk_display_pointer_ungrab (gtk_widget_get_display (shell->canvas),
                              gdk_event_get_time (event));

  shell->pointer_grabbed = FALSE;
}

gboolean
gimp_display_shell_keyboard_grab (GimpDisplayShell *shell,
                                  const GdkEvent   *event)
{
  GdkGrabStatus status;

  g_return_val_if_fail (GIMP_IS_DISPLAY_SHELL (shell), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);
  g_return_val_if_fail (shell->keyboard_grabbed == FALSE, FALSE);

  status = gdk_keyboard_grab (gtk_widget_get_window (shell->canvas),
                              FALSE,
                              gdk_event_get_time (event));

  if (status == GDK_GRAB_SUCCESS)
    {
      shell->keyboard_grabbed = TRUE;

      return TRUE;
    }

  g_printerr ("%s: gdk_keyboard_grab failed with status %d\n",
              G_STRFUNC, status);

  return FALSE;
}

void
gimp_display_shell_keyboard_ungrab (GimpDisplayShell *shell,
                                    const GdkEvent   *event)
{
  g_return_if_fail (GIMP_IS_DISPLAY_SHELL (shell));
  g_return_if_fail (event != NULL);
  g_return_if_fail (shell->keyboard_grabbed == TRUE);

  gdk_display_keyboard_ungrab (gtk_widget_get_display (shell->canvas),
                               gdk_event_get_time (event));

  shell->keyboard_grabbed = FALSE;
}
