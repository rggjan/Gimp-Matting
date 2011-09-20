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

/*
 * EXIF-handling code for the jpeg plugin.  May eventually be better
 * to move this stuff into libgimpbase or a new libgimpmetadata and
 * make it available for other plugins.
 */

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <setjmp.h>

#include <jpeglib.h>
#include <jerror.h>

#include <libexif/exif-content.h>
#include <libexif/exif-data.h>
#include <libexif/exif-utils.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

#include "jpeg-exif.h"
#include "gimpexif.h"
#include "jpeg.h"
#include "jpeg-settings.h"

#include "libgimp/stdplugins-intl.h"


#define THUMBNAIL_SIZE             128
#define JPEG_EXIF_ROTATE_PARASITE  "exif-orientation-rotate"


static gboolean  jpeg_exif_rotate_query_dialog (gint32 image_ID);


/*  Replacement for exif_data_new_from_file() to work around
 *  filename encoding problems (see bug #335391).
 */
ExifData *
jpeg_exif_data_new_from_file (const gchar  *filename,
                              GError      **error)
{
  ExifData    *data;
  GMappedFile *file;

  g_return_val_if_fail (filename != NULL, NULL);
  g_return_val_if_fail (error == NULL || *error == NULL, NULL);

  file = g_mapped_file_new (filename, FALSE, error);
  if (! file)
    return NULL;

  data = exif_data_new_from_data ((guchar *) g_mapped_file_get_contents (file),
                                  g_mapped_file_get_length (file));

  g_mapped_file_unref (file);

  return data;
}


gint
jpeg_exif_get_orientation (ExifData *exif_data)
{
  ExifEntry *entry;
  gint       byte_order = exif_data_get_byte_order (exif_data);

  /* get orientation and rotate image accordingly if necessary */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_ORIENTATION)))
    {
      return exif_get_short (entry->data, byte_order);
    }

  return 0;
}


gboolean
jpeg_exif_get_resolution (ExifData *exif_data,
                          gdouble  *xresolution,
                          gdouble  *yresolution,
                          gint     *unit)
{
  gboolean       success;
  ExifEntry     *entry;
  gint           byte_order;
  gdouble        xres;
  gdouble        yres;
  gint           ruint;
  ExifRational   r;

  success = FALSE;
  byte_order = exif_data_get_byte_order (exif_data);

  do
    {
      entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                      EXIF_TAG_X_RESOLUTION);
      if (!entry)
        break;

      r = exif_get_rational (entry->data, byte_order);
      if (r.denominator == 0.0)
        break;
      xres = r.numerator / r.denominator;

      entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                      EXIF_TAG_Y_RESOLUTION);
      if (!entry)
        break;

      r = exif_get_rational (entry->data, byte_order);
      if (r.denominator == 0.0)
        break;
      yres = r.numerator / r.denominator;

      entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                      EXIF_TAG_RESOLUTION_UNIT);
      if (!entry)
        break;

      ruint = exif_get_short (entry->data, byte_order);
      if ((ruint != 2) && /* inches */
          (ruint != 3))   /* centimetres */
        break;

      success = TRUE;
    }
  while (0);

  if (success)
    {
      *xresolution = xres;
      *yresolution = yres;
      *unit = ruint;
    }

  return success;
}


void
jpeg_setup_exif_for_save (ExifData     *exif_data,
                          const gint32  image_ID)
{
  ExifRational  r;
  gdouble       xres, yres;
  ExifEntry    *entry;
  gint          byte_order = exif_data_get_byte_order (exif_data);

  /* set orientation to top - left */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_ORIENTATION)))
    {
      exif_set_short (entry->data, byte_order, (ExifShort) 1);
    }

  /* set x and y resolution */
  gimp_image_get_resolution (image_ID, &xres, &yres);
  r.numerator =   xres;
  r.denominator = 1;
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_X_RESOLUTION)))
    {
      exif_set_rational (entry->data, byte_order, r);
    }
  r.numerator = yres;
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_Y_RESOLUTION)))
    {
      exif_set_rational (entry->data, byte_order, r);
    }

  /* set resolution unit, always inches */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_RESOLUTION_UNIT)))
    {
      exif_set_short (entry->data, byte_order, (ExifShort) 2);
    }

  /* set software to "GIMP" and include the version number */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_SOFTWARE)))
    {
      const gchar *name = "GIMP " GIMP_VERSION;

      entry->data = (guchar *) g_strdup (name);
      entry->size = strlen (name) + 1;
      entry->components = entry->size;
    }

  /* set the width and height */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_EXIF],
                                       EXIF_TAG_PIXEL_X_DIMENSION)))
    {
      exif_set_long (entry->data, byte_order,
                     (ExifLong) gimp_image_width (image_ID));
    }
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_EXIF],
                                       EXIF_TAG_PIXEL_Y_DIMENSION)))
    {
      exif_set_long (entry->data, byte_order,
                     (ExifLong) gimp_image_height (image_ID));
    }

  /*
   * set the date & time image was saved
   * note, date & time of original photo is stored elsewwhere, we
   * aren't losing it.
   */
  if ((entry = exif_content_get_entry (exif_data->ifd[EXIF_IFD_0],
                                       EXIF_TAG_DATE_TIME)))
    {
      /* small memory leak here */
      entry->data = NULL;
      exif_entry_initialize (entry, EXIF_TAG_DATE_TIME);
    }

  /* should set components configuration, don't know how */

  /*
   * remove entries that don't apply to jpeg
   * (may have come from tiff or raw)
   */
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_COMPRESSION);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_IMAGE_WIDTH);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_IMAGE_LENGTH);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_BITS_PER_SAMPLE);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_SAMPLES_PER_PIXEL);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_PHOTOMETRIC_INTERPRETATION);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_STRIP_OFFSETS);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_PLANAR_CONFIGURATION);
  gimp_exif_data_remove_entry (exif_data, EXIF_IFD_0, EXIF_TAG_YCBCR_SUB_SAMPLING);

  /* should set thumbnail attributes */
}

void
jpeg_exif_rotate (gint32 image_ID,
                  gint   orientation)
{
  switch (orientation)
    {
    case 1:  /* standard orientation, do nothing */
      break;

    case 2:  /* flipped right-left               */
      gimp_image_flip (image_ID, GIMP_ORIENTATION_HORIZONTAL);
      break;

    case 3:  /* rotated 180                      */
      gimp_image_rotate (image_ID, GIMP_ROTATE_180);
      break;

    case 4:  /* flipped top-bottom               */
      gimp_image_flip (image_ID, GIMP_ORIENTATION_VERTICAL);
      break;

    case 5:  /* flipped diagonally around '\'    */
      gimp_image_rotate (image_ID, GIMP_ROTATE_90);
      jpeg_swap_original_settings (image_ID);
      gimp_image_flip (image_ID, GIMP_ORIENTATION_HORIZONTAL);
      break;

    case 6:  /* 90 CW                            */
      gimp_image_rotate (image_ID, GIMP_ROTATE_90);
      jpeg_swap_original_settings (image_ID);
      break;

    case 7:  /* flipped diagonally around '/'    */
      gimp_image_rotate (image_ID, GIMP_ROTATE_90);
      jpeg_swap_original_settings (image_ID);
      gimp_image_flip (image_ID, GIMP_ORIENTATION_VERTICAL);
      break;

    case 8:  /* 90 CCW                           */
      gimp_image_rotate (image_ID, GIMP_ROTATE_270);
      jpeg_swap_original_settings (image_ID);
      break;

    default: /* shouldn't happen                 */
      break;
    }
}

void
jpeg_exif_rotate_query (gint32 image_ID,
                        gint   orientation)
{
  GimpParasite *parasite;
  gboolean      query = load_interactive;

  if (orientation < 2 || orientation > 8)
    return;

  parasite = gimp_get_parasite (JPEG_EXIF_ROTATE_PARASITE);

  if (parasite)
    {
      if (strncmp (gimp_parasite_data (parasite), "yes",
                   gimp_parasite_data_size (parasite)) == 0)
        {
          query = FALSE;
        }
      else if (strncmp (gimp_parasite_data (parasite), "no",
                        gimp_parasite_data_size (parasite)) == 0)
        {
          gimp_parasite_free (parasite);
          return;
        }

      gimp_parasite_free (parasite);
    }

  if (query && ! jpeg_exif_rotate_query_dialog (image_ID))
    return;

  jpeg_exif_rotate (image_ID, orientation);
}

static gboolean
jpeg_exif_rotate_query_dialog (gint32 image_ID)
{
  GtkWidget *dialog;
  GtkWidget *hbox;
  GtkWidget *vbox;
  GtkWidget *label;
  GtkWidget *toggle;
  GdkPixbuf *pixbuf;
  gint       response;

  dialog = gimp_dialog_new (_("Rotate Image?"), PLUG_IN_ROLE,
                            NULL, 0, NULL, NULL,

                            _("_Keep Orientation"), GTK_RESPONSE_CANCEL,
                            GIMP_STOCK_TOOL_ROTATE, GTK_RESPONSE_OK,

                            NULL);

  gtk_dialog_set_alternative_button_order (GTK_DIALOG (dialog),
                                           GTK_RESPONSE_OK,
                                           GTK_RESPONSE_CANCEL,
                                           -1);

  gtk_window_set_resizable (GTK_WINDOW (dialog), FALSE);
  gimp_window_set_transient (GTK_WINDOW (dialog));

  hbox = gtk_hbox_new (FALSE, 12);
  gtk_container_set_border_width (GTK_CONTAINER (hbox), 12);
  gtk_box_pack_start (GTK_BOX (gtk_dialog_get_content_area (GTK_DIALOG (dialog))),
                      hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  vbox = gtk_vbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (hbox), vbox, FALSE, FALSE, 0);
  gtk_widget_show (vbox);

  pixbuf = gimp_image_get_thumbnail (image_ID,
                                     THUMBNAIL_SIZE, THUMBNAIL_SIZE,
                                     GIMP_PIXBUF_SMALL_CHECKS);

  if (pixbuf)
    {
      GtkWidget *image;
      gchar     *name;

      image = gtk_image_new_from_pixbuf (pixbuf);
      g_object_unref (pixbuf);

      gtk_box_pack_start (GTK_BOX (vbox), image, FALSE, FALSE, 0);
      gtk_widget_show (image);

      name = gimp_image_get_name (image_ID);

      label = gtk_label_new (name);
      gtk_label_set_ellipsize (GTK_LABEL (label), PANGO_ELLIPSIZE_MIDDLE);
      gimp_label_set_attributes (GTK_LABEL (label),
                                 PANGO_ATTR_STYLE,  PANGO_STYLE_ITALIC,
                                 -1);
      gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 0);
      gtk_widget_show (label);

      g_free (name);
    }

  vbox = gtk_vbox_new (FALSE, 12);
  gtk_box_pack_start (GTK_BOX (hbox), vbox, TRUE, TRUE, 0);
  gtk_widget_show (vbox);

  label = g_object_new (GTK_TYPE_LABEL,
                        "label",   _("According to the EXIF data, "
                                     "this image is rotated."),
                        "wrap",    TRUE,
                        "justify", GTK_JUSTIFY_LEFT,
                        "xalign",  0.0,
                        "yalign",  0.5,
                        NULL);
  gimp_label_set_attributes (GTK_LABEL (label),
                             PANGO_ATTR_SCALE,  PANGO_SCALE_LARGE,
                             PANGO_ATTR_WEIGHT, PANGO_WEIGHT_BOLD,
                             -1);
  gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  label = g_object_new (GTK_TYPE_LABEL,
                        "label",   _("Would you like GIMP to rotate it "
                                     "into the standard orientation?"),
                        "wrap",    TRUE,
                        "justify", GTK_JUSTIFY_LEFT,
                        "xalign",  0.0,
                        "yalign",  0.5,
                        NULL);
  gtk_box_pack_start (GTK_BOX (vbox), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  toggle = gtk_check_button_new_with_mnemonic (_("_Don't ask me again"));
  gtk_box_pack_end (GTK_BOX (vbox), toggle, FALSE, FALSE, 0);
  gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON (toggle), FALSE);
  gtk_widget_show (toggle);

  response = gimp_dialog_run (GIMP_DIALOG (dialog));

  if (gtk_toggle_button_get_active (GTK_TOGGLE_BUTTON (toggle)))
    {
      GimpParasite *parasite;
      const gchar  *str = (response == GTK_RESPONSE_OK) ? "yes" : "no";

      parasite = gimp_parasite_new (JPEG_EXIF_ROTATE_PARASITE,
                                    GIMP_PARASITE_PERSISTENT,
                                    strlen (str), str);
      gimp_attach_parasite (parasite);
      gimp_parasite_free (parasite);
    }

  gtk_widget_destroy (dialog);

  return (response == GTK_RESPONSE_OK);
}
