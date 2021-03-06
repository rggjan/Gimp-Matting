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

#include "libgimpbase/gimpbase.h"
#include "libgimpconfig/gimpconfig.h"
#include "libgimpwidgets/gimpwidgets.h"

#include "tools-types.h"

#include "core/gimpmarshal.h"

#include "gimpalignoptions.h"
#include "gimptooloptions-gui.h"

#include "gimp-intl.h"


enum
{
  ALIGN_BUTTON_CLICKED,
  LAST_SIGNAL
};

enum
{
  PROP_0,
  PROP_ALIGN_REFERENCE,
  PROP_OFFSET_X,
  PROP_OFFSET_Y
};


static void   gimp_align_options_set_property (GObject      *object,
                                               guint         property_id,
                                               const GValue *value,
                                               GParamSpec   *pspec);
static void   gimp_align_options_get_property (GObject      *object,
                                               guint         property_id,
                                               GValue       *value,
                                               GParamSpec   *pspec);


G_DEFINE_TYPE (GimpAlignOptions, gimp_align_options, GIMP_TYPE_TOOL_OPTIONS)

#define parent_class gimp_selection_options_parent_class

static guint align_options_signals[LAST_SIGNAL] = { 0 };


static void
gimp_align_options_class_init (GimpAlignOptionsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property  = gimp_align_options_set_property;
  object_class->get_property  = gimp_align_options_get_property;

  klass->align_button_clicked = NULL;

  align_options_signals[ALIGN_BUTTON_CLICKED] =
    g_signal_new ("align-button-clicked",
                  G_TYPE_FROM_CLASS (klass),
                  G_SIGNAL_RUN_FIRST,
                  G_STRUCT_OFFSET (GimpAlignOptionsClass,
                                   align_button_clicked),
                  NULL, NULL,
                  gimp_marshal_VOID__ENUM,
                  G_TYPE_NONE, 1,
                  GIMP_TYPE_ALIGNMENT_TYPE);

  GIMP_CONFIG_INSTALL_PROP_ENUM (object_class, PROP_ALIGN_REFERENCE,
                                 "align-reference",
                                 N_("Reference image object a layer will be aligned on"),
                                 GIMP_TYPE_ALIGN_REFERENCE_TYPE,
                                 GIMP_ALIGN_REFERENCE_FIRST,
                                 GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_DOUBLE (object_class, PROP_OFFSET_X,
                                   "offset-x",
                                   N_("Horizontal offset for distribution"),
                                   -GIMP_MAX_IMAGE_SIZE, GIMP_MAX_IMAGE_SIZE, 0,
                                   GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_DOUBLE (object_class, PROP_OFFSET_Y,
                                   "offset-y",
                                   N_("Vertical offset for distribution"),
                                   -GIMP_MAX_IMAGE_SIZE, GIMP_MAX_IMAGE_SIZE, 0,
                                   GIMP_PARAM_STATIC_STRINGS);
}

static void
gimp_align_options_init (GimpAlignOptions *options)
{
}

static void
gimp_align_options_set_property (GObject      *object,
                                 guint         property_id,
                                 const GValue *value,
                                 GParamSpec   *pspec)
{
  GimpAlignOptions *options = GIMP_ALIGN_OPTIONS (object);

  switch (property_id)
    {
    case PROP_ALIGN_REFERENCE:
      options->align_reference = g_value_get_enum (value);
      break;

    case PROP_OFFSET_X:
      options->offset_x = g_value_get_double (value);
      break;

    case PROP_OFFSET_Y:
      options->offset_y = g_value_get_double (value);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_align_options_get_property (GObject    *object,
                                 guint       property_id,
                                 GValue     *value,
                                 GParamSpec *pspec)
{
  GimpAlignOptions *options = GIMP_ALIGN_OPTIONS (object);

  switch (property_id)
    {
    case PROP_ALIGN_REFERENCE:
      g_value_set_enum (value, options->align_reference);
      break;

    case PROP_OFFSET_X:
      g_value_set_double (value, options->offset_x);
      break;

    case PROP_OFFSET_Y:
      g_value_set_double (value, options->offset_y);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_align_options_button_clicked (GtkButton        *button,
                                   GimpAlignOptions *options)
{
  GimpAlignmentType action;

  action = GPOINTER_TO_INT (g_object_get_data (G_OBJECT (button),
                                               "align-action"));

  g_signal_emit (options, align_options_signals[ALIGN_BUTTON_CLICKED], 0,
                 action);
}

static GtkWidget *
gimp_align_options_button_new (GimpAlignOptions  *options,
                               GimpAlignmentType  action,
                               GtkWidget         *parent,
                               const gchar       *tooltip)
{
  GtkWidget   *button;
  GtkWidget   *image;
  const gchar *stock_id = NULL;

  switch (action)
    {
    case GIMP_ALIGN_LEFT:
      stock_id = GIMP_STOCK_GRAVITY_WEST;
      break;
    case GIMP_ALIGN_HCENTER:
      stock_id = GIMP_STOCK_HCENTER;
      break;
    case GIMP_ALIGN_RIGHT:
      stock_id = GIMP_STOCK_GRAVITY_EAST;
      break;
    case GIMP_ALIGN_TOP:
      stock_id = GIMP_STOCK_GRAVITY_NORTH;
      break;
    case GIMP_ALIGN_VCENTER:
      stock_id = GIMP_STOCK_VCENTER;
      break;
    case GIMP_ALIGN_BOTTOM:
      stock_id = GIMP_STOCK_GRAVITY_SOUTH;
      break;
    case GIMP_ARRANGE_LEFT:
      stock_id = GIMP_STOCK_GRAVITY_WEST;
      break;
    case GIMP_ARRANGE_HCENTER:
      stock_id = GIMP_STOCK_HCENTER;
      break;
    case GIMP_ARRANGE_RIGHT:
      stock_id = GIMP_STOCK_GRAVITY_EAST;
      break;
    case GIMP_ARRANGE_TOP:
      stock_id = GIMP_STOCK_GRAVITY_NORTH;
      break;
    case GIMP_ARRANGE_VCENTER:
      stock_id = GIMP_STOCK_VCENTER;
      break;
    case GIMP_ARRANGE_BOTTOM:
      stock_id = GIMP_STOCK_GRAVITY_SOUTH;
      break;
    default:
      g_return_val_if_reached (NULL);
      break;
    }

  button = gtk_button_new ();
  gtk_widget_set_sensitive (button, FALSE);
  gtk_widget_show (button);

  image = gtk_image_new_from_stock (stock_id, GTK_ICON_SIZE_BUTTON);
  gtk_container_add (GTK_CONTAINER (button), image);
  gtk_widget_show (image);

  gtk_box_pack_start (GTK_BOX (parent), button, FALSE, FALSE, 0);
  gtk_widget_show (button);

  gimp_help_set_help_data (button, _("Align left edge of target"), NULL);

  g_object_set_data (G_OBJECT (button), "align-action",
                     GINT_TO_POINTER (action));
  g_signal_connect (button, "clicked",
                    G_CALLBACK (gimp_align_options_button_clicked),
                    options);

  return button;
}

GtkWidget *
gimp_align_options_gui (GimpToolOptions *tool_options)
{
  GObject          *config  = G_OBJECT (tool_options);
  GimpAlignOptions *options = GIMP_ALIGN_OPTIONS (tool_options);
  GtkWidget        *vbox    = gimp_tool_options_gui (tool_options);
  GtkWidget        *align_vbox;
  GtkWidget        *hbox;
  GtkWidget        *frame;
  GtkWidget        *label;
  GtkWidget        *spinbutton;
  GtkWidget        *combo;
  gint              n = 0;

  frame = gimp_frame_new (_("Align"));
  gtk_box_pack_start (GTK_BOX (vbox), frame, FALSE, FALSE, 0);
  gtk_widget_show (frame);

  align_vbox = gtk_vbox_new (FALSE, 6);
  gtk_container_add (GTK_CONTAINER (frame), align_vbox);
  gtk_widget_show (align_vbox);

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  frame = gimp_frame_new (_("Relative to:"));
  gtk_box_pack_start (GTK_BOX (align_vbox), frame, FALSE, FALSE, 0);
  gtk_widget_show (frame);

  combo = gimp_prop_enum_combo_box_new (config, "align-reference", 0, 0);
  gtk_container_add (GTK_CONTAINER (frame), combo);
  gtk_widget_show (combo);

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_LEFT, hbox,
                                   _("Align left edge of target"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_HCENTER, hbox,
                                   _("Align center of target"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_RIGHT, hbox,
                                   _("Align right edge of target"));

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_TOP, hbox,
                                   _("Align top edge of target"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_VCENTER, hbox,
                                   _("Align middle of target"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ALIGN_BOTTOM, hbox,
                                   _("Align bottom of target"));

  frame = gimp_frame_new (_("Distribute"));
  gtk_box_pack_start (GTK_BOX (vbox), frame, FALSE, FALSE, 0);
  gtk_widget_show (frame);

  align_vbox = gtk_vbox_new (FALSE, 6);
  gtk_container_add (GTK_CONTAINER (frame), align_vbox);
  gtk_widget_show (align_vbox);

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_LEFT, hbox,
                                   _("Distribute left edges of targets"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_HCENTER, hbox,
                                   _("Distribute horizontal centers of targets"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_RIGHT, hbox,
                                   _("Distribute right edges of targets"));

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_TOP, hbox,
                                   _("Distribute top edges of targets"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_VCENTER, hbox,
                                   _("Distribute vertical centers of targets"));

  options->button[n++] =
    gimp_align_options_button_new (options, GIMP_ARRANGE_BOTTOM, hbox,
                                   _("Distribute bottoms of targets"));

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (align_vbox), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  label = gtk_label_new (_("Offset:"));
  gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  spinbutton = gimp_prop_spin_button_new (config, "offset-x",
                                          1, 20, 0);
  gtk_box_pack_start (GTK_BOX (hbox), spinbutton, FALSE, FALSE, 0);
  gtk_widget_show (spinbutton);

  return vbox;
}
