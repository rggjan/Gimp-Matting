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

#include "libgimpcolor/gimpcolor.h"
#include "libgimpconfig/gimpconfig.h"
#include "libgimpwidgets/gimpwidgets.h"

#include "tools-types.h"

//#include "base/matting.h"

#include "widgets/gimpwidgets-utils.h"

#include "gimpforegroundselectoptions.h"
#include "gimptooloptions-gui.h"

#include "gimp-intl.h"

enum
{
  PROP_0,
  PROP_CONTINUOUS,
  PROP_DRAW_MODE,
  PROP_STROKE_WIDTH,
  PROP_START_PERCENTAGE,
  PROP_MASK_COLOR,
};


static void   gimp_foreground_select_options_set_property (GObject      *object,
                                                           guint         property_id,
                                                           const GValue *value,
                                                           GParamSpec   *pspec);
static void   gimp_foreground_select_options_get_property (GObject      *object,
                                                           guint         property_id,
                                                           GValue       *value,
                                                           GParamSpec   *pspec);


G_DEFINE_TYPE (GimpForegroundSelectOptions, gimp_foreground_select_options,
               GIMP_TYPE_SELECTION_OPTIONS)


static void
gimp_foreground_select_options_class_init (GimpForegroundSelectOptionsClass *klass)
{
  GObjectClass *object_class = G_OBJECT_CLASS (klass);

  object_class->set_property = gimp_foreground_select_options_set_property;
  object_class->get_property = gimp_foreground_select_options_get_property;

  GIMP_CONFIG_INSTALL_PROP_BOOLEAN (object_class, PROP_CONTINUOUS,
                                    "continuous",
                                    N_("Do continuous updating of the mask"),
                                    TRUE,
                                    GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_ENUM (object_class, PROP_DRAW_MODE,
                                 "draw-mode",
                                 N_("Paint over areas to mark color values for "
                                    "inclusion or exclusion from selection"),
                                 GIMP_TYPE_MATTING_DRAW_MODE,
                                 GIMP_MATTING_DRAW_MODE_FOREGROUND,
                                 GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_INT (object_class, PROP_STROKE_WIDTH,
                                "stroke-width",
                                N_("Size of the brush used for refinements"),
                                1, 80, 18,
                                GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_INT (object_class, PROP_START_PERCENTAGE,
                                "start-percentage",
                                N_("What is the minimum percent of knokwn "
                                   "pixels to start with the algorithm"),
                                0, 100, 70,
                                GIMP_PARAM_STATIC_STRINGS);

  GIMP_CONFIG_INSTALL_PROP_ENUM (object_class, PROP_MASK_COLOR,
                                 "mask-color",
                                 N_("Color of selection preview mask"),
                                 GIMP_TYPE_CHANNEL_TYPE,
                                 GIMP_BLUE_CHANNEL,
                                 GIMP_PARAM_STATIC_STRINGS);
}

static void
gimp_foreground_select_options_init (GimpForegroundSelectOptions *options)
{
}

static void
gimp_foreground_select_options_set_property (GObject      *object,
                                             guint         property_id,
                                             const GValue *value,
                                             GParamSpec   *pspec)
{
  GimpForegroundSelectOptions *options = GIMP_FOREGROUND_SELECT_OPTIONS (object);

  switch (property_id)
    {
    case PROP_CONTINUOUS:
      options->continuous = g_value_get_boolean (value);
      break;

    case PROP_DRAW_MODE:
      options->draw_mode = g_value_get_enum (value);
      break;

    case PROP_STROKE_WIDTH:
      options->stroke_width = g_value_get_int (value);
      break;

    case PROP_START_PERCENTAGE:
      options->start_percentage = g_value_get_int (value);
      break;

    case PROP_MASK_COLOR:
      options->mask_color = g_value_get_enum (value);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

static void
gimp_foreground_select_options_get_property (GObject    *object,
                                             guint       property_id,
                                             GValue     *value,
                                             GParamSpec *pspec)
{
  GimpForegroundSelectOptions *options = GIMP_FOREGROUND_SELECT_OPTIONS (object);

  switch (property_id)
    {
    case PROP_CONTINUOUS:
      g_value_set_boolean (value, options->continuous);
      break;

    case PROP_DRAW_MODE:
      g_value_set_enum (value, options->draw_mode);
      break;

    case PROP_STROKE_WIDTH:
      g_value_set_int (value, options->stroke_width);
      break;

    case PROP_START_PERCENTAGE:
      g_value_set_int (value, options->start_percentage);
      break;

    case PROP_MASK_COLOR:
      g_value_set_enum (value, options->mask_color);
      break;

    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
      break;
    }
}

GtkWidget *
gimp_foreground_select_options_gui (GimpToolOptions *tool_options)
{
  GObject   *config = G_OBJECT (tool_options);
  GtkWidget *vbox   = gimp_tool_options_gui (tool_options);
  GtkWidget *hbox;
  GtkWidget *button;
  GtkWidget *frame;
  GtkWidget *scale;
  GtkWidget *label;
  GtkWidget *menu;
  GtkWidget *inner_frame;
  GtkWidget *table;
  gchar     *title;

  /*  single / multiple objects  */
  button = gimp_prop_check_button_new (config, "continuous", _("Continuous"));
  gtk_box_pack_start (GTK_BOX (vbox), button, FALSE, FALSE, 0);
  gtk_widget_show (button);

  /*  foreground / background  */
  title = g_strdup_printf (_("Interactive refinement  (%s)"),
                           gimp_get_mod_string (GDK_CONTROL_MASK));

  frame = gimp_prop_enum_radio_frame_new (config, "draw-mode", title, -1, -1);
  g_free (title);

  gtk_box_pack_start (GTK_BOX (vbox), frame, FALSE, FALSE, 0);
  gtk_widget_show (frame);

  /*  stroke width  */
  inner_frame = gtk_vbox_new (FALSE, 0);
  gtk_box_pack_start (GTK_BOX (gtk_bin_get_child (GTK_BIN (frame))),
                      inner_frame, FALSE, FALSE, 2);
  gtk_widget_show (inner_frame);

  hbox = gtk_hbox_new (FALSE, 6);
  gtk_box_pack_start (GTK_BOX (inner_frame), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  label = gtk_label_new (_("Small brush"));
  gimp_label_set_attributes (GTK_LABEL (label),
                             PANGO_ATTR_STYLE, PANGO_STYLE_ITALIC,
                             PANGO_ATTR_SCALE,  PANGO_SCALE_SMALL,
                             -1);
  gtk_box_pack_start (GTK_BOX (hbox), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  label = gtk_label_new (_("Large brush"));
  gimp_label_set_attributes (GTK_LABEL (label),
                             PANGO_ATTR_STYLE, PANGO_STYLE_ITALIC,
                             PANGO_ATTR_SCALE,  PANGO_SCALE_SMALL,
                             -1);
  gtk_box_pack_end (GTK_BOX (hbox), label, FALSE, FALSE, 0);
  gtk_widget_show (label);

  scale = gimp_prop_hscale_new (config, "stroke-width", 1.0, 5.0, 0);
  gtk_scale_set_draw_value (GTK_SCALE (scale), FALSE);
  gtk_box_pack_start (GTK_BOX (inner_frame), scale, FALSE, FALSE, 0);
  gtk_widget_show (scale);

  /*  start_percentage  */
  table = gtk_table_new (2, 3, FALSE);
  gtk_table_set_row_spacings (GTK_TABLE (table), 2);
  gtk_table_set_col_spacings (GTK_TABLE (table), 2);
  gtk_box_pack_start (GTK_BOX (vbox), table, FALSE, FALSE, 0);
  gtk_widget_show (table);

  scale = gimp_prop_hscale_new (config, "start-percentage", 0.1, 1.0, 0);
  gtk_scale_set_value_pos (GTK_SCALE (scale), GTK_POS_RIGHT);
  gimp_table_attach_aligned (GTK_TABLE (table), 0, 0,
                             _("Start Percentage:"), 0.0, 0.5, scale, 2, FALSE);

  /*  mask color */
  menu = gimp_prop_enum_combo_box_new (config, "mask-color",
                                       GIMP_RED_CHANNEL, GIMP_BLUE_CHANNEL);
  gimp_table_attach_aligned (GTK_TABLE (table), 0, 1,
                             _("Preview color:"), 0.0, 0.5, menu, 2, FALSE);

  return vbox;
}

void
gimp_foreground_select_options_get_mask_color (GimpForegroundSelectOptions *options,
                                               GimpRGB                     *color)
{
  g_return_if_fail (GIMP_IS_FOREGROUND_SELECT_OPTIONS (options));
  g_return_if_fail (color != NULL);

  switch (options->mask_color)
    {
    case GIMP_RED_CHANNEL:
      gimp_rgba_set (color, 1, 0, 0, 0.5);
      break;

    case GIMP_GREEN_CHANNEL:
      gimp_rgba_set (color, 0, 1, 0, 0.5);
      break;

    case GIMP_BLUE_CHANNEL:
      gimp_rgba_set (color, 0, 0, 1, 0.5);
      break;

    default:
      g_warn_if_reached ();
      break;
    }
}
