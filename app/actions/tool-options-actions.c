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

#include "libgimpwidgets/gimpwidgets.h"

#include "actions-types.h"

#include "core/gimp.h"
#include "core/gimpcontext.h"
#include "core/gimplist.h"
#include "core/gimptoolinfo.h"

#include "widgets/gimpactiongroup.h"
#include "widgets/gimphelp-ids.h"

#include "tool-options-actions.h"
#include "tool-options-commands.h"

#include "gimp-intl.h"


/*  local function prototypes  */

static void tool_options_actions_update_presets (GimpActionGroup *group,
                                                 const gchar     *action_prefix,
                                                 GCallback        callback,
                                                 const gchar     *help_id,
                                                 GimpContainer   *presets);


/*  global variables  */

static const GimpActionEntry tool_options_actions[] =
{
  { "tool-options-popup", GIMP_STOCK_TOOL_OPTIONS,
    NC_("tool-options-action", "Tool Options Menu"), NULL, NULL, NULL,
    GIMP_HELP_TOOL_OPTIONS_DIALOG },

  { "tool-options-save-preset-menu", GTK_STOCK_SAVE,
    NC_("tool-options-action", "_Save Tool Preset"), "", NULL, NULL,
    GIMP_HELP_TOOL_OPTIONS_SAVE },

  { "tool-options-restore-preset-menu", GTK_STOCK_REVERT_TO_SAVED,
    NC_("tool-options-action", "_Restore Tool Preset"), "", NULL, NULL,
    GIMP_HELP_TOOL_OPTIONS_RESTORE },

  { "tool-options-edit-preset-menu", GTK_STOCK_EDIT,
    NC_("tool-options-action", "E_dit Tool Preset"), NULL, NULL, NULL,
    GIMP_HELP_TOOL_OPTIONS_EDIT },

  { "tool-options-delete-preset-menu", GTK_STOCK_DELETE,
    NC_("tool-options-action", "_Delete Tool Preset"), "", NULL, NULL,
    GIMP_HELP_TOOL_OPTIONS_DELETE },

  { "tool-options-save-new-preset", GTK_STOCK_NEW,
    NC_("tool-options-action", "_New Tool Preset..."), "", NULL,
    G_CALLBACK (tool_options_save_new_preset_cmd_callback),
    GIMP_HELP_TOOL_OPTIONS_SAVE },

  { "tool-options-reset", GIMP_STOCK_RESET,
    NC_("tool-options-action", "R_eset Tool Options"), "",
    NC_("tool-options-action", "Reset to default values"),
    G_CALLBACK (tool_options_reset_cmd_callback),
    GIMP_HELP_TOOL_OPTIONS_RESET },

  { "tool-options-reset-all", GIMP_STOCK_RESET,
    NC_("tool-options-action", "Reset _all Tool Options"), "",
    NC_("tool-options-action", "Reset all tool options"),
    G_CALLBACK (tool_options_reset_all_cmd_callback),
    GIMP_HELP_TOOL_OPTIONS_RESET }
};


/*  public functions  */

#define SET_VISIBLE(action,condition) \
        gimp_action_group_set_action_visible (group, action, (condition) != 0)
#define SET_HIDE_EMPTY(action,condition) \
        gimp_action_group_set_action_hide_empty (group, action, (condition) != 0)

void
tool_options_actions_setup (GimpActionGroup *group)
{
  gimp_action_group_add_actions (group, "tool-options-action",
                                 tool_options_actions,
                                 G_N_ELEMENTS (tool_options_actions));

  SET_HIDE_EMPTY ("tool-options-restore-preset-menu", FALSE);
  SET_HIDE_EMPTY ("tool-options-edit-preset-menu",    FALSE);
  SET_HIDE_EMPTY ("tool-options-delete-preset-menu",  FALSE);
}

void
tool_options_actions_update (GimpActionGroup *group,
                             gpointer         data)
{
  GimpContext  *context   = gimp_get_user_context (group->gimp);
  GimpToolInfo *tool_info = gimp_context_get_tool (context);

  SET_VISIBLE ("tool-options-save-preset-menu",    tool_info->presets);
  SET_VISIBLE ("tool-options-restore-preset-menu", tool_info->presets);
  SET_VISIBLE ("tool-options-edit-preset-menu",    tool_info->presets);
  SET_VISIBLE ("tool-options-delete-preset-menu",  tool_info->presets);

  tool_options_actions_update_presets (group, "tool-options-save-preset",
                                       G_CALLBACK (tool_options_save_preset_cmd_callback),
                                       GIMP_HELP_TOOL_OPTIONS_SAVE,
                                       tool_info->presets);

  tool_options_actions_update_presets (group, "tool-options-restore-preset",
                                       G_CALLBACK (tool_options_restore_preset_cmd_callback),
                                       GIMP_HELP_TOOL_OPTIONS_RESTORE,
                                       tool_info->presets);

  tool_options_actions_update_presets (group, "tool-options-edit-preset",
                                       G_CALLBACK (tool_options_edit_preset_cmd_callback),
                                       GIMP_HELP_TOOL_OPTIONS_EDIT,
                                       tool_info->presets);

  tool_options_actions_update_presets (group, "tool-options-delete-preset",
                                       G_CALLBACK (tool_options_delete_preset_cmd_callback),
                                       GIMP_HELP_TOOL_OPTIONS_DELETE,
                                       tool_info->presets);
}


/*  private function  */

static void
tool_options_actions_update_presets (GimpActionGroup *group,
                                     const gchar     *action_prefix,
                                     GCallback        callback,
                                     const gchar     *help_id,
                                     GimpContainer   *presets)
{
  GList *list;
  gint   n_children = 0;
  gint   i;

  for (i = 0; ; i++)
    {
      gchar     *action_name;
      GtkAction *action;

      action_name = g_strdup_printf ("%s-%03d", action_prefix, i);
      action = gtk_action_group_get_action (GTK_ACTION_GROUP (group),
                                            action_name);
      g_free (action_name);

      if (! action)
        break;

      gtk_action_group_remove_action (GTK_ACTION_GROUP (group), action);
    }

  if (presets)
    n_children = gimp_container_get_n_children (presets);

  if (n_children > 0)
    {
      GimpEnumActionEntry entry;

      entry.name           = NULL;
      entry.label          = NULL;
      entry.accelerator    = "";
      entry.tooltip        = NULL;
      entry.value          = 0;
      entry.value_variable = FALSE;
      entry.help_id        = help_id;

      for (list = GIMP_LIST (presets)->list, i = 0;
           list;
           list = g_list_next (list), i++)
        {
          GimpObject *options = list->data;

          entry.name     = g_strdup_printf ("%s-%03d", action_prefix, i);
          entry.label    = gimp_object_get_name (options);
          entry.stock_id = gimp_viewable_get_stock_id (GIMP_VIEWABLE (options));
          entry.value    = i;

          gimp_action_group_add_enum_actions (group, NULL, &entry, 1, callback);

          g_free ((gchar *) entry.name);
        }
    }
}

#undef SET_VISIBLE
#undef SET_HIDE_EMPTY
