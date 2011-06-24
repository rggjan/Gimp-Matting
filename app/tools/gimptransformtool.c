/* GIMP - The GNU Image Manipulation Program
 * Copyright (C) 1995-2001 Spencer Kimball, Peter Mattis, and others
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

#include <stdlib.h>

#include <gegl.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>

#include "libgimpmath/gimpmath.h"
#include "libgimpconfig/gimpconfig.h"
#include "libgimpwidgets/gimpwidgets.h"

#include "tools-types.h"

#include "base/boundary.h"
#include "base/tile-manager.h"

#include "core/gimp.h"
#include "core/gimpcontext.h"
#include "core/gimpchannel.h"
#include "core/gimpdrawable-transform.h"
#include "core/gimperror.h"
#include "core/gimpimage.h"
#include "core/gimpimage-undo.h"
#include "core/gimpimage-undo-push.h"
#include "core/gimpitem-linked.h"
#include "core/gimplayer.h"
#include "core/gimpprogress.h"
#include "core/gimptoolinfo.h"
#include "core/gimp-transform-utils.h"
#include "core/gimp-utils.h"

#include "vectors/gimpvectors.h"
#include "vectors/gimpstroke.h"

#include "widgets/gimpdialogfactory.h"

#include "display/gimpcanvasgroup.h"
#include "display/gimpdisplay.h"
#include "display/gimpdisplayshell.h"
#include "display/gimpdisplayshell-transform.h"
#include "display/gimptooldialog.h"

#include "gimptoolcontrol.h"
#include "gimpperspectivetool.h"
#include "gimptransformoptions.h"
#include "gimptransformtool.h"
#include "gimptransformtoolundo.h"

#include "gimp-intl.h"


#define RESPONSE_RESET  1
#define MIN_HANDLE_SIZE 6


static void      gimp_transform_tool_finalize               (GObject               *object);

static gboolean  gimp_transform_tool_initialize             (GimpTool              *tool,
                                                             GimpDisplay           *display,
                                                             GError               **error);
static void      gimp_transform_tool_control                (GimpTool              *tool,
                                                             GimpToolAction         action,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_button_press           (GimpTool              *tool,
                                                             const GimpCoords      *coords,
                                                             guint32                time,
                                                             GdkModifierType        state,
                                                             GimpButtonPressType    press_type,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_button_release         (GimpTool              *tool,
                                                             const GimpCoords      *coords,
                                                             guint32                time,
                                                             GdkModifierType        state,
                                                             GimpButtonReleaseType  release_type,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_motion                 (GimpTool              *tool,
                                                             const GimpCoords      *coords,
                                                             guint32                time,
                                                             GdkModifierType        state,
                                                             GimpDisplay           *display);
static gboolean  gimp_transform_tool_key_press              (GimpTool              *tool,
                                                             GdkEventKey           *kevent,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_modifier_key           (GimpTool              *tool,
                                                             GdkModifierType        key,
                                                             gboolean               press,
                                                             GdkModifierType        state,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_oper_update            (GimpTool              *tool,
                                                             const GimpCoords      *coords,
                                                             GdkModifierType        state,
                                                             gboolean               proximity,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_cursor_update          (GimpTool              *tool,
                                                             const GimpCoords      *coords,
                                                             GdkModifierType        state,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_options_notify         (GimpTool              *tool,
                                                             GimpToolOptions       *options,
                                                             const GParamSpec      *pspec);

static void      gimp_transform_tool_draw                   (GimpDrawTool          *draw_tool);

static void      gimp_transform_tool_dialog_update          (GimpTransformTool     *tr_tool);

static TileManager *
                 gimp_transform_tool_real_transform         (GimpTransformTool     *tr_tool,
                                                             GimpItem              *item,
                                                             TileManager           *orig_tiles,
                                                             gint                   orig_offset_x,
                                                             gint                   orig_offset_y,
                                                             gint                  *new_offset_x,
                                                             gint                  *new_offset_y);

static void      gimp_transform_tool_set_function           (GimpTransformTool *tr_tool,
                                                             TransformAction    function);
static void      gimp_transform_tool_bounds                 (GimpTransformTool     *tr_tool,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_dialog                 (GimpTransformTool     *tr_tool);
static void      gimp_transform_tool_prepare                (GimpTransformTool     *tr_tool,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_transform              (GimpTransformTool     *tr_tool,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_transform_bounding_box (GimpTransformTool     *tr_tool);
static void      gimp_transform_tool_grid_recalc            (GimpTransformTool     *tr_tool);

static void      gimp_transform_tool_handles_recalc         (GimpTransformTool     *tr_tool,
                                                             GimpDisplay           *display);
static void      gimp_transform_tool_response               (GtkWidget             *widget,
                                                             gint                   response_id,
                                                             GimpTransformTool     *tr_tool);


G_DEFINE_TYPE (GimpTransformTool, gimp_transform_tool, GIMP_TYPE_DRAW_TOOL)

#define parent_class gimp_transform_tool_parent_class


static void
gimp_transform_tool_class_init (GimpTransformToolClass *klass)
{
  GObjectClass      *object_class = G_OBJECT_CLASS (klass);
  GimpToolClass     *tool_class   = GIMP_TOOL_CLASS (klass);
  GimpDrawToolClass *draw_class   = GIMP_DRAW_TOOL_CLASS (klass);

  object_class->finalize          = gimp_transform_tool_finalize;

  tool_class->initialize          = gimp_transform_tool_initialize;
  tool_class->control             = gimp_transform_tool_control;
  tool_class->button_press        = gimp_transform_tool_button_press;
  tool_class->button_release      = gimp_transform_tool_button_release;
  tool_class->motion              = gimp_transform_tool_motion;
  tool_class->key_press           = gimp_transform_tool_key_press;
  tool_class->modifier_key        = gimp_transform_tool_modifier_key;
  tool_class->active_modifier_key = gimp_transform_tool_modifier_key;
  tool_class->oper_update         = gimp_transform_tool_oper_update;
  tool_class->cursor_update       = gimp_transform_tool_cursor_update;
  tool_class->options_notify      = gimp_transform_tool_options_notify;

  draw_class->draw                = gimp_transform_tool_draw;

  klass->dialog                   = NULL;
  klass->dialog_update            = NULL;
  klass->prepare                  = NULL;
  klass->motion                   = NULL;
  klass->recalc_matrix            = NULL;
  klass->transform                = gimp_transform_tool_real_transform;
}

static void
gimp_transform_tool_init (GimpTransformTool *tr_tool)
{
  GimpTool *tool = GIMP_TOOL (tr_tool);

  gimp_tool_control_set_action_value_1 (tool->control,
                                        "tools/tools-transform-preview-opacity-set");

  gimp_tool_control_set_scroll_lock (tool->control, TRUE);
  gimp_tool_control_set_preserve    (tool->control, FALSE);
  gimp_tool_control_set_dirty_mask  (tool->control,
                                     GIMP_DIRTY_IMAGE_SIZE |
                                     GIMP_DIRTY_DRAWABLE   |
                                     GIMP_DIRTY_SELECTION);
  gimp_tool_control_set_precision   (tool->control,
                                     GIMP_CURSOR_PRECISION_SUBPIXEL);

  tr_tool->function      = TRANSFORM_CREATING;

  gimp_matrix3_identity (&tr_tool->transform);

  tr_tool->handle_w      = GIMP_TOOL_HANDLE_SIZE_LARGE;
  tr_tool->handle_h      = GIMP_TOOL_HANDLE_SIZE_LARGE;

  tr_tool->progress_text = _("Transforming");
}

static void
gimp_transform_tool_finalize (GObject *object)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (object);

  if (tr_tool->dialog)
    {
      gtk_widget_destroy (tr_tool->dialog);
      tr_tool->dialog = NULL;
    }

  if (tr_tool->grid_coords)
    {
      g_free (tr_tool->grid_coords);
      tr_tool->grid_coords = NULL;
    }

  G_OBJECT_CLASS (parent_class)->finalize (object);
}

static gboolean
gimp_transform_tool_initialize (GimpTool     *tool,
                                GimpDisplay  *display,
                                GError      **error)
{
  GimpTransformTool *tr_tool  = GIMP_TRANSFORM_TOOL (tool);
  GimpImage         *image    = gimp_display_get_image (display);
  GimpDrawable      *drawable = gimp_image_get_active_drawable (image);

  if (! GIMP_TOOL_CLASS (parent_class)->initialize (tool, display, error))
    {
      return FALSE;
    }

  if (gimp_item_is_content_locked (GIMP_ITEM (drawable)))
    {
      g_set_error_literal (error, GIMP_ERROR, GIMP_FAILED,
			   _("The active layer's pixels are locked."));
      return FALSE;
    }

  if (display != tool->display)
    {
      gint i;

      /*  Set the pointer to the active display  */
      tool->display  = display;
      tool->drawable = drawable;

      /*  Initialize the transform tool dialog */
      if (! tr_tool->dialog)
        gimp_transform_tool_dialog (tr_tool);

      /*  Find the transform bounds for some tools (like scale,
       *  perspective) that actually need the bounds for initializing
       */
      gimp_transform_tool_bounds (tr_tool, display);

      /*  Inizialize the tool-specific trans_info, and adjust the
       *  tool dialog
       */
      gimp_transform_tool_prepare (tr_tool, display);

      /*  Recalculate the transform tool  */
      gimp_transform_tool_recalc_matrix (tr_tool);

      /*  start drawing the bounding box and handles...  */
      gimp_draw_tool_start (GIMP_DRAW_TOOL (tool), display);

      tr_tool->function = TRANSFORM_CREATING;

      /*  Save the current transformation info  */
      for (i = 0; i < TRANS_INFO_SIZE; i++)
        {
          tr_tool->old_trans_info[i]  = tr_tool->trans_info[i];
          tr_tool->prev_trans_info[i] = tr_tool->trans_info[i];
        }
    }

  return TRUE;
}

static void
gimp_transform_tool_control (GimpTool       *tool,
                             GimpToolAction  action,
                             GimpDisplay    *display)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (tool);

  switch (action)
    {
    case GIMP_TOOL_ACTION_PAUSE:
      break;

    case GIMP_TOOL_ACTION_RESUME:
      gimp_transform_tool_bounds (tr_tool, display);
      gimp_transform_tool_recalc_matrix (tr_tool);
      break;

    case GIMP_TOOL_ACTION_HALT:
      tr_tool->function = TRANSFORM_CREATING;

      if (tr_tool->dialog)
        gimp_dialog_factory_hide_dialog (tr_tool->dialog);

      tool->drawable = NULL;
      break;
    }

  GIMP_TOOL_CLASS (parent_class)->control (tool, action, display);
}

static void
gimp_transform_tool_button_press (GimpTool            *tool,
                                  const GimpCoords    *coords,
                                  guint32              time,
                                  GdkModifierType      state,
                                  GimpButtonPressType  press_type,
                                  GimpDisplay         *display)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (tool);

  if (tr_tool->function == TRANSFORM_CREATING)
    gimp_transform_tool_oper_update (tool, coords, state, TRUE, display);

  tr_tool->lastx = coords->x;
  tr_tool->lasty = coords->y;

  gimp_tool_control_activate (tool->control);
}

static void
gimp_transform_tool_button_release (GimpTool              *tool,
                                    const GimpCoords      *coords,
                                    guint32                time,
                                    GdkModifierType        state,
                                    GimpButtonReleaseType  release_type,
                                    GimpDisplay           *display)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (tool);
  gint               i;

  gimp_tool_control_halt (tool->control);

  /*  if we are creating, there is nothing to be done...exit  */
  if (tr_tool->function == TRANSFORM_CREATING && tr_tool->use_grid)
    return;

  if (release_type != GIMP_BUTTON_RELEASE_CANCEL)
    {
      /* Shift-clicking is another way to approve the transform  */
      if ((state & GDK_SHIFT_MASK) || ! tr_tool->use_grid)
        {
          gimp_transform_tool_response (NULL, GTK_RESPONSE_OK, tr_tool);
        }

      /*  Store current trans_info  */
      for (i = 0; i < TRANS_INFO_SIZE; i++)
        tr_tool->prev_trans_info[i] = tr_tool->trans_info[i];
    }
  else
    {
      gimp_draw_tool_pause (GIMP_DRAW_TOOL (tool));

      /*  Restore the previous trans_info  */
      for (i = 0; i < TRANS_INFO_SIZE; i++)
        tr_tool->trans_info[i] = tr_tool->prev_trans_info[i];

      /*  reget the selection bounds  */
      gimp_transform_tool_bounds (tr_tool, display);

      /*  recalculate the tool's transformation matrix  */
      gimp_transform_tool_recalc_matrix (tr_tool);

      gimp_draw_tool_resume (GIMP_DRAW_TOOL (tool));
    }
}

static void
gimp_transform_tool_motion (GimpTool         *tool,
                            const GimpCoords *coords,
                            guint32           time,
                            GdkModifierType   state,
                            GimpDisplay      *display)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (tool);

  /*  if we are creating, there is nothing to be done so exit.  */
  if (tr_tool->function == TRANSFORM_CREATING || ! tr_tool->use_grid)
    return;

  gimp_draw_tool_pause (GIMP_DRAW_TOOL (tool));

  tr_tool->curx = coords->x;
  tr_tool->cury = coords->y;

  /*  recalculate the tool's transformation matrix  */
  if (GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->motion)
    {
      GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->motion (tr_tool);

      gimp_transform_tool_recalc_matrix (tr_tool);
    }

  tr_tool->lastx = tr_tool->curx;
  tr_tool->lasty = tr_tool->cury;

  gimp_draw_tool_resume (GIMP_DRAW_TOOL (tool));
}

static gboolean
gimp_transform_tool_key_press (GimpTool    *tool,
                               GdkEventKey *kevent,
                               GimpDisplay *display)
{
  GimpTransformTool *trans_tool = GIMP_TRANSFORM_TOOL (tool);
  GimpDrawTool      *draw_tool  = GIMP_DRAW_TOOL (tool);

  if (display == draw_tool->display)
    {
      switch (kevent->keyval)
        {
        case GDK_KEY_Return:
        case GDK_KEY_KP_Enter:
        case GDK_KEY_ISO_Enter:
          gimp_transform_tool_response (NULL, GTK_RESPONSE_OK, trans_tool);
          return TRUE;

        case GDK_KEY_BackSpace:
          gimp_transform_tool_response (NULL, RESPONSE_RESET, trans_tool);
          return TRUE;

        case GDK_KEY_Escape:
          gimp_transform_tool_response (NULL, GTK_RESPONSE_CANCEL, trans_tool);
          return TRUE;
        }
    }

  return FALSE;
}

static void
gimp_transform_tool_modifier_key (GimpTool        *tool,
                                  GdkModifierType  key,
                                  gboolean         press,
                                  GdkModifierType  state,
                                  GimpDisplay     *display)
{
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tool);

  if (key == GDK_CONTROL_MASK)
    g_object_set (options,
                  "constrain", ! options->constrain,
                  NULL);
}

static void
gimp_transform_tool_oper_update (GimpTool         *tool,
                                 const GimpCoords *coords,
                                 GdkModifierType   state,
                                 gboolean          proximity,
                                 GimpDisplay      *display)
{
  GimpTransformTool *tr_tool   = GIMP_TRANSFORM_TOOL (tool);
  GimpDrawTool      *draw_tool = GIMP_DRAW_TOOL (tool);
  TransformAction    function  = TRANSFORM_HANDLE_NONE;

  if (display != tool->display)
    {
      gimp_transform_tool_set_function (tr_tool, function);
      return;
    }

  if (tr_tool->use_handles)
    {
      gdouble closest_dist;
      gdouble dist;

      dist = gimp_draw_tool_calc_distance_square (draw_tool, display,
                                                  coords->x, coords->y,
                                                  tr_tool->tx1, tr_tool->ty1);
      closest_dist = dist;
      function = TRANSFORM_HANDLE_NW;

      dist = gimp_draw_tool_calc_distance_square (draw_tool, display,
                                                  coords->x, coords->y,
                                                  tr_tool->tx2, tr_tool->ty2);
      if (dist < closest_dist)
        {
          closest_dist = dist;
          function = TRANSFORM_HANDLE_NE;
        }

      dist = gimp_draw_tool_calc_distance_square (draw_tool, display,
                                                  coords->x, coords->y,
                                                  tr_tool->tx3, tr_tool->ty3);
      if (dist < closest_dist)
        {
          closest_dist = dist;
          function = TRANSFORM_HANDLE_SW;
        }

      dist = gimp_draw_tool_calc_distance_square (draw_tool, display,
                                                  coords->x, coords->y,
                                                  tr_tool->tx4, tr_tool->ty4);
      if (dist < closest_dist)
        {
          closest_dist = dist;
          function = TRANSFORM_HANDLE_SE;
        }

      if (tr_tool->use_mid_handles)
        {
          if (gimp_canvas_item_hit (tr_tool->handles[TRANSFORM_HANDLE_N],
                                    coords->x, coords->y))
            {
              function = TRANSFORM_HANDLE_N;
            }
          else if (gimp_canvas_item_hit (tr_tool->handles[TRANSFORM_HANDLE_E],
                                         coords->x, coords->y))
            {
              function = TRANSFORM_HANDLE_E;
            }
          else if (gimp_canvas_item_hit (tr_tool->handles[TRANSFORM_HANDLE_S],
                                         coords->x, coords->y))
            {
              function = TRANSFORM_HANDLE_S;
            }
          else if (gimp_canvas_item_hit (tr_tool->handles[TRANSFORM_HANDLE_W],
                                         coords->x, coords->y))
            {
              function = TRANSFORM_HANDLE_W;
            }
        }
    }

  if (tr_tool->use_center &&
      gimp_canvas_item_hit (tr_tool->handles[TRANSFORM_HANDLE_CENTER],
                            coords->x, coords->y))
    {
      function = TRANSFORM_HANDLE_CENTER;
    }

  gimp_transform_tool_set_function (tr_tool, function);
}

static void
gimp_transform_tool_cursor_update (GimpTool         *tool,
                                   const GimpCoords *coords,
                                   GdkModifierType   state,
                                   GimpDisplay      *display)
{
  GimpTransformTool    *tr_tool = GIMP_TRANSFORM_TOOL (tool);
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tool);
  GimpCursorType        cursor;
  GimpCursorModifier    modifier = GIMP_CURSOR_MODIFIER_NONE;
  GimpImage            *image    = gimp_display_get_image (display);

  cursor = gimp_tool_control_get_cursor (tool->control);

  if (tr_tool->use_handles)
    {
      switch (tr_tool->function)
        {
        case TRANSFORM_HANDLE_NW:
          cursor = GIMP_CURSOR_CORNER_TOP_LEFT;
          break;

        case TRANSFORM_HANDLE_NE:
          cursor = GIMP_CURSOR_CORNER_TOP_RIGHT;
          break;

        case TRANSFORM_HANDLE_SW:
          cursor = GIMP_CURSOR_CORNER_BOTTOM_LEFT;
          break;

        case TRANSFORM_HANDLE_SE:
          cursor = GIMP_CURSOR_CORNER_BOTTOM_RIGHT;
          break;

        case TRANSFORM_HANDLE_N:
          cursor = GIMP_CURSOR_SIDE_TOP;
          break;

        case TRANSFORM_HANDLE_E:
          cursor = GIMP_CURSOR_SIDE_RIGHT;
          break;

        case TRANSFORM_HANDLE_S:
          cursor = GIMP_CURSOR_SIDE_BOTTOM;
          break;

        case TRANSFORM_HANDLE_W:
          cursor = GIMP_CURSOR_SIDE_LEFT;
          break;

        default:
          cursor = GIMP_CURSOR_CROSSHAIR_SMALL;
          break;
        }
    }

  if (tr_tool->use_center && tr_tool->function == TRANSFORM_HANDLE_CENTER)
    {
      modifier = GIMP_CURSOR_MODIFIER_MOVE;
    }

  switch (options->type)
    {
      GimpDrawable *drawable;

    case GIMP_TRANSFORM_TYPE_LAYER:
      drawable = gimp_image_get_active_drawable (image);
      if (gimp_item_is_content_locked (GIMP_ITEM (drawable)))
        modifier = GIMP_CURSOR_MODIFIER_BAD;
      break;

    case GIMP_TRANSFORM_TYPE_SELECTION:
      break;

    case GIMP_TRANSFORM_TYPE_PATH:
      if (! gimp_image_get_active_vectors (image))
        modifier = GIMP_CURSOR_MODIFIER_BAD;
      break;
    }

  gimp_tool_control_set_cursor          (tool->control, cursor);
  gimp_tool_control_set_cursor_modifier (tool->control, modifier);

  GIMP_TOOL_CLASS (parent_class)->cursor_update (tool, coords, state, display);
}

static void
gimp_transform_tool_options_notify (GimpTool         *tool,
                                    GimpToolOptions  *options,
                                    const GParamSpec *pspec)
{
  GimpTransformTool *tr_tool = GIMP_TRANSFORM_TOOL (tool);

  GIMP_TOOL_CLASS (parent_class)->options_notify (tool, options, pspec);

  if (tr_tool->use_grid)
    {
      gimp_draw_tool_pause (GIMP_DRAW_TOOL (tr_tool));

      if (! strcmp (pspec->name, "type") ||
          ! strcmp (pspec->name, "direction"))
        {
          if (tr_tool->function != TRANSFORM_CREATING)
            {
              if (tool->display)
                {
                  /*  reget the selection bounds  */
                  gimp_transform_tool_bounds (tr_tool, tool->display);

                  /*  recalculate the tool's transformation matrix  */
                  gimp_transform_tool_recalc_matrix (tr_tool);
                }
            }
        }

      if (tr_tool->function != TRANSFORM_CREATING)
        {
          gimp_transform_tool_grid_recalc (tr_tool);
          gimp_transform_tool_transform_bounding_box (tr_tool);
        }

      gimp_draw_tool_resume (GIMP_DRAW_TOOL (tr_tool));
    }

  if (! strcmp (pspec->name, "constrain"))
    {
      gimp_transform_tool_dialog_update (tr_tool);
    }
}

static void
gimp_transform_tool_draw (GimpDrawTool *draw_tool)
{
  GimpTool             *tool    = GIMP_TOOL (draw_tool);
  GimpTransformTool    *tr_tool = GIMP_TRANSFORM_TOOL (draw_tool);
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tool);
  GimpImage            *image   = gimp_display_get_image (tool->display);
  gint                  i;

  for (i = 0; i < G_N_ELEMENTS (tr_tool->handles); i++)
    tr_tool->handles[i] = NULL;

  if (tr_tool->use_grid)
    {
      GimpCanvasGroup *stroke_group;

      if (gimp_transform_options_show_preview (options))
        {
          gimp_draw_tool_add_transform_preview (draw_tool,
                                                tool->drawable,
                                                &tr_tool->transform,
                                                tr_tool->x1,
                                                tr_tool->y1,
                                                tr_tool->x2,
                                                tr_tool->y2,
                                                GIMP_IS_PERSPECTIVE_TOOL (tr_tool),
                                                options->preview_opacity);
        }

      stroke_group = gimp_draw_tool_add_stroke_group (draw_tool);

      gimp_draw_tool_push_group (draw_tool, stroke_group);

      /*  draw the grid  */
      if (tr_tool->grid_coords &&
          gimp_transform_polygon_is_convex (tr_tool->tx1, tr_tool->ty1,
                                            tr_tool->tx2, tr_tool->ty2,
                                            tr_tool->tx3, tr_tool->ty3,
                                            tr_tool->tx4, tr_tool->ty4))
        {
          gint k = tr_tool->ngx + tr_tool->ngy;
          gint i, gci;

          for (i = 0, gci = 0; i < k; i++, gci += 4)
            {
              gdouble x1, y1, x2, y2;

              gimp_matrix3_transform_point (&tr_tool->transform,
                                            tr_tool->grid_coords[gci],
                                            tr_tool->grid_coords[gci + 1],
                                            &x1, &y1);
              gimp_matrix3_transform_point (&tr_tool->transform,
                                            tr_tool->grid_coords[gci + 2],
                                            tr_tool->grid_coords[gci + 3],
                                            &x2, &y2);

              gimp_draw_tool_add_line (draw_tool, x1, y1, x2, y2);
            }
        }

      /*  draw the bounding box  */
      gimp_draw_tool_add_line (draw_tool,
                               tr_tool->tx1, tr_tool->ty1,
                               tr_tool->tx2, tr_tool->ty2);
      gimp_draw_tool_add_line (draw_tool,
                               tr_tool->tx2, tr_tool->ty2,
                               tr_tool->tx4, tr_tool->ty4);
      gimp_draw_tool_add_line (draw_tool,
                               tr_tool->tx3, tr_tool->ty3,
                               tr_tool->tx4, tr_tool->ty4);
      gimp_draw_tool_add_line (draw_tool,
                               tr_tool->tx3, tr_tool->ty3,
                               tr_tool->tx1, tr_tool->ty1);

      gimp_draw_tool_pop_group (draw_tool);
    }

  gimp_transform_tool_handles_recalc (tr_tool, tool->display);

  if (tr_tool->use_handles)
    {
      /*  draw the tool handles  */
      tr_tool->handles[TRANSFORM_HANDLE_NW] =
        gimp_draw_tool_add_handle (draw_tool,
                                   GIMP_HANDLE_SQUARE,
                                   tr_tool->tx1, tr_tool->ty1,
                                   tr_tool->handle_w, tr_tool->handle_h,
                                   GIMP_HANDLE_ANCHOR_CENTER);

      tr_tool->handles[TRANSFORM_HANDLE_NE] =
        gimp_draw_tool_add_handle (draw_tool,
                                   GIMP_HANDLE_SQUARE,
                                   tr_tool->tx2, tr_tool->ty2,
                                   tr_tool->handle_w, tr_tool->handle_h,
                                   GIMP_HANDLE_ANCHOR_CENTER);

      tr_tool->handles[TRANSFORM_HANDLE_SW] =
        gimp_draw_tool_add_handle (draw_tool,
                                   GIMP_HANDLE_SQUARE,
                                   tr_tool->tx3, tr_tool->ty3,
                                   tr_tool->handle_w, tr_tool->handle_h,
                                   GIMP_HANDLE_ANCHOR_CENTER);

      tr_tool->handles[TRANSFORM_HANDLE_SE] =
        gimp_draw_tool_add_handle (draw_tool,
                                   GIMP_HANDLE_SQUARE,
                                   tr_tool->tx4, tr_tool->ty4,
                                   tr_tool->handle_w, tr_tool->handle_h,
                                   GIMP_HANDLE_ANCHOR_CENTER);

      if (tr_tool->use_mid_handles)
        {
          gdouble x, y;

          x = (tr_tool->tx1 + tr_tool->tx2) / 2.0;
          y = (tr_tool->ty1 + tr_tool->ty2) / 2.0;

          tr_tool->handles[TRANSFORM_HANDLE_N] =
            gimp_draw_tool_add_handle (draw_tool,
                                       GIMP_HANDLE_SQUARE,
                                       x, y,
                                       tr_tool->handle_w, tr_tool->handle_h,
                                       GIMP_HANDLE_ANCHOR_CENTER);

          x = (tr_tool->tx2 + tr_tool->tx4) / 2.0;
          y = (tr_tool->ty2 + tr_tool->ty4) / 2.0;

          tr_tool->handles[TRANSFORM_HANDLE_E] =
            gimp_draw_tool_add_handle (draw_tool,
                                       GIMP_HANDLE_SQUARE,
                                       x, y,
                                       tr_tool->handle_w, tr_tool->handle_h,
                                       GIMP_HANDLE_ANCHOR_CENTER);

          x = (tr_tool->tx3 + tr_tool->tx4) / 2.0;
          y = (tr_tool->ty3 + tr_tool->ty4) / 2.0;

          tr_tool->handles[TRANSFORM_HANDLE_S] =
            gimp_draw_tool_add_handle (draw_tool,
                                       GIMP_HANDLE_SQUARE,
                                       x, y,
                                       tr_tool->handle_w, tr_tool->handle_h,
                                       GIMP_HANDLE_ANCHOR_CENTER);

          x = (tr_tool->tx3 + tr_tool->tx1) / 2.0;
          y = (tr_tool->ty3 + tr_tool->ty1) / 2.0;

          tr_tool->handles[TRANSFORM_HANDLE_W] =
            gimp_draw_tool_add_handle (draw_tool,
                                       GIMP_HANDLE_SQUARE,
                                       x, y,
                                       tr_tool->handle_w, tr_tool->handle_h,
                                       GIMP_HANDLE_ANCHOR_CENTER);
        }
    }

  /*  draw the center  */
  if (tr_tool->use_center)
    {
      GimpCanvasGroup *stroke_group;
      gint             d = MIN (tr_tool->handle_w, tr_tool->handle_h);

      stroke_group = gimp_draw_tool_add_stroke_group (draw_tool);

      tr_tool->handles[TRANSFORM_HANDLE_CENTER] = GIMP_CANVAS_ITEM (stroke_group);

      gimp_draw_tool_push_group (draw_tool, stroke_group);

      gimp_draw_tool_add_handle (draw_tool,
                                 GIMP_HANDLE_CIRCLE,
                                 tr_tool->tcx, tr_tool->tcy,
                                 d, d,
                                 GIMP_HANDLE_ANCHOR_CENTER);
      gimp_draw_tool_add_handle (draw_tool,
                                 GIMP_HANDLE_CROSS,
                                 tr_tool->tcx, tr_tool->tcy,
                                 d, d,
                                 GIMP_HANDLE_ANCHOR_CENTER);

      gimp_draw_tool_pop_group (draw_tool);
    }

  if (tr_tool->handles[tr_tool->function])
    {
      gimp_canvas_item_set_highlight (tr_tool->handles[tr_tool->function],
                                      TRUE);
    }

  if (options->type == GIMP_TRANSFORM_TYPE_SELECTION)
    {
      GimpMatrix3     matrix = tr_tool->transform;
      const BoundSeg *orig_in;
      const BoundSeg *orig_out;
      BoundSeg       *segs_in;
      BoundSeg       *segs_out;
      gint            num_segs_in;
      gint            num_segs_out;

      gimp_channel_boundary (gimp_image_get_mask (image),
                             &orig_in, &orig_out,
                             &num_segs_in, &num_segs_out,
                             0, 0, 0, 0);

      segs_in  = g_memdup (orig_in,  num_segs_in  * sizeof (BoundSeg));
      segs_out = g_memdup (orig_out, num_segs_out * sizeof (BoundSeg));

      if (segs_in)
        {
          for (i = 0; i < num_segs_in; i++)
            {
              gdouble tx, ty;

              gimp_matrix3_transform_point (&matrix,
                                            segs_in[i].x1, segs_in[i].y1,
                                            &tx, &ty);
              segs_in[i].x1 = RINT (tx);
              segs_in[i].y1 = RINT (ty);

              gimp_matrix3_transform_point (&matrix,
                                            segs_in[i].x2, segs_in[i].y2,
                                            &tx, &ty);
              segs_in[i].x2 = RINT (tx);
              segs_in[i].y2 = RINT (ty);
            }

          gimp_draw_tool_add_boundary (draw_tool,
                                       segs_in, num_segs_in,
                                       NULL,
                                       0, 0);
          g_free (segs_in);
        }

      if (segs_out)
        {
          for (i = 0; i < num_segs_out; i++)
            {
              gdouble tx, ty;

              gimp_matrix3_transform_point (&matrix,
                                            segs_out[i].x1, segs_out[i].y1,
                                            &tx, &ty);
              segs_out[i].x1 = RINT (tx);
              segs_out[i].y1 = RINT (ty);

              gimp_matrix3_transform_point (&matrix,
                                            segs_out[i].x2, segs_out[i].y2,
                                            &tx, &ty);
              segs_out[i].x2 = RINT (tx);
              segs_out[i].y2 = RINT (ty);
            }

          gimp_draw_tool_add_boundary (draw_tool,
                                       segs_out, num_segs_out,
                                       NULL,
                                       0, 0);
          g_free (segs_out);
        }
    }
  else if (options->type == GIMP_TRANSFORM_TYPE_PATH)
    {
      GimpVectors *vectors;
      GimpStroke  *stroke = NULL;
      GimpMatrix3  matrix = tr_tool->transform;

      vectors = gimp_image_get_active_vectors (image);

      if (vectors)
        {
          if (options->direction == GIMP_TRANSFORM_BACKWARD)
            gimp_matrix3_invert (&matrix);

          while ((stroke = gimp_vectors_stroke_get_next (vectors, stroke)))
            {
              GArray   *coords;
              gboolean  closed;

              coords = gimp_stroke_interpolate (stroke, 1.0, &closed);

              if (coords && coords->len)
                {
                  gint i;

                  for (i = 0; i < coords->len; i++)
                    {
                      GimpCoords *curr = &g_array_index (coords, GimpCoords, i);

                      gimp_matrix3_transform_point (&matrix,
                                                    curr->x, curr->y,
                                                    &curr->x, &curr->y);
                    }

                  gimp_draw_tool_add_strokes (draw_tool,
                                              &g_array_index (coords,
                                                              GimpCoords, 0),
                                              coords->len, FALSE);
                }

              if (coords)
                g_array_free (coords, TRUE);
            }
        }
    }
}

static void
gimp_transform_tool_dialog_update (GimpTransformTool *tr_tool)
{
  if (tr_tool->dialog &&
      GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->dialog_update)
    {
      GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->dialog_update (tr_tool);
    }
}

static TileManager *
gimp_transform_tool_real_transform (GimpTransformTool *tr_tool,
                                    GimpItem          *active_item,
                                    TileManager       *orig_tiles,
                                    gint               orig_offset_x,
                                    gint               orig_offset_y,
                                    gint              *new_offset_x,
                                    gint              *new_offset_y)
{
  GimpTool             *tool    = GIMP_TOOL (tr_tool);
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tool);
  GimpContext          *context = GIMP_CONTEXT (options);
  GimpProgress         *progress;
  TileManager          *ret  = NULL;

  progress = gimp_progress_start (GIMP_PROGRESS (tool),
                                  tr_tool->progress_text, FALSE);

  if (gimp_item_get_linked (active_item))
    gimp_item_linked_transform (active_item, context,
                                &tr_tool->transform,
                                options->direction,
                                options->interpolation,
                                options->recursion_level,
                                options->clip,
                                progress);

  if (orig_tiles)
    {
      /*  this happens when transforming a selection cut out of a
       *  normal drawable, or the selection
       */

      GimpTransformResize clip_result = options->clip;

      /*  always clip the selction and unfloated channels
       *  so they keep their size
       */
      if (GIMP_IS_CHANNEL (active_item) &&
          tile_manager_bpp (orig_tiles) == 1)
        clip_result = GIMP_TRANSFORM_RESIZE_CLIP;

      ret = gimp_drawable_transform_tiles_affine (GIMP_DRAWABLE (active_item),
                                                  context,
                                                  orig_tiles,
                                                  orig_offset_x,
                                                  orig_offset_y,
                                                  &tr_tool->transform,
                                                  options->direction,
                                                  options->interpolation,
                                                  options->recursion_level,
                                                  clip_result,
                                                  new_offset_x,
                                                  new_offset_y,
                                                  progress);
    }
  else
    {
      /*  this happens for entire drawables, paths and layer groups  */

      gimp_item_transform (active_item, context,
                           &tr_tool->transform,
                           options->direction,
                           options->interpolation,
                           options->recursion_level,
                           options->clip,
                           progress);
    }

  if (progress)
    gimp_progress_end (progress);

  return ret;
}

static void
gimp_transform_tool_transform (GimpTransformTool *tr_tool,
                               GimpDisplay       *display)
{
  GimpTool             *tool           = GIMP_TOOL (tr_tool);
  GimpTransformOptions *options        = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tool);
  GimpContext          *context        = GIMP_CONTEXT (options);
  GimpImage            *image          = gimp_display_get_image (display);
  GimpItem             *active_item    = NULL;
  TileManager          *orig_tiles     = NULL;
  gint                  orig_offset_x;
  gint                  orig_offset_y;
  TileManager          *new_tiles;
  gint                  new_offset_x;
  gint                  new_offset_y;
  const gchar          *null_message   = NULL;
  const gchar          *locked_message = NULL;
  gboolean              new_layer;

  switch (options->type)
    {
    case GIMP_TRANSFORM_TYPE_LAYER:
      active_item = GIMP_ITEM (gimp_image_get_active_drawable (image));
      null_message   = _("There is no layer to transform.");
      locked_message = _("The active layer's pixels are locked.");
      break;

    case GIMP_TRANSFORM_TYPE_SELECTION:
      active_item = GIMP_ITEM (gimp_image_get_mask (image));
      /* cannot happen, so don't translate these messages */
      null_message   = "There is no selection to transform.";
      locked_message = "The selection's pixels are locked.";
      break;

    case GIMP_TRANSFORM_TYPE_PATH:
      active_item = GIMP_ITEM (gimp_image_get_active_vectors (image));
      null_message   = _("There is no path to transform.");
      locked_message = _("The active path's strokes are locked.");
      break;
    }

  if (! active_item)
    {
      gimp_tool_message_literal (tool, display, null_message);
      gimp_tool_control (tool, GIMP_TOOL_ACTION_HALT, display);
      return;
    }

  if (gimp_item_is_content_locked (active_item))
    {
      gimp_tool_message_literal (tool, display, locked_message);
      gimp_tool_control (tool, GIMP_TOOL_ACTION_HALT, display);
      return;
    }

  if (tr_tool->dialog)
    gimp_dialog_factory_hide_dialog (tr_tool->dialog);

  gimp_set_busy (display->gimp);

  /* undraw the tool before we muck around with the transform matrix */
  gimp_draw_tool_stop (GIMP_DRAW_TOOL (tr_tool));

  /*  We're going to dirty this image, but we want to keep the tool around  */
  gimp_tool_control_set_preserve (tool->control, TRUE);

  gimp_image_undo_group_start (image, GIMP_UNDO_GROUP_TRANSFORM,
                               tr_tool->undo_desc);

  switch (options->type)
    {
    case GIMP_TRANSFORM_TYPE_LAYER:
      if (! gimp_viewable_get_children (GIMP_VIEWABLE (tool->drawable)) &&
          ! gimp_channel_is_empty (gimp_image_get_mask (image)))
        {
          orig_tiles = gimp_drawable_transform_cut (tool->drawable,
                                                    context,
                                                    &orig_offset_x,
                                                    &orig_offset_y,
                                                    &new_layer);
        }
      break;

    case GIMP_TRANSFORM_TYPE_SELECTION:
      orig_tiles = tile_manager_ref (gimp_drawable_get_tiles (GIMP_DRAWABLE (active_item)));
      orig_offset_x = 0;
      orig_offset_y = 0;
      break;

    case GIMP_TRANSFORM_TYPE_PATH:
      break;
    }

  /*  Send the request for the transformation to the tool...
   */
  new_tiles = GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->transform (tr_tool,
                                                                  active_item,
                                                                  orig_tiles,
                                                                  orig_offset_x,
                                                                  orig_offset_y,
                                                                  &new_offset_x,
                                                                  &new_offset_y);

  if (orig_tiles)
    tile_manager_unref (orig_tiles);

  switch (options->type)
    {
    case GIMP_TRANSFORM_TYPE_LAYER:
      if (new_tiles)
        {
          /*  paste the new transformed image to the image...also implement
           *  undo...
           */
          gimp_drawable_transform_paste (tool->drawable,
                                         new_tiles,
                                         new_offset_x, new_offset_y,
                                         new_layer);
          tile_manager_unref (new_tiles);
        }
      break;

     case GIMP_TRANSFORM_TYPE_SELECTION:
      if (new_tiles)
        {
          gimp_channel_push_undo (GIMP_CHANNEL (active_item), NULL);

          gimp_drawable_set_tiles (GIMP_DRAWABLE (active_item),
                                   FALSE, NULL, new_tiles,
                                   gimp_drawable_type (GIMP_DRAWABLE (active_item)));
          tile_manager_unref (new_tiles);
        }
      break;

    case GIMP_TRANSFORM_TYPE_PATH:
      /*  Nothing to be done  */
      break;
    }

  gimp_image_undo_push (image, GIMP_TYPE_TRANSFORM_TOOL_UNDO,
                        GIMP_UNDO_TRANSFORM, NULL,
                        0,
                        "transform-tool", tr_tool,
                        NULL);

  gimp_image_undo_group_end (image);

  /*  We're done dirtying the image, and would like to be restarted if
   *  the image gets dirty while the tool exists
   */
  gimp_tool_control_set_preserve (tool->control, FALSE);

  gimp_tool_control (tool, GIMP_TOOL_ACTION_HALT, display);

  gimp_unset_busy (display->gimp);

  gimp_image_flush (image);
}

static void
gimp_transform_tool_set_function (GimpTransformTool *tr_tool,
                                  TransformAction    function)
{
  if (function != tr_tool->function)
    {
      if (tr_tool->handles[tr_tool->function] &&
          gimp_draw_tool_is_active (GIMP_DRAW_TOOL (tr_tool)))
        {
          gimp_canvas_item_set_highlight (tr_tool->handles[tr_tool->function],
                                          FALSE);
        }

      tr_tool->function = function;

      if (tr_tool->handles[tr_tool->function] &&
          gimp_draw_tool_is_active (GIMP_DRAW_TOOL (tr_tool)))
        {
          gimp_canvas_item_set_highlight (tr_tool->handles[tr_tool->function],
                                          TRUE);
        }
    }
}

static void
gimp_transform_tool_transform_bounding_box (GimpTransformTool *tr_tool)
{
  g_return_if_fail (GIMP_IS_TRANSFORM_TOOL (tr_tool));

  gimp_matrix3_transform_point (&tr_tool->transform,
                                tr_tool->x1, tr_tool->y1,
                                &tr_tool->tx1, &tr_tool->ty1);
  gimp_matrix3_transform_point (&tr_tool->transform,
                                tr_tool->x2, tr_tool->y1,
                                &tr_tool->tx2, &tr_tool->ty2);
  gimp_matrix3_transform_point (&tr_tool->transform,
                                tr_tool->x1, tr_tool->y2,
                                &tr_tool->tx3, &tr_tool->ty3);
  gimp_matrix3_transform_point (&tr_tool->transform,
                                tr_tool->x2, tr_tool->y2,
                                &tr_tool->tx4, &tr_tool->ty4);

  gimp_matrix3_transform_point (&tr_tool->transform,
                                tr_tool->cx, tr_tool->cy,
                                &tr_tool->tcx, &tr_tool->tcy);
}

static void
gimp_transform_tool_bounds (GimpTransformTool *tr_tool,
                            GimpDisplay       *display)
{
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tr_tool);
  GimpImage            *image   = gimp_display_get_image (display);

  switch (options->type)
    {
    case GIMP_TRANSFORM_TYPE_LAYER:
      {
        GimpDrawable *drawable;
        gint          offset_x;
        gint          offset_y;

        drawable = gimp_image_get_active_drawable (image);

        gimp_item_get_offset (GIMP_ITEM (drawable), &offset_x, &offset_y);

        gimp_item_mask_bounds (GIMP_ITEM (drawable),
                               &tr_tool->x1, &tr_tool->y1,
                               &tr_tool->x2, &tr_tool->y2);
        tr_tool->x1 += offset_x;
        tr_tool->y1 += offset_y;
        tr_tool->x2 += offset_x;
        tr_tool->y2 += offset_y;
      }
      break;

    case GIMP_TRANSFORM_TYPE_SELECTION:
    case GIMP_TRANSFORM_TYPE_PATH:
      gimp_channel_bounds (gimp_image_get_mask (image),
                           &tr_tool->x1, &tr_tool->y1,
                           &tr_tool->x2, &tr_tool->y2);
      break;
    }

  tr_tool->cx = (gdouble) (tr_tool->x1 + tr_tool->x2) / 2.0;
  tr_tool->cy = (gdouble) (tr_tool->y1 + tr_tool->y2) / 2.0;

  tr_tool->aspect = ((gdouble) (tr_tool->x2 - tr_tool->x1) /
                     (gdouble) (tr_tool->y2 - tr_tool->y1));

  /*  changing the bounds invalidates any grid we may have  */
  if (tr_tool->use_grid)
    gimp_transform_tool_grid_recalc (tr_tool);
}

static void
gimp_transform_tool_grid_recalc (GimpTransformTool *tr_tool)
{
  GimpTransformOptions *options = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tr_tool);

  if (tr_tool->grid_coords != NULL)
    {
      g_free (tr_tool->grid_coords);
      tr_tool->grid_coords = NULL;
    }

  if (options->preview_type != GIMP_TRANSFORM_PREVIEW_TYPE_GRID &&
      options->preview_type != GIMP_TRANSFORM_PREVIEW_TYPE_IMAGE_GRID)
    return;

  switch (options->grid_type)
    {
    case GIMP_TRANSFORM_GRID_TYPE_N_LINES:
    case GIMP_TRANSFORM_GRID_TYPE_SPACING:
      {
        gint     i, gci;
        gdouble *coords;
        gint     width, height;

        width  = MAX (1, tr_tool->x2 - tr_tool->x1);
        height = MAX (1, tr_tool->y2 - tr_tool->y1);

        if (options->grid_type == GIMP_TRANSFORM_GRID_TYPE_N_LINES)
          {
            if (width <= height)
              {
                tr_tool->ngx = options->grid_size;
                tr_tool->ngy = tr_tool->ngx * MAX (1, height / width);
              }
            else
              {
                tr_tool->ngy = options->grid_size;
                tr_tool->ngx = tr_tool->ngy * MAX (1, width / height);
              }
          }
        else /* GIMP_TRANSFORM_GRID_TYPE_SPACING */
          {
            gint grid_size = MAX (2, options->grid_size);

            tr_tool->ngx = width  / grid_size;
            tr_tool->ngy = height / grid_size;
          }

        tr_tool->grid_coords = coords =
          g_new (gdouble, (tr_tool->ngx + tr_tool->ngy) * 4);

        gci = 0;

        for (i = 1; i <= tr_tool->ngx; i++)
          {
            coords[gci]     = tr_tool->x1 + (((gdouble) i) / (tr_tool->ngx + 1) *
                                             (tr_tool->x2 - tr_tool->x1));
            coords[gci + 1] = tr_tool->y1;
            coords[gci + 2] = coords[gci];
            coords[gci + 3] = tr_tool->y2;

            gci += 4;
          }

        for (i = 1; i <= tr_tool->ngy; i++)
          {
            coords[gci]     = tr_tool->x1;
            coords[gci + 1] = tr_tool->y1 + (((gdouble) i) / (tr_tool->ngy + 1) *
                                             (tr_tool->y2 - tr_tool->y1));
            coords[gci + 2] = tr_tool->x2;
            coords[gci + 3] = coords[gci + 1];

            gci += 4;
          }
      }

    default:
      break;
    }
}

static void
gimp_transform_tool_handles_recalc (GimpTransformTool *tr_tool,
                                    GimpDisplay       *display)
{
  gint dx1, dy1;
  gint dx2, dy2;
  gint dx3, dy3;
  gint dx4, dy4;
  gint x1, y1;
  gint x2, y2;

  gimp_display_shell_transform_xy (gimp_display_get_shell (display),
                                   tr_tool->tx1, tr_tool->ty1,
                                   &dx1, &dy1);
  gimp_display_shell_transform_xy (gimp_display_get_shell (display),
                                   tr_tool->tx2, tr_tool->ty2,
                                   &dx2, &dy2);
  gimp_display_shell_transform_xy (gimp_display_get_shell (display),
                                   tr_tool->tx3, tr_tool->ty3,
                                   &dx3, &dy3);
  gimp_display_shell_transform_xy (gimp_display_get_shell (display),
                                   tr_tool->tx4, tr_tool->ty4,
                                   &dx4, &dy4);

  x1 = MIN4 (dx1, dx2, dx3, dx4);
  y1 = MIN4 (dy1, dy2, dy3, dy4);
  x2 = MAX4 (dx1, dx2, dx3, dx4);
  y2 = MAX4 (dy1, dy2, dy3, dy4);

  tr_tool->handle_w = CLAMP ((x2 - x1) / 3,
                             MIN_HANDLE_SIZE, GIMP_TOOL_HANDLE_SIZE_LARGE);
  tr_tool->handle_h = CLAMP ((y2 - y1) / 3,
                             MIN_HANDLE_SIZE, GIMP_TOOL_HANDLE_SIZE_LARGE);
}

static void
gimp_transform_tool_dialog (GimpTransformTool *tr_tool)
{
  GimpTool     *tool      = GIMP_TOOL (tr_tool);
  GimpToolInfo *tool_info = tool->tool_info;
  const gchar  *stock_id;

  if (! GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->dialog)
    return;

  stock_id = gimp_viewable_get_stock_id (GIMP_VIEWABLE (tool_info));

  tr_tool->dialog = gimp_tool_dialog_new (tool_info,
                                          gimp_display_get_shell (tool->display),
                                          tool_info->blurb,
                                          GIMP_STOCK_RESET, RESPONSE_RESET,
                                          GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                          stock_id,         GTK_RESPONSE_OK,
                                          NULL);
  gtk_dialog_set_default_response (GTK_DIALOG (tr_tool->dialog),
                                   GTK_RESPONSE_OK);
  gtk_dialog_set_alternative_button_order (GTK_DIALOG (tr_tool->dialog),
                                           RESPONSE_RESET,
                                           GTK_RESPONSE_OK,
                                           GTK_RESPONSE_CANCEL,
                                           -1);

  g_signal_connect (tr_tool->dialog, "response",
                    G_CALLBACK (gimp_transform_tool_response),
                    tr_tool);

  GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->dialog (tr_tool);
}

static void
gimp_transform_tool_prepare (GimpTransformTool *tr_tool,
                             GimpDisplay       *display)
{
  if (tr_tool->dialog)
    {
      GimpTransformOptions *options  = GIMP_TRANSFORM_TOOL_GET_OPTIONS (tr_tool);
      GimpImage            *image    = gimp_display_get_image (display);
      GimpDrawable         *drawable = gimp_image_get_active_drawable (image);

      gimp_viewable_dialog_set_viewable (GIMP_VIEWABLE_DIALOG (tr_tool->dialog),
                                         GIMP_VIEWABLE (drawable),
                                         GIMP_CONTEXT (options));
      gimp_tool_dialog_set_shell (GIMP_TOOL_DIALOG (tr_tool->dialog),
                                  gimp_display_get_shell (display));
    }

  if (GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->prepare)
    GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->prepare (tr_tool);
}

void
gimp_transform_tool_recalc_matrix (GimpTransformTool *tr_tool)
{
  g_return_if_fail (GIMP_IS_TRANSFORM_TOOL (tr_tool));

  if (GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->recalc_matrix)
    GIMP_TRANSFORM_TOOL_GET_CLASS (tr_tool)->recalc_matrix (tr_tool);

  gimp_transform_tool_transform_bounding_box (tr_tool);

  gimp_transform_tool_dialog_update (tr_tool);

  if (tr_tool->dialog)
    gtk_widget_show (tr_tool->dialog);
}

static void
gimp_transform_tool_response (GtkWidget         *widget,
                              gint               response_id,
                              GimpTransformTool *tr_tool)
{
  GimpTool *tool = GIMP_TOOL (tr_tool);

  switch (response_id)
    {
    case RESPONSE_RESET:
      {
        gint i;

        gimp_draw_tool_pause (GIMP_DRAW_TOOL (tool));

        /*  Restore the previous transformation info  */
        for (i = 0; i < TRANS_INFO_SIZE; i++)
          tr_tool->trans_info[i] = tr_tool->old_trans_info[i];

        /*  reget the selection bounds  */
        gimp_transform_tool_bounds (tr_tool, tool->display);

        /*  recalculate the tool's transformation matrix  */
        gimp_transform_tool_recalc_matrix (tr_tool);

        gimp_draw_tool_resume (GIMP_DRAW_TOOL (tool));
      }
      break;

    case GTK_RESPONSE_OK:
      g_return_if_fail (tool->display != NULL);
      gimp_transform_tool_transform (tr_tool, tool->display);
      break;

    default:
      gimp_tool_control (tool, GIMP_TOOL_ACTION_HALT, tool->display);
      break;
    }
}
