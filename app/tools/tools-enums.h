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

#ifndef __TOOLS_ENUMS_H__
#define __TOOLS_ENUMS_H__

/*
 * these enums are registered with the type system
 */

#define GIMP_TYPE_BUTTON_PRESS_TYPE (gimp_button_press_type_get_type ())

GType gimp_button_press_type_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_BUTTON_PRESS_NORMAL,
  GIMP_BUTTON_PRESS_DOUBLE,
  GIMP_BUTTON_PRESS_TRIPLE
} GimpButtonPressType;


#define GIMP_TYPE_BUTTON_RELEASE_TYPE (gimp_button_release_type_get_type ())

GType gimp_button_release_type_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_BUTTON_RELEASE_NORMAL,
  GIMP_BUTTON_RELEASE_CANCEL,
  GIMP_BUTTON_RELEASE_CLICK,
  GIMP_BUTTON_RELEASE_NO_MOTION
} GimpButtonReleaseType;


#define GIMP_TYPE_RECTANGLE_CONSTRAINT (gimp_rectangle_constraint_get_type ())

GType gimp_rectangle_constraint_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_RECTANGLE_CONSTRAIN_NONE,
  GIMP_RECTANGLE_CONSTRAIN_IMAGE,
  GIMP_RECTANGLE_CONSTRAIN_DRAWABLE
} GimpRectangleConstraint;


#define GIMP_TYPE_RECTANGLE_PRECISION (gimp_rectangle_precision_get_type ())

GType gimp_rectangle_precision_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_RECTANGLE_PRECISION_INT,
  GIMP_RECTANGLE_PRECISION_DOUBLE,
} GimpRectanglePrecision;


#define GIMP_TYPE_RECTANGLE_TOOL_FIXED_RULE (gimp_rectangle_tool_fixed_rule_get_type ())

GType gimp_rectangle_tool_fixed_rule_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_RECTANGLE_TOOL_FIXED_ASPECT, /*< desc="Aspect ratio" >*/
  GIMP_RECTANGLE_TOOL_FIXED_WIDTH,  /*< desc="Width"        >*/
  GIMP_RECTANGLE_TOOL_FIXED_HEIGHT, /*< desc="Height"       >*/
  GIMP_RECTANGLE_TOOL_FIXED_SIZE,   /*< desc="Size"         >*/
} GimpRectangleToolFixedRule;


#define GIMP_TYPE_RECT_SELECT_MODE (gimp_rect_select_mode_get_type ())

GType gimp_rect_select_mode_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_RECT_SELECT_MODE_FREE,        /*< desc="Free select"        >*/
  GIMP_RECT_SELECT_MODE_FIXED_SIZE,  /*< desc="Fixed size"         >*/
  GIMP_RECT_SELECT_MODE_FIXED_RATIO  /*< desc="Fixed aspect ratio" >*/
} GimpRectSelectMode;


#define GIMP_TYPE_TRANSFORM_TYPE (gimp_transform_type_get_type ())

GType gimp_transform_type_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_TRANSFORM_TYPE_LAYER,     /*< desc="Layer"     >*/
  GIMP_TRANSFORM_TYPE_SELECTION, /*< desc="Selection" >*/
  GIMP_TRANSFORM_TYPE_PATH       /*< desc="Path"      >*/
} GimpTransformType;


#define GIMP_TYPE_TRANSFORM_PREVIEW_TYPE (gimp_transform_preview_type_get_type ())

GType gimp_transform_preview_type_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_TRANSFORM_PREVIEW_TYPE_OUTLINE,     /*< desc="Outline"      >*/
  GIMP_TRANSFORM_PREVIEW_TYPE_GRID,        /*< desc="Grid"         >*/
  GIMP_TRANSFORM_PREVIEW_TYPE_IMAGE,       /*< desc="Image"        >*/
  GIMP_TRANSFORM_PREVIEW_TYPE_IMAGE_GRID   /*< desc="Image + Grid" >*/
} GimpTransformPreviewType;


#define GIMP_TYPE_TRANSFORM_GRID_TYPE (gimp_transform_grid_type_get_type ())

GType gimp_transform_grid_type_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_TRANSFORM_GRID_TYPE_N_LINES,  /*< desc="Number of grid lines" >*/
  GIMP_TRANSFORM_GRID_TYPE_SPACING   /*< desc="Grid line spacing"    >*/
} GimpTransformGridType;


#define GIMP_TYPE_VECTOR_MODE (gimp_vector_mode_get_type ())

GType gimp_vector_mode_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_VECTOR_MODE_DESIGN,      /*< desc="Design" >*/
  GIMP_VECTOR_MODE_EDIT,        /*< desc="Edit"   >*/
  GIMP_VECTOR_MODE_MOVE         /*< desc="Move"   >*/
} GimpVectorMode;


#define GIMP_TYPE_TOOL_ACTION (gimp_tool_action_get_type ())

GType gimp_tool_action_get_type (void) G_GNUC_CONST;

typedef enum
{
  GIMP_TOOL_ACTION_PAUSE,
  GIMP_TOOL_ACTION_RESUME,
  GIMP_TOOL_ACTION_HALT
} GimpToolAction;


/*
 * non-registered enums; register them if needed
 */

typedef enum /*< skip >*/
{
  SELECTION_SELECT,
  SELECTION_MOVE_MASK,
  SELECTION_MOVE,
  SELECTION_MOVE_COPY,
  SELECTION_ANCHOR
} SelectFunction;

/*  Modes of GimpEditSelectionTool  */
typedef enum /*< skip >*/
{
  GIMP_TRANSLATE_MODE_VECTORS,
  GIMP_TRANSLATE_MODE_CHANNEL,
  GIMP_TRANSLATE_MODE_LAYER_MASK,
  GIMP_TRANSLATE_MODE_MASK,
  GIMP_TRANSLATE_MODE_MASK_TO_LAYER,
  GIMP_TRANSLATE_MODE_MASK_COPY_TO_LAYER,
  GIMP_TRANSLATE_MODE_LAYER,
  GIMP_TRANSLATE_MODE_FLOATING_SEL
} GimpTranslateMode;

/*  Motion event report modes  */
typedef enum /*< skip >*/
{
  GIMP_MOTION_MODE_EXACT,
  GIMP_MOTION_MODE_HINT,
  GIMP_MOTION_MODE_COMPRESS
} GimpMotionMode;


#endif /* __TOOLS_ENUMS_H__ */
