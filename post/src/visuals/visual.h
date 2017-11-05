/*****************************************************************************
 *
 *  Elmer, A Finite Element Software for Multiphysical Problems
 *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
 * 
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program (in file fem/GPL-2); if not, write to the 
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 *  Boston, MA 02110-1301, USA.
 *
 *****************************************************************************/

/*******************************************************************************
 *
 * Visual class definitions.
 *
 *******************************************************************************
 *
 *                     Author:       Juha Ruokolainen
 *
 *                    Address: CSC - IT Center for Science Ltd.
 *                                Keilaranta 14, P.O. BOX 405
 *                                  02101 Espoo, Finland
 *                                  Tel. +358 0 457 2723
 *                                Telefax: +358 0 457 2302
 *                              EMail: Juha.Ruokolainen@csc.fi
 *
 *                       Date: 27 Sep 1995
 *
 *
 * Modification history:
 *
 * 28 Sep 1995, modified visual_t and visual_type_t structures to make a list
 *              rather than an array of visual types.
 *
 ******************************************************************************/

#ifdef EXT
#undef EXT
#endif

#ifdef MODULE_VISUALS
#define EXT 
#else
#define EXT extern
#endif

#ifdef MODULE_VISUALS
    double TooLong1 = 0.30,TooLong2 = 0.5;
#else
    extern double TooLong1,TooLong2;
#endif

#define VIS_VISUAL_PARAM_LOGICAL 1
#define VIS_VISUAL_PARAM_INT     2
#define VIS_VISUAL_PARAM_FLOAT   4
#define VIS_VISUAL_PARAM_POINTER 8

typedef struct visual_param_s
{
   /*
    * Name string for the parameter
    */
   char *Name;

   /*
    * Format string for scanf when scanning parameter value
    */
   char *Format;
  
   /*
    * Parameter offset from structure beginning
    */
   long int Offset;

   /*
    *  type of the parameter (one of VIS_VISUAL_PARAM_*)
    */
   int   ParamType;

   /*
    * Initial values for the parameters
    */
   int IntValue;
   double FloatValue;
   void *PointerValue; 
} visual_param_t;

typedef struct visual_type_s
{
    /*
     * Next visual type in the list
    */
    struct visual_type_s *Next;

    /*
     * One line description of the visual
     */
    char *VisualName;

    /*
     *  Alloc memory for parameters
     */
    void *(*AllocParams)();

    /*
     * Dispose a visual
     */
    void (*DeleteParams)(void *);

    /*
     * Realize a visual instance
     */
    int (*RealizeVisual)(geometry_t *,element_model_t *,void *,double);

    /*
     * Structure holding names and offsets of particular
     * visual types parameters
     */
    visual_param_t *VisualParams;
} visual_type_t;

/*
 *  visual type list (remove this)
 */
typedef struct visual_defs_s
{
    int NumberOfVisualTypes;
    visual_type_t *VisualTypes;
} visual_defs_t;

/*
 * Global variable holding start of the list of visual types
 */
EXT visual_defs_t VisualDefs;

typedef struct visual_s
{
    /*
     *  Next visual in list
     */
    struct visual_s *Next;

    /*
     * Name of the visual instance
     */
    char *Name;

   /*
    * this  is pointer to arrow_t,mesh_t etc. structures, see below
    */
    void *VisualParams;

    /*
     * Pointer to visual type structure holding the function pointers
     * acting on this type of visual 
     */

    visual_type_t *VisualType;
} visual_t;

typedef enum
{
    mesh_style_none,mesh_style_line,mesh_style_surf,mesh_style_line_and_surf
} mesh_style_t;

typedef enum
{
    arrow_style_stick,arrow_style_arrow
} arrow_style_t;

typedef enum
{
    particle_style_vector,particle_style_sphere
} particle_style_t;

typedef enum
{
    particle_integ_euler, particle_integ_runge_kutta
} particle_integ_method_t;

typedef enum
{
    particle_policy_fixed, particle_policy_adaptive
} particle_integ_policy_t;


int vis_initialize_arrow_visual();
int vis_initialize_colscale_visual();
int vis_initialize_contour_line_visual();
int vis_initialize_isosurface_visual();
void vis_triangle ( triangle_t *t,vertex_t *v,double *color,double CScl,double CAdd);
int vis_initialize_mesh_visual();
int vis_initialize_particle_visual();
int vis_initialize_sphere_visual();
int vis_add_visual_type( visual_type_t *VisualDef );
visual_type_t *vis_get_visual_type(char *name);
char *vis_get_visual_type_name(visual_t *visual);
char *vis_get_visual_name(visual_t *visual);
int vis_set_param( visual_t *visual, char *name,int intvalue,double doublevalue,void *ptrvalue );
void vis_default_params( visual_t *visual );
visual_t *vis_new_visual(char *name);
visual_t *vis_link_visual( visual_t *list,visual_t *new );
visual_t *vis_add_visual(visual_t *visual,char *name);
void vis_delete_visual( visual_t *visual );
int vis_display_visual( geometry_t *geometry, element_model_t *model, visual_t *VL,double t );
int vis_display_list( geometry_t *geometry, element_model_t *model, visual_t *VL,double t );
int vis_initialize_visual_types();
