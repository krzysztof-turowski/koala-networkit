/* dyn_tree.h == Dynamic tree data structure type definitions 

   Based on Sleator, Tarjan : Self-Adjusting Binary Search Trees
   JACM Vol. 32., No. 3, July 1985 pp.652-686 
   
   Slight difference: The Real roots have only middle children,
   so that infinite value on the root can be handled */
/*
   Implemented by 
   Tamas Badics, 1991, 
   Rutgers University, RUTCOR
   P.O.Box 5062
   New Brunswick, NJ, 08901
 
  e-mail: badics@rutcor.rutgers.edu
*/

#ifndef DYN_T
#define DYN_T

/*********************************************************
  Usage of the dynamic tree routines:

  The main purpose of using this data structure is to deal with
  trees of any type of objects with a value associated with them. 
  See the literature! (linking-cutting trees)
  
  The user need only the dyn_tree.h header file and the object
  files which can be created with the enclosed makefile.
  All the other sourcefiles could be referred to as a black box.

  So the only requirement from the objects is to be linkable
  into a dynamic tree. Hence the objects (which are naturally 
  structures) have to be a member which is a void pointer.

  In fact the object is connected through this void pointer to
  the dynamic_tree node (which has another pointer pointing
  back to the dyn_item, and created by dyn_make_tree).
  This allowes to convert a pointer to the object into a pointer
  to dyn_item. 

  For example:
  typedef structure OBJECT {
     anything                 \
	 .                         >  with length offset 
	 .                        /
	 void * dyn_point;
	 .
	 .
	 .
  }object;
  
  object * objptr;

  It is wiser to place our "void * dyn_point" into the first field
  of the structure, but not necessary.
  And the dynamic_tree functions can be called as follows:

  dyn_make_tree((dyn_item *)((int)(objptr) + offset), value);   

  and so on.

  If offset == 0 then our objptr can be converted directly.
  dyn_make_tree((dyn_item *)(objptr), value);  

  And similarly careful returnvalue handling is necessary.
  If offset == 0 then objptr can be identified with a dyn_item * 
  type pointer from the dynamic tree routins' point of view.

  The DOUBLE type is a macro, which can be changed to int or float or
  any number type under the compilation from the command line or
  from the makefile.

*******************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

#ifndef DOUBLE
#define DOUBLE int
#endif

typedef struct DYN_ITEM{
	void * back;
}dyn_item;

void dyn_make_tree(dyn_item * item, int value);
                                        /* Put item to a new 1-node
										   tree with value */

void dyn_link(dyn_item * oldroot, dyn_item * newfather, DOUBLE new_value);
                                        /* Hang the tree rooted in oldroot
										   to newfather with new_value. 
										   Must be in different trees! */

void dyn_cut(dyn_item * cuthere);       /* Cut the tree between cuthere 
										   and its father */

void dyn_add_value(dyn_item * from, DOUBLE value);    
                                        /* Add value to each node 
										   on the path from 'from' 
										   to the root */

DOUBLE dyn_find_value(dyn_item * item); /* Find the value of item */

dyn_item * dyn_find_bottleneck(dyn_item * from, DOUBLE neck);
                                        /* Find the bottleneck on the path 
										   from 'from' to its root. That is
										   return the nearest ancestor of
										   'from', whose value is <= neck. */

dyn_item * dyn_find_root(dyn_item * item);
                                        /* Find the root of item's tree */

dyn_item * dyn_find_father(dyn_item * item);
                                        /* Find the father in the tree */

#endif

#ifdef __cplusplus
}
#endif
