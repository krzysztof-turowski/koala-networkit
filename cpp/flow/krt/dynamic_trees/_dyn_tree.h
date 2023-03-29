/* _dyn_tree.h == local type and function definitions*/
/*
   Implemented by 
   Tamas Badics, 1991, 
   Rutgers University, RUTCOR
   P.O.Box 5062
   New Brunswick, NJ, 08901
 
  e-mail: badics@rutcor.rutgers.edu
*/

#ifndef _DYN_T
#define _DYN_T

#include "dyn_tree.h"

#define dmg  (g->dmin)
#define dmf  (f->dmin)
#define dmh (ch->dmin)
#define dmr	 (r->dmin)
#define dmm	 (m->dmin)
#define dma	 (a->dmin)
#define dmb	 (b->dmin)
#define dmc	 (c->dmin)
#define dmd	 (d->dmin)
#define dml	 (l->dmin)

#define dvg  (g->dval)
#define dvf  (f->dval)
#define dvh (ch->dval)
#define dvr	 (r->dval)
#define dvm	 (m->dval)
#define dva	 (a->dval)
#define dvb	 (b->dval)
#define dvc	 (c->dval)
#define dvd	 (d->dval)
#define dvl	 (l->dval)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DYN_NODE{
	dyn_item		* item;
	struct DYN_NODE	* left, * right, * father;
	DOUBLE            dmin;
	DOUBLE            dval;
}dyn_node;


dyn_node * dyn_splay(dyn_node * child ); /* Brings up the child in its virtual 
											tree. After this procedure child will
											be a middle-child of The Root.
											Returns the pointer to child.*/

void dyn_splay_solid(dyn_node * child );   /* Brings up the child in its 
											  solid subtree. After this 
											  procedure child will be the 
											  root of the solid tree*/

void splice(dyn_node * m);   /* m is a solid root, middle-child 
								of its f, who is also a solid root.
								Change m to be the left-child. 
								(the old left-child will be a middle)*/

#endif

#ifdef __cplusplus
}
#endif