/* dyn_tree.c == Dynamic tree routines, see dyn_tree.h */
/*
   Implemented by
   Tamas Badics, 1991,
   Rutgers University, RUTCOR
   P.O.Box 5062
   New Brunswick, NJ, 08901

  e-mail: badics@rutcor.rutgers.edu
*/

#include "_dyn_tree.h"
#include "dyn_tree.h"
#include "macros.h"
#include <stdio.h>
#include <stdlib.h>

/*==================================================================*/
void dyn_make_tree(dyn_item *item, int value)
/* Put item to a new 1-node
   tree with value. */
{
    dyn_node *node;

    node = (dyn_node *) malloc(sizeof(dyn_node));

    item->back = node;
    node->item = item;

    node->father = node->left = node->right = NULL;

    node->dval = value;
    node->dmin = 0;
}

/*==================================================================*/
void dyn_link(dyn_item *oldroot, dyn_item *newfather, DOUBLE new_value)
/* Hang the tree rooted in oldroot
   to newfather. Must be in
   different trees! */
{
    dyn_node *old_r = oldroot->back;
    dyn_node *nf = newfather->back;

    if (old_r->father)
        return;

    dyn_splay(nf);

    if (nf->father == old_r) {
        printf("Error in dyn_link: Both in the same tree!\n");
        return;

        exit(1);
    }
    old_r->father = nf;
    old_r->dval = new_value;
    old_r->dmin = 0;
}

/*==================================================================*/
void dyn_cut(dyn_item *cuthere)   /* Cut the tree between cuthere
                    and its father */

{
    dyn_node *r, *c, *l;

    if ((c = cuthere->back)->father == NULL)
        return;

    dyn_splay(c);

    if ((r = c->right)) {
        r->father = c->father;
        c->right = NULL;
        dvr += dvc;
    }
    if ((l = c->left)) {
        c->left = NULL;
        dvl += dvc;
    }
    c->father = NULL;
    dmc = 0;
}

/*==================================================================*/
dyn_item *dyn_find_root(dyn_item *item)
/* Find the root of item's tree */

{
    dyn_node *r = item->back;

    return (r->father ? dyn_splay(r)->father->item : item);
}

/*==================================================================*/
void dyn_add_value(dyn_item *from, DOUBLE value)
/* Add value to each node
   on the path from 'from'
   to the root */
{
    dyn_node *m, *r, *f = from->back;

    dyn_splay(f);

    dmf = MAX2(0, (r = f->right) ? dmr - dvr : 0);
    if ((m = f->left)) {
        dvf += value;
        dvm -= value;
        dmf = MAX2(dmf, dmm - dvm);
    } else
        dvf += value;
}

/*==========================================================*/
DOUBLE dyn_find_value(dyn_item *item)   /* Find the value of item */

{
    dyn_node *n = item->back;

    dyn_splay(n);
    return (n->dval);
}

/*==========================================================*/
dyn_item *dyn_find_bottleneck(dyn_item *from, DOUBLE neck)
/* Find the bottleneck on the path
   from 'from' to its root. That is
   return the nearest ancestor of
   'from', whose value is <= neck.
   Otherwise return the root. */
{
    int i;
    dyn_node *f, *c, *l;  /* f == from
                c == candidate
                l == c->left
                */

    DOUBLE valc;

    f = from->back;

    if (f->father == NULL)              /* from is The Root. */
        return from;

    dyn_splay(f);

    if (dvf - dmf > neck)
        return (f->father->item);          /* No bottleneck */

    if (dvf <= neck)
        return (from);                     /* from is the bottleneck */

    if ((c = f->right) == NULL || (valc = dvf + dvc) - dmc > neck)
        return (f->father->item);          /* No bottleneck */

    /* There is a bottleneck in the right subtree of f
       (rooted with c) */

    l = c->left;
    while ((i = (l ? valc + dvl - dml <= neck : 0)) || valc > neck) {
        if (i)
            c = l;
        else
            c = c->right;

        l = c->left;
        valc += dvc;
    }

    return (dyn_splay(c)->item);
}

/*==================================================================*/
dyn_item *dyn_find_father(dyn_item *item) {
    dyn_node *n = item->back;
    dyn_node *r;

    if (n->father == NULL)
        return NULL;

    dyn_splay_solid(n);

    if (r = n->right) {
        while (r->left)
            r = r->left;
        return r->item;
    } else
        return (n->father->item);
}
