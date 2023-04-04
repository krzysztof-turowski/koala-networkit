/* dyn_splay.c ==== rotate subroutines for the splay step in dynamic tree
   data structure */
/*
   Implemented by
   Tamas Badics, 1991,
   Rutgers University, RUTCOR
   P.O.Box 5062
   New Brunswick, NJ, 08901

  e-mail: badics@rutcor.rutgers.edu
*/

#include <stdio.h>
#include "_dyn_tree.h"
#include "macros.h"

/*==========================================================================*/
static void __dyn_step_r( dyn_node * f, dyn_node * ch)
{
  DOUBLE tmp;
  dyn_node * b, * c, *chf;

  b = ch->left;
  c = f->left;

  if( NULL != (f->right = b) )
    b->father = f;
  ch->left = f;

  chf = ch->father = f->father;
  if( chf->right == f)
    chf->right = ch;
  else if(chf->left == f)
    chf->left = ch;

  f->father = ch;

  /* Value update: */

  tmp = dvh;
  if (b)
    dvb += tmp;
  dmh = tmp + dmf;
  dvh += dvf;
  dvf = -tmp;
  dmf = MAX3( 0, (b ? dmb - dvb : 0), (c ? dmc - dvc : 0));
}

/*==========================================================================*/
static void __dyn_step_l( dyn_node  * f, dyn_node * ch)
{
  DOUBLE tmp;
  dyn_node * b, * c, * chf;

  b = ch->right;
  c = f->right;

  if(NULL != (f->left = b) )
    b->father = f;
  ch->right = f;

  chf = ch->father = f->father;
  if( chf->left == f)
    chf->left = ch;
  else if(chf->right == f)
    chf->right = ch;

  f->father = ch;

  /* Value update: */

  tmp = dvh;
  if (b)
    dvb += tmp;
  dmh = tmp + dmf;
  dvh += dvf;
  dvf = -tmp;
  dmf = MAX3( 0, (b ? dmb - dvb : 0), (c ? dmc - dvc : 0));
}

/*==========================================================================*/
static void __dyn_step_rl(dyn_node * g, dyn_node * f, dyn_node * ch)
{
  DOUBLE tmp;
  dyn_node * a, * b, * c, * d, * chf;

  a = f->right;
  b = ch->right;
  c = ch->left;
  d = g->left;

  if( NULL != (g->right = c))
    g->right->father = g;

  if (NULL != (f->left = b))
    f->left->father = f;
  ch->left   = g;
  ch->right   = f;
  chf = ch->father = g->father;

  if( chf->right == g)
    chf->right = ch;
  else if(chf->left == g)
    chf->left = ch;

  g->father = f->father = ch;

  /*======== Value update: ========*/

  tmp = dvh + dvf;
  if (b)
    dvb += dvh;
  if (c)
    dvc += tmp;
  dvf  = -dvh;
  dvh = tmp + dvg;
  dvg  = -tmp;
  dmh = dmg + tmp;
  dmg = MAX3( 0, (d ? dmd - dvd : 0), (c ? dmc - dvc : 0));
  dmf = MAX3( 0, (a ? dma - dva : 0), (b ? dmb - dvb : 0));
}

/*==========================================================================*/
static void __dyn_step_lr(dyn_node * g, dyn_node * f, dyn_node * ch)
{
  DOUBLE tmp;
  dyn_node * a, * b, * c, * d, * chf;

  a = f->left;
  b = ch->left;
  c = ch->right;
  d = g->right;

  if (NULL != (g->left = c) )
    g->left->father = g;
  if( NULL != (f->right = b) )
    f->right->father = f;
  ch->right   = g;
  ch->left   = f;
  chf = ch->father = g->father;

  if( chf->right == g)
    chf->right = ch;
  else if(chf->left == g)
    chf->left = ch;

  g->father = f->father = ch;

  /*======== Value update: ========*/

  tmp = dvh + dvf;
  if (b)
    dvb += dvh;
  if (c)
    dvc += tmp;
  dvf  = -dvh;
  dvh = tmp + dvg;
  dvg  = -tmp;
  dmh = dmg + tmp;
  dmg = MAX3( 0, (d ? dmd - dvd : 0), (c ? dmc - dvc : 0));
  dmf = MAX3( 0, (a ? dma - dva : 0), (b ? dmb - dvb : 0));

}

/*==========================================================================*/
static void __dyn_step_rr( dyn_node * g, dyn_node * f, dyn_node * ch)
{
  DOUBLE tmp, tmp1;
  dyn_node * b, * c, * d, * chf;

  b = ch->left;
  c = f->left;
  d = g->left;

  if( NULL != (g->right = c) )
    g->right->father = g;
  f->left   = g;
  if( NULL != (f->right = b))
    f->right->father = f;
  ch->left = f;
  chf = ch->father = g->father;

  if( chf->right == g)
    chf->right = ch;
  else if(chf->left == g)
    chf->left = ch;
  ch->father   = g->father;

  g->father = f;
  f->father = ch;

  /*======== Value update: ========*/

  if (b)
    dvb += dvh;
  if (c)
    dvc += dvf;
  tmp = dvh + dvf;
  dmh = tmp + dmg;
  dmg = MAX3( 0, (c ? dmc - dvc : 0), (d ? dmd - dvd : 0));
  dmf = MAX3( 0, (b ? dmb - dvb : 0), dvf + dmg);
  tmp1 = -dvh;
  dvh = tmp + dvg;
  dvg = -dvf;
  dvf = tmp1;
}

/*==========================================================================*/
static void __dyn_step_ll( dyn_node * g, dyn_node * f, dyn_node * ch)
{
  DOUBLE tmp, tmp1;
  dyn_node * b, * c, * d, * chf;

  b = ch->right;
  c = f->right;
  d = g->right;

  if( NULL != (g->left = c) )
    g->left->father = g;
  f->right   = g;
  if (NULL != (f->left = b) )
    f->left->father = f;
  ch->right   = f;
  chf = ch->father = g->father;

  if( chf->right == g)
    chf->right = ch;
  else if(chf->left == g)
    chf->left = ch;
  ch->father   = g->father;

  g->father = f;
  f->father = ch;

  /*======== Value update: ========*/

  if (b)
    dvb += dvh;
  if (c)
    dvc += dvf;
  tmp = dvh + dvf;
  dmh = tmp + dmg;
  dmg = MAX3( 0, (c ? dmc - dvc : 0), ( d ? dmd - dvd : 0));
  dmf = MAX3( 0, (b ? dmb - dvb : 0), dvf + dmg);
  tmp1 = -dvh;
  dvh = tmp + dvg;
  dvg = -dvf;
  dvf = tmp1;
}

/*==========================================================================*/
void splice(dyn_node * m)   /* m is a solid root, middle-child
                 of its f, who is also a solid root.
                 Change m to be the left-child.
                 (the old left-child will be a middle)*/
{
  dyn_node * f = m->father;
  dyn_node * l = f->left;
  dyn_node * r = f->right;

  f->left = m;  /* Changing the solid edge f-l to f-m */

  /* Here comes the value update */

  dvm -= dvf;
  if (l)
    dvl += dvf;
  dmf = MAX3( 0,(r ? dmr - dvr : 0), dmm - dvm);
}

/*==========================================================================*/
void dyn_splay_solid(dyn_node * ch )    /* Brings up the ch in its
                       solid subtree. After this
                       procedure ch will be the
                       root of the solid tree.

                       ch cannot be The Root!*/

{
  dyn_node * g ;
  dyn_node * f ;
  for(;;){
    f = ch->father;
    g = f->father;

    if (f->left == ch){
      if (g->left == f)
        __dyn_step_ll(g, f, ch);
      else if (g->right == f)
        __dyn_step_rl(g, f, ch);
      else
        __dyn_step_l(f, ch);
    }else if (f->right == ch){
      if (g->left == f)
        __dyn_step_lr(g, f, ch);
      else if (g->right == f)
        __dyn_step_rr(g, f, ch);
      else
        __dyn_step_r(f, ch);
    }else
      break;
  }
}
/*==========================================================================*/

dyn_node * dyn_splay(dyn_node * ch )   /* Brings up the ch in its virtual
                      tree. After this procedure, ch
                      will be a mid-child of the root
                      of the virtual tree*/
{
  dyn_node * f;


  if (ch->father == NULL)   /* If ch is The Root */
    return ch;

  dyn_splay_solid(ch);

  while ((f = ch->father)->father){  /* While f is not The Root */

    dyn_splay_solid(f);    /* After this, ch is a middle-child
                 of its f, and f is
                 a solid root */
    splice(ch);
    __dyn_step_l(f, ch);
  }
  return (ch);
}
