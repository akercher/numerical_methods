/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: quadtreePR.c                           */
/* Description: Uses quadtreePR.h to create a      */
/*   point region quadtree data structure.         */
/*-------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quadtreePR.h"

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* corresponds to quadrants */
enum{NE = 1,NW = 2,SW = 3,SE = 4};

/* Global Variables */
tree mytree;
int id = -1;

/*------------*/
/* treeCreate */
/*------------*/
void treeCreate(int npoin,double* x,double* y)
{
  int i;
  double* bbox = NULL;
  double eps = 0.005;
  node* inode;

  /* compute bounding box   */
  bbox = nodeBoundBox(npoin,x,y);

  /* initialize tree by creating root node   */
  inode = malloc(sizeof(node));  
  inode = nodeCreate(bbox[0]-eps,bbox[1]+eps,bbox[2]-eps,bbox[3]+eps);
  id = id + 1;

  free(bbox);

  mytree.root = inode;

  /* create the leaves of the root node */
  nodeCreateLeaf(mytree.root);

  /* insert points into tree */
  for(i=0;i<npoin;i++)
    {
      pointInsert(mytree.root,&x[i],&y[i]);
      id = id + 1;
    }

  /* define point and radius to base search off of */
  double xfind = 0.34;
  double yfind = 0.1;
  double rfind2 = xfind*xfind + yfind*yfind;
  double rad = 0.15;

  /* find points within radius (defined above) point that is defined above */
  nodeSeperation(mytree.root,xfind,yfind,rfind2,rad);
}

/*------------*/
/* nodeCreate */
/*------------*/
node* nodeCreate(double u0,double u1,double v0,double v1)
{
  double xm,ym;
  node* v;
  v = (node*) malloc(sizeof(node));
  
  xm = (u0+u1)/2.0;
  ym = (v0+v1)/2.0;

  v -> nid = id;
  v -> pnode = NULL;
  v -> npts = 0;

  v -> px = NULL;
  v -> py = NULL;

  v -> xmin = u0;
  v -> xmax = u1;
  v -> ymin = v0;
  v -> ymax = v1;

  v -> qNE = NULL;
  v -> qNW = NULL;
  v -> qSW = NULL;
  v -> qSE = NULL;

  v -> norm2 = xm*xm + ym*ym;
  v -> inRadius = 0;

  
  return v;
}

/*----------------*/
/* nodeCreateLeaf */
/*----------------*/
void nodeCreateLeaf(node* p)
{
  double xm, ym;
  node* inode;

  inode = malloc(sizeof(node));

  /* midpoint of quadrant   */
  xm = (p->xmax + p->xmin)/2.0;
  ym = (p->ymax + p->ymin)/2.0;

  /* northeast leaf */  
  inode = nodeCreate(xm,p->xmax,ym,p->ymax);
  inode->pnode = p;
  p->qNE = inode;

  /* northwest leaf */
  inode = nodeCreate(p->xmin,xm,ym,p->ymax);
  inode->pnode = p;
  p->qNW = inode;

  /* southwest leaf */
  inode = nodeCreate(p->xmin,xm,p->ymin,ym);
  inode->pnode = p;
  p->qSW = inode;

  /* southeast leaf */
  inode = nodeCreate(xm,p->xmax,p->ymin,ym);
  inode->pnode = p;
  p->qSE = inode;
  
}

/*-------------*/
/* pointInsert */
/*-------------*/
void pointInsert(node* p,double* xptr,double* yptr)
{

  int qp;

  /* parent */
  if(p->nid < 0){
    /* find quadrant point is in w.r.t. parent */
    qp = pointQuadrant(p,xptr,yptr);

    if(qp == NE){
      pointInsert(p->qNE,xptr,yptr);
    }
    if(qp == NW){
      pointInsert(p->qNW,xptr,yptr);
    }
    if(qp == SW){
      pointInsert(p->qSW,xptr,yptr);
    }
    if(qp == SE){
      pointInsert(p->qSE,xptr,yptr);
    }
  }
  /* leaf */
  else{
    /* leaf not occupied */
    if((p->px == NULL) && (p->py == NULL)){
      p->nid = id;
      p->px = xptr;
      p->py = yptr;
      p->norm2 = (*xptr)*(*xptr) + (*yptr)*(*yptr);
      p->npts = 1;
    }
    /* leaf is occupied */
    else{
      /* devide leaf into quadrants and reinsert point associated with leaf */
      nodeQuadDevide(p);

      /* now insert point in appropriate quadrant */
      qp = pointQuadrant(p,xptr,yptr);
      if(qp == NE){
	pointInsert(p->qNE,xptr,yptr);
      }
      if(qp == NW){
	pointInsert(p->qNW,xptr,yptr);
      }
      if(qp == SW){
	pointInsert(p->qSW,xptr,yptr);
      }
      if(qp == SE){
	pointInsert(p->qSE,xptr,yptr);
      }
    }
  }
}

/*---------*/
/* nodeGet */
/*---------*/
void nodeGet(node* v,node** nget,double xfind,double yfind)
{
  if(!v){
    return;
  }

  /* leaf contains point --> found node */
  if((v->px != NULL) && (v->py != NULL)){
    if((*(v->px) == xfind) && (*(v->py) == yfind)){
      *nget = v;
      return;
    }
  }

  nodeGet(v->qNE,nget,xfind,yfind);
  nodeGet(v->qNW,nget,xfind,yfind);
  nodeGet(v->qSW,nget,xfind,yfind);
  nodeGet(v->qSE,nget,xfind,yfind);
}

/*------------*/
/* nodeNchild */
/*------------*/
int nodeNchild(node* v)
{
  int cnt = 0;

  /* only count children that contain points */
  if((v->qNE->px != NULL) && (v->qNE->py !=NULL)){
    cnt = cnt + 1;
  }
  if((v->qNW->px != NULL) && (v->qNW->py !=NULL)){
    cnt = cnt + 1;
  }
  if((v->qSW->px != NULL) && (v->qSW->py !=NULL)){
    cnt = cnt + 1;
  }
  if((v->qSE->px != NULL) && (v->qSE->py !=NULL)){
    cnt = cnt + 1;
  }
  
  return cnt;
}

/*--------------*/
/* nodeQuadrant */
/*--------------*/
int nodeQuadrant(node* v)
{
  /* quandrant w.r.t. parent */
  if(v->pnode->qNE == v){
    return NE;
  }
  else if(v->pnode->qNW == v){
    return NW;
  }
  else if(v->pnode->qSW == v){
    return SW;
  }
  else{
    return SE;
  }
}

/*---------------*/
/* pointQuadrant */
/*---------------*/
int pointQuadrant(node* v,double* xptr,double* yptr)
{

  double xm;
  double ym;

  xm = (v->xmax + v->xmin)/2.0;
  ym = (v->ymax + v->ymin)/2.0;

  /* quandrant w.r.t. parent */
  if(*xptr > xm){
    if(*yptr > ym){
      return NE;
    }
    else{
      return SE;
      }
  }
  else{
    if(*yptr > ym){
      return NW;
    }
    else{
      return SW;
    }
  }
}

/*----------------*/
/* nodeQuadDevide */
/*----------------*/
void nodeQuadDevide(node* p)
{

  int idold;
  int qn;
  int qp;
  double xm;
  double ym;
  double xold;
  double yold;
  double* xptr;
  double* yptr;

  /* get information that corresponds to leaf */
  xptr = p->px;
  yptr = p->py;
  xold = *(p->px);
  yold = *(p->py);

  idold = p->nid;

  /* erase leaf's information */
  p->px = NULL;
  p->py = NULL;
  xm = (p->xmax + p->xmin)/2.0;
  ym = (p->ymax + p->ymin)/2.0;
  p->norm2 = xm*xm + ym*ym;
  p->npts = 0;

  /* define id based on quadrant w.r.t. parent */
  qn = nodeQuadrant(p);

  if(qn == NE){
    p->nid = (p->pnode->nid)*4 + 2;
  }
  if(qn == NW){
    p->nid = (p->pnode->nid)*4 + 1;
  }
  if(qn == SW){
    p->nid = (p->pnode->nid)*4;
  }
  if(qn == SE){
    p->nid = (p->pnode->nid)*4 - 1;
  }

  /* create node's quadrants */
  nodeCreateLeaf(p);

  /* place point that was previously in leaf into proper quadrant */
  qp = pointQuadrant(p,&xold,&yold);

  if(qp == NE){
    p->qNE->nid = idold;
    p->qNE->px = xptr;
    p->qNE->py = yptr;
    p->qNE->norm2 = xold*xold + yold*yold;
    p->qNE->npts = 1;    
  }
  if(qp == NW){
    p->qNW->nid = idold;
    p->qNW->px = xptr;
    p->qNW->py = yptr;
    p->qNW->norm2 = xold*xold + yold*yold;
    p->qNW->npts = 1;    
  }
  if(qp == SW){
    p->qSW->nid = idold;
    p->qSW->px = xptr;
    p->qSW->py = yptr;
    p->qSW->norm2 = xold*xold + yold*yold;
    p->qSW->npts = 1;    
  }
  if(qp == SE){
    p->qSE->nid = idold;
    p->qSE->px = xptr;
    p->qSE->py = yptr;
    p->qSE->norm2 = xold*xold + yold*yold;
    p->qSE->npts = 1;    
  }

}

/*------------*/
/* nodeRemove */
/*------------*/
void nodeRemove(node* v)
{
  int nchild;
  int nQuad; 
  node* parent;

  /* if node does not exist, exit */
  if(!v){
    return;
  }

  nQuad = nodeQuadrant(v);
  parent = v->pnode;

  /* sever my tie with my parent */
  switch(nQuad)
    {
    case NE:
      parent->qNE = NULL;
      break;
    case NW:
      parent->qNW = NULL;
      break;
    case SW:
      parent->qSW = NULL;
      break;
    case SE:
      parent->qSE = NULL;
      break;
    }

  nchild = nodeNchild(v->pnode);

  /* if only one child left --> reinsert */
  if(nchild == 2){
    if(parent->qNE != NULL){
      parent->px = parent->qNE->px;
      parent->py = parent->qNE->py;
      parent->nid = parent->qNE->nid;      
      parent->qNE = NULL;
    }
    if(v->qNW != NULL){
      parent->px = parent->qNW->px;
      parent->py = parent->qNW->py;
      parent->nid = parent->qNW->nid;      
      parent->qNW = NULL;
    }
    if(v->qSW != NULL){
      parent->px = parent->qSW->px;
      parent->py = parent->qSW->py;
      parent->nid = parent->qSW->nid;      
      parent->qSW = NULL;
    }
    if(v->qSE != NULL){
      parent->px = parent->qSE->px;
      parent->py = parent->qSE->py;
      parent->nid = parent->qSE->nid;      
      parent->qSE = NULL;
    }
  }
  
  /* deallocate memory */
  free(v);
}

/*----------------*/
/* nodeSeperation */
/*----------------*/
void nodeSeperation(node* v,double xp,double yp,double rp,double r)
{

  if(!v){
    return;
  }

  double xm;
  double ym;
  double d;

  /* only check if node is leaf */
  if((v->px != NULL) && (v->py != NULL)){
    d = v->norm2 + rp - 2*((*(v->px)*xp) + (*(v->py)*yp));
    if(d <= r*r){
      v->inRadius = 1;
    }
  }
  else{
    xm = (v->xmax + v->xmin)/2.0;
    ym = (v->ymax + v->ymin)/2.0;
    if((xm < xp+r) && (ym < yp+r)){
      nodeSeperation(v->qNE,xp,yp,rp,r);
    }
    if((xm > xp-r) && (ym < yp+r)){
      nodeSeperation(v->qNW,xp,yp,rp,r);
    }
    if((xm > xp-r) && (ym > yp-r)){
      nodeSeperation(v->qSW,xp,yp,rp,r);
    }
    if((xm < xp+r) && (ym > yp-r)){
      nodeSeperation(v->qSE,xp,yp,rp,r);
    }
  }
}

/*--------------*/
/* nodeBoundBox */
/*--------------*/
double* nodeBoundBox(int npoin,double* x,double* y)
{

  int i;
  double* rect = NULL;

  rect = malloc(sizeof(double)*4);

  rect[0] = x[0];
  rect[1] = x[1];
  rect[2] = y[0];
  rect[3] = y[1];

  for(i=1;i<npoin;i++){
    rect[0] = min(rect[0],x[i]);
    rect[1] = max(rect[1],x[i]);
    rect[2] = min(rect[2],y[i]);
    rect[3] = max(rect[3],y[i]);
  }

  return rect;
}

/*-----------*/
/* nodePrint */
/*-----------*/
void nodePrint(node* v)
{
  if(!v){
    return;
  }

  double* xp;
  double* yp;

  if(v->nid == -1){
    printf("node(%d) --> root\n",v->nid);
  }
  else{
    int nquad = nodeQuadrant(v);
    if(v->nid < 0){
      printf("node(%d) --> Q = %d, P = %d\n",v->nid,nquad,v->pnode->nid);
    }
    if(v->px != NULL){
      xp = v->px;
      yp = v->py;
      printf("node(%d): x = %lf y = %lf --> Q = %d, P = %d\n",v->nid,*xp,*yp,nquad,v->pnode->nid);
    }
  }

  nodePrint(v->qNE);
  nodePrint(v->qNW);
  nodePrint(v->qSW);
  nodePrint(v->qSE);

}
