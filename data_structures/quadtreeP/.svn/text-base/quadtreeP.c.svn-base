/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: quadtreeP.c                            */
/* Description: Uess quadtreeP.h to create a       */
/*   quadtree data structure.                      */
/*-------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "quadtreeP.h"

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* corresponds to quadrants */
enum{NE = 1,NW = 2,SW = 3,SE = 4};

/* Global Variables */
int id = 0;
tree mytree;

double x0,x1;
double z0,z1;

/*------------*/
/* treeCreate */
/*------------*/
void treeCreate(int npoin,double x[],double y[],double rad)
{
  int i,j;
  double xfind;
  double yfind;
  double rfind2;
  node* inode;

  inode = nodeCreate(x[0],y[0]);
  mytree.root = inode;
  id = id + 1;

  /* put remaining points into tree */
  for(i=1;i<npoin;i++)
    {
      inode = nodeCreate(x[i],y[i]);

      nodeInsert(mytree.root,inode);
      id = id + 1;
    }
	
  /* search for a point */
  xfind = 0.34;
  yfind = 0.1;
  node* nget = NULL;
  nodeGet(mytree.root,&nget,xfind,yfind);

  if(nget != NULL){

    xfind = nget->x;
    yfind = nget->y;
    rfind2 = nget->norm2;
  }
  else{
    rfind2 = xfind*xfind + yfind*yfind;
  }

  nodeSeperation(mytree.root,xfind,yfind,rfind2,rad);
}

/*------------*/
/* nodeCreate */
/*------------*/
node* nodeCreate(double x,double y)
{
  node* v;
  v = (node*) malloc(sizeof(node));
  
  v -> nid=id;
  v -> pnode = NULL;

  v -> x = x;
  v -> y = y;

  v -> norm2 = x*x + y*y;
  v -> inRadius = 0;

  v -> cnodeNE = NULL;
  v -> cnodeNW = NULL;
  v -> cnodeSW = NULL;
  v -> cnodeSE = NULL;
  
  return v;
}

void nodeInsert(node* p, node* v)  
{	
  if(v->x - p->x > 0){
    /* east */
    if(v->y-p->y > 0){
      /* northeast */
      if(p->cnodeNE != NULL){
	nodeInsert(p->cnodeNE,v);
      }
      else{
	p->cnodeNE = v;
	p->cnodeNE->pnode = p;
      }
    }
    else{ 
      /* southeast */
      if(p->cnodeSE != NULL){
	nodeInsert(p->cnodeSE,v);
      }
      else{
	p->cnodeSE = v;
	p->cnodeSE->pnode = p;
      }
    }
  }
  else{ 
    /* west */
    if(v->y-p->y > 0){
      /* northwest */
      if(p->cnodeNW != NULL){
	nodeInsert(p->cnodeNW,v);
      }
      else{
	p->cnodeNW = v;
	p->cnodeNW->pnode = p;
      }
    }
    else{ 
      /* southwest */
      if(p->cnodeSW != NULL){
	nodeInsert(p->cnodeSW,v);
      }
      else{
	p->cnodeSW = v;
	p->cnodeSW->pnode = p;
      }
    }
  }
}

/*-----------*/
/* nodePrint */
/*-----------*/
void nodePrint(node* v)
{
  if(!v){
    return;
  }

  if(v->nid == 0){
    printf("node(%d) --> root x = %f y = %f\n",v->nid,v->x,v->y);
  }
  else{
  int nquad = nodeQuadrant(v);
  printf("node(%d): x = %f y = %f --> Q = %d, P = %d\n",v->nid,v->x,v->y,nquad,v->pnode->nid);
  }

  nodePrint(v->cnodeNE);
  nodePrint(v->cnodeNW);
  nodePrint(v->cnodeSW);
  nodePrint(v->cnodeSE);	 
}

void nodeGet(node* v,node** nget,double x,double y)
{
  if(!v){
    return;
  }
  
  if(v->x == x && v->y == y){
    *nget = v;
    return;
  }
  
  nodeGet(v->cnodeNE,nget,x,y);
  nodeGet(v->cnodeNW,nget,x,y);
  nodeGet(v->cnodeSW,nget,x,y);
  nodeGet(v->cnodeSE,nget,x,y);		
}

/*------------*/
/* nodeNchild */
/*------------*/
int nodeNchild(node* v)
{
  int cnt = 0;
  
  if(v->cnodeNE != NULL){
    cnt = cnt + 1;
  }
  if(v->cnodeNW != NULL){
    cnt = cnt + 1;    
  }  
  if(v->cnodeSW != NULL){
    cnt = cnt + 1;
  }  
  if(v->cnodeSE != NULL){
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
  if(v->pnode->cnodeNE == v){
    return NE;
  } 
  else if(v->pnode->cnodeNW == v){
    return NW;
  }
  else if(v->pnode->cnodeSW == v){
    return SW;
  }
  else{
    return SE;
  }
}

/*------------*/
/* nodeRemove */
/*------------*/
void nodeRemove(node* v)
{
  int nchild = nodeNchild(v);
  int nquad = nodeQuadrant(v);
  node* parent = v->pnode;			

  /* if node does not exist, exit */
  if(!v){
    return;
  }  
  
  /* sever my tie with my parent */
  switch(nquad)
    {
    case NE:
      parent->cnodeNE = NULL;
      break;
    case NW:
      parent->cnodeNW = NULL;
      break;
    case SW:
      parent->cnodeSW = NULL;
      break;
    case SE:
      parent->cnodeSE = NULL;
      break;
    }
  
  /* if created orphans --> reinsert them */
  if(nchild>0){
    if(v->cnodeNE != NULL){
      nodeInsert(parent,v->cnodeNE);
    }
    if(v->cnodeNW != NULL){
      nodeInsert(parent,v->cnodeNW);
    }
    if(v->cnodeSW != NULL){
      nodeInsert(parent,v->cnodeSW);
    }
    if(v->cnodeSE != NULL){
      nodeInsert(parent,v->cnodeSE);
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
  double d;

  if(!v) return;

  d = v->norm2 + rp - 2*((v->x*xp) + (v->y*yp));
  if(d <= r*r){
    v->inRadius = 1;
  }
  
  if((v->x < xp+r) && (v->y < yp+r)){
    nodeSeperation(v->cnodeNE,xp,yp,rp,r);
  }
  if((v->x > xp-r) && (v->y < yp+r)){
    nodeSeperation(v->cnodeNW,xp,yp,rp,r);
  }
  if((v->x > xp-r) && (v->y > yp-r)){
    nodeSeperation(v->cnodeSW,xp,yp,rp,r);
  }
  if((v->x < xp+r) && (v->y > yp-r)){
    nodeSeperation(v->cnodeSE,xp,yp,rp,r);
  }
}

/*------------*/
/* nodeBounds */
/*------------*/
double nodeBounds(node* v,double rmax)
{
  double r;

  if(!v){
    r = rmax;
    return r;
  }

  if(v->norm2 > rmax){
    r = v->norm2;
/*     printf("node(%d) within %f of (%f,%f)\n",v->nid,r,xp,yp); */
  }
  else{
    r = rmax;
  }

  if(v->cnodeNE != NULL){
    r = nodeBounds(v->cnodeNE,r);
  }
  if(v->cnodeNW != NULL){  
    r = nodeBounds(v->cnodeNW,r);
  }
  if(v->cnodeSW != NULL){  
    r = nodeBounds(v->cnodeSW,r);
  }
  if(v->cnodeSE != NULL){  
    r = nodeBounds(v->cnodeSE,r);
  }

  return r;  
}
/*---------------*/
/* nodeSubBounds */
/*---------------*/
void nodeSubBounds(node* v,double x,double y)
{

  if(v->pnode != NULL){
    if(v->pnode->x < x){
      x0 = max(x0,v->pnode->x);
    }
    if(v->pnode->x > x){
      x1 = min(x1,v->pnode->x);
    }
    if(v->pnode->y < y){
      z0 = max(z0,v->pnode->y);
    }
    if(v->pnode->y > y){
      z1 = min(z1,v->pnode->y);
    }
    nodeSubBounds(v->pnode,x,y);
  }
  return;
}
