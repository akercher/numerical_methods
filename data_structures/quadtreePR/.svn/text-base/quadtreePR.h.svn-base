/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: quadtreePR.h                           */
/* Description: Defines data structure and lists   */
/*   functions related to a point region data      */
/*   structure.                                    */ 
/*-------------------------------------------------*/

/* defines attributes of node in quadtree */
typedef struct _node node;

struct _node
{
  /* node id */
  int nid;

  /* parent node */
  node* pnode;  

  /* number of points in quadrant */
  int npts;  
  
  /* point coordinates */
  double* px;
  double* py;

  /* box coordinates */
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  
  /* quadrants */
  node* qNE;
  node* qNW;
  node* qSW;
  node* qSE;

  /* Euclidean 2-norm */
  double norm2;
  int inRadius;
	
};

typedef struct _tree tree;
struct _tree 
{
  node* root;

};

#ifndef ADD_H_GUARD
#define ADD_H_GUARD
/* Prototypes */
node* nodeCreate(double u0,double u1,double v0,double v1);
void nodeCreateLeaf(node* p);
void pointInsert(node* p,double* xptr,double* yptr);
void nodePrint(node* v);
void nodeGet(node* v,node** nget,double xfind,double yfind);
int nodeNchild(node* v);
int nodeQuadrant(node* v);
int pointQuadrant(node* v,double* xptr,double* yptr);
void nodeQuadDevide(node* p);
void nodeRemove(node* v);
void nodeSeperation(node* v,double xp,double yp,double rp,double r);
double* nodeBoundBox(int npoin,double* x,double* y);
void treeCreate(int npoin,double* x,double* y);
#endif
