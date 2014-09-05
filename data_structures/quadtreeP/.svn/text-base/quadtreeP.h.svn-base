/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: quadtreeP.h                            */
/* Description: Defines data structure and         */
/*   operations related to a quadtree data         */
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
  
  /* point coordinates */
  double x;
  double y;

  /* Euclidean 2-norm */
  double norm2;
  int inRadius;
  
  /* children nodes */
  node* cnodeNE;
  node* cnodeNW;
  node* cnodeSW;
  node* cnodeSE;
	
};

typedef struct _tree tree;
struct _tree 
{
  node* root;

};

#ifndef ADD_H_GUARD
#define ADD_H_GUARD
/* Prototypes */
node* nodeCreate(double x,double y);
void nodeInsert(node* p,node* v);
void nodePrint(node* v);
void nodeGet(node* v,node** nget,double x,double y);
int nodeNchild(node* v);
int nodeQuadrant(node* v);
void nodeRemove(node* v);
void nodeSeperation(node* v,double xp,double yp,double rp,double r);
double nodeBounds(node* v,double rmax);
void nodeSubBounds(node* v,double x,double y);
void treeCreate(int npoin,double* x,double* y,double r);
#endif

