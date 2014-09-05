/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/

This folder contains material relavent to CSI 703 
assignment #6.  The following is a list of the
included fi1es: 

  1. readme.txt : this file
  2. Makefile : used to compile code on LINUX
  3. searchCompare.c : described below
  4. quadtreePR.h : defines data structure
  5. quadtreePR.c : described below
  6. akercher_assignment06.c : described below

The code is complied and run in the terminal.

  $ make
  $ ./run

When promt, enter the path or file name.

  $ Enter Path/File name: myfile.zfem

/*-----------------*/
/* searchCompare.c */
/*-----------------*/

Prints, to standard output, the time it took to 
find all points within a radius of 0.15 of the 
point (0.34,0.1) using the quadtree and using
a brute force algorithm.  After the quadtree is
created in akercher_assignment06.c, this function
is called.  

/*--------------*/
/* quadtreePR.c */
/*--------------*/

Uses structure defined in quadtreePR.h to create a 
point region quadtree based on the points defined in 
the .zfem file.  The functions within the file are:

  - void treeCreate(int npoin,double* x,double* y);
    ==> main function.  Excepts x and y array defining
        the point coordinates and the number of points
	in the arrays.

  - double* nodeBoundBox(int npoin,double* x,double* y);
    ==> returns bounding box of point coordinates.
	
  - node* nodeCreate(double u0,double u1,double v0,double v1);
    ==> creates node defined by its bounding box.

  - void pointInsert(node* p,double* xptr,double* yptr);	
    ==> Inserts point into tree.

  - void nodeCreateLeaf(node* p);
    ==> If node needs to be divided into quadrants, this 
        creates the node's children.

  - void nodeQuadDevide(node* p);
    ==> Called when attempting to place node in quadrant
        that already contains the maximum number of points.
        Calls nodeCreateLeaf();

  - void nodeGet(node* v,node** nget,double xfind,double yfind);
    ==> Returns node corresponding to point coordinates defined
        by xfind and yfind.  NULL is returned if node does not
	exist.

  - void nodeRemove(node* v);
    ==> Removes node from tree.  Node to be removed is found
        and returned by nodeGet() first;

  - void nodeSeperation(node* v,double xp,double yp,double rp,double r);
    ==> Calculates the distance between nodes of tree and point
        coordinates defined by xp and yp.

  - int pointQuadrant(node* v,double* xptr,double* yptr);
    ==> Returns quadrant of point coordinate w.r.t. particular
        node.

  - int nodeQuadrant(node* v);
    ==> Returns quadrant of node w.r.t. to parent node.

  - int nodeNchild(node* v);
    ==> Returns number of non-empty children corresponding
        to particular node.

  - void nodePrint(node* v);
    ==> Prints all nodes of tree.

/*-------------------------*/
/* akercher_assignment06.c */
/*-------------------------*/

Reads in data from .zfem file and calls treeCreate(), the
main function in quadtreePR.c.  Then uses the returned quadtree
to visualize point coordinates, the building of the tree and
the results of the search, which were points within a radius 
of 0.15 of the point (0.34,0.1).  These point are shown in 
yellow --> green depending on there distance from the point
(0.34,0.1) and are circled in the visualization.  
Options options for minipulating the visualization, and 
corresponding keys are the following:
  1. Points ==> press 'p'.
  2. Points and Quadtree ==> press 't'.
  3. Original View ==> press 'o'
When the mouse is over the active window, the
program can preform the following:
  1. Rotation ==> hold 'left mouse'
  2. Scaling ==> scroll 'mouse wheel'
To exit the program, either press 'q' or 'ESC' on the
keyboard or close the application window.








