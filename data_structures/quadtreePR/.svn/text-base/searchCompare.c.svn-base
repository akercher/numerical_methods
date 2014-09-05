/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: searchCompare.c                        */
/* Description: Uses root node of point region     */
/*   quadtree data structure to compare the time   */ 
/*   it takes to find all points within a radius   */
/*   of 0.15 of the point (0.34,01) using the      */ 
/*   quadtree and searching by brute force.        */
/*-------------------------------------------------*/

#define ndim 3

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/types.h>
#include <math.h>
#include "quadtreePR.h"


/* Prototypes */
void timeSearch(node* v,int npoin, double* x,double* y);
void bfSeperation(int npoin,double* x,double* y);

/*------------*/
/* timeSearch */
/*------------*/
void timeSearch(node* v,int npoin,double* x,double* y)
{

  int N = 1000000;
  int M = 50;
  int i,j;

  time_t qt0,qt1,bt0,bt1;
  clock_t qc0,qc1,bc0,bc1;
  double qdt,bdt,qdc,bdc;

  double xfind = 0.34;
  double yfind = 0.1;
  double rfind2 = xfind*xfind + yfind*yfind;
  double rt = 0.15;

  printf("\n\tAVERAGE EXECUTION TIME OF SEARCH:\n");
  qt0 = time(NULL);
  qc0 = clock();
  for(j=0;j<N;j++){
    nodeSeperation(v,xfind,yfind,rfind2,rt);
  }
  qt1 = time(NULL);
  qc1 = clock();
  qdt = qt1-qt0;
  qdc = (qc1 - qc0)/CLOCKS_PER_SEC;
  printf("\t\t(1) QUADTREE SEARCH: %E\n",qdt/(double) N);

  bt0 = time(NULL);
  bc0 = clock();
  for(j=0;j<M;j++){
    bfSeperation(npoin,x,y);
  }
  bt1 = time(NULL);
  bc1 = clock();
  bdt = bt1 - bt0;
  bdc = (bc1 - bc0)/CLOCKS_PER_SEC;
  printf("\t\t(2) BRUTE FORCE SEARCH: %E\n",bdt/(double) M);

}

/*--------------*/
/* bfSeperation */
/*--------------*/
void bfSeperation(int npoin,double* x,double* y)
{
  int i,j,k;
  double xp = 0.34;
  double yp = 0.1;
  double rp = xp*xp + yp*yp;
  double r = 0.15;
  double u,v,d;
  int* inRadius = NULL;

  inRadius = malloc(sizeof(int)*npoin*npoin);

  k = 0;
  for(i=0;i<npoin;i++){
    u = x[i];
    for(j=0;j<npoin;j++){
      inRadius[k] = 0;
      v = y[j];
      d = (u*u + v*v) + rp - 2*((u*xp) + (v*yp));
      if(d <= r*r){
	inRadius[i] = 1;
      }
      k = k + 1;
    }
  }
  free(inRadius);
}

/*------*/
/* Read */
/*------*/
/* void getdata(void){ */

/*   int i; */
/*   char buffer[240]; */
/*   char buffer2[240]; */
/*   char buffer3[240]; */
/*   char elem_type[240]; */

/*   for(i=0;i<3;i++){ */
/*     fgets(buffer,240,fid); */
/*   } */

/*   /\* read number of points *\/ */
/*   fscanf(fid,"%d",&npoin); */

/*   /\* allocate memory for array of point coordinates *\/ */
/*   xyz = (double*) malloc(sizeof(double)*ndim*npoin); */

/*   /\* read in point coordinates *\/ */
/*   for(i=0;i<npoin;i++){ */
/*     fscanf(fid,"%lf %lf %lf",&xyz[ndim*i],&xyz[ndim*i+1],&xyz[ndim*i+2]); */
/*   } */

/*   /\* disregard \n and next line, 'END' *\/ */
/*   fgets(buffer2,240,fid); */
/*   fgets(buffer2,240,fid); */

/*   /\* read element type, 'TRIANGLE | TETRA' *\/ */
/*   fgets(elem_type,10,fid); */

/*   /\* disregard next line, 'title' *\/ */
/*   fgets(buffer3,240,fid); */

/*   /\* read number of elements *\/ */
/*   fscanf(fid,"%d",&nelem); */

/*   /\* determine number of nodes for elements *\/ */
/*   if (strcmp(elem_type,"TRIANGLE\n") == 0){ */
/*     nnode = 3; */
/*   } */
/*   else{ */
/*     nnode = 4; */
/*   }   */
/*   /\* allocate memory for array of field values *\/ */
/*   lnode = (int*) malloc(sizeof(int)*nnode*nelem); */

/*   if (nnode == 3){ */
/*     for(i=0;i<nelem;i++){ */
/*       fscanf(fid,"%d %d %d",&lnode[nnode*i],&lnode[nnode*i+1],&lnode[nnode*i+2]); */
/*     } */
/*   } */
/*   else{ */
/*     for(i=0;i<nelem;i++){ */
/*       fscanf(fid,"%d %d %d %d",&lnode[nnode*i],&lnode[nnode*i+1],&lnode[nnode*i+2],&lnode[nnode*i+3]); */
/*     } */
/*   } */

/*   return; */
/* } */
