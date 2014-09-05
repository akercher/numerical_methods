/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 722: Advancing front grid generation        */
/* Last Updated: 10/03/2011                        */
/*-------------------------------------------------*/
/* Program: faces.h                                */
/* Description: Defines data structure related to  */
/*   faces and grid.                               */ 
/*-------------------------------------------------*/

/* defines attributes of face*/
typedef struct _face face;

struct _face
{
  /* face id */
  int fid;

  /* is face part of active front? 0:no; 1:yes */
  int state;  
  
  /* point coordinates */
  int* fnode;
};

/* defines attributes of grid*/
typedef struct _grid grid;

struct _grid
{
  /* grid id */
  int gid;

  /* number of faces */
  int nface;

  /* is face part of active front? 0:no; 1:yes */
  int* state;  
  
  /* point coordinates */
  int* fnode;
};
