/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: akercher_assignment06.c                */
/* Description: The program reads a .zfem file     */ 
/*   containing 3D point coordinates  and element  */
/*   connectivity.  Points within a radius of 0.15 */
/*   of the point (0.34,0.1) are shown in green.   */
/*   Options, and corresponding key are the        */
/*   following:                                    */
/*     1. Points ==> press 'p'.                    */
/*     2. Points and Quadtree ==> press 't'.       */
/*     3. Original View ==> press 'o'              */
/*   When the mouse is over the active window, the */
/*   program can preform the following:            */
/*     1. Rotation ==> hold 'left mouse'           */
/*     2. Scaling ==> scroll 'mouse wheel'         */
/*   See readme file for more information.         */
/*   To exit the program, either press 'q' on the  */
/*   keyboard or close the application window.     */
/*   Finally, the name of the file used to         */
/*   visualize the set must be less than 100       */
/*   characters in length.                         */
/* Example: The program is compiled and run using  */ 
/*   the following terminal commands:              */
/*     $ make                                      */
/*     $ ./run                                     */
/*   The user will be prompt to enter the name of  */
/*   the file they wish to visualize.              */
/*     $ Enter Path/File name: myfile.zfem         */
/*-------------------------------------------------*/


/*----------*/
/* Includes */
/*----------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include "quadtreePR.h"

/* define the number of dimensions for grid (1D,2D,3D,...,ND) */
#define ndim 3
/* define number of models */
#define nkey 2

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* corresponds to mouse buttons */
enum{UP = 1,DOWN = 0};

/* Prototypes */
void init(void);
void state_reset(void);
void getdata(void);
void drawBounds(void);
void drawQuadrants(node*);
void drawSearch(void);
void drawElements(void);
void display(void);
void reshape(int,int);

void keyboard(unsigned char,int,int);
void motion(int,int);
void mouse(int,int,int,int);
void mousewheel(int,int,int,int);

/* Global Variables */
tree mytree;

FILE *fid;
int npoin;
int nnode;
int nelem;
int *lnode;
double *xyz;

/* aspect ratio for window */
GLfloat aspect_ratio = 1.0;

/* how model is viewed; points or points and bounds */
unsigned int model = 1;
unsigned int key_state[nkey];

/* window width and height */
int winW = 625;
int winH = 625;

/* previous mouse position */
int xold     = 0;
int yold     = 0;
int zold     = 0;

/* state of left/right mouse buttons */
int mouseL_state = UP;
int mouseR_state = UP;

/* initial translation/scaling via mouse */
GLfloat t0[] = {0.0,0.0,0.0};
GLfloat s0[] = {0.9,0.9,0.9};
GLfloat mtran[] = {0.0,0.0,0.0};
GLfloat mscal[] = {0.9,0.9,0.9};

/*---------------*/
/* Main Function */
/*---------------*/
int main(int argc,char *argv[])
{

  int i,j;

  if(argc < 2){
    char fname[100];

    /* read file name */
    printf("Enter Path/File Name: ");
    scanf("%s",fname);

    /* open file */
    fid = fopen(fname,"r");
  }
  else{
    fid = fopen(argv[1],"r");
  }

  /* read data */
  getdata();
  /* close file */
  fclose(fid);

  /* create the quadtree*/
  double* x;
  double* y;

  x = malloc(sizeof(double)*npoin);
  y = malloc(sizeof(double)*npoin);

  /* define x,y point coordinates for quadtree */
  for(i=0;i<npoin;i++){
    x[i] = xyz[ndim*i];
    y[i] = xyz[ndim*i+1];
  }

  /* create quadtree   */
  treeCreate(npoin,x,y);

  /* compare searh methods */
  timeSearch(mytree.root,npoin,x,y);

  /* GLUT management */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(winW,winH);
  glutInitWindowPosition(0,0);
  glutCreateWindow("CSI 703: Assignment #6");

  /* keyboard input function */
  glutKeyboardFunc(keyboard);

  /* drawing function */
  glutDisplayFunc(display);

  /* properly resizing window function */
  glutReshapeFunc(reshape);

  /* motion function for mouse movements */
  glutMotionFunc(motion);

  /* mouse function for translations and rotations */
  glutMouseFunc(mouse);

  /* mouse wheel function for scaling */
  glutMouseWheelFunc(mousewheel);

  /* call the initialization function */
  init();

  /* get in the infinite loop */
  glutMainLoop();

  return 0;
}

/*-------------------------*/
/* Initialization Function */
/*-------------------------*/
void init(void)
{

  /* set background color RGB plus alpha */
  glClearColor(0.1,0.1,0.1,1.0);

  glShadeModel(GL_SMOOTH);
  
  glEnable(GL_DEPTH_TEST);

  state_reset();

  return;
}

/*-------------*/
/* Draw Bounds */
/*-------------*/
void drawBounds(void)
{
  node* inode;

  inode = mytree.root;

  glColor3f(0.0, 0.0, 1.0);
  glBegin(GL_LINES);

  glVertex3f(inode->xmax,inode->ymin,0.0);
  glVertex3f(inode->xmax,inode->ymax,0.0);

  glVertex3f(inode->xmin,inode->ymin,0.0);
  glVertex3f(inode->xmin,inode->ymax,0.0);

  glVertex3f(inode->xmin,inode->ymax,0.0);
  glVertex3f(inode->xmax,inode->ymax,0.0);

  glVertex3f(inode->xmin,inode->ymin,0.0);
  glVertex3f(inode->xmax,inode->ymin,0.0);

  glEnd();
  
  return;
}

/*----------------*/
/* Draw Quadrants */
/*----------------*/
void drawQuadrants(node* v)
{
  int i;

  /* parent */
  if(v->nid < 0){

    /* midpoint     */
    double xm;
    double ym;

    xm = (v->xmax + v->xmin)/2.0;
    ym = (v->ymax + v->ymin)/2.0;

    glVertex3f(xm,v->ymin,0.0);
    glVertex3f(xm,v->ymax,0.0);
    
    glVertex3f(v->xmin,ym,0.0);
    glVertex3f(v->xmax,ym,0.0);

    drawQuadrants(v->qNE);
    drawQuadrants(v->qNW);
    drawQuadrants(v->qSW);
    drawQuadrants(v->qSE);
  }
    
  return;
}

/*---------------------*/
/* Draw Search Results */
/*---------------------*/
void drawSearch(void)
{
  int i;
  double rs = 0.15;
  double phi;

  glPushMatrix();
  glTranslatef(0.34,0.1,0.0);
  glColor3f(0.0,1.0,1.0);
  glBegin(GL_LINE_STRIP);
  for(i=0;i<=2*360;i++){
    phi = (i/2.0)*M_PI/180.0;
    glVertex2f(rs*cos(phi),rs*sin(phi));
  }
  glEnd();
  glPopMatrix();
  
  return;
}

/*---------------*/
/* Draw Elements */
/*---------------*/
void drawElements(void)
{
  int i;
  double  rp = 0.34*0.34 +0.1*0.1;
  double d;

  glPointSize(1.25);

  /* Visualize points/points and quadrants */
  switch(model){

  /* draw points */
  case 0:
    glBegin(GL_POINTS);
    for(i=0;i<npoin;i++){

      node* nget = NULL;
      nodeGet(mytree.root,&nget,xyz[ndim*i],xyz[ndim*i+1]);
      
      glColor3f(1.0,0.0,0.0);
      if(nget->inRadius == 1){
	d = (nget->norm2) + rp - 2*(*(nget->px)*0.34 + *(nget->py)*0.1);
	d = sqrt(d)/0.15;
	glColor3f(1.0 - d,1.0,0.0);		
      }
      if((*(nget->px) == 0.34) && (*(nget->py) == 0.1)){
	glColor3f(1.0,0.0,1.0);		
      }

      /* colors and vertices */
      glVertex2f(*(nget->px),*(nget->py));
    }
    glEnd();

    break;

  /* visualize quadtree*/
  case 1:    
    /* draw points */
    glBegin(GL_POINTS);
    for(i=0;i<npoin;i++){

      node* nget = NULL;
      nodeGet(mytree.root,&nget,xyz[ndim*i],xyz[ndim*i+1]);
      
      glColor3f(1.0,0.0,0.0);
      if(nget->inRadius == 1){
	d = (nget->norm2) + rp - 2*(*(nget->px)*0.34 + *(nget->py)*0.1);
	d = sqrt(d)/0.15;
	glColor3f(1.0 - d,1.0,0.0);	
      }
      if((*(nget->px) == 0.34) && (*(nget->py) == 0.1)){
	glColor3f(1.0,0.0,1.0);		
      }

      /* colors and vertices */
      glVertex2f(*(nget->px),*(nget->py));
    }
    glEnd();

    /* draw quadrants */
    glColor3f(0.0, 0.0, 1.0);
    glBegin(GL_LINES);
    drawQuadrants(mytree.root);
    glEnd();

    /* draw search area */
    drawSearch();

    break;
  }
  return;
}

/*------------------*/
/* Display Function */
/*------------------*/
void display(void)
{

  /* clean all pixels */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* Switch to ModelView */
  glMatrixMode(GL_MODELVIEW);
  /* loads the identity */
  glLoadIdentity();

  /* use mouse to translate */
  glTranslatef(mtran[0],mtran[1],mtran[2]);

  /* use mouse to scale */
  glScalef(mscal[0],mscal[1],mscal[2]);

  drawBounds();
  drawElements();

  glutSwapBuffers();

  return;
}

/*----------*/
/* keyboard */
/*----------*/
void keyboard(unsigned char key,int x,int y){

  int i;
  switch(key){
   /* press 'ESC' or 'q' to exit*/
  case 27:  
  case 'q':
    exit(1);
    break;
  case 'p':
    model=0;
    if(key_state[model]==0){
      printf("==> showing points\n");
      state_reset();
      key_state[model]=1;
    }
    break;
  case 't':
    model=1;
    if(key_state[model]==0){
      printf("==> visualizing quadtree\n");
      state_reset();
      key_state[model]=1;
    }
    break;
  case 'o':
    printf("==> returning to original position\n");
    for(i=0;i<3;i++){
      mtran[i] = t0[i];
      mscal[i] = s0[i];
    }
    break;
  }
  glutPostRedisplay();   
  return;
}

/*---------*/
/* Reshape */
/*---------*/
void reshape(int w,int h){

   /* prevent problems */
  if(h == 0){ 
    h = 1;
  }

  winW = w;
  winH = h;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* reset vieport to window dimensions */
  glViewport(0,0,winW,winH);

  aspect_ratio = winW/winH;

  /* Reset coordinate system */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(0.0,aspect_ratio,0.1,100.0);
     
  glutPostRedisplay();
  return;
}

/*-------------------*/
/* set state to zero */
/*-------------------*/
void state_reset(void){
   int i;
   for(i=0;i<nkey;i++){
     key_state[i]=0;
   }  
   return;
}

/*----------------------------------------------*/
/* motion tiggered when mouse moves into window */
/*----------------------------------------------*/
void motion(int x,int y){

  /* left mouse down, then preform translation */
  if(mouseL_state == DOWN){
    mtran[0] -= (double)(xold - x)/200;
    mtran[1] += (double)(yold - y)/200;
    glutPostRedisplay();
  }

  xold = x;
  yold = y;
  zold = y;
  return;
}

/*---------------------------------------------------*/
/* mouse interactions for translations and rotations */
/*---------------------------------------------------*/
void mouse(int button,int state,int x,int y){

  switch(button){
  case GLUT_LEFT_BUTTON:
    if(state == GLUT_DOWN){
      mouseL_state = DOWN;
      mouseR_state = UP;
      xold = x;
      yold = y;
    }
    break;
  case GLUT_RIGHT_BUTTON:
      mouseL_state = UP;
      mouseR_state = UP;
    break;
  }
  return;
}

/*-------------------------*/
/* mousewheel interactions */
/*-------------------------*/
void mousewheel(int button,int dir,int x,int y) 
{
  /* zoom in */
  if(dir > 0){
    mscal[0] *= (double) 1.1;
    mscal[1] *= (double) 1.1;
    mtran[2] *= (double) 1.1;
    glutPostRedisplay();
  }
  /* zoom out */
  else if(dir < 0){
    mscal[0] *= (double) 0.90;
    mscal[1] *= (double) 0.90;
    mtran[2] *= (double) 0.90;
    glutPostRedisplay();
  }
  return;
}

/*------*/
/* Read */
/*------*/
void getdata(void){

  int i;
  char buffer[240];
  char buffer2[240];
  char buffer3[240];
  char elem_type[240];

  for(i=0;i<3;i++){
    fgets(buffer,240,fid);
  }

  /* read number of points */
  fscanf(fid,"%d",&npoin);

  /* allocate memory for array of point coordinates */
  xyz = (double*) malloc(sizeof(double)*ndim*npoin);

  /* read in point coordinates */
  for(i=0;i<npoin;i++){
    fscanf(fid,"%lf %lf %lf",&xyz[ndim*i],&xyz[ndim*i+1],&xyz[ndim*i+2]);
  }

  /* disregard \n and next line, 'END' */
  fgets(buffer2,240,fid);
  fgets(buffer2,240,fid);

  /* read element type, 'TRIANGLE | TETRA' */
  fgets(elem_type,10,fid);

  /* disregard next line, 'title' */
  fgets(buffer3,240,fid);

  /* read number of elements */
  fscanf(fid,"%d",&nelem);

  /* determine number of nodes for elements */
  if (strcmp(elem_type,"TRIANGLE\n") == 0){
    nnode = 3;
  }
  else{
    nnode = 4;
  }  
  /* allocate memory for array of field values */
  lnode = (int*) malloc(sizeof(int)*nnode*nelem);

  if (nnode == 3){
    for(i=0;i<nelem;i++){
      fscanf(fid,"%d %d %d",&lnode[nnode*i],&lnode[nnode*i+1],&lnode[nnode*i+2]);
    }
  }
  else{
    for(i=0;i<nelem;i++){
      fscanf(fid,"%d %d %d %d",&lnode[nnode*i],&lnode[nnode*i+1],&lnode[nnode*i+2],&lnode[nnode*i+3]);
    }
  }

  return;
}
