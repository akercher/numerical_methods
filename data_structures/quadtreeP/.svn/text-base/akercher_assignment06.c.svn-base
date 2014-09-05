/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 703: Assignment #6                          */
/* 03/30/2010                                      */
/*-------------------------------------------------*/
/* Program: akercher_assignment06.c                */
/* Description: The program reads a .zfem file     */ 
/*   containing 3D point coordinates  and element  */
/*   connectivity.  Points within a radius of 0.15 */
/*   of the point (0.34,0.1) are shown in yellow.  */
/*   Options, and corresponding key are the        */
/*   following:                                    */
/*     1. Points ==> press 'p'.                    */
/*     2. Points and Quadtree ==> press 't'.       */
/*     3. Points ==> press 'p'                     */
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

#include "quadtreeP.h"

/* Define the Number of Dimensions for Grid (1D,2D,3D,...,ND)*/
#define ndim 3
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
double rmax;

double x0,x1;
double z0,z1;

/* aspect ratio for window */
GLfloat aspect_ratio = 1.0;

/* how model is viewed; points, filled triangles, or wireframe */
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
GLfloat r0[] = {0.0,0.0,0.0};
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
  double rt = 0.15;
  double* x;
  double* y;

  x = malloc(sizeof(double)*npoin);
  y = malloc(sizeof(double)*npoin);

  /* define x,y point coordinates for quadtree */
  for(i=0;i<npoin;i++){
    x[i] = xyz[ndim*i];
    y[i] = xyz[ndim*i+1];
  }

  treeCreate(npoin,x,y,rt);

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
  int i;
  double phi;

  rmax = nodeBounds(mytree.root,0.0);
  rmax = sqrt(rmax);
  rmax = rmax + rmax/50;

  glColor3f(1.0, 0.0, 0.0);
  glBegin(GL_LINE_STRIP);
  for(i=0;i<=360;i++){
    phi = i*M_PI/180.0;
    glVertex2f(rmax*cos(phi),rmax*sin(phi));
  }
  glEnd();
  
  return;
}

/*---------------*/
/* Draw Elements */
/*---------------*/
void drawElements(void)
{
  int i;
  int cnt;
  double phi,theta,w,u;

  glPointSize(2.0);

  /* Visualize points/elements and field */
  switch(model){

  /* draw points */
  case 0:
    glBegin(GL_POINTS);
    for(i=0;i<npoin;i++){

      node* nget = NULL;
      nodeGet(mytree.root,&nget,xyz[ndim*i],xyz[ndim*i+1]);
      
      glColor3f(0.0,0.0,1.0);
      if(nget->inRadius == 1){
	glColor3f(1.0,1.0,0.0);	
      }

      /* colors and vertices */
      glVertex2f(nget->x,nget->y);
    }
    glEnd();

    break;

  /* visualize quadtree*/
  case 1:    
    glBegin(GL_POINTS);
    for(i=0;i<npoin;i++){

      node* nget = NULL;
      nodeGet(mytree.root,&nget,xyz[ndim*i],xyz[ndim*i+1]);
      
      glColor3f(0.0,0.0,1.0);
      if(nget->inRadius == 1){
	glColor3f(1.0,1.0,0.0);	
      }

      /* colors and vertices */
      glVertex2f(nget->x,nget->y);
    }
    glEnd();


    glColor3f(1.0, 0.0, 0.0);
    glBegin(GL_LINES);
    for(i=0;i<npoin;i++){

      node* nget = NULL;
      nodeGet(mytree.root,&nget,xyz[ndim*i],xyz[ndim*i+1]);

      cnt = 0;
      cnt = nodeNchild(nget);
      if(cnt > 0){
	theta = acos(nget->x/rmax);
	phi = asin(nget->y/rmax);

	u = rmax*cos(phi);
	w = rmax*sin(theta);

	x0 = min(-u,u);
	x1 = max(-u,u);
	z0 = min(-w,w);
	z1 = max(-w,w);

	nodeSubBounds(nget,xyz[ndim*i],xyz[ndim*i+1]);

	glColor3f(1.0,0.0,0.0);
	glVertex3f(x0,nget->y,0.0);
	glVertex3f(x1,nget->y,0.0);

	glColor3f(1.0,0.0,0.0);
	glVertex3f(nget->x,z0,0.0);
	glVertex3f(nget->x,z1,0.0);
      }
    }
    glEnd ();

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
