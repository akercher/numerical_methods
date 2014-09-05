/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 722: Advancing front grid generation        */
/* Last Updated: 10/22/2011                        */
/*-------------------------------------------------*/
/* Program: csi722_gridVisual.c                    */
/* Description: see readme.txt                     */ 
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
#include "faces.h"

/* Define the Number of Dimensions for Grid (1D,2D,3D,...,ND)*/
#define ndim 2
#define nnofa 2 //nodes per face
#define ngrid 148 //number of grids

#ifndef max
#define max(a,b) ((a)>(b)?(a):(b))
#endif
#ifndef min
#define min(a,b) ((a)<(b)?(a):(b))
#endif

/* correspnd to mouse buttons */
enum{UP = 1,DOWN = 0};

/* Prototypes */
void init(void);
void getdata(void);

void clamp(float []);
void colorMap(GLfloat fin,GLfloat* fcol);
void drawFaces(void);
void display(void);
void reshape(int,int);

void keyboard(unsigned char,int,int);
void motion(int,int);
void mouse(int,int,int,int);

/* Global Variables */
FILE *grid_in;
FILE *face_in;
FILE *pts_in;
char* filesIN[16];
int viewS;
int igrid;
int npoin;
int nnode;
int n;
int nelem;
int nface;
face* faces;
grid* grids;
int *lnode;
float *xyz;			/* point coordinates */
float pi = 3.1415926535898;

/* aspect ratio for window */
GLfloat aspect_ratio = 1.0;

/* how model is viewed */
unsigned int model = 4;

/* window width and height */
int winW = 625;
int winH = 625;

/* previous mouse position */
int xold     = 0;
int yold     = 0;
int zold     = 0;

/* state of left/right/center mouse buttons */
int mouseL_state = UP;
int mouseR_state = UP;

/* initial translation/rotation/scaling via mouse */
GLfloat r0[] = {0.0,0.0,0.0};
GLfloat t0[] = {0.0,0.0,0.0};
GLfloat s0[] = {0.45,0.45,0.45};

GLfloat mrot[] = {0.0,0.0,0.0};
GLfloat mtran[] = {-.97,-.80,0.0};
GLfloat mscal[] = {0.55,0.55,0.55};

/* object to draw cylinder */
GLUquadricObj *qobj;

/*---------------*/
/* Main Function */
/*---------------*/
int main(int argc,char* argv[])
{

  int i;

  igrid = 0;

  for(i=0;i<argc;i++){
    filesIN[i] = argv[i];  
  }

  grid_in = fopen(argv[1],"r");
  /* face_in = fopen(argv[1],"r"); */
  pts_in = fopen(argv[2],"r");

  /* read data */
  getdata();
  /* close file */
  fclose(grid_in);
  fclose(pts_in);

  /* GLUT management */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(winW,winH);
  glutInitWindowPosition(0,0);
  glutCreateWindow("CSI 722: Fall 2011");

  glShadeModel(GL_SMOOTH);  
  glEnable(GL_DEPTH_TEST);

  /* drawing function */
  glutDisplayFunc(display);

  /* properly resizing window function */
  glutReshapeFunc(reshape);

  /* motion function for mouse movements */
  glutMotionFunc(motion);

  /* keyboard input function */
  glutKeyboardFunc(keyboard);

  /* mouse function for translations and rotations */
  glutMouseFunc(mouse);

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

  /* reset vieport to window dimensions */
  glViewport(0,0,winW,winH);

  /* Reset coordinate system */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(0.0,aspect_ratio,1.0,100.0);

  glutPostRedisplay();

  printf("  ==> To cycle forward through grid generation press: w\n");
  printf("  ==> To cycle backward through grid generation press: s\n");
  printf("  ==> To exit press: ESC or q\n");

  return;
}


/*----------*/
/* colorMap */
/*----------*/
void colorMap(GLfloat fin,GLfloat* fcol)
{

  GLfloat dx = 0.05;
  GLfloat ga = (6.0 - 2.0*dx)*fin + dx;

  fcol[0] = max(0,(3-fabs(ga-4)-fabs(ga-5))/2);
  fcol[1] = max(0,(4-fabs(ga-2)-fabs(ga-4))/2);
  fcol[2] = max(0,(3-fabs(ga-1)-fabs(ga-2))/2);

}

/*------------*/
/* Draw Faces */
/*------------*/
void drawFaces(void)
{

  int i,j,k;
  GLfloat* fcol;

  fcol = malloc(sizeof(GLfloat)*nnode);

  /* draw faces */
  nface = grids[igrid].nface;

  for (i=0;i<nface;i++){

    GLfloat gstate = (float) grids[igrid].state[i];
    GLfloat line_width = 1.0 + 2.0*gstate;

    GLfloat c1[3] = {1.0,1.0 - gstate,0.0};      

    glLineWidth(line_width);
    glBegin(GL_LINE_STRIP);
    for(j=0;j<nnofa;j++){
      glColor3fv(c1);
      glVertex2f(xyz[ndim*grids[igrid].fnode[nnofa*i+j]],xyz[ndim*grids[igrid].fnode[nnofa*i+j] + 1]);
    }
    glEnd();
  }

  /* draw points */
  glPointSize(8.0);
  glBegin(GL_POINTS);
  for(i=0;i<nface;i++){
    
    GLfloat c2[3] = {0.0,0.0,1.0};

    glColor3fv(c2);    
    glVertex2f(xyz[ndim*grids[igrid].fnode[nnofa*i+j]],xyz[ndim*grids[igrid].fnode[nnofa*i+j] + 1]);
  }   
  glEnd();

  free(fcol);

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

  /* use mouse to rotate */
  glRotatef(mrot[0],1.0,0.0,0.0);
  glRotatef(mrot[1],0.0,1.0,0.0);

  drawFaces();

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
     exit(1);
     break;
   case 'q':  
     exit(1);
     break;

   /* press spacebar to cycle through grid generation */  
   /* case 32: */
   /*   igrid = igrid + 1; */
   /*   if(igrid > ngrid-1){ */
   /*     igrid = 0; */
   /*   } */     
     /* break; */
   /* cycle forward through grid generation */  
   case 'w':
     igrid = igrid + 1;
     if(igrid > ngrid-1){
       igrid = 0;
     }
     
     break;
   /* cycle backward through grid generation */  
   case 's':
     igrid = igrid - 1;
     if(igrid < 0){
       igrid = ngrid-1;
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
  gluPerspective(0.0,aspect_ratio,1,1000.0);
     
  glutPostRedisplay();
  return;
}

/*----------------------------------------*/
/* clamps range of rotation to [-360,360] */
/*----------------------------------------*/
void clamp(float theta[])
{
   int i;

   for(i=0;i<3;i++){
     if(theta[i] > 360 || theta[i] < -360){
       theta[i] = 0.0f;
     }
   }
   return;
}

/*----------------------------------------------*/
/* motion tiggered when mouse moves into window */
/*----------------------------------------------*/
void motion(int x,int y){

  /* left mouse down, then preform translation */
  if(mouseL_state == DOWN){
    mtran[0] -= (float)(xold - x)/200;
    mtran[1] += (float)(yold - y)/200;
    glutPostRedisplay();
  }
  /* right mouse down, then preform rotation */
  else if(mouseR_state == DOWN){
    mrot[0] -= (float)((yold - y)*180.0)/200.0;
    mrot[1] -= (float)((xold - x)*180.0)/200.0;
    clamp(mrot);
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
    if(state == GLUT_DOWN){
      mouseR_state = DOWN;
      mouseL_state = UP;

      xold = x;
      yold = y;
    }
    break;

  case 3:
    /*zoom in*/  
    mscal[0] *= (float) 1.025;
    mscal[1] *= (float) 1.025;
    mtran[2] *= (float) 1.025;
    glutPostRedisplay();
    break;
  case 4:
    /*zoom out*/  
    mscal[0] *= (float) 0.975;
    mscal[1] *= (float) 0.975;
    mtran[2] *= (float) 0.975;
    glutPostRedisplay();
    break;

  }
  return;
}


/*------*/
/* Read */
/*------*/
void getdata(void){

  int i,j;

  r0[0] = 0.0;
  mrot[0] = r0[0];

  grids = malloc(sizeof(grid)*ngrid);  

  /* read in grids */
  for(i=0;i<ngrid;i++){

    /* read number of faces */
    fscanf(grid_in,"%d",&nface);

    grids[i].nface = nface;

    grids[i].state = (int*) malloc(sizeof(int)*nface);
    grids[i].fnode = (int*) malloc(sizeof(int)*nnofa*nface);
    grids[i].gid = i;

    for(j=0;j<nface;j++){
      fscanf(grid_in," %d %d %d",&grids[i].state[j],&grids[i].fnode[nnofa*j],&grids[i].fnode[nnofa*j+1]);
    }
  }  

  /* /\* read number of lines *\/ */
  /* fscanf(face_in,"%d",&nface); */
  /* faces = malloc(sizeof(face)*nface); */

  /* printf("%d\n",nface); */

  /* /\* read in nodes defining faces *\/ */
  /* for(i=0;i<nface;i++){ */
  /*   /\* read number of faces *\/ */
  /*   faces[i].fid = i; */
  /*   faces[i].state = 0; */
  /*   faces[i].fnode = (int*) malloc(sizeof(int)*nnofa); */
  /*   /\* for(j=0;j<nnofa;j++){ *\/ */
  /*     fscanf(face_in,"%d %d",&faces[i].fnode[0],&faces[i].fnode[1]); */
  /*   /\* } *\/ */
  /* } */

  /* read number of points */
  fscanf(pts_in,"%d",&npoin);

  /* allocate memory for array of point coordinates */
  xyz = (float*) malloc(sizeof(float)*ndim*npoin);
    
  /* read in point coordinates */
  for(i=0;i<npoin;i++){
    fscanf(pts_in,"%f %f",&xyz[ndim*i],&xyz[ndim*i+1]);
  }    

  return;
}
