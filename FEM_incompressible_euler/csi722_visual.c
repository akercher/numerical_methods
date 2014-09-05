/*-------------------------------------------------*/
/* Andrew Kercher                                  */
/* CSI 721: Project 1, Spring 2011                 */
/* Last Updated: 04/01/2011                        */
/*-------------------------------------------------*/
/* Program: csi721_visual.c                        */
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

/* Define the Number of Dimensions for Grid (1D,2D,3D,...,ND)*/
#define ndim 3
#define nkey 5

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
void state_reset(void);
void getdata(void);

void clamp(float []);
void drawAxis(void);
/* void drawBounds(void); */
/* void drawProblem(void); */
/* void printDescription(void); */
void colorMap_rainbow(GLfloat fin,GLfloat* fcol);
void colorMap_jet(GLfloat fin,GLfloat* fcol);
void drawElements(void);
void vector_fld(GLfloat [],GLfloat []);
void display(void);
void reshape(int,int);

void keyboard(unsigned char,int,int);
void motion(int,int);
void mouse(int,int,int,int);
void mousewheel(int,int,int,int);

/* Global Variables */
FILE *fid;
char* filesIN[2];
int npoin;
int nnode;
int n;
int nelem;
int *lnode;
float *xyz;			/* point coordinates */
float *fxyz;			/* scaler field */
float *Vxyz;			/* vector field */
float pi = 3.1415926535898;

/* aspect ratio for window */
GLfloat aspect_ratio = 1.0;

/* how model is viewed */
unsigned int model = 0;
unsigned int key_state[nkey];

/* window width and height */
int winW = 1200;
int winH = 690;

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
GLfloat s0[] = {0.15,0.15,0.15};
GLfloat s1[] = {0.5,0.5,0.5};

GLfloat mrot[] = {0.0,0.0,0.0};
GLfloat mtran[] = {0.0,0.0,0.0};
GLfloat mscal[] = {0.75,0.75,0.75};

/* object to draw cylinder */
GLUquadricObj *qobj;

/*---------------*/
/* Main Function */
/*---------------*/
int main(int argc,char* argv[])
{

  int i;

  for(i=0;i<argc;i++){
    filesIN[i] = argv[i];  
  }

  for(i=0;i<argc-1;i++){
    fid = fopen(argv[i+1],"r");
    /* read data */
    getdata();
    /* close file */
    fclose(fid);
  }

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

  /* mouse wheel function for scaling */
  /* glutMouseWheelFunc(mousewheel); */

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
  /* glClearColor(0.1,0.1,0.1,1.0); /\* dark gray *\/ */
  glClearColor(0.3,0.3,0.3,1.0); /* light gray */

  /* reset vieport to window dimensions */
  glViewport(0,0,winW,winH);

  /* Reset coordinate system */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(0.0,aspect_ratio,1.0,100.0);

  state_reset();

  glutPostRedisplay();

  return;
}

/*-----------*/
/* Draw Axis */
/*-----------*/
void drawAxis(void)
{
   glColor3f(0.5,0.5,0.5);
   glBegin(GL_LINES);
      glColor3f(0.5,0.0,0.0);
      glVertex3f(-50.0,0.0,0.0);
      glVertex3f(50.0,0.0,0.0);

      glColor3f(0.0,0.5,0.0);
      glVertex3f(0.0,50.0,0.0);
      glVertex3f(0.0,-50.0,0.0);

      glColor3f(0.0,0.0,0.5);
      glVertex3f(0.0,0.0,-50.0);
      glVertex3f(0.0,0.0,50.0);
   glEnd ();
   return;
}

/*------------------*/
/* colorMap_rainbow */
/*------------------*/
void colorMap_rainbow(GLfloat fin,GLfloat* fcol)
{

  int i;
  GLfloat xmin,xmax;
  GLfloat dx = 0.10;
  GLfloat ga = (6.0 - 2.0*dx)*fin + dx;

  fcol[0] = max(0,(3-fabs(ga-4)-fabs(ga-5))/2);
  fcol[1] = max(0,(4-fabs(ga-2)-fabs(ga-4))/2);
  fcol[2] = max(0,(3-fabs(ga-1)-fabs(ga-2))/2);

  for(i=0;i<3;i++){
    xmin = max(fcol[i],0.0);
    fcol[i] = min(xmin,1.0);
  }

}

/*--------------*/
/* colorMap_jet */
/*--------------*/
void colorMap_jet(GLfloat fin,GLfloat* fcol)
{

  int i;
  GLfloat xmin,xmax;
  GLfloat four_fin = 4.0*fin;

  fcol[0] = min(four_fin - 0.5,-four_fin + 4.5);
  fcol[1] = min(four_fin - 0.5,-four_fin + 3.5);
  fcol[2] = min(four_fin + 0.5,-four_fin + 2.5);

  for(i=0;i<3;i++){
    xmin = max(fcol[i],0.0);
    fcol[i] = min(xmin,1.0);
  }
}

/*---------------*/
/* Draw Elements */
/*---------------*/
void drawElements(void)
{

  int i,j,k;
  int npl;
  GLfloat fmax;
  GLfloat* fcol;

  fcol = malloc(sizeof(GLfloat)*nnode);

  /* Visualize points/elements and field */
  switch(model){

  /* draw shaded triangles */
  case 0:

    fmax = 0.0;
    for(i=0;i<npoin;i++){
      fmax = max(fmax,fxyz[i]);
    }

    glBegin(GL_TRIANGLES);
    for(i=0;i<nelem;i++){

      /* nodes of element i */
      GLint a = lnode[i*nnode] - 1;
      GLint b = lnode[i*nnode+1] - 1;
      GLint c = lnode[i*nnode+2] - 1;

      /* GLfloat fa = (1.0 + fxyz[a*n])/2; */
      /* GLfloat fb = (1.0 + fxyz[b*n])/2; */
      /* GLfloat fc = (1.0 + fxyz[c*n])/2; */
      GLfloat fa = fxyz[a]/fmax;
      GLfloat fb = fxyz[b]/fmax;      
      GLfloat fc = fxyz[c]/fmax;

      colorMap_rainbow(fa,fcol);
      colorMap_jet(fa,fcol);
      GLfloat c1[3] = {fcol[0],fcol[1],fcol[2]};

      colorMap_rainbow(fb,fcol);
      colorMap_jet(fb,fcol);
      GLfloat c2[3] = {fcol[0],fcol[1],fcol[2]};

      colorMap_rainbow(fc,fcol);
      colorMap_jet(fc,fcol);
      GLfloat c3[3] = {fcol[0],fcol[1],fcol[2]};

      /* define vertices */
      GLfloat v1[3] = {xyz[a*ndim],xyz[a*ndim+1],fxyz[a]/fmax};
      GLfloat v2[3] = {xyz[b*ndim],xyz[b*ndim+1],fxyz[b]/fmax};
      GLfloat v3[3] = {xyz[c*ndim],xyz[c*ndim+1],fxyz[c]/fmax};

      /* colors and vertices */
      glColor3fv(c1);
      glVertex3fv(v1);
      glColor3fv(c2);
      glVertex3fv(v2);
      glColor3fv(c3);
      glVertex3fv(v3);
    }
    glEnd();
    break;
    
  /* draw points */
  case 1:

    fmax = 0.0;
    for(i=0;i<npoin;i++){
      fmax = max(fmax,fxyz[i]);
    }

    glPointSize(3.0);
    glBegin(GL_POINTS);
    for(i=0;i<npoin;i++){

      /* colors and vertices */
      /* GLfloat fa = (1.0 + fxyz[i*n])/2; */
      GLfloat fa = fxyz[i]/fmax;

      /* colorMap_rainbow(fa,fcol); */
      colorMap_jet(fa,fcol);
      GLfloat c1[3] = {fcol[0],fcol[1],fcol[2]};
      /* GLfloat c1[3] = {1,0,0}; */

      glColor3fv(c1);    
      glVertex3f(xyz[ndim*i],xyz[ndim*i+1],fxyz[i]/fmax);

    }   
    glEnd();
    break;

  /* draw wireframe */
  case 2:

    fmax = 0.0;
    for(i=0;i<npoin;i++){
      fmax = max(fmax,fxyz[i]);
    }

    /* for(i=0;i<n*nelem;i++){ */
    for(i=0;i<nelem;i++){

      /* nodes of element i */
      GLint a = lnode[i*nnode] - 1;
      GLint b = lnode[i*nnode+1] - 1;
      GLint c = lnode[i*nnode+2] - 1;

      /* GLfloat fa = (1.0 + fxyz[a])/2; */
      /* GLfloat fb = (1.0 + fxyz[b])/2; */
      /* GLfloat fc = (1.0 + fxyz[c])/2; */
      GLfloat fa = fxyz[a]/fmax;
      GLfloat fb = fxyz[b]/fmax;
      GLfloat fc = fxyz[c]/fmax;

      /* colorMap_rainbow(fa,fcol); */
      colorMap_jet(fa,fcol);
      GLfloat c1[3] = {fcol[0],fcol[1],fcol[2]};

      /* colorMap_rainbow(fb,fcol); */
      colorMap_jet(fb,fcol);
      GLfloat c2[3] = {fcol[0],fcol[1],fcol[2]};

      /* colorMap_rainbow(fc,fcol); */
      colorMap_jet(fc,fcol);
      GLfloat c3[3] = {fcol[0],fcol[1],fcol[2]};

      /* GLfloat c1[3] = {1.0,0.0,0.0}; */
      /* GLfloat c2[3] = {1.0,0.0,0.0}; */
      /* GLfloat c3[3] = {1.0,0.0,0.0}; */

      /* define vertices */
      GLfloat v1[3] = {xyz[a*ndim],xyz[a*ndim+1],fxyz[a]/fmax};
      GLfloat v2[3] = {xyz[b*ndim],xyz[b*ndim+1],fxyz[b]/fmax};
      GLfloat v3[3] = {xyz[c*ndim],xyz[c*ndim+1],fxyz[c]/fmax};

      /* colors and vertices */      
      glBegin(GL_LINE_STRIP);
      glColor3fv(c1);
      glVertex3fv(v1);
      glColor3fv(c2);
      glVertex3fv(v2);
      glColor3fv(c3);
      glVertex3fv(v3);
      glColor3fv(c1);
      glVertex3fv(v1);
      glEnd();
    }
    break;
    
  /* draw vector field */
  case 3:    
    for(i=0;i<npoin;i++){
      
      GLfloat x0[3] = {0.0,0.0,0.0};
      GLfloat x1[3] = {0.0,0.0,0.0};
      
      GLfloat fmag = sqrt((Vxyz[i*n]*Vxyz[i*n])+(Vxyz[i*n+1]*Vxyz[i*n+1]));
      /* GLfloat c1[3] = {fmag,0.0,1-fmag}; */

      GLfloat fa = fmag;
      colorMap_rainbow(fa,fcol);
      /* colorMap_jet(fa,fcol); */
      GLfloat c1[3] = {fcol[0],fcol[1],fcol[2]};

      for(j=0;j<3;j++){
      	  x0[j] = xyz[ndim*i + j] - (0.05)*Vxyz[n*i + j]/(2.0*fmag);
      	  x1[j] = xyz[ndim*i + j] + (0.05)*Vxyz[n*i + j]/(2.0*fmag);
      }

      /* colors and vertices */
      glLineWidth(1.5);
      glColor3fv(c1);
      vector_fld(x0,x1);
    }
    break;
  }
  free(fcol);
  return;
}

/*-----------------*/
/* Vector Function */
/*-----------------*/
void vector_fld(GLfloat xi[],GLfloat xf[]){

  int i,axis;
  float r,d;
  GLfloat x,y,z;
  GLfloat u[3],v[3],w[3];
  static GLfloat xaxis[3] = {1.0,0.0,0.0};
  static GLfloat yaxis[3] = {0.0,1.0,0.0};
  static GLfloat zaxis[3] = {0.0,0.0,1.0};

  /* define length of vector components */
  w[0] = xf[0] - xi[0];
  w[1] = xf[1] - xi[1];
  w[2] = xf[2] - xi[2];
  
  /* x-axis */
  axis = 1;
  float mag = fabs(w[0]);
  if(fabs(w[1])  > mag){
    /* y-axis */
    axis = 2;
    mag = fabs(w[1]);
  }
  if(fabs(w[2])  > mag ){
    /* z-axis */
    axis = 3;
    mag = fabs(w[2]);
  }
  
  r = sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
  if(r > 0.0){
    w[0] = w[0]/r;
    w[1] = w[1]/r;
    w[2] = w[2]/r;
  }

  /* define length of sides of the arrow */
  d = (0.20)*r;  

  /* draw main part of vector */
  glBegin(GL_LINE_STRIP);
  glVertex3fv(xi);
  glVertex3fv(xf);
  glEnd();

  if(axis != 1){  
    /* cross seperation vector with x-axis */
    v[0] = w[1]*xaxis[2] - xaxis[1]*w[2];
    v[1] = w[2]*xaxis[0] - xaxis[2]*w[0];
    v[2] = w[0]*xaxis[1] - xaxis[0]*w[1];

    /* make a unit vector */
    r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(r > 0.0){
      v[0] = v[0]/r;
      v[1] = v[1]/r;
      v[2] = v[2]/r;
    }

    u[0] = v[1]*w[2] - w[1]*v[2];
    u[1] = v[2]*w[0] - w[2]*v[0];
    u[2] = v[0]*w[1] - w[0]*v[1];

    x = xf[0] + d*(u[0] - w[0]);
    y = xf[1] + d*(u[1] - w[1]);
    z = xf[2] + d*(u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();

    x = xf[0] + d*(-u[0] - w[0]);
    y = xf[1] + d*(-u[1] - w[1]);
    z = xf[2] + d*(-u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();      
  }

  if(axis != 2){
    /* cross seperation vector with y-axis */
    v[0] = w[1]*yaxis[2] - yaxis[1]*w[2];
    v[1] = w[2]*yaxis[0] - yaxis[2]*w[0];
    v[2] = w[0]*yaxis[1] - yaxis[0]*w[1];

    /* make a unit vector */
    r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(r > 0.0){
      v[0] = v[0]/r;
      v[1] = v[1]/r;
      v[2] = v[2]/r;
    }

    u[0] = v[1]*w[2] - w[1]*v[2];
    u[1] = v[2]*w[0] - w[2]*v[0];
    u[2] = v[0]*w[1] - w[0]*v[1];

    x = xf[0] + d*(u[0] - w[0]);
    y = xf[1] + d*(u[1] - w[1]);
    z = xf[2] + d*(u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();

    x = xf[0] + d*(-u[0] - w[0]);
    y = xf[1] + d*(-u[1] - w[1]);
    z = xf[2] + d*(-u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();      
  }

  if(axis != 3){
    /* cross seperation vector with z-axis */
    v[0] = w[1]*zaxis[2] - zaxis[1]*w[2];
    v[1] = w[2]*zaxis[0] - zaxis[2]*w[0];
    v[2] = w[0]*zaxis[1] - zaxis[0]*w[1];

    /* make a unit vector */
    r = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    if(r > 0.0){
      v[0] = v[0]/r;
      v[1] = v[1]/r;
      v[2] = v[2]/r;
    }

    u[0] = v[1]*w[2] - w[1]*v[2];
    u[1] = v[2]*w[0] - w[2]*v[0];
    u[2] = v[0]*w[1] - w[0]*v[1];

    x = xf[0] + d*(u[0] - w[0]);
    y = xf[1] + d*(u[1] - w[1]);
    z = xf[2] + d*(u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();

    x = xf[0] + d*(-u[0] - w[0]);
    y = xf[1] + d*(-u[1] - w[1]);
    z = xf[2] + d*(-u[2] - w[2]);

    /* draw sides of arrow */
    glBegin(GL_LINE_STRIP);
    glVertex3fv(xf);
    glVertex3f(x,y,z);
    glEnd();      
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

  /* use mouse to rotate */
  glRotatef(mrot[0],1.0,0.0,0.0);
  glRotatef(mrot[1],0.0,1.0,0.0);

/*   drawAxis(); */
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
   case 'o':
     for(i=0;i<3;i++){
       mrot[i] = r0[i];
       mtran[i] = t0[i];
       mscal[i] = s0[i];
     }
     break;
   case 't':
     /* n = 1; */
     model=0;
     if(key_state[model] == 0){
       printf("==> showing triangular elements\n");
       state_reset();
       key_state[model]=1;
     }
     break;
   case 'p':
     /* n = 1; */
     model=1;
     if(key_state[model] == 0){
       printf("==> showing points\n");
       state_reset();
       key_state[model]=1;
     }
     break;
   case 'w':
     /* n = 1; */
     model=2;
     if(key_state[model] == 0){
       printf("==> showing triangular wireframe\n");
       state_reset();
       key_state[model]=1;
     }
     break;
   case 'v':
     n = 3;
     model=3;
     if(key_state[model] == 0){
       printf("==> showing vector field\n");
       state_reset();
       key_state[model]=1;
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

/*-------------------------*/
/* mousewheel interactions */
/*-------------------------*/
void mousewheel(int button,int dir,int x,int y) 
{
  
  /* zoom in */
  if(dir > 0){
    mscal[0] *= (float) 1.025;
    mscal[1] *= (float) 1.025;
    mtran[2] *= (float) 1.025;
    glutPostRedisplay();
  }
  /* zoom out */
  else if(dir < 0){
    mscal[0] *= (float) 0.975;
    mscal[1] *= (float) 0.975;
    mtran[2] *= (float) 0.975;
    glutPostRedisplay();
  }

  mouseL_state = UP;
  mouseR_state = UP;

  return;
}

/*------*/
/* Read */
/*------*/
void getdata(void){

  int i,j;
  int npl;
  char buffer[240];
  char domainType[240];
  char domainTitle[240];
  char keyword[240];
  char fldType[240];
  char fldTitle[240];
  char elemType[240];
  char elemTitle[240];

  fgets(domainType,240,fid);
  fgets(domainTitle,240,fid);
  fgets(keyword,240,fid);

  /* read number of points */
  fscanf(fid,"%d",&npoin);

  /* allocate memory for array of point coordinates */
  xyz = (float*) malloc(sizeof(float)*ndim*npoin);

  /* read in point coordinates */
  for(i=0;i<npoin;i++){
    fscanf(fid,"%f %f %f",&xyz[ndim*i],&xyz[ndim*i+1],&xyz[ndim*i+2]);
  }
  
  /* disregard '\n', newline */
  fgets(buffer,240,fid);
  
  /* read field type, 'REAL | INTEGER' */
  fgets(fldType,9,fid);
  
  /* title of field*/
  fgets(fldTitle,240,fid);
  
  /* read number of dimensions for field */
  fscanf(fid,"%d",&n);
  
  /* scaler field */
  if(n == 1){
    /* allocate memory for array of field values */
    fxyz = (float*) malloc(sizeof(float)*n*npoin);
    for(i=0;i<n*npoin;i++){
      fscanf(fid,"%f",&fxyz[n*i]);
    }
  }
  /* vector field */
  else if(n == 3){

    /* allocate memory for array of field values */
    Vxyz = (float*) malloc(sizeof(float)*n*npoin);
    for(i=0;i<n*npoin;i++){
      fscanf(fid,"%f %f %f",&Vxyz[n*i],&Vxyz[n*i+1],&Vxyz[n*i+2]);
    }
  }
  
  /* read next keyword */
  fgets(keyword,13,fid);
  
  if(strcmp(keyword," END_FIELDS\n") != 0){
    fgets(buffer,240,fid);
  }
  
  /* read element type, 'TRIANGLE | STRUCTURED' */
  fgets(elemType,11,fid);
  
  /* title of elements */
  fgets(elemTitle,240,fid);
  
  /* determine number of nodes for elements */
  if (strcmp(elemType," TRIANGLE\n") == 0){
    nnode = 3;
    
    /* read number of elements */
    fscanf(fid,"%d",&nelem);
  }
  else{
    nnode = 4;
    nelem = 1;
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
