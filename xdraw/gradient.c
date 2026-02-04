/******************************************************************************
**  NAME	GRADIENT.C
**  AUTHOR	Sheryl M. Glasser
**
**  DESCRIPTION
**     Supporting functions for xwindows
**
**  Copyright (c) CounterPoint Graphics 1993.  All rights reserved.
******************************************************************************/
/*
Last  correction of bugs: 30.11.2000 : aspect+zoom, zoom 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*#include <fcntl.h>*/

#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#define zprintf if (zdebug) printf

extern int dialogwindow;
extern Window dialog_win;

FILE *gfile = NULL;
int tell_grad = 0;
int finding_contour_value = 0;

extern int debug_m;
extern float *buf;
//extern FloatSc xscale, xoffset, yscale, yoffset;
extern int wx0,wy0;


#include "gendefs.h"
#include "curves.h"		/* (xcontour.h uses CURVE_SET) */
#include "xtools.h"
#include "xdraw.h"
#include "xcontour.h"
#include "xinit.h"
#include "xedit.h"
#include "setcolor.h"

extern void postscript(char *);	/* near heap problem, try elim .h's*/
extern LOOP loop[];
extern int nnode_r;

#ifdef USE_MENU
#include "device.h"
#include "menuwin.h"
#include "menu.h"		/* for set_selected_iwin */
#include "glx.h"
static int zdebug=1;
#else
static int zdebug=0;
#endif				/* ...of USE_MENU */

/*=============================================================================
**            DEFINITIONS AND VARIABLES
**===========================================================================*/
extern Display *mydisplay;
extern int myscreen;
extern Colormap cmap;
extern Window root;
extern Font font;
extern XFontStruct *font_struct;

extern int nwindow;
extern VIEW view[];
extern int font_height;
extern unsigned long myforeground, mybackground;
extern CURVE_SET curveset[];
extern int exitflag, redrawflag, titledrawn;
extern unsigned int xscreen, yscreen;
extern int ftype;
extern int ntime;
extern float data_limits;
extern void get_world_coordinates(int, int, double *, double *, char *);


extern Window zoomwin;

typedef struct {
  double dist;
  float x, y, psi;
  int   ir, iz;
  long section;
} NEAREST;

static NEAREST nearList[50];

static int n_near;
static long n_segment, n_segment_z;
static double min_dist, max_dist;
static float psi_of_min, psi_ignore, psi_calc;
static int min_ixy, min_iblock;

extern int coord_on;
extern int n_nearest_seg;
extern NEAR_SEG seg[];

float xworld, yworld;

extern float psi_on_seg;
extern float x_on_seg, y_on_seg;
extern float r_on_seg, z_on_seg;

/*-----------------------------------------------------------------------------
|	get_world
-----------------------------------------------------------------------------*/
float get_world(int xcurs, int ycurs, XRectangle clipbox, 
		float xmin, float ymin, float xmax, float ymax,
		float *px, float *py)
{
  float dx, dy, result;

  dx = (float)(xcurs - clipbox.x) / (float)clipbox.width;
  dy = (float)(ycurs - clipbox.y) / (float)clipbox.height;

  *px = xmin + dx * (xmax - xmin);
  *py = ymin + dy * (ymax - ymin);
}

/*-----------------------------------------------------------------------------
|	matrix_tell
-----------------------------------------------------------------------------*/
void matrix_tell(FILE *fd, double *cMat, int nrow, int ncol, char *caption)
{
  int irow, icol, i;
  if (nrow > 1 ) fprintf(fd, "Matrix %s\n", caption);
  else fprintf(fd, "Vector %s\n", caption);
  for(irow=0; irow<nrow; irow++) {
    for (icol=0; icol<ncol; icol++) {
      i = irow*ncol + icol;
      fprintf(fd, "%11.4lg ", *(cMat+i));
      //fprintf(fd, "%6.2lg ", *(cMat+i));
    }
    fprintf(fd, "\n");
  }
}

/*-----------------------------------------------------------------------------
|	matrix_vector_mult
-----------------------------------------------------------------------------*/
double matrix_vector_mult(double *xyMat, double *c, int i, int n,
			  FILE *fd,  double *df)
{
  double dpsi;
  int j;
  j = n * i;
  dpsi = xyMat[j+0]*c[0] + xyMat[j+1]*c[1] + 
	 xyMat[j+2]*c[2] + xyMat[j+3]*c[3] + xyMat[j+4]*c[4];
  if (tell_grad) fprintf(fd, "Check %d: %lg vs %lg\n", i, dpsi, df[i]);
  return dpsi;
}

#ifdef NOT_YET
/*-----------------------------------------------------------------------------
|	get_angle
-----------------------------------------------------------------------------*/
double get_angle(ir0, iz0, dr, dz, nrz, xoff, yoff)
{
  double dx, dy;
  dx0 = *(buf + xoff + ir0);
  dy0 = *(buf + yoff + ir0 + nrz);
  dx1 = *(buf + xoff + ir0 + dr);
  dy1 = *(buf
}

/*-----------------------------------------------------------------------------
|	point_is_between
|	* theta's are in radians 0 to 2*pi (atan2 gives -pi to pi)
-----------------------------------------------------------------------------*/
void axis_point_is_between(float dy, float dx, float ir0, float iz0)
{
  theta1 = get_angle(ir0, iz0,  1,  0);
  theta2 = get_angle(ir0, iz0,  0,  1, theta1);
  theta3 = get_angle(ir0, iz0, -1,  0, theta1);
  theta4 = get_angle(ir0, iz0,  0, -1, theta1);
  theta = atan2(dy, dx) - theta1;
  theta = angle_range(theta);

  if      (theta>=theta1 && theta<theta2) { *ir5 =  1; *iz5 =  1; }
  else if (theta>=theta2 && theta<theta3) { *ir5 = -1; *iz5 =  1; }
  else if (theta>=theta3 && theta<theta4) { *ir5 = -1; *iz5 = -1; }
  else 					  { *ir5 =  1; *iz5 = -1; }
}


/*-----------------------------------------------------------------------------
|	interpolate_x
|	* (x-x1) = b(r-r1) + a(r-r1)**2
|	* (x2-x1) =  b + a
|	  (x0-x1) = -b + a
|	* dx/dr = b + 2*a*(r-r1)
|	* x1, r1 are the 2nd of 3 points
|	* we want the value of dx/dr at dr (offset from center)
-----------------------------------------------------------------------------*/
float interpolate_x(long xoffi, int di, int ir0, int ir2, int nr, 
		   char *caption)
{
  double a, b, r1, r2, c, ra, rb;
  float x0, x1, x2, dr, *p, xresult;
  long i0;

  if (ir0<0)
    { i0 = 0; dr = -1.; }
  else if (ir2>=nr)
    { i0 = -2*di; dr = 1.; }
  else
    { i0 = -di; dr = 0; }

  //printf("interpolate_x %s: off=%ld, ir0=%d, ir2=%d, nr=%d, i0=%d di=%d\n",
  // caption, xoffi, ir0, ir2, nr, i0, di);

  p = buf + xoffi + i0; x0 = *p;	/* problem when eg x = *(buf+xoffi+i0)?*/
  p += di; x1 = *p;
  p += di; x2 = *p;

  //printf("interpolate_x %s: i0=%d, xoffi=%ld, values %g %g %g\n",
  // caption, i0, xoffi, x0, x1, x2);

  a = (double)((x2 - x0) / 2.);
  b = (double)( x2 -x1) - a;
  //printf("   a=%lg, b=%lg\n", a, b);

  if (a==0.0 && b==0.0) 	/* x0 = x1 = x2 */
    xresult = 0;

  else if (a==0.0) 		/* dx = b * dr */
    xresult = (float)b;

  else {
    xresult = (float)(b + (2.0 * a * dr));
  }

  printf("%s: x's %g %g %g, a,b = %lg %lg, --> %g\n", caption, x0, x1, x2, a, b, xresult);
  return xresult;
}

/*-----------------------------------------------------------------------------
|	test_this_block
|	* nearest_M_value calls redraw_mb calls here
|	* nr, nz is # nodes, e.g. 101
-----------------------------------------------------------------------------*/
void test_this_block(int iblock, int nr, int nz, 
		     long xoff, long yoff, long zoff)
{
  long irz, nrz, ir, iz, min_irz, ir1, iz1;
  int ir0, ir2, iz0, iz2;
  float x, y, x1, y1,  z;
  double dx, dy, dist;
  float dxdr, dydr, dxdz, dydz;
  float drdx, drdy, dzdx, dzdy, det;
  float dfdr, dfdz, dfdx, dfdy, fxy;

  nrz = nr * nz;
  min_iblock = iblock;

  /*------ Find nearest point in the block */

  for (irz=0; irz<nrz; irz++) {
    x = *(buf + xoff + irz);
    y = *(buf + yoff + irz);
    z = *(buf + zoff + irz);
    dx = (double)(x - xworld);
    dy = (double)(y - yworld);
    dist = sqrt(dx*dx + dy*dy);
    if (irz==0 || dist<min_dist) {
      min_dist = dist;
      min_irz = irz; 
      x1 = x;
      y1 = y;
      psi_of_min = z;
    }
  }

  irz = min_irz;
  iz1 = irz / nr;
  ir1 = irz - (iz1*nr);
  printf("Min d=%lg, value %g, in block %d at (%ld,%ld)\n", 
	 min_dist, psi_of_min, iblock, ir1, iz1, nr, nz);

  /*------ Find nearest neighbors in ir, iz space */

  ir0 = ir1 - 1; ir2 = ir1 + 1;
  iz0 = iz1 - 1; iz2 = iz1 + 1;

  /*------ Find dx/dr, dy/dr, dx/dz, dy/dz */

  xoff += irz;
  yoff += irz;
  zoff += irz;

  dxdr = interpolate_x(xoff, 1,  ir0, ir2, nr, "dx/dr");
  dydr = interpolate_x(yoff, 1,  ir0, ir2, nr, "dy/dr");
  dxdz = interpolate_x(xoff, nr, iz0, iz2, nz, "dx/dz");
  dydz = interpolate_x(yoff, nr, iz0, iz2, nz, "dy/dz");

  /*------ Find dr/dx, dr/dy, dz/dx, dz/dy */

  det = dxdr * dydz - dxdz * dydr;
  drdx =  dydr / det;
  drdy = -dxdz / det;
  dzdx = -dydr / det;
  dzdy =  dxdr / det;

  /*------ Find df/dr, df/dz */

  dfdr = interpolate_x(zoff, 1,  ir0, ir2, nr, "df/dr");
  dfdz = interpolate_x(zoff, nr, iz0, iz2, nz, "df/dz");

  /*------ Find df/dx, df/dy */

  dfdx = dfdr*drdx + dfdz*dzdx;
  dfdy = dfdr*drdy + dfdz*dzdy;

  /*------ Find value */

  printf("dx=%g, dy=%g\n", xworld-x1, yworld - y1);
  psi_calc = psi_of_min + dfdx * (xworld - x1) + dfdy * (yworld - y1);
}

#endif

/*-----------------------------------------------------------------------------
|	add_this_point
|	* test if current point may be added to nearList[]
|	* qn is next available in nearList
-----------------------------------------------------------------------------*/
int add_this_point(int npoint, double dist, 
		   float x, float y, float z,
		   int irz, int mr, int nmax)
{
  NEAREST *p, *pn, *pn2, *qn;

  if (npoint == 0) {
    p = nearList; npoint++; goto ATP_INSERT; 
  }
  
  for(p=nearList, qn=p+npoint; p<qn; p++) {
    if (x==p->x && y==p->y) return npoint;	/* e.g. "center" of polar */
  }

  /*------ Test if Insert (one may drop out) */

  for(p=nearList; p<qn; p++) {
    if (dist <= p->dist) goto ATP_OPEN;
  }

  /*------ This value is greater than anything in list so far */
  
  if (npoint == nmax) return npoint;
  p = qn;
  npoint++;
  goto ATP_INSERT;

  /*------ open list if needed */

 ATP_OPEN:
  if (npoint<nmax) { pn2 = nearList + npoint; npoint++; }
  else pn2 = nearList + npoint - 1;
  for(pn = pn2; pn > p; pn--) 
    *pn = *(pn-1);

  /*------ insert here */

 ATP_INSERT:

  //printf("add at %ld; npoints=%d\n", p - nearList, npoint);
  p->dist = dist;
  p->psi = z;
  p->x = x;
  p->y = y;
  p->ir = irz % mr;;
  p->iz = irz / mr;
  return npoint;
}

/*----------------------------------------------------------------------------------
|	determinant
----------------------------------------------------------------------------------*/
double determinant(double *aMat, int n, FILE *fd)
{
  double bMat[25], det, det1, det2, *p;
  int hcol, irow, icol, ij;

  fprintf(fd, "Det %d", n);
  if (n==2) {
    det1 = aMat[0] * aMat[3];
    det2 = aMat[1] * aMat[2];
    det = det1 - det2;
    fprintf(fd, ": %lg %lg %lg %lg = %lg\n", 
	    aMat[0], aMat[1], aMat[2], aMat[3], det);
    return det;
  }

  fprintf(fd, "\n");
  det = 0.0;
  for(hcol=0; hcol<n; hcol++) {
    for(irow=1, p=bMat, ij=n; irow<n; irow++) {
      for(icol=0; icol<n; icol++,ij++) {
	if (icol != hcol) *p++ = aMat[ij];
      }
    }
    fprintf(fd, "Det %d, %d\n", n, hcol);
    det1 = determinant(bMat, n-1, fd);
    fprintf(fd, "Det %d, %d: %g * %g\n", n, hcol, det1, aMat[hcol]);
    det += (det1 * aMat[hcol]);
  }
  return det;
}

/*----------------------------------------------------------------------------------
|	myLU
|	* LU decomposition
----------------------------------------------------------------------------------*/
int myLU(double *M, double *a, double *b)
{
  int i, ij, indx;

   for(i=0; i<25;i++) {
     a[i]=b[i]=0.0;
     if (i==0 || i==6 || i==12 || i==18 || i==24) a[i] = 1.0;
   }

  for(i=0; i<5; i++) b[i]    = M[i];
  indx = 0;
  if (b[0] == 0.0) goto LU_ABORT;

  for(i=1, ij=i*5; i<5; i++, ij+=5)	a[ij]   =  M[ij] / b[0];

  for(i=1; i<5; i++)			b[5+i]  = M[5+i]   - a[5]  * b[i];
  indx = 1;
  if (b[6] == 0.0) goto LU_ABORT;
  
  for(i=2, ij=5*i; i<5; i++, ij+=5)	a[ij+1] = (M[ij+1] - a[ij] * b[1]) / b[6];

  for(i=2; i<5; i++) 			b[10+i] = M[10+i]  - a[10] * b[i] - a[11] * b[5+i];
  indx = 2;
  if (b[10+2] == 0.0) goto LU_ABORT;

  for(i=3, ij=5*i; i<5; i++, ij+=5)	a[ij+2] = (M[ij+2] - a[ij] * b[2] - a[ij+1] * b[5+2]) / b[10+2];

  for(i=3; i<5; i++)			b[15+i] = M[15+i]  - a[15] * b[i] - a[16] * b[5+i] - a[17] * b[10+i];
  indx = 3;
  if (b[15+3] == 0.0) goto LU_ABORT;

  for(i=4, ij=5*i; i<5; i++, ij+=5)	a[ij+3] = (M[ij+3] - a[ij] * b[3] - a[ij+1] * b[5+3] - a[ij+2] * b[10+3]) / b[15+3];
  for(i=4; i<5; i++)			b[20+i] = M[20+i]  - a[20] * b[i] - a[21] * b[5+i] - a[22] * b[10+i] - a[23] * b[15+i];

  return 0;

 LU_ABORT:
   printf("Matrix cannot be decomposed (%d)\n", indx);
   return 1;
}

/*----------------------------------------------------------------------------------
|	LUsimple
|	* LU decomposition for 5 x 5 matrix
----------------------------------------------------------------------------------*/
int LUsimple(double *Mat, double *a, double *b)
{
   int i, indx;
   for(i=0; i<25;i++) {
     a[i]=b[i]=0.0;
     if (i==0 || i==6 || i==12 || i==18 || i==24) a[i] = 1.0;
   }

   indx = 0;
   *(b+0) = *(Mat+0);
   *(b+1) = *(Mat+1);
   *(b+2) = *(Mat+2);
   *(b+3) = *(Mat+3);
   *(b+4) = *(Mat+4);
   if (*(b+0)==0.0) goto ABORT;

   *(a+5)  =  *(Mat+5)  / *(b+0);
   *(a+10) =  *(Mat+10) / *(b+0);
   *(a+15) =  *(Mat+15) / *(b+0);
   *(a+20) =  *(Mat+20) / *(b+0);

   indx = 6;
   *(b+6)  =  *(Mat+6) - (*(a+5) * *(b+1));
   *(b+7)  =  *(Mat+7) - (*(a+5) * *(b+2));
   *(b+8)  =  *(Mat+8) - (*(a+5) * *(b+3));
   *(b+9)  =  *(Mat+9) - (*(a+5) * *(b+4));
   if (*(b+6)==0.0) goto ABORT;

   *(a+11) = (*(Mat+11) - *(a+10) * *(b+1)) / *(b+6);
   *(a+16) = (*(Mat+16) - *(a+15) * *(b+1)) / *(b+6);
   *(a+21) = (*(Mat+21) - *(a+20) * *(b+1)) / *(b+6);

   indx = 12;
   *(b+12) =  *(Mat+12) - *(a+10) * *(b+2) - *(a+11) * *(b+7);
   *(b+13) =  *(Mat+13) - *(a+10) * *(b+3) - *(a+11) * *(b+8);
   *(b+14) =  *(Mat+14) - *(a+10) * *(b+4) - *(a+11) * *(b+9);
   if (*(b+12)==0.0) goto ABORT;

   /* ai2  =   (mi2       - ai0 * b02        - ai1 * b12)       / b22, i=3,4 */
   *(a+17) = (*(Mat+17) - *(a+15) * *(b+2) - *(a+16) * *(b+7)) / *(b+12);
   *(a+22) = (*(Mat+22) - *(a+20) * *(b+2) - *(a+21) * *(b+7)) / *(b+12);

   /* b3i       m3i         a30 * b0i        - a31 * b1i      - a32 * b2i, i=3,4 */
   indx = 18;
   *(b+18) =  *(Mat+18) - *(a+15) * *(b+3) - *(a+16) * *(b+8) - *(a+17) * *(b+13);
   *(b+19) =  *(Mat+19) - *(a+15) * *(b+4) - *(a+16) * *(b+9) - *(a+17) * *(b+14);
   if (*(b+18)==0.0) goto ABORT;

   /* ai3  =   (mi3       - ai0 * b03        - ai1 * b13        - ai2 * b23)       / b33, i=4 */
   *(a+23) = (*(Mat+23) - *(a+20) * *(b+3) - *(a+21) * *(b+8) - *(a+22) * *(b+13)) / *(b+18);
   
   /* b4i       m4i         a40 * b0i        - a41 * b1i      - a42 * b2i         - a43 * b3i, i=4 */
   *(b+24) =  *(Mat+24) - *(a+20) * *(b+4) - *(a+21) * *(b+9) - *(a+22) * *(b+14) - *(a+23) * *(b+19);
   return 0;

ABORT:
   printf("Matrix cannot be decomposed (%d)\n", indx);
   return 1;
}

/*-----------------------------------------------------------------------------
|	nearest_M_value_xx
|	* xwin, ywin are world coordinates
|	* return value is psi at nearest point, psixy is calc value
|	* method: f = a*dx + b*dy + c*dx*dx + d*dy*dy + e*dx*dy
|	  find 5 "nearest" f,dx,dy; use LU decomposition to solve
|	  status: correct results for sample data, real data wrong sign
|	  not handling special cases incl. case where det=0
|	* can't use LU on e.g. 1,0; -1,0; 0,1; 0,-1; 1 1 because
|	  determinant is zero!  These values would apply to any
|	  Cartesian system
-----------------------------------------------------------------------------*/
int nearest_M_value_xx(float xwin, float ywin, CURVE_SET *cp,
		       float *psixy, float *gradx, float *grady, 
		       float *psi_min_dist, float *min_dist)
{
  float value, x, y, z, ysc, dist;
  double dx, dy, dx0, dy0;
  int mr, mz, ir, iz, ir0, iz0, imin, npoints;
  int i, j, k, ik, kj;
  long irz, mrz;
  CVAR xp, yp, zp;
  VIEW *v;
  FILE *fd;
  NEAREST *qn, *qn0, *p, *p0;
  double xyMat[25], A[25], B[25], xyMat2[25], det;
  double c[5], g[5], psi0, x0, y0, dpsi, sum;
  int matSize, nMat, iError, section, irow, icol;
  double df[6];

  //xwin = 0.105;
  //ywin = 0.72;
   
  xworld = xwin;
  yworld = ywin;

  //next line crashes because function is in a different file!?!
  //xwin_to_integer(xwin, ywin, &ixworld, &iyworld);

  printf("mrz=(%d,%d), nnode_r=%d, i_block=%d\n",
	 loop[2].count, loop[1].count,  nnode_r, cp->i_block);

  mr = loop[2].count;
  mz = (cp->i_block >= nnode_r) ? 1 : loop[1].count;
  mrz = mr*mz;

  /*------ output gradient.out for debug */

  if (tell_grad) {
    fd = fopen("gradient.out", "wt");
    fprintf(fd, "nearest_M_value, find coord in '%s', gtype %c, mr=%d, mz=%d\n",
	    cp->title, cp->gtype, mr, mz);
    fprintf(fd, "for point at %g %g\n", xwin, ywin);
    printf("Writing to NEAREST.OUT\n");
    //fprintf(fd, "point at (%g, %g), or (%d, %d)\n",
    //xwin, ywin, ixworld, iyworld);
  }

  /*------ Add points sorted by distance from xwin, ywin */

  xp = cp->ix;
  yp = cp->iy;
  zp = cp->iz;

  for(irz=npoints=0; irz<mrz; irz++) {
    x = *(buf + xp.off0 +irz);
    y = *(buf + yp.off0 +irz);
    z = *(buf + zp.off0 +irz);
    dx = (double)(x-xwin);
    dy = (double)(y-ywin);
    dist = sqrt(dx*dx + dy*dy);

    if (tell_grad)
      fprintf(fd, "%6d %15.5g, %15.5g: %15.5g, d=%lg\n", irz, x, y, z, dist);

    //printf("%ld of %ld\n", irz, mrz);
    npoints = add_this_point(npoints, dist, x, y, z, irz, mr, 50);
  }

  /*------ Abort if not ok */

  if (npoints < 6) {
    printf("Insufficient points for extrapolation (%d, need 6).\n", npoints);
    return 1;
  }

  /*------ Assign sections: 1,0=0  1,1=1  2,0=2  2,1=3  2,2=4 */

  for (i=0, p=nearList; i<npoints; i++, p++) {
    if (i==0) { p->section = 0; ir0 = p->ir; iz0 = p->iz; }
    else {
      ir = abs(p->ir - ir0);
      iz = abs(p->iz - iz0);
      if (iz > ir) { j = iz; iz = ir; ir = j; }		/* ir>=iz and >0 */
      if (ir > 2 || iz > 2 ) p->section = 9;
      else p->section = 2 * (ir-1) + iz;
    }
  }

  /*------ Report (debug) */

  if (tell_grad) {
    for(i=0, p=nearList; i<npoints; i++, p++) {
      if (i==0) { dx0 = p->x; dy0 = p->y; }
      dx = p->x - dx0;
      dy = p->y - dy0;
      fprintf(fd, "%d. (%d,%d:%d) = (%lg,%lg) %.5g %.5g psi=%.5g dist=%.5g\n",
	      i, p->ir, p->iz, p->section, dx, dy,
	      p->x, p->y, p->psi, p->dist);
    }
  }

  /*------ SAMPLE Load Matrix: df = dx + 2dy + 3dx2 - dy2 - 2dxdy */

  //#ifdef DEAD_SAMPLE
  p = nearList;
  x0 = 3.; y0 = 5.;

  p->psi =  0.0;  p->x = 3.0;     p->y = 5.0;     p++;
  p->psi =  4.0;  p->x = x0 + 1.; p->y = y0 + 0.; p++;
  p->psi =  2.0;  p->x = x0 - 1.; p->y = y0 + 0.; p++;
  p->psi =  1.0;  p->x = x0 + 0.; p->y = y0 + 1.; p++;
  p->psi = -3.0;  p->x = x0 + 0.; p->y = y0 - 1.; p++;
  p->psi =  3.0;  p->x = x0 + 1.; p->y = y0 + 1.; 
  npoints = 6;
  //#endif

  /*------ Load Matrix - note nearList[0] contains the really nearest point */

  x0 =   (double)nearList[0].x;
  y0 =   (double)nearList[0].y;
  psi0 = (double)nearList[0].psi;

  matSize = 5;

  for(section=nMat=0; section<=4; section++) {
    for(i=1, p=nearList+1; i<npoints; i++,p++) {
      if (p->section == section) {
	j = nMat*5;

	dx = (double)p->x - x0;
	dy = (double)p->y - y0;
	df[nMat] = (double)p->psi - psi0;

	if (tell_grad) fprintf(fd, "Using i=%d (section %d), dx, dy= %g, %g; dpsi = %g\n", 
			       i, p->section, dx, dy, df[nMat]);
	xyMat[j+0] = dx;
	xyMat[j+1] = dy;
	xyMat[j+2] = dx * dx;
	xyMat[j+3] = dy * dy;
	xyMat[j+4] = dx * dy;
	nMat++;
	if (nMat >= 5) goto NMV_READY;
      }
    }
  }

 NMV_READY:
  if (tell_grad) {
    matrix_tell(fd, xyMat, matSize, matSize, "xyMat");
    det = determinant(xyMat, 5, fd);
    fprintf(fd, "Det = %lg\n\n", det);
  }

  /*------ LU Decomposition: df = xyMat*c = L*U*c == A*B*c == A*g */

  //iError = LUsimple(xyMat, A, B);
  iError =myLU(xyMat, A, B);

  if (tell_grad) {
    matrix_tell(fd, A, matSize, matSize, "A (lower dcomposition, replace diags with 1.0)");
    matrix_tell(fd, B, matSize, matSize, "B (upper decomposition)");
    matrix_tell(fd, df, 1,      matSize, "df[]");

    /*---- Calculate xyMat = A * b */

    for(irow=0; irow<5; irow++) {
      for (icol=0; icol<5; icol++) {
	fprintf(fd, "(%2d,%2d) ", irow, icol);
	for(k=0, sum=0.0; k<5; k++) {
	  ik = irow*5 + k;
	  kj = k*5 + icol;
	  //fprintf(fd, "(%d,%d) %.4lg * %.4lg + ", ik, kj,  A[ik], B[kj]);
	  fprintf(fd, "%11.4lg * %11.4lg + ", A[ik], B[kj]);
	  sum += A[ik] * B[kj];
	}
	fprintf(fd, "= %11.4lg\n", sum);
      }
      fprintf(fd, "\n");
    }
  }

  if (iError) {
    if (tell_grad)
      fprintf(fd, "Matrix cannot be decomposed\n");
    goto NMV_DONE;
  }

  /*------ Solve 'df = A*g' for y */

  g[0] = df[0];
  g[1] = df[1] - A[5 ]*g[0];
  g[2] = df[2] - A[10]*g[0] - A[11]*g[1];
  g[3] = df[3] - A[15]*g[0] - A[16]*g[1] - A[17]*g[2];
  g[4] = df[4] - A[20]*g[0] - A[21]*g[1] - A[22]*g[2] - A[23]*g[3];

  /*------ Solve 'g = B*c' for c */

  c[4] =  g[4] / B[24];
  c[3] = (g[3] - B[19]*c[4]) / B[18];
  c[2] = (g[2] - B[14]*c[4] - B[13]*c[3]) / B[12];
  c[1] = (g[1] - B[9 ]*c[4] - B[8 ]*c[3] - B[7]*c[2]) / B[6];
  c[0] = (g[0] - B[4 ]*c[4] - B[3 ]*c[3] - B[2]*c[2] - B[1]*c[1]) / B[0];

  if (tell_grad)
     fprintf(fd, "\nCoefficients %lg %lg %lg %lg %lg \n\n",
	c[0],c[1],c[2],c[3],c[4]);

  /*------ Test (debug) */

  dpsi = matrix_vector_mult(xyMat, c, 0, 5,  fd, df);
  dpsi = matrix_vector_mult(xyMat, c, 1, 5,  fd, df);
  dpsi = matrix_vector_mult(xyMat, c, 2, 5,  fd, df);
  dpsi = matrix_vector_mult(xyMat, c, 3, 5,  fd, df);
  dpsi = matrix_vector_mult(xyMat, c, 4, 5,  fd, df);

  /*------ Calculate Psi, Gradx, Grady for actual offset */

  *psixy   = c[0]*dx + c[1]*dy + c[2]*dx*dx + c[3]*dy*dy + c[4]*dx*dy;
  *gradx   = c[0] + 2.*c[2]*dx + c[4]*dy;
  *grady   = c[1] + 2.*c[3]*dy + c[4]*dx;
  if (tell_grad) {
    fprintf(fd, "Psi at min dist = %g\n", nearList->psi);
    fprintf(fd, "Calculate psi at %lg, %lg\n", dx, dy);
    fprintf(fd, "Result = %lg\n", *psixy);
  }

 NMV_DONE:
  if (tell_grad) fclose(fd);
  *psi_min_dist = nearList[0].psi;
  *min_dist     = nearList[0].dist;
  return 0;
}			/* ... done nearest_M_value_xx */

/*-----------------------------------------------------------------------------
|	cAngle
-----------------------------------------------------------------------------*/
double cAngle(float xc, float yc, int ir, int iz, int mr,
	      float *px, float *py, double twopi,
	      float *x2, float *y2, FILE *fd)
{
  long irz;
  double dx, dy, theta;
  irz  = iz  * mr + ir;

  *x2 = *(px + irz);
  *y2 = *(py + irz);

  dx = (double)(*x2 - xc);
  dy = (double)(*y2 - yc);
  theta = atan2(dy, dx);
  if (theta < 0.0) theta += twopi;

  if (tell_grad) fprintf(fd, "cAngle, x=%13.4g, %13.4g dx=%13.4lg,%13.4lg, theta=%lg\n", 
			 *x2, *y2, dx, dy, theta);
  return theta;
}

/*-----------------------------------------------------------------------------
|	isBetween
-----------------------------------------------------------------------------*/
int isBetween(double theta, double theta1, double theta2, double twopi)
{
  if (theta1 < theta2) {
    return (theta>=theta1 && theta<theta2) ? 1: 0;
  }

  if (theta>=theta1 && theta<=twopi) return 1;
  else if (theta>=0.0 && theta<theta2) return 1;
  return 0;
}

/*-----------------------------------------------------------------------------
  -----------------------------------------------------------------------------*/
void tell_calc_value(FILE *fd, char *caption, char *s, 
		     float *px, double *a, float idr, float idz, int mr)
{
  double x, dr, dz;
  long irz;
  dr = (double)idr;
  dz = (double)idz;
  x = *px + a[0]*dr + a[1]*dz + a[2]*dr*dr + a[3]*dz*dz + a[4]*dr*dz;
  fprintf(fd, "%s %s: %10.4lg\n", caption, s, x);
}

/*-----------------------------------------------------------------------------
|	myCoefficients
|	* dx = a*dr + b*dz + c*dr*dr + d*dz*dz + e*dr*dz
-----------------------------------------------------------------------------*/
void myCoefficients(float *px, int mr, float r5, float z5, float x5,
		    int ir5, int iz5, int irc, int izc, 
		    double *ax, FILE *fd, char *caption)
{
  double x0, dx30, dx10, dx03, dx01, dx22, fr, fz;
  float *pxi;
  int irz, use_ir5;

  use_ir5 = 0;

  px += izc*mr + irc;

  x0 = *px;

  dx30 = *(px + 1) - x0;
  dx10 = *(px - 1) - x0;

  dx03 = *(px + mr) - x0;
  dx01 = *(px - mr) - x0;

  dx22 = x5 - x0;

  if (dx22 == 0.0) {
    dx22 = *(px + (ir5-irc) + mr*(iz5-izc)) - x0;
    r5 = ir5;
    z5 = iz5;
    use_ir5 = 1;
  }

  ax[0] = (dx30 - dx10) / 2.0;
  ax[2] = (dx30 + dx10) / 2.0;

  ax[1] = (dx03 - dx01) / 2.0;
  ax[3] = (dx03 + dx01) / 2.0;

  fr = r5 - (float)irc;
  fz = z5 - (float)izc;
  //fr = (float)(ir5-irc);
  //fz = (float)(iz5-izc);
  //printf("myCoeff %s, dr5, dz5 = %lg, %lg\n", caption, fr, fz);

  /* If fr or fz is zero, we have a problem! */

  ax[4] = (dx22 - ax[0] * fr - ax[1] * fz - 
	   ax[2] * fr*fr - ax[3] * fz*fz) / (fr*fz);

  if (tell_grad) {
    fprintf(fd, "\n");
    fprintf(fd, "     d%s: %13.4lg, %13.4lg, %13.4lg, %13.4lg, %13.4lg\n",
	    caption, dx10, dx30, dx01, dx03, dx22);
    fprintf(fd, "Coeff %s: %13.4lg, %13.4lg, %13.4lg, %13.4lg, %13.4lg\n",
	    caption, ax[0], ax[1], ax[2], ax[3], ax[4]);

    tell_calc_value(fd, caption, "Right", px, ax,  1,  0, mr);
    tell_calc_value(fd, caption, "Self ", px, ax,  0,  0, mr);
    tell_calc_value(fd, caption, "Left ", px, ax, -1,  0, mr);
    tell_calc_value(fd, caption, "Above", px, ax,  0,  1, mr);
    tell_calc_value(fd, caption, "Below", px, ax,  0, -1, mr);
    if (use_ir5)
      tell_calc_value(fd, caption, "Fifth", px, ax, r5-irc, z5-izc, mr);
    else
      fprintf(fd, "%s %s: %10.4lg at (%g,%g)\n", caption, "Fifth", x5, r5, z5);
  }
}

/*-----------------------------------------------------------------------------
|	solve_derivs
-----------------------------------------------------------------------------*/
int solve_derivs(double *ax, double *ay, double drc, double dzc,
		 double *drdx, double *dzdx, double *drdy, double *dzdy,
		 FILE *fd)
{
  double dxdr, dxdz, dydr, dydz, Det;

  dxdr = ax[0] + 2.0 * ax[2] * drc + ax[4] * dzc;
  dxdz = ax[1] + 2.0 * ax[3] * dzc + ax[4] * drc;

  dydr = ay[0] + 2.0 * ay[2] * drc + ay[4] * dzc;
  dydz = ay[1] + 2.0 * ay[3] * dzc + ay[4] * drc;

  Det = dxdr * dydz - dydr * dxdz;

  if (tell_grad) {
    //fprintf(fd, "dfdr, dfdz = %14.4lg, %14.4lg\n", dfdr, dfdz);
    fprintf(fd, "dxdr, dxdz = %14.4lg, %14.4lg\n", dxdr, dxdz);
    fprintf(fd, "dydr, dydz = %14.4lg, %14.4lg\n", dydr, dydz);
    fprintf(fd, "Det = %lg\n", Det);
  }

  if (Det == 0.0) return 1;

  *drdx =  dydz / Det;
  *drdy = -dxdz / Det;
  *dzdx = -dydr / Det;
  *dzdy =  dxdr / Det;

  return 0;
}

/*-----------------------------------------------------------------------------
  -----------------------------------------------------------------------------*/
//int get_contour_value(float *psixy, float *psi_seg, float *gradx, float *grady)
int get_contour_value(int *pi0, int *pi1,
		      double *pdx, double *pdy, double *pa, double *pb)
{
  double distmin, Det, a, b;
  double dx10, dx20, dy10, dy20, df10, df20, df;
  double x0, y0, x1, y1, x2, y2, dx, dy;
  int i, i0, i1, i2;

  if (n_nearest_seg == 0) {
    printf("No nearest segment%c\n", '\007');
    return 1;
  }

  if (tell_grad) printf("Get_contour_value for %g, %g\n", xworld, yworld);

  for (i=0; i<n_nearest_seg; i++) {
    //printf("%d. psi=%g, dist=%lg at %lg,%lg\n", 
    //i, seg[i].psi, seg[i].distmin, seg[i].xNearest, seg[i].yNearest);
    if (i==0 || seg[i].distmin < distmin) {
      i0 = i;
      distmin = seg[i].distmin;
    }
  }

  for (i=0, i1=-1; i<n_nearest_seg; i++) {
    if (i==i0) ;
    else if (i1==-1 || seg[i].distmin < distmin) {
      i1 = i;
      distmin = seg[i].distmin;
    }
  }

  for (i=0, i2=-1; i<n_nearest_seg; i++) {
    if (i==i0 || i==i1) ;
    else if (i2==-1 || seg[i].distmin < distmin) {
      i2 = i;
      distmin = seg[i].distmin;
    }
  }

  get_world_coordinates(i0, 0, &x0, &y0, "i0");
  get_world_coordinates(i1, 0, &x1, &y1, "i1");
  get_world_coordinates(i2, 0, &x2, &y2, "i2");

#ifdef DEAD_CODE
  dx10 = x1 - x0;
  dy10 = y1 - y0;

  dx20 = x2 - x0;
  dy20 = y2 - y0;
  
  df10 = seg[i1].psi - seg[i0].psi;
  df20 = seg[i2].psi - seg[i0].psi;

  Det = (dx10 * dy20) - (dx20 * dy10);
  printf("Det = %lg from (%lg * %lg - %lg * %lg)\n", 
	 Det, dx20, dy10, dx10, dy20);
  printf("df10=%lg, df20=%lg\n", df10, df20);

  a = (df10 * dy20 - df20 * dy10) / Det;
  b = (df10 * dx20 - df20 * dx10) / Det;

  dx = (double)xworld - x0;
  dy = (double)yworld - y0;

  df = a * dx + b * dy + seg[i0].psi;
  if (tell_grad)
    printf("dx,dy = %lg, %lg\na,b = %lg, %lg\ndf=%lg\n\n", dx, dy, a, b, df);

  *pdx = dx;
  *pdy = dy;
  *pa = a;
  *pb = b;
#endif

  *pi0 = i0;
  *pi1 = i1;

  return 0;

#ifdef DEAD_CODE
  *psixy = (float) (seg[i1].psi + a * dx + b * dy);
  *psi_seg = seg[i1].psi;
  *gradx = (float) a;
  *grady = (float) b;
#endif

  return 0;
}

/*-----------------------------------------------------------------------------
|	nearest_M_value
|	* xwin, ywin are world coordinates
|	* return iError
|	* method: df = a*dr + b*dz + c*dr*dr + d*dz*dz + e*dr*dz
|	  similar for dx and dy; must find dr, dz 
|	* NOTE.  If included e.g. arg 'float *psixy', I got a segmentation
|	  fault when I loaded a perfectly ok value.
-----------------------------------------------------------------------------*/
int nearest_M_value(float xwin, float ywin, CURVE_SET *cp, char *vtext, double *psi_calc)
{
  int mr, mz, mrz, i;
  long irz, irz5;
  int irc, izc;			/* r,z of "center" for calculation */
  int ir0, iz0;			/* r,z of point nearest xwin, ywin */
  int ir5, iz5;			/* r,z of 5th point to use in calculation */
  int i0, i1;			/* nearest, 2nd nearest values in seg[] */

  int ir50, iz50, jr5, jz5, jrz, jr, jz;
  double dr0c, dz0c;		/* dist of nearest from center: ir0-irc */
  double dxwc, dywc;		/* dist of xwin, ywin from xc, yc */

  double theta1, theta2, theta3, theta4, thetaw, twopi;
  float xc, yc, zc, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;
  float x5b, y5b, r5, z5;

  float x, y, z, x0, y0, dxc, dyc, psiMin, psi;

  double dx, dy, dr, df, a, b, x0n, y0n, x1n, y1n;
  double dxw0, dyw0, drw0, dzw0, dw, dxw, dyw, costh, sinth;
  double dist, distMin, xDiff, yDiff, xSum, ySum;
  double ax[5], ay[5], af[5];
  double dfdr, dfdz; // dxdr, dxdz, dydr, dydz;
  double drdx, drdy, dzdx, dzdy;
  double xcalc, ycalc, rcalc, zcalc, xtest, ytest, drc, dzc, xtarget, ytarget;
  double xnew, ynew, rnew, znew;
  double psi0, dpsi, psi0_new, psi0_orig;
  int iter, psi0_inited, ipass, use_limit, show_extra;
  double gradx, grady, ratio;
  char text1[200], text2[200], text3[200], text4[300];

  CVAR xp, yp, zp;
  VIEW *v;
  FILE *fd;
  float *px, *py, *pz;
  float dfdr_on_seg, dfdz_on_seg;
  int use_sample, iError, any;

  /*printf("mrz=(%d,%d), nnode_r=%d, i_block=%d\n",
    loop[2].count, loop[1].count,  nnode_r, cp->i_block);*/

  mr = loop[2].count;
  mz = (cp->i_block >= nnode_r) ? 1 : loop[1].count;
  mrz = mr*mz;

  use_sample = 0;
  if (use_sample) {
    xwin = 0.25; ywin = 0.75;
    xwin = ywin = 0.5;
    xwin = .01162; ywin = .8211;
    xwin = 0.00380581; ywin =  0.00122227;
    printf("USING SAMPLE x,y\n");
  }

  xworld = xwin;
  yworld = ywin;
  iError = 0;

  /*printf("Nearest_M_Value, xworld, yworld = %g, %g\n", xworld, yworld);*/

  /*------ output gradient.out for debug */

  if (tell_grad) {
    printf("\nWriting to GRADIENT.OUT\n");
    fd = fopen("gradient.out", "wt");
    fprintf(fd, "nearest_M_value, find coord in '%s', gtype %c, mr=%d, mz=%d\n",
	    cp->title, cp->gtype, mr, mz);
    fprintf(fd, "Selected point: %g %g\n", xwin, ywin);
  }

  /*------ Get offsets */

  xp = cp->ix;  px = buf + xp.off0;
  yp = cp->iy;  py = buf + yp.off0;
  zp = cp->iz;  pz = buf + zp.off0;

  /*------ Scan Window, find Nearest Segment */

  if (cp->gtype == 'M') {
    i = cp - curveset;
    v = view + i;

    finding_contour_value = 1;
    n_nearest_seg = 0;
    psi0_inited = 0;
    ipass = 0;
    use_limit = 0;
    enable_contour_test(1);

    /*------ Scan all current contour lines (without drawing) */
 
    redraw_mb(cp, v->xmin, v->xmax, v->ymin, v->ymax );
    iError = get_contour_value(&i0, &i1, &dx, &dy, &a, &b);

    if (!psi0_inited) psi0_orig = seg[i0].psi;
    psi0_inited = 1;
    psi0 = seg[i0].psi;
    dpsi = (seg[i0].psi - seg[i1].psi) / 2.0;

    /*---- Scan 2 additional contour lines */

  SCAN_AGAIN:
    if (fabs(dpsi) < 1.0e-8) goto SCAN_DONE;
    if (tell_grad) 
      printf("Scan 2 more based on psi=%.10lg, dpsi=%.10lg\n", 
	     psi0, dpsi);

    psi = (float)(psi0 + dpsi);  set_extra_psi(2, psi);
    redraw_mb(cp, v->xmin, v->xmax, v->ymin, v->ymax );
    iError = get_contour_value(&i0, &i1, &dx, &dy, &a, &b);

    psi = (float)(psi0 - dpsi); set_extra_psi(2, psi);
    redraw_mb(cp, v->xmin, v->xmax, v->ymin, v->ymax );
    iError = get_contour_value(&i0, &i1, &dx, &dy, &a, &b);

    ratio = (double)(seg[i0].distmin / seg[i1].distmin);
    get_world_coordinates(i0, 0, &x0n, &y0n, "");
    get_world_coordinates(i1, 0, &x1n, &y1n, "");
    dx = x1n - x0n;
    dy = y1n - y0n;

    /*------ Scan an "intermediate" contour line */

    if (ratio > 0.75 && ipass < 50) {
      if (tell_grad) 
	printf("Distmin ratio %lg, dx=%lg, dy=%lg\n", ratio, dx, dy);
      if (fabs(x1n-x0n) < 1.0e-8 && fabs(y1n-y0n) < 1.0e-8)
	goto SCAN_DONE;
      /*if (ipass == 49) gfile = fopen("gfile.out", "wt");*/
      psi  = psi0 = (seg[i0].psi + seg[i1].psi) / 2.0;  set_extra_psi(2, psi);
      dpsi = (seg[i0].psi - seg[i1].psi) / 2.0;
      if (tell_grad) 
	printf("Scan 1 at psi=%.10lg, pass=%d\n", psi, ipass);
      redraw_mb(cp, v->xmin, v->xmax, v->ymin, v->ymax );
      iError = get_contour_value(&i0, &i1, &dx, &dy, &a, &b);
      ipass++;
      if (gfile != NULL) { fclose(gfile); gfile = NULL; }
      goto SCAN_AGAIN;
    }

  SCAN_DONE:
    if (tell_grad || ipass>=49)
      printf("Done, pass=%d, Distmin ratio %lg, dx=%lg, dy=%lg\n", ipass, ratio, dx, dy);

    set_extra_psi(0, 0.0);
    finding_contour_value = 0;
    enable_contour_test(0);

    /*---- Use the final contour lines */

    //*psixy = *psi_seg = 999.0; gradx = grady = 0;
    //iError = get_contour_value(psixy, psi_seg, gradx, grady);
    //psiMin = (float) (seg[i0].psi + a * dx + b * dy);
    //*gradx = (float) a;
    //*grady = (float) b;

    get_world_coordinates(i0, 0, &x0n, &y0n, "");
    get_world_coordinates(i1, 0, &x1n, &y1n, "");
    //printf("world %lg,%lg and %lg,%lg\n", x0n, y0n, x1n, y1n);

    df = seg[i1].psi - seg[i0].psi;
    dx = x1n - x0n;
    dy = y1n - y0n;
    dr = sqrt(dx*dx + dy*dy);
    costh = dx / dr;
    sinth = dy / dr;

    dxw0 = (double)xworld - x0n;
    dyw0 = (double)yworld - y0n;

    dxw = dxw0 * costh + dyw0 * sinth;
    dyw = dyw0 * costh - dxw0 * sinth;

    dw = sqrt(dxw0*dxw0 + dyw0*dyw0);
    if (tell_grad) printf("df=%lg, dx=%lg, dy=%lg, dr = %lg,\ndw = %lg, dxw = %lg\n", 
	   df, dx, dy, dr, dw, dxw);

    psi0_new = seg[i0].psi;		/* segmentation error when I try to print this! */
    psiMin = (float)(psi0_new + dxw * df / dr);
    if (psiMin < cp->iz.min) { psiMin = cp->iz.min; use_limit = 1; }
    if (psiMin > cp->iz.max) { psiMin = cp->iz.max; use_limit = 1; }

    show_extra = 0;
    if (psi0_new != psi0_orig && show_extra) {
	set_extra_psi(2, psi0_new);
	redraw_mb(cp, v->xmin, v->xmax, v->ymin, v->ymax );
	set_extra_psi(0, 0.00);
    }

    //if (coord_on != 4 && iError == 0) {

    if (psi0_new == psi0_orig) 
      sprintf(vtext, "     Value = %.10g (nearest contour=%.10g)",
	      psiMin, psi0_new); //psi0_orig);
    else
      sprintf(vtext, "     Value = %g (nearest contour: vis=%lg, added=%lg)",
	      psiMin, psi0_orig, psi0_new);
    //}

    *psi_calc = (double)psiMin;

#ifdef DEAD_CODE
    if (coord_on != 4) 
      sprintf(text, "%s", text1);

    else {
      dpsi = seg[i1].psi - seg[i0].psi;
      dx   = x1n - x0n;
      dy   = y1n - y0n;
      
      gradx = dpsi / dx;
      grady = dpsi / dy;

      sprintf(text4, "based on dpsi = (%lg - %lg) = %lg, dx=%lg, dy=%lg",
	      seg[i1].psi, seg[i0].psi, dpsi, dx, dy);

      if (gradx == -gradx) sprintf(text2, "inf");
      else sprintf(text2, "%lg", gradx);

      if (grady == -grady) sprintf(text3, "inf");
      else sprintf(text3, "%lg", grady);

      sprintf(text, "%s\n     Grad  = %s, %s", text1, text2, text3);
      if (tell_grad) printf("Gradient is %s\n", text4);
    }
#endif
  }

#ifdef DEAD_CODE_CALC
  /*------ Find nearest point */

  for(irz=0; irz<mrz; irz++) {
    x = *(px +irz);		/* float */
    y = *(py +irz);
    z = *(pz +irz);
    dx = (double)(xwin - x);
    dy = (double)(ywin - y);
    dist = sqrt(dx*dx + dy*dy);
    if (tell_grad) fprintf(fd, "%2d,%2d: %11.4g, %11.4g, %11.4g: d=%lg\n",
			   irz%mr, irz/mr, x, y, z, dist);

    if (irz == 0 || dist < distMin) { 
      ir0 = irz % mr;			/* column */
      iz0 = irz / mr;			/* row */
      x0 = x; dxw0 = dx;
      y0 = y; dyw0 = dy;
      distMin = dist;
      psiMin  = z;
    }
  }

  if (coord_on == 4 && iError)
    sprintf(text, "\nNearest grid value = %g (unable to interpolate)%c", psiMin, '\007');

  //*psi_min_dist = psiMin;
  //*min_dist     = distMin;

  goto NMV2_DONE;

  /*------ Find appropriate irc, izc as "center" for calculation */

  irc = ir0;
  izc = iz0;

  if (irc == 0) irc++;
  else if (irc == mr-1) irc--;
  if (izc == 0) izc++;
  else if (izc == mz-1) izc--;

  dr0c = ir0 - irc;		/* these values both zero except at edges */
  dz0c = iz0 - izc;

  /*------ Find dx, dy from "center" */

  irz = izc * mr + irc;
  xc = *(buf + xp.off0 + irz);
  yc = *(buf + yp.off0 + irz);
  zc = *(buf + zp.off0 + irz);

  dxwc = (double)(xwin - xc);
  dywc = (double)(ywin - yc);

  if (tell_grad) {
    fprintf(fd, "\nSelected point: %g, %g\n", xwin, ywin);
    fprintf(fd,   "Nearest point: %2d,%2d (%g, %g) %g\n", 
	    ir0, iz0, x0, y0, psiMin);
    fprintf(fd,   "Center point:' %2d,%2d (%g, %g) %g\n\n",
	    irc, izc, xc, yc, zc);
  }

  /*------ Find angle of each "arm" from center, find which "arms" xwin, ywin is between */

  twopi = 2.0 * atan2(0.0, -1.0);

  theta1 = cAngle(xc, yc, irc+1, izc, mr, px, py, twopi, &x1, &y1, fd);
  theta2 = cAngle(xc, yc, irc, izc+1, mr, px, py, twopi, &x2, &y2, fd);
  theta3 = cAngle(xc, yc, irc-1, izc, mr, px, py, twopi, &x3, &y3, fd);
  theta4 = cAngle(xc, yc, irc, izc-1, mr, px, py, twopi, &x4, &y4, fd);

  if (dxwc==0.0 && dywc==0.0) thetaw = theta1;
  else {
    thetaw = atan2(dywc, dxwc);
    if (thetaw < 0.0) thetaw += twopi;
  }

  if (isBetween(thetaw, theta1, theta2, twopi))      { ir5 = 1;  iz5 = 1; }
  else if (isBetween(thetaw, theta2, theta3, twopi)) { ir5 = -1; iz5 = 1; }
  else if (isBetween(thetaw, theta3, theta4, twopi)) { ir5 = -1; iz5 = -1; }
  else { ir5 = 1; iz5 = -1; }

  ir5 += irc;
  iz5 += izc;

  irz5 = iz5 * mr + ir5;
  x5 = *(px + irz5);
  y5 = *(py + irz5);

  if (tell_grad) {
    fprintf(fd, "\n5th = (%d,%d) for angle %lg\n\n", 
	    ir5-irc, iz5-izc, thetaw);
    fprintf(fd, "Right: %15.4g, %15.4g: %15.4g\n", *(px+irz+1), *(py+irz+1), *(pz+irz+1));
    fprintf(fd, "Self:  %15.4g, %15.4g: %15.4g\n", *(px+irz), *(py+irz), *(pz+irz) );
    fprintf(fd, "Left:  %15.4g, %15.4g: %15.4g\n", *(px+irz-1), *(py+irz-1), *(pz+irz-1) );
    fprintf(fd, "Above: %15.4g, %15.4g: %15.4g\n", *(px+irz+mr), *(py+irz+mr), *(pz+irz+mr));
    fprintf(fd, "Below: %15.4g, %15.4g: %15.4g\n", *(px+irz-mr), *(py+irz-mr), *(pz+irz-mr));
    fprintf(fd, "Fifth: %15.4g, %15.4g: %15.4g\n", x5, y5, *(pz+irz+irz5));
  }

  /*------ If x,y at ir5, iz5 is same as x,y for another "arm" */

  if ((x5==x1 && y5==y1) || (x5==x2 && y5==y2) || (x5==x3 && y5==y3) ||
      (x5==y4 && y5==y4)) {
    ir50 = ir5;
    iz50 = iz5;
    printf("ir50, iz50=%d,%d, irc,izc=%d,%d\n", ir50, iz50, irc,izc);
    for (i=any=0; i<4; i++) {
      if (i==0) { jr5 =  1; jz5 = 1; }
      if (i==1) { jr5 = -1; jz5 = 1; }
      if (i==2) { jr5 = -1; jz5 = -1; }
      if (i==3) { jr5 =  1; jz5 = -1; }
      jrz = (jz5+izc) * mr + (jr5+irc);
      x5b = *(px + jrz);
      y5b = *(py + jrz);
      if ((x5b==x1 && y5b==y1) || (x5b==x2 && y5b==y2) || (x5b==x3 && y5b==y3) ||
	  (x5b==y4 && y5b==y4)) ;
      else {
	dist = (float)sqrt((x5b-x0)*(x5b-x0) + (y5b-y0)*(y5b-y0));
	if (i==0 || dist<distMin) {
	  distMin = dist;
	  ir5 = jr5;
	  iz5 = jz5;
	  x5 = x5b;
	  y5 = y5b;
	  irz5 = jrz;
	  any++;
	}
      }
    }

    if (any==0) printf("NO FIFTH POINT!!%c\n", '\007');
    if (tell_grad) fprintf(fd, "Fifth: %15.4g, %15.4g: %15.4g (at %d,%d)\n", 
			   x5, y5, *(pz+irz5), ir5, iz5);
    ir5 += irc;
    iz5 += izc;
  }

  /*------ Solve for coefficients from eg dx = ax[0]*dr + ax[1]*dz + ax[2]*dr*dr etc */

  //myCoefficients(px, mr, ir5, iz5, irc, izc, ax, fd, "X");
  //myCoefficients(py, mr, ir5, iz5, irc, izc, ay, fd, "Y");
  //myCoefficients(pz, mr, ir5, iz5, irc, izc, af, fd, "F");

  r5 = r_on_seg;
  z5 = z_on_seg;

  myCoefficients(px, mr, r5, z5, x_on_seg,   ir5, iz5, irc, izc, ax, fd, "X");
  myCoefficients(py, mr, r5, z5, y_on_seg,   ir5, iz5, irc, izc, ay, fd, "Y");
  myCoefficients(pz, mr, r5, z5, psi_on_seg, ir5, iz5, irc, izc, af, fd, "F");

  iter = 0;
  xtarget = (double)xwin;
  ytarget = (double)ywin;
  rcalc = (double)ir0; 
  zcalc = (double)iz0;
  xcalc = (double)x0;
  ycalc = (double)y0;
  if (tell_grad) fprintf(fd, "\nBegin Iteration\n");
  goto CMP_ITER;

  /*------ Solve for dr, dz */

 NEXT_ITER:
  iter++;
  drc = rcalc - (double)irc;
  dzc = zcalc - (double)izc;
  if (tell_grad) 
    fprintf(fd, "\nIteration %d: rcalc,zcalc = (%g, %g), irc,izc = (%g, %g)\n",
	    iter, rcalc, zcalc, (double)irc, (double)izc);

  /*------ Get drdx etc at rcalc, zcalc */

  iError = solve_derivs(ax, ay, drc, dzc, &drdx, &dzdx, &drdy, &dzdy, fd);
  if (iError) goto NMV2_DONE;

  /*------ Get rnew = rcalc + drdx*dx + drdy * dy, znew = ditto */
  /*       Get xnew, ynew, then a new rcalc, zcalc, xcalc, ycalc */
  /*       On successive iterations xcalc, ycalc move towards xtarget, ytarget */

  rnew = rcalc + drdx * (xtarget - xcalc) + drdy * (ytarget - ycalc);
  znew = zcalc + dzdx * (xtarget - xcalc) + dzdy * (ytarget - ycalc);

  drc = rnew - (double)irc;
  dzc = znew - (double)izc;

  xnew = xc + ax[0]*drc + ax[1]*dzc + ax[2]*drc*drc + ax[3]*dzc*dzc + ax[4]*drc*dzc;
  ynew = yc + ay[0]*drc + ay[1]*dzc + ay[2]*drc*drc + ay[3]*dzc*dzc + ay[4]*drc*dzc;

  rcalc = rnew; 
  zcalc = znew;
  xcalc = xnew; 
  ycalc = ynew;

  /*------ Test if Done */

 CMP_ITER:
  if (fabs(xtarget) >= 1.0) xtest = xcalc;
  else if (fabs(xtarget) == 0.0) xtest = xcalc;
  else xtest = (xcalc - xtarget) / xtarget;
  xtest = fabs(xtest);

  if (fabs(ytarget) >= 1.0) ytest = ycalc;
  else if (fabs(ytarget) == 0.0) ytest = ycalc;
  else ytest = (ycalc - ytarget) / ytarget;
  ytest = fabs(ytest);

  if (tell_grad) {
    fprintf(fd, "New r,z calc = (%lg, %lg), x,y calc = (%lg, %lg) vs (%lg,%lg)\n",
	    rcalc, zcalc, xcalc, ycalc, xtarget, ytarget);
    fprintf(fd, "x, y test = %lg, %lg\n", xtest, ytest);
  }

  if (xtest <= 1.0e-10 && ytest <= 1.0e-10) goto DONE_ITER;
  if (iter > 200) goto DONE_ITER;

  goto NEXT_ITER;

  /*------ Done,  we have ax[] and bx[], can obtain r, z for any point */
  /*       First do calculations relative to rc, zc */

DONE_ITER:

  drc = rcalc - (double)irc;
  dzc = zcalc - (double)izc;

#endif	/* DEAD_CODE_CALC */

#ifdef DEAD_CODE
  *psixy   = zc + af[0]*drc + af[1]*dzc + 
    af[2]*drc*drc + af[3]*dzc*dzc + af[4]*drc*dzc;

  dfdr   = af[0] + 2.*af[2]*drc + af[4]*dzc;
  dfdz   = af[1] + 2.*af[3]*dzc + af[4]*drc;

  *gradx = dfdr * drdx + dfdz * dzdx;
  *grady = dfdr * drdy + dfdz * dzdy;

  *psi_seg = psi_on_seg;

  if (tell_grad) {
    printf("\nNearest (r,z) = %d,%d, (x,y) = %g,%g, value %g\n",
	   ir0, iz0, x0, y0, psiMin);
    printf("Calc    (r,z) = %lg,%lg, (dr,dz) from nearest = %lg,%lg\n",
	   rcalc, zcalc, rcalc-(double)ir0, zcalc-(double)iz0);
    printf("Coordinate = %.10lg, %.10lg (calculated)\n", xcalc, ycalc);
    fprintf(fd, "\ndr = %lg, dz = %lg\n", drc, dzc);
    fprintf(fd, "Psi at nearest grid point = %g\n", psiMin);
    fprintf(fd, "Psi at nearest segment    = %g\n", psi_on_seg);
    fprintf(fd, "Psi calculated            = %lg\n", *psixy);
  }

  /*------ Better to calculate psi from x_on_seg etc */

  drc = r_on_seg - (double)irc;
  dzc = z_on_seg - (double)izc;

  dfdr_on_seg   = af[0] + 2.*af[2]*drc + af[4]*dzc;
  dfdz_on_seg   = af[1] + 2.*af[3]*dzc + af[4]*drc;

  *psixy   = psi_on_seg + 
    dfdr_on_seg * (rcalc - r_on_seg) + dfdz_on_seg * (zcalc - z_on_seg);

  if (tell_grad) {
    printf("\nNearest (r,z) = %d,%d, (x,y) = %g,%g, value %g\n",
	   ir0, iz0, x0, y0, psiMin);
    printf("Calc    (r,z) = %lg,%lg, (dr,dz) from nearest = %lg,%lg\n",
	   rcalc, zcalc, rcalc-(double)ir0, zcalc-(double)iz0);
    printf("Coordinate = %.10lg, %.10lg (calculated)\n", xcalc, ycalc);
    fprintf(fd, "\ndr = %lg, dz = %lg\n", drc, dzc);
    fprintf(fd, "Psi at nearest grid point = %g\n", psiMin);
    fprintf(fd, "Psi at nearest segment    = %g\n", psi_on_seg);
    fprintf(fd, "Psi calculated            = %lg\n", *psixy);
  }
#endif

 NMV2_DONE:
  if (tell_grad) fclose(fd);
  return iError;
}			/* ... done nearest_M_value */

