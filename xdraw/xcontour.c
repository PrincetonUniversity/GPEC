/******************************************************************************
**  NAME      XCONTOUR.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Handles contour plots for xdraw
**      Notation: "ir" = x-axis, "iz" = y-axis
**
**  Copyright (c) GlassWare 1992.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#if !defined(__MACH__)
#include <malloc.h>
#endif
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

#define Char short
#define Float double

#include "gendefs.h"

#define XCONTOUR		/* used in xcontour.h */

#define ptFlt float
typedef struct
{				/* (locate here because hxdraw.h uses) */
  byte ix, iy;			/* actually ir,iz, but "cartesian" helps me think */
  byte where, used;
  ptFlt x, y;			/* if not M, x = ir*rscale + roffset, NOT in physical space */
} POINTS;			/* if M, x = (x-x_min)*xscale + wx0, YES physical space */

typedef struct
  {
    byte ix, iy;
    byte where;
  } FLAG;

typedef struct
  {
    int level;
    Char rsign, asign;
#ifdef OBSOLETE
    float ax, bx, cx;	/* unused? */
    float ay, by, cy;
#endif
  } BIN;

typedef struct
  {
    float a, b, c, g;	// g is now unused!
    int flag;
  } QCOEFF;

#define POINTS_ POINTS
#define BIN_ BIN
#define NCOLOR ncurve
#define Malloc(i,j) malloc(i*j)

#include "curves.h"
#include "xcontour.h"
#include "xdraw.h"		/* for freadf() */
#include "xtools.h"
#include "setcolor.h"
#include "ps.h"

#define FZ (float)0

static int gtype, polar;	/* WARNING! formerly external, now not set */
static int conlabels;
extern int ps_modeon;
extern int redrawflag;
extern float *buf;
extern LOOP loop[];
extern int nloop;
extern CURVE_SET curveset[];
extern int ftype;
extern BLOCK_XY *xylist;
extern BLOCK_F  *flist;
extern int nnode,nqty,nnode_r;
extern Colormap cmap;
extern int read_entire_topology, multi_topology;

extern float xworld, yworld;
extern int finding_contour_value;
static int need_xywin, need_distmin_inited;
static double distmin;

float psi_on_seg;
float r_on_seg, z_on_seg;
float x_on_seg, y_on_seg;

static BLOCK_XY *block_xy_of_distmin;
static int irmin1, izmin1, irmin2, izmin2;
static byte where_min1, where_min2;
static double xmin1, xmin2, ymin1, ymin2, xNearest, yNearest;
//extern int contour_value_method;

NEAR_SEG seg[5];
int n_nearest_seg = 0;

extern void getlimits_qty(CURVE_SET *cp,int first_block,int last_block,
			  float *min, float *max,
			  int iqt, int itime);
static float *p00;
static float *p10,*p20,*pm0;
static float *p01,*p02,*p0m;

static CURVE_SET *lastcp, *current_cp;
static long off0, xstep, ystep;
static int mr, mz;
static int overflow;

int debug_index = -1;		/* -1 - no debug.  ii=149, i=150 */
static int debug1 = 0;		/* write to debug00.dat .. debugii.dat via savelevel() */
				/* can be auto set via debug_index>=0 */

static int debug2 = 0;		/* write to grid.dat via call savegrid()*/
static int debug3 = 0;		/* write to psi.dat, save details from has_psi() */
static int debug4 = 0;		/* currently just reports call to redraw1 from redraw_mb */
int   show_values = 0;
extern int special_debug;	/* toggle on/off via '!', 'd' */
extern int debug_M2_detail;
extern int tell_grad;
int tell_save_extremum;

extern int debug_m, debug_mb;
/* see alse tellMe in xdraw.c, use of ncurve_dbg etc in redraw1 */

static FILE *file3;

static FILE *dbg_file=NULL;	/* used before drawcontour() if debugging */
static int dbgx1, dbgx2, dbgy1, dbgy2;
static BLOCK_XY *current_block_xy;
static int current_iblock;

static double  twopi, dtheta;
static double  pi = 3.1415926535897931;  
static float arrow_length = 10.0;
static double ar_angle= 0.31415926535897931; /* 18 degrees */ 
static LOOP *rlp;
int flags;

/*=============================================================================
**                    SUPPORTING FUNCTIONS
**===========================================================================*/
extern Display *mydisplay;
extern Window redraw_win;
extern GC redraw_gc;
extern FloatSc xscale, xoffset, yscale, yoffset;
static FloatSc rscale, roffset, zscale, zoffset;
extern int wx0,wy0;
extern int myscreen;

static int isep = -1;
static int ncurve=32, last_ncurve, hzcurve=0;
extern int default_ncurve;
static int nrz, nrza;
static float psi0, dpsi, dpsi0;
static float force;
static int forced = 0;
static int use_extra_psi = 0;
static float extra_psi_value;
static int quad_always = 1;
static int i0 = 0, kz0, kzn, kr0, krn;
/*NN*/
static int xoff,yoff;
static float x_min,x_max;	/* via get_box_info, then arg of redraw_mb eg. */
static float y_min,y_max;
/*NN-end*/
/*-NN - temporary*/
static float hr,hz;
static  Float *rspl,*zspl;

#define PS_AT           0x1	/* bits for ps->where */
#define PS_RIGHT        0x2
#define PS_ABOVE        0x4
#define PS_FROM         0x1	/* bits for ps->used */
#define PS_TO           0x2
#define PS_USED        (PS_FROM | PS_TO)

static POINTS_ *list;
static BIN_ *bin;

/*-----------------------------------------------------------------------------
|	set_extra_psi
  -----------------------------------------------------------------------------*/
void set_extra_psi(int method, float psi)
{
  if (method == 1) {
    printf("Enter psi, or 999 to disable: ");
    scanf("%g", &extra_psi_value);
    printf("draw extra contour, psi = %g\n", extra_psi_value);
    if (extra_psi_value == 999.0) use_extra_psi = 0;
    else use_extra_psi = 1;
    redrawflag = 1;
  }

  if (method == 2) {
    use_extra_psi = 2;
    extra_psi_value = psi;
    //printf("Set psi = %g\n", psi);
  }

  if (method == 0) use_extra_psi = 0;
}

/*-----------------------------------------------------------------------------
|	get_dpsi
-----------------------------------------------------------------------------*/
static void get_dpsi(CURVE_SET *cp, float *zmin, float *zmax)
{
  float min, max, f, delta, df;
  int i, icrit;
  float dist, distmin, psi;
  char text[80];
  
  min = cp->iz.min;
  max = cp->iz.max;
  delta = max - min;
  //printf("get_dpsi, min/max = %g,%g, delta =  %g\n", min, max, delta);

  if (ncurve == 1)
    {
      dpsi = delta;
      psi0 = forced ? cp->force : (min + max) / 2.;
      goto DONE;
    }

  if (!forced)
    {
      df = (ncurve & 1) ? delta / (float)(ncurve+1) : delta / (float)ncurve;
     
      min += df / 2.;
      max -= df / 2.;
      delta = max - min;
      dpsi = delta / (float)(ncurve - 1); psi0 = min;
      goto DONE;
    }

  df = cp->forcemin * delta;
  if (min < force && min+df <= force) min += df;
  df = cp->forcemax * delta;
  if (max > force && max-df >= force) max -= df;

  /*printf("get_dpsi, forcemin/max %g, %g\nmin/max %g, %g\n", 
    cp->forcemin, cp->forcemax, min, max);*/
  
  f = cp->force;			/* ....Forced! */
  if (f < min) f = min;
  if (f > max) f = max;
  
  dpsi = delta / (float)(ncurve-1);		/* trial dpsi */
  icrit=-1; distmin=delta;			/* find which i is closest...*/
  /*printf("get_dpsi, forced value = %g, initial dpsi=%g\n", f, dpsi);*/

  for(i=0; i<ncurve; i++)			/* ... to f */
    {
      psi = min + i * dpsi;
      dist = fabs(psi-f);
      if (i==0 || dist < distmin) { 
	/*printf(" %d. psi=%g, dist=%g\n", i, psi, dist);*/
	distmin=dist; icrit=i;
      }
    }
  /*printf("get_dpsi, i=%d, distmin=%g, i*dpsi=%g\n", icrit, distmin, (float)(icrit*dpsi));*/

#ifdef DEAD_CODE
  if ((icrit == 0 || f > icrit*dpsi) && icrit != ncurve-1)
    dpsi = (max-f)/(float)(ncurve-1-icrit);
  else
    dpsi = (f-min)/(float)icrit;
#endif

  /*------ Get psi0 and psi2 based on current dpsi */

  psi0 = f - (float)icrit * dpsi;	/* may be < or > than min */
  //psi2 = psi0 + (ncurve-1)*dpsi;

  /*------ Adjust */

  if (psi0 < min || psi0 > max) ncurve++;

  /*printf("get_dpsi, new dpsi=%g, new psi0=%g, hzcurve=%d\n", dpsi, psi0, hzcurve);*/

  psi0 = f - (float)icrit * dpsi;

  if (hzcurve)
    {
      dpsi0 = dpsi / (hzcurve + 1);
      hzcurve = 2 * hzcurve + 1;
    }
DONE:
  *zmin = min;
  *zmax = max;
  if (show_values)
    {
       sprintf(text,"Min, max of data: %g,%g\n", cp->iz.min, cp->iz.max);
       xprintf(text);
       sprintf(text,"Min, max used:    %g %g\n", min, max);
       xprintf(text);
       sprintf(text,"%d levels, delta = %g\n", ncurve, dpsi);
       xprintf(text);

       if (ftype==6) {
	 sprintf(text,"Time step = %d, qty number =%d \n", cp->itime, cp->iqty);
	 xprintf(text);
	 if(cp->f_block)
	   sprintf(text,"First block = %d, Last block =%d \n", cp->f_block,
		   cp->l_block);
	 else
	   sprintf(text," All blocks ( %d)\n ", nnode);
	 xprintf(text);
       }
    }
}

/*-----------------------------------------------------------------------------
|	mb_dots.  draw dots for multiblock
|	xmin..ymax are limits according to current zoom
-----------------------------------------------------------------------------*/
void mb_dots(CURVE_SET *cp, float xmin, float xmax, float ymin, float ymax)
{
  float min, max;
  printf("mb_dots, mr=%d, mz=%d, values %f to %f\n",
	 mr, mz, cp->iz.min, cp->iz.max);
}

/*-----------------------------------------------------------------------------
|	redraw1.  draw contours
|	xmin..ymax are limits according to current zoom
|	cp for multiBlock: iblock = current block,
|	ix, iy, iz are type CVAR, have off0, min, max
|	get xp, yp, zp as pointers for current block
-----------------------------------------------------------------------------*/
void redraw1(CURVE_SET *cp, float xmin, float xmax, float ymin, float ymax)
{
  float *p, *pr, *pa, *p0;
  float min, max, psi;
  float slope;
  float dr, dz, fr, fz, rdelt, zdelt, r1a, z1a;
  int count[50];
  CVAR *xp, *yp;
  int i, i0, ii, j, ir, iz, r1, z1, r2, z2, nc;
  int rect_type;
  static int inited = 0;
  double theta, rad;
  float x,y,dx;
  BIN_ *b;
  BIN_ *b2;
  char text[80];
  int ncurve_dbg, icurve1_dbg, icurve2_dbg;
  int tell_bin;

  int excl_flag;
  FILE *fd, *tbFile;

  ncurve_dbg = -1;		// uncommented=no debugging
  //ncurve_dbg = 1;		// # contours (via h or d) for debug section
  icurve1_dbg = 0; 
  icurve2_dbg = 0;

  if (debug2) printf("Write to GRID.DAT\n");
  if (debug3 || ncurve_dbg>=0) {
    printf("Write to PSI.DAT, curves %d to %d\n", icurve1_dbg, icurve2_dbg);
    file3 = fopen("psi.dat", "wt");
    fprintf(file3, "PSI.DAT.  A line is drawn through each section whose last column is not blank.\n");
    fprintf(file3, "interpolation: '-' for linear, '>' for extremum\n");
    fprintf(file3, "A: *p00<psi<=*p10, B: *p00<psi<*p01, C: *p00>psi>=*p10, D: *p00>psi>*p01\n");
    fprintf(file3, "xmin..ymax: %f to %f,  %f to %f\n\n", x_min, x_max, y_min, y_max);
  }

  if (show_values) printf("**redraw1, debug_mb=%d\n", debug_mb);
  if (debug_mb) {
    fd=fopen("debugmb.dat", "a+t");
    fprintf(fd, "**redraw1\n");
  }

  if (cp->iz.min == cp->iz.max)
    {
      /* xprintf("Can't draw contour plot fmin=fmax \n");*/
      return ;
    }

  x_min = xmin; x_max = xmax;
  y_min = ymax; y_max = ymin;
  
  pi = 3.1415926535897931;   /* Some  values in case polar */
  twopi = 2. * pi;
  
  current_cp = cp;
  off0 = cp->iz.off0;
  xstep = cp->iz.dfaml;
  ystep = cp->iz.dstep;
  if (debug_M2_detail) 
    printf("redraw1, block %d, off0=%ld, xstep=%d, ystep=%d\n", 
	   current_iblock, off0, xstep, ystep);

  xp = &cp->ix;
  yp = &cp->iy;

  /*------ Obtain mr, mz */

  if (ftype < 5)
    {
      i = xp->index & INDBITS;  mr = loop[i].count;
      i = yp->index & INDBITS;  mz = loop[i].count;
    }
  else if (ftype == 5 )
    {
      mr = loop[1].count;
      mz = loop[0].count;
      /*printf(" mr=%d mz=%d \n",mr,mz);*/

      hr = 1.0/(float)(mr-1);
      hz = twopi/(float)(mz-1);
      rspl=(Float *) malloc(mr*sizeof(Float));
      zspl=(Float *) malloc(mz*sizeof(Float));
    }

  else if (ftype ==6)		/* mr,mz depends on if drawing 1 block or all */
    { 
	mr = loop[2].count;
        mz = (cp->i_block>= nnode_r)?1:loop[1].count;
    }

  nrz = mr * mz;
  /*if (debug4) xprintf("redraw1, ftype=%d, mr=%d, mz=%d\n", ftype, mr, mz);*/

  /*------ Obtain xoff, yoff */

  if  ((ftype ==5) || (ftype == 6) )
    {
      xoff = xp->off0;
      yoff = yp->off0;
    }

  /*------ Obtain fr, fz */

  rlp = &loop[i];
  dx = (float).1;			/* for solid; should reflect mz! */
  dtheta = twopi / (double)mr;

  i = (cp->gtype == 'S') ? 0 : 1;
  if (ftype != 6)
    {
      if (!(mr-i) || !(mz-i)) return;    /* NN */
      fr = (float)(mr-i);
      fz = (float)(mz-i);
    }
  else
  {
    if (!(nrz-i)) return;
    fr = fz = (float)(nrz-i);    
  }

  /*------ Obtain dr, dz etc */

  dr = (xp->max - xp->min) / fr;
  dz = (yp->max - yp->min) / fz;

  rdelt = rscale = dr * xscale;		/* Use these ==> ir, iz = type X */
  zdelt = zscale = dz * yscale;
  roffset = xp->min * xscale + xoffset;
  zoffset = yp->min * yscale + yoffset;
  //printf("rdelt, zdelt=%g,%g, xscale,yscale=%lg,%lg\n", rdelt, zdelt, xscale, yscale);
  
  gtype = cp->gtype;
  polar = cp->flags & POLAR;
  conlabels = cp->flags & LABELF;
  
  ncurve = cp->ncurve;
  //printf("ncurve=%d\n", ncurve);
  if (ncurve==0) ncurve = cp->ncurve = default_ncurve;
  hzcurve = cp->hzcurve;
  //printf("redraw1, ncurve=%d, hzcurve=%d\n", ncurve, hzcurve);

  forced = cp->flags & FORCED;
  force = cp->force;		/* if not forced, this was set to avg */
  /*flags = cp->flags;*/
  excl_flag = (abs(cp->exclmax - cp->exclmin)<= 0.000001)? 0 : 1;

  /*-------- find kr0..kzn - index limits for tracing contour */

  if( ftype < 5) {
    kr0 = (int)((xmin - xp->min) * fr / (xp->max - xp->min));  kr0 -= 2;
    krn = (int)((xmax - xp->min) * fr / (xp->max - xp->min));  krn += 2;
    kz0 = (int)((ymax - yp->min) * fz / (yp->max - yp->min));  kz0 -= 2;
    kzn = (int)((ymin - yp->min) * fz / (yp->max - yp->min));  kzn += 2;

    if (kr0 < 0) kr0=0;
    if (kz0 < 0) kz0=0;
    if (krn >mr) krn=mr;
    if (kzn >mz) kzn=mz;
  }

  else {
    kr0 = kz0 = 0;
    krn = mr;
    kzn = mz;
  }

  if (polar) { kzn = mz; krn = mr; }

  /*--------- Get dpsi (contour values) */

  get_dpsi(cp, &min, &max);
						     
  /*--------- Allocate list[], bin[] */

  if (!inited)				/* Allocate list[] and bin[] */
    {
      nrza = nrz;
      list = (POINTS_ *)Malloc((long)nrz,sizeof(POINTS));
      bin = (BIN_ *)Malloc((long)nrz,sizeof(BIN));
      if (list == NULL || bin == NULL)
        {
	  printf("Unable to allocate auxiliary array%s\n",
		  list ? "." : "s.");
	  exit(1);
        }
    }
  else if (nrza < nrz)
    {
      nrza = nrz;
      list = (POINTS_ *)realloc(list, (long)nrz * sizeof(POINTS));
      bin = (BIN_ *)realloc(bin, (long)nrz * sizeof(BIN));
    }

  /*--------- Load bin[], one entry for each ir, iz: level, rsign, asign */

  /*if (debug4) 
    xprintf("Ready to load bin[]: ncurve=%d, isep=%d\n", ncurve, isep);*//*isep=-1 normally*/

  if (!inited || ncurve!=last_ncurve ||		/* Load bin[] for points */
      cp != lastcp) {
    lastcp = cp;
    p0 = buf + off0;
    tell_bin = 0;
    if (tell_bin) { tbFile = fopen("bin.out", "wt"); printf("Write BIN.OUT\n"); }

    for (iz=0,b=bin; iz<mz; iz++,p0+=ystep)
      for (ir=0,p=p0; ir<mr; ir++,p+=xstep,b++) {

	/*---- level: determines color */

	if  (*p < min)		b->level = 0;
	else if (*p > max)	b->level = (int) ((max - min) / dpsi);
	else			b->level = (int) ((*p - min) / dpsi);

	b->rsign = b->asign = 0;
	pr = p+xstep;		/* for value at pt to right of current pt */
	pa = p+ystep;		/* for value at pt above current pt */

	/*---- rsign: 1=value of pt to right is gt than current pt, 0=same, -1=less  */
	  
	if (ir == mr-1) slope = FZ;
	else {
	  slope = *pr - *p;
	  if      (slope < FZ) b->rsign = -1;
	  else if (slope > FZ) b->rsign = 1;
	}

	/*---- asign: 1=value of pt above is gt than current pt, 0=same, -1=less*/

	if (iz == mz-1) slope = FZ;
	else {
	  slope = *pa - *p;
	  if      (slope < FZ) b->asign = -1;
	  else if (slope > FZ) b->asign = 1;
	}
	if (tell_bin) fprintf(tbFile, "%3d,%3d: %3d, %2d, %2d based on self=%.10f, rit=%.10f, abv=%.10f\n",
			      ir, iz, b->level, b->rsign, b->asign, *p, *pr, *pa);
      }
  }
  if (tell_bin) fclose(tbFile);

  b2 = bin + nrz;
  last_ncurve = ncurve;
  /*inited = 1; */		/* ???? */

  /*------ Draw solid area - uses filled rectangles */
  /*	   xmin, xmax are for window, xmin_data and xmax_data may be different */

  if (gtype == 'S' || cp->multiDots==1)
    {
      if (ps_modeon)
	return;
      for (i = 0; i <= ncurve; i++)
	count[i] = 0;

      ir = mr - 1;
      iz = -1;

#ifdef DEAD_CODE
      /* txc.bin had x from 0.0 to 1.0 but xmin from -0.248 to 1.248 */

      printf("xp: %g to %g vs %g to %g\n", xp->min, xp->max, xmin, xmax);
      printf("yp: %g to %g vs %g to %g\n", yp->min, yp->max, ymin, ymax);

      r1a = xmin * xscale + xoffset;		/* r1a, z1a are float */
      z1a = ymax * yscale + yoffset;		/* but r1..z2 are int */
#endif

      r1a = xp->min * xscale + xoffset;		/* r1a, z1a are float */
      z1a = yp->min * yscale + yoffset;		/* but r1..z2 are int */

      if (r1a != r1a || z1a != z1a)
	printf("nan for r1a or z1a: multidots\n");

      if (debug_mb) {
	fprintf(fd, "\nSolid contour, x=%g to %g, y=%g to %g\n",
		xmin, xmax, ymax, ymin);
	fprintf(fd, "r1a=%g, z1a=%g\n\n",r1a, z1a);
      }

      for (b = bin; b < b2; b++)		/* (i=0; i<nrz; i++) */
	{
	  ir++;
	  if (ir >= mr)				/* next row */
	    {
	      ir = 0;
	      iz++;
	      r1 = (int)r1a;
	      fz = z1a + (float)iz * zdelt;
	      z1 = (int)fz;
	      z2 = (int)(fz + zdelt);
	    }
	  else					/* same row */
	    r1 = r2;

	  r2 = (int)(r1a + (float)(ir+1) * rdelt);

	  if (debug_mb) {
	    x = (r1 - xoffset)/xscale - xmin;
	    y = (z1 - yoffset)/yscale - ymin;
	    fprintf(fd, "%4d. %4d %4d = %15.4g %15.4g\n", ir, r1,z1,x,y);
	  }

	  if (polar)
	    {
	      /*if (iz==0 && ir>0) continue;*/
	      theta = (double)(ir * 2. * 3.1415926535897931 / mr);
	      rad = (double)(iz * cp->iz.max / (mz-1));
	      x = (float)(rad * cos(theta));
	      y = (float)(rad * sin(theta));
	      r1 = (int)((x-dx) * xscale + xoffset);
	      r2 = (int)((x+dx) * xscale + xoffset);
	      z1 = (int)((y-dx) * yscale + yoffset);
	      z2 = (int)((y+dx) * yscale + yoffset);
	    }

	  i = b->level;
	  count[i]++;
	  setcolor(i, ncurve, isep);		/* color for solid */
	  rect_type = 2;

	  if (ftype == 6) {
	    /*printf("%d, %d: %d\n", ir, iz, i);*/
	    i = iz*mr + ir;
	    //x = buf[i];
	    //y = buf[i+mr*mz];
	    x = buf[xp->off0 + i];
	    y = buf[yp->off0 + i];
	    r1 = (x*xscale + xoffset);
	    z2 = (y*yscale + yoffset);
	    if (r1 != r1 || z1 != z1) 
	      printf("nan encountered at Draw Dot\n");
	    rect_type = 3;
	  }

	  /* Draw Dot: 0=point, 1=big pixels, 2=half size pixeld, 3=pixels fixed delta */

	  if (debug_mb)
	    fprintf(fd, "%4d. %4d %4d = %15.4g %15.4g\n", ir, r1,z1,x,y);

	  if (rect_type == 0)
	    XDrawPoint(mydisplay, redraw_win, redraw_gc, r1, z2);

	  if (rect_type == 1 )
	    XFillRectangle(mydisplay, redraw_win, redraw_gc,
			   r1, z2, (r2 - r1), (z1 - z2));

	  if (rect_type == 2)
	    XFillRectangle(mydisplay, redraw_win, redraw_gc,
			   (r1 + r2)/2, (z1 + z2)/2, (r2 - r1)/2, (z1 - z2)/2);

	  if (rect_type == 3)
	    XFillRectangle(mydisplay, redraw_win, redraw_gc,
			   r1, z2, 4, 4);
	}
    }

  /*------ Draw curves (contours) */

  else
    {
      rand_label(0);

      nc = (use_extra_psi==0) ?  ncurve : (ncurve+1);
      i0 = (use_extra_psi==2) ? (nc-1): 0;
      //printf("nc=%d, i0=%d\n", nc, i0);

      for (i = i0; i < nc; i++) {			/* For each contour.... */
	  psi = psi0 + i * dpsi;
	  if ((psi >= cp->exclmin) && (psi <= cp->exclmax ) ) continue;
	  if (i == ncurve) { 
	    psi = extra_psi_value; 
	    //if (tell_grad) printf("USING psi=%f\n", psi); 
	  }
	  ii = (int) ((psi - min) / dpsi);

	  if ((cp->flags & SINGLE))  {
	      if ( (i !=cp->i_single) && (i != cp->j_single))
	      continue;
	  }
 
	  setcolor(i, ncurve, isep);			/* color for contour */
	  if (debug_index == ii) { 
	    debug1 = 1;
	    xprintf("See debug.dat: savelevel() for psi=%f i=%d, ii=%d\n", psi, i, ii); 
	  }

	  /*---- Draw contour for isep not equal -1 */

	  if (i == isep && hzcurve) {
	    for (psi -= dpsi - dpsi0, j = 0; j < hzcurve; j++, psi += dpsi0) {
	      drawcontour(psi, ii);
	    }
	    psi -= dpsi;
	    continue;
	  }

	  if (show_values) {
	    /*sprintf(text, "Contour %d: %g\n", ii, psi); */
	    sprintf(text, "Contour %d/%d: %g\n", i, ncurve-1, psi);
	    xprintf(text);
	  }

	  /*---- Draw contour for type 6 case */

	  if ((ftype == 6) && (cp->i_block >=nnode_r)) {
	    if (finding_contour_value) need_distmin_inited = 1;
	    //printf("DrawContour_t for psi=%.10g\n", psi);
	    drawcontour_t(psi, ii, xoff, off0, mr, loop[1].count);
	    if (finding_contour_value && !need_distmin_inited)
	      transfer_distmin();
	  }

	  /*---- Draw contour - special debug case */

	  else if (ncurve == ncurve_dbg) {		/* debug drawing: ncurve_dbg is normally -1 */
	    if (i>=icurve1_dbg && i <=icurve2_dbg)  {
	      setcolor(i-icurve1_dbg, icurve2_dbg - icurve1_dbg +1, isep);
	      //psi = 0.00027394;
	      //xprintf("%d. psi=%.10f FORCED\n", i, psi);

	      debug3 = 1;
	      dbg_file = fopen("lines.dat", "wt");
	      fprintf(dbg_file, "Draw contour for psi=%.10f\n", psi);
	      dbgx1=dbgy1= 10000;
	      dbgx2=dbgy2=-10000;

	      drawcontour(psi, ii);

	      fprintf(dbg_file, "limits x=%d to %d, y=%d to %d\n", 
		      dbgx1, dbgx2, dbgy1, dbgy2);
	      fclose(dbg_file);
	      dbg_file = NULL;
	      debug3 = 0;
	    }
	  }

	  /*---- Draw contour - all other cases*/

	  else {
	    if (finding_contour_value) need_distmin_inited = 1;
	    //printf("DrawContour for psi=%.10g, block %d\n", psi, current_iblock);
	    drawcontour(psi, ii);
	    if (finding_contour_value && !need_distmin_inited) 
	      transfer_distmin();
	  }
	  //debug1 = 0;
      }
    }

  show_values = 0;
  if (debug_mb) { fprintf(fd, "**done redraw1\n"); fclose(fd); }
  if (debug3 || ncurve_dbg>=0) fclose(file3);
}

/*-----------------------------------------------------------------------------
|	enable_contour_test
-----------------------------------------------------------------------------*/
void enable_contour_test(int enable)
{
    need_xywin = enable;
    need_distmin_inited = enable;
}

/*-----------------------------------------------------------------------------
|	calculate_xywin
-----------------------------------------------------------------------------*/
void calculate_xywin(float xworld, float yworld, float *px, float *py)
{
  FloatSc ysc, x, y;
  int ix, iy;
  ysc = yscale;
  if (yscale<0) ysc = -yscale;

  x = (xworld - x_min) * xscale;
  *px = x + wx0;

  y = (y_max - yworld) * ysc;
  *py = y + wy0;

  if (x != x || y != y)
    printf("nan encountered in calculate_xywin\n");
}

/*-----------------------------------------------------------------------------
|	test_this_segment
|	* x2 etc are from eg pnow->x (struct POINTS), from list[]
|	* list[] was loaded in has_psi
|	  not-M: ir*rscale + roffset
|	  M:     (x-x_min)*xscale + wx0
-----------------------------------------------------------------------------*/
extern FILE *gfile;
void test_this_segment(float psi, ptFlt x2, ptFlt y2, ptFlt x1, ptFlt y1,
		       int ir2, int iz2, byte where2,
		       int ir1, int iz1, byte where1)
{
  static float xwin, ywin;
  double dist, dr, dx, dy, costh, sinth, x, y;
  double dist1, dist2;
  int abortMe;

  abortMe = 0;

  if (need_xywin) {
    calculate_xywin(xworld, yworld, &xwin, &ywin);
    need_xywin = 0;
    //printf("test_this_segment %g, %g = %d, %d\n", 
    //xworld, yworld, xwin, ywin);
  }

  if (x1==x2 && y1==y2) {
    dx = (double)(x2 - xwin);
    dy = (double)(y2 - ywin);
    dist = sqrt(dx*dx + dy*dy);
  }

  else {
    dx = (double)(x2 - x1);
    dy = (double)(y2 - y1);
    dr = sqrt(dx*dx + dy*dy);
    costh = dx / dr;
    sinth = dy / dr;

    dx = (double)(xwin - x1);
    dy = (double)(ywin - y1);
    x = dx * costh + dy * sinth;
    y = dy * costh - dx * sinth;
    if (x < 0.0 || x > dr) abortMe = 1;
    dist = fabs(y);
  }

  if (gfile != NULL) fprintf(gfile, "psi=%g, dist=%lg, need=%d, block=%d, abort=%d: (%g,%g) to (%g,%g) vs (%g,%g)\n",
			     psi, dist, need_distmin_inited, current_iblock, abortMe,
			     x1, y1, x2, y2, xwin, ywin);
  if (abortMe) {
    dist1 = sqrt(dx*dx + dy*dy);
    dx = (double)(xwin - x2);
    dy = (double)(ywin - y2);
    dist2 = sqrt(dx*dx + dy*dy);
    dist = (dist1 < dist2) ? dist1 : dist2;
    //return;
  }

  if (need_distmin_inited != 0 || dist < distmin) {
    distmin = dist;
    psi_on_seg = psi;

    irmin1 = ir1; izmin1 = iz1; where_min1 = where1;
    irmin2 = ir2; izmin2 = iz2; where_min2 = where2;
    block_xy_of_distmin = current_block_xy;

    xmin1 = x1; ymin1 = y1;	/* these are the (int) end points of the ...*/
    xmin2 = x2; ymin2 = y2;	/* ... contour segment */

    if (gfile != NULL) fprintf(gfile, "entered\n");

    //printf("test_this_segment, psi=%g: (%g,%g) to (%g,%g), dist=%lg: (%d,%d,%d) (%d,%d,%d)\n", 
    //psi, x1, y1, x2, y2, distmin, irmin1, izmin1, where_min1, irmin2, izmin2, where_min2);

    if (x1==x2 && y1==y2) {
      xNearest = x1;
      yNearest = y1;
    }

    else {
      xNearest = x * costh + x1;
      yNearest = x * sinth + y1;
    }

    need_distmin_inited = 0;
    //if (distmin == 0.0) finding_contour_value = -1;

    //printf("test_this_segment, psi=%g: (%d,%d) to (%d,%d), dist=%lg\n", 
    //psi, x1, y1, x2, y2, distmin);
  }
}

/*-----------------------------------------------------------------------------
|	transfer_distmin
-----------------------------------------------------------------------------*/
void transfer_distmin()
{
  int iseg, i;
  double dist_max;
  if (gfile != NULL)
    fprintf(gfile, "transfer_distmin, psi=%g, distmin=%lg, block=%d\n", 
	    psi_on_seg, distmin, current_iblock);

  /*------ If not yet 5 items in seg[], iseg = next */

  if (n_nearest_seg < 5) {
    for(i=0, iseg=-1; i<n_nearest_seg; i++)
      if (seg[i].psi == psi_on_seg) iseg = i;

    if (iseg == -1) iseg = n_nearest_seg++;
  }

  /*------ Find maximum "distmin" if current 5 seg[] */

  else {
    for(i=0; i<n_nearest_seg; i++) {
      if (i==0 || seg[i].distmin > dist_max) {
	dist_max = seg[i].distmin;
	iseg = i;
      }
    }
    if (dist_max <= distmin) return;
    for(i=0; i<n_nearest_seg; i++) {
      if (seg[i].psi == psi_on_seg) iseg = i;
    }
  }

  //printf("transfer_distmin, psi=%g, dist=%lg at %d, block %d\n", 
  //psi_on_seg, distmin, iseg, current_iblock);

  seg[iseg].psi		= psi_on_seg;
  seg[iseg].distmin	= distmin;
  seg[iseg].xNearest	= xNearest;
  seg[iseg].yNearest	= yNearest;
  seg[iseg].irmin1	= irmin1;
  seg[iseg].irmin2	= irmin2;
  seg[iseg].izmin1	= izmin1;
  seg[iseg].izmin2	= izmin2;
  seg[iseg].where_min1	= where_min1;
  seg[iseg].where_min2	= where_min2;
  seg[iseg].xmin1	= xmin1;
  seg[iseg].xmin2	= xmin2;
  seg[iseg].ymin1	= ymin1;
  seg[iseg].ymin2	= ymin2;
  seg[iseg].block_xy	= block_xy_of_distmin;
}

/*-----------------------------------------------------------------------------
|	calc_frfz
|	* x1, y1, ir1, iz1 = values at the vertex
|	* where1: 1=at, 2=right, 4=above
|	* x,y = actual coordinate on ir1, iz1, where1
-----------------------------------------------------------------------------*/
void calc_frfz(float x1, float y1, int ir1, int iz1, byte where1,
	       float *px, float *py, int mr, float *fr1, float *fz1)
{
  double dx, dy, dr, df;
  float x1v, y1v, x1v2, y1v2;
  long irz;

  /*------ x, y of the first vertex point of line containing x1, y1 */

  irz = block_xy_of_distmin->offset + iz1*mr + ir1;
  x1v = x1v2 = *(px + irz);
  y1v = y1v2 = *(py + irz);
  *fr1 = ir1;
  *fz1 = iz1;

  /*------ x, y of the second vertex point */

  if (where1 == 2) {
    x1v2 = *(px + irz + 1);
    y1v2 = *(py + irz + 1);
  }
  if (where1 == 4) {
    x1v2 = *(px + irz + mr);
    y1v2 = *(py + irz + mr);
  }

  //printf("Seek %g,%g in %g,%g to %g,%g\n", x1, y1, x1v, y1v, x1v2, y1v2);

  /*------ Distance between the vertex points */

  dx = x1v2 - x1v;
  dy = y1v2 - y1v;
  dr = sqrt(dx*dx + dy*dy);

  /*------ Fractional distance of x1, y1 along that segment */

  dx = x1 - x1v;
  dy = y1 - y1v;
  df = sqrt(dx*dx + dy*dy) / dr;

  if (where1 == 2) *fr1 += df;
  if (where1 == 4) *fz1 += df;

  //printf("dist between vertices %lg, fraction %lg: r,z=%g, %g\n",
  //dr, df, *fr1, *fz1);

  printf("x,y of seg vertex =      %g,%g, r,z,where=%d,%d,%d: %g,%g\n",
	 x1, y1, ir1, iz1, (int)where1, *fr1, *fz1);

}

/*-----------------------------------------------------------------------------
|	get_world_coordinates
-----------------------------------------------------------------------------*/
void get_world_coordinates(int i, int which,  double *x, double *y, char *caption)
{
  FloatSc ysc;
  float x1, y1;
  ysc = yscale;
  if( yscale<0.0 ) ysc=-yscale;

  if (which == 0) { x1 = seg[i].xNearest; y1 = seg[i].yNearest; }
  if (which == 1) { x1 = seg[i].xmin1;    y1 = seg[i].ymin1; }
  if (which == 2) { x1 = seg[i].xmin2;    y1 = seg[i].ymin2; }

  *x =  (x1 - wx0) / xscale + x_min;
  *y = -(y1 - wy0) / ysc    + y_max;

  if (strlen(caption) && tell_grad)
    printf("%s = %d, psi = %g at %lg, %lg, dist %lg\n", 
	 caption, i, seg[i].psi, *x, *y, seg[i].distmin);
}

/*-----------------------------------------------------------------------------
|	tell_contour_value
-----------------------------------------------------------------------------*/
void tell_contour_value(float *px, float *py, float *pz, int mr,
			float *psixy, float *psi_seg,  
			float *gradx, float *grady)
{
  float x1, y1, x2, y2, x3, y3;
  float fr1, fz1, fr2, fz2, fr3, fz3;
  double dx, dy, dr, f;
  FloatSc ysc;
  int i;

  if (n_nearest_seg == 0) printf("No nearest segment%c\n", '\007');
  else for (i=0; i<n_nearest_seg; i++) {
    printf("%d. psi=%g, dist=%lg\n", i, seg[i].psi, seg[i].distmin);
  }
  return;

  ysc = yscale;
  if( yscale<0.0 ) ysc=-yscale;
  printf("\n");

  /*------ End points of the contour segment */

  x1 = (xmin1 - wx0) / xscale + x_min;
  x2 = (xmin2 - wx0) / xscale + x_min;
  y1 = y_max - (ymin1 - wy0) / ysc;
  y2 = y_max - (ymin2 - wy0) / ysc;

  //printf("x,y of nearest segment = %g,%g to %g,%g\n",
  //x1, y1, x2, y2);

  //printf("x,y of nearest pt = %lg, %lg on %lg,%lg to %lg,%lg\n",
  //xNearest, yNearest, xmin1, ymin1, xmin2, ymin2);

  /*------ Find r,z of the end points */

  calc_frfz(x1, y1, irmin1, izmin1, where_min1,
	    px, py, mr, &fr1, &fz1);

  calc_frfz(x2, y2, irmin2, izmin2, where_min2,
	    px, py, mr, &fr2, &fz2);

  dx = x2 - x1;
  dy = y2 - y1;
  dr = sqrt(dx*dx + dy*dy);
  
  x3 = (xNearest - wx0) / xscale + x_min;
  y3 = y_max - (yNearest - wy0) / ysc;

  dx = x3 - x1;
  dy = y3 - y1;
  f = sqrt(dx*dx + dy*dy) / dr;

  r_on_seg = fr1 + (fr2 - fr1) * f;
  z_on_seg = fz1 + (fz2 - fz1) * f;
  
  x_on_seg = x3;
  y_on_seg = y3;

  printf("x,y of point on segment= %lg,%lg, r,z=%g,%g\n", x3, y3, 
	 r_on_seg, z_on_seg);

  printf("Psi of nearest segment %g: (%d,%d) and (%d,%d), block offset %lx\n",
	 psi_on_seg, irmin1, izmin1, irmin2, izmin2, 
	 block_xy_of_distmin->offset);
}

/*-----------------------------------------------------------------------------
|	drawcontour
-----------------------------------------------------------------------------*/
void drawcontour(float psi, int ii)
{
  POINTS *p, *p1, *p2, *pnow, *pfrom1, *pfrom2, *q, *qln, *qlf, *qlft;
  POINTS *qlf0, *qlft0, *psave, *qsave1, *qsave2;
  int ir, iz, jr, jz, first, pass;
  int n, nvec, ivec, labeli, lx, ly;
  int ind, ind2, di, dj, dii;
  int count, dc_tell;
  byte where, used;
  char text[80];
  FLAG ilist[8], *il, *il2;
  BIN_ *b;
  BIN_ *b0;

  float *xbuf0,*xbuf,*ybuf0,*ybuf;	/*NN*/
  float ysc;				/*NN*/

  /*if (show_values) printf("Level %d: %g\n", ii, psi);*/
  b0 = bin + kz0 * mr;	/* 'bin' is an allocated buffer type BIN */
  xbuf0 = buf + xoff;
  ybuf0 = buf + yoff;
  xbuf = xbuf0;
  ybuf = ybuf0;

  if (debug3==1) {
    fprintf(file3, "drawcontour, Psi=%.8f\n", psi);
    fprintf(file3, "off0=%ld, xstep=%ld, ystep=%ld\n", off0, xstep, ystep);
    //fprintf(file3, "   ir    iz  pm0      p00      p10      p20       p0m      p00      p01      p02\n");
    //fprintf(file3, "   --    --  ---      ---      ---      ---       ---      ---      ---      ---\n");
    fprintf(file3, "         x          y    ir    iz  pm0      p00      p10      p20       p0m      p00      p01      p02\n");
    fprintf(file3, "         -          -    --    --  ---      ---      ---      ---       ---      ---      ---      ---\n");
  }

  /*printf("krn, kzn = %d,%d\n", krn, kzn);*/
  dc_tell = 0;
  if (dc_tell) 
    printf("Psi=%g, limits %g,%g to %g,%g\n", psi, x_min, y_min, x_max, y_max);

  /*------ Load array list[] with all ir, iz which might contain psi */

  overflow = 0;
  for (iz=kz0,p=list; iz<kzn; iz++,b0+= mr,xbuf0+=mr,ybuf0+=mr)
    {						/* ir,iz containing psi */
      if (dc_tell) printf("iz=%d:", iz);
      for (ir = kr0; ir < krn; ir++)
	{
	  if (dc_tell) printf(" %d", ir);

	  /*---- if multiblock, ignore points not in current block */
	  /*     WARN.  x_min etc are limits of current zoom window */

          if (( ftype == 5) || (ftype == 6))
          {
           xbuf = xbuf0 + ir;
           ybuf = ybuf0 + ir;
           if (*xbuf < current_block_xy->xmin || *xbuf > current_block_xy->xmax ||
	       *ybuf < current_block_xy->ymin || *ybuf > current_block_xy->ymax )
	     continue;
           //if ( *xbuf < x_min || *xbuf > x_max ||
	   // *ybuf < y_min || *ybuf > y_max ) continue;
	  }

	  b = b0 + ir;

	  if (debug2) savegrid(ir, iz);		/* write to grid.dat */

	  /*---- Add to list[] only if "nearby" */

	  dii = b->level - ii;
	  if (abs(dii) <= 2) goto MAYBE;
	  if (ir < mr-1 && ((b+1)->level-ii)*dii <= 0)
	    goto MAYBE;
	  if (iz < mz-1 && ((b+mr)->level-ii)*dii <= 0)
	    goto MAYBE;
	  continue;

	  /*---- Load into list[] if "cell" ir,iz contains psi */

	MAYBE:
	  if (dc_tell) printf("!");
	  p = has_psi(psi, ir, iz, *xbuf, *ybuf, p1=p);
	}
      if (dc_tell) printf("\n");
    }

  //if (overflow) printf("Overflow for psi=%g\n", psi);	/* if overflow, program doesn't add point */
  //if (debug3==1) fclose(file3);

  p2 = p;
  nvec = p2 - list;
  ivec = -1;
  labeli = rand_label(nvec);
  count = 0;
  
  if (debug1) {
    if (ii<=9) sprintf(text, "0%d", ii);
    else sprintf(text, "%d", ii);
    savelevel(1, p2 - list, text);	/* Open file, print entire list */
  }

  /*------ 'Dot' option - draw points in list[], don't try to connect */

  if (gtype == 'D')
    {
      for (p = list; p < p2; p++)
	{
	  jr = (int)p->x;
	  jz = (int)p->y;
	  if (!ps_modeon)
	    XDrawLine(mydisplay, redraw_win, redraw_gc, jr, jz, jr, jz);
	  else
	    {
	      ps_moveto(jr, jz);
	      ps_lineto(jr, jz);
	    }
	}
      if (ps_modeon)
	ps_stroke();
      return;
    }

/*----------- list[] = all gridpoints bordering psi ----*/
/*            next line starts with q: a totally unused point */

  if (debug1) {
    savelevel(2, 0, "\nScan list[], seek connections\n");
    savelevel(2, 0, "list[] is struct POINTS: ix,iy, where, used, x,y\n");
    savelevel(2, 0, "offsets in this section link to offsets in previous section\n");
  }

  quad_always = (gtype != 'L');
  quad_always = 0;
  p1 = list - 1;
  pass = -1;

START:					/* Get 1st point on a contour */
  for (q=p1+1; q<p2 && q->used; q++);	/* = totally unused point */
  if (q >= p2)				/* If none found, done! */
    {
      if (debug1)
	savelevel(0, p2 - list, NULL);	/* print new list, close */

      //if (conlabels) {
      if (conlabels && finding_contour_value==0) {
	sprintf(text, "%g", psi);
	XDrawString(mydisplay, redraw_win, redraw_gc, lx, ly,
		    text, strlen(text));
      }
      return;
    }

  if (debug1)
    savelevel(2, pass, "New line, pass=%d\n");

  qlf0 = qlft0 = 0;
  pass = 0;
  first = -1;

  /*--------------- Seek next point q ------*/

  for (p1=q,pnow=pfrom1=0;;)
    {
      qln = qlf = qlft = 0;
      pfrom2 = pfrom1;
      pfrom1 = pnow;
      pnow = q;						/* pnow=last pt drawn */
      nextpoint((pass==3)?0:first, pnow, pfrom1);	/* init angle test */

      /*---- Load ilist[] with rel indices of points surrounding pnow */

      il2 = surround(ilist, pnow);
      if (debug1)
	savelevel(2, pnow - list, "Line from %d ");
      n = il2 - ilist;
      if (ps_modeon && first)
	ps_moveto((int)pnow->x, (int)pnow->y);

      q = (first || pass==3) ? list : p1;

      /*---- Search for all q's that are reasonable nearby (di,dj ok) */

      for (; q < p2; q++)
	{
	  if (q==pnow || q==pfrom1 || q==pfrom2)
	    continue;
	  dj = q->iy - pnow->iy;	/* get q = likely candidate */
	  if (dj < -1) continue;	/* try only i,j in -1..1 */
	  if (dj > 1)  break;
	  di = q->ix - pnow->ix;
	  if (polar && (q->ix==0 || pnow->ix==0) &&
	      (di==mr-1 || di==1-mr));
	  else if (di < -1 || di > 1)
	    continue;
	  if (debug1)
	    savelevel(2, q - list, " Try %d");

	  ind = q->ix + 256 * q->iy;
	  where = q->where;
	  used = q->used;
	  for (il = ilist; il < il2; il++)	/* find it in ilist[]... */
	    {
	      ind2 = il->ix + 256 * il->iy;
	      if (ind != ind2)
		continue;
	      if (!(where & il->where))
		continue;
	      if (!used)			/* ...and set qln etc */
		qln = nextpoint(1, q, qln);	/* ..based on "empty"..*/
	      else if (used != PS_USED)		/* ..or min angle */
		qlf = nextpoint(2, q, qlf);
	      else
		qlft = nextpoint(3, q, qlft);
	    }
	}

      if (first)
	{
	  qlf0 = qlf;
	  qlft0 = qlft;
	}
      if (qln)
	q = qln;		/* connect to unused pt */
      else if (qlf)
	q = qlf;		/* connect- close contour */
      else if (qlft)
	q = qlft;		/* connect to prior contour */
      else			/* no connection: check if bkwrd from... */
	{			/* beginning, if not break out */
	  if (debug1)
	    savelevel(2, 0, "\n**Break\n");
	  goto ONE_MORE;
	}
      if (debug1)
	savelevel(2, q-list, " Using %d\n");
      pnow->used |= PS_FROM;
      q->used |= PS_TO;

      /*------ If finding nearest contour line, test here */

      if (finding_contour_value > 0) {	/* note pnow->x etc are integers (struct POINTS, array list )! */
	test_this_segment(psi, pnow->x, pnow->y, q->x, q->y,
			  pnow->ix, pnow->iy, pnow->where,
			  q->ix, q->iy, q->where);
	//if (finding_contour_value < 0) return;	/* set neg if dist = 0 found */
      }

      /*------ Draw contour line segment on screen here */

      else if (!ps_modeon) {
	  XDrawLine(mydisplay, redraw_win, redraw_gc,
		    (int)pnow->x, (int)pnow->y, (int)q->x, (int)q->y);
	  if (dbg_file != NULL) {
	     fprintf(dbg_file, "%4d. %4d %4d to %4d %4d; %4d %4d to %4d %4d\n", 
		     ++count, (int)pnow->ix, (int)pnow->iy, (int)q->ix, (int)q->iy,
		     pnow->x, pnow->y, q->x, q->y);
	     if ((int)pnow->x < dbgx1) dbgx1 = (int)pnow->x;
	     if ((int)pnow->x > dbgx2) dbgx2 = (int)pnow->x;
	     if ((int)pnow->y < dbgy1) dbgy1 = (int)pnow->y;
	     if ((int)pnow->y > dbgy2) dbgy2 = (int)pnow->y;
	  }
	  ivec++;
	  if (conlabels && (ivec==labeli || ivec==0)) { lx = (int)q->x; ly = (int)q->y; }
        }

      /*------ Print here */

      else
	ps_lineto((int)q->x, (int)q->y);

      /*------ Prepare for next line segment */

      if (pass == 0)			/* Save info re:1st pt, reverse dir */
	{
	  psave = pnow;
	  qsave1 = q;
	  qsave2 = 0;
	  pass = 1;
	}
      else if (pass == 1)
	{
	  pass = 2;
	  qsave2 = q;
	}

      first = 0;
      if (!qln)			/* done if connect was already used */
	{
	ONE_MORE:		/* but first see if start of line connects too */
	  if (pass == 3) break;
	  if (first)
	    {
	      XDrawPoint(mydisplay, redraw_win, redraw_gc, (int)pnow->x, (int)pnow->y);
	      break;
	    }
	  pass = 3;		/* pass=3 ==> reverse of original direction */
	  q = psave;
	  pnow = qsave1;
	  pfrom1 = qsave2;
	  first = -1;
	}
    }

  if (ps_modeon)
    ps_stroke();
  goto START;
}                              /*... end of drawcontour */

/*-----------------------------------------------------------------------------
|	nextpoint
|	* return q of possible next point with minimum angle
|	* k=-1 or 0 ==> initialize ("first" true or false)
|	* k=1,2,3 ==> e.g. closing contour (qln,qlf,qlft)  
-----------------------------------------------------------------------------*/
POINTS_ *nextpoint(int k, POINTS_ * q, POINTS_ * qln)
{
  static double diff[3], angle1[3], angle2[3];
  double angle, newdiff;
  static int inited[3];
  static POINTS_ *pnow;
  static POINTS_ *pfrom;
  static int first;
  int i;
 
#define xangle(f,p) atan2((double)(f->y - p->y),(double)(f->x - p->x))
#define tangle(f,p) fabs(f->y - p->y)<0.0000001 &&fabs(f->x - p->x )<0.00000001
  /*
#define tangle(f,p) f->y==p->y && f->x==p->x

*/
#define tprint(f,p,sf,sp) printf("%s=%d, %s=%d\n", sf, f-list, sp, p-list)

  if (k <= 0)		/* first or next */
    {
      first = k;
      pnow = q;
      pfrom = qln;
      for (i = 0; i < 3; i++)
	inited[i] = 0;
      return (pnow);
    }

  i = k - 1;
  if (qln == NULL)
    qln = q;			/* If don't have one yet, take it */
  else if (!first)		/* Else if already have one, and.. */
    {				/* point is not first on line.. */
      if (!inited[i])
	{
	  /*if (tangle(pnow, pfrom))
	    tprint(pnow,pfrom,"pnow","pfrom");
	  if (tangle(qln, pnow))
	    tprint(qln, pnow, "qln","pnow");*/
  
          if (tangle(pnow, pfrom))
	      angle1[i] = 0;
           else
           angle1[i] = xangle(pnow, pfrom);	/* atan returns -pi to pi */
          if (tangle(qln, pnow))
              angle2[i] = 0;
          else
	      angle2[i] = xangle(qln, pnow);
	  diff[i] = anglecmp(angle1[i], angle2[i]);
	  inited[i]++;
	}
          if (tangle(q, pnow))
            angle = 0;
          else
            angle = xangle(q, pnow);
      /*if (tangle(q,pnow))
	tprint(q,pnow,"q","pnow");*/
      newdiff = anglecmp(angle1[i], angle);
      if (newdiff < diff[i])
	{
	  diff[i] = newdiff;
	  angle2[i] = angle;
	  qln = q;
	}
    }
  return (qln);
}                                                /*... end of nextpoint */

/*-----------------------------------------------------------------------------
|	anglecmp
-----------------------------------------------------------------------------*/
double anglecmp(double a1, double a2)
{
  double d;
  d = fabs(a2 - a1);		/* value 0 to 2*pi */
  if (d > pi)
    d = twopi - d;
  return (d);
}

/*-----------------------------------------------------------------------------
|    surround -- load array 'ilist' with indices of points surround 'ind'
-----------------------------------------------------------------------------*/
FLAG *surround(FLAG * ilist, POINTS * p)
{
  FLAG *ls;
  int iz,ir, i,j, i1,i2, j1,j2, ii,jj, ni;
  int rpolar, apolar;
  byte *q;
  static int inited = 0;
  static byte at[9];
  static byte right[6];
  static byte above[6];
#define X    PS_AT
#define R    PS_RIGHT
#define A    PS_ABOVE

  if (!inited)				/* load these arrays once only */
    {
      q = at;
      *q++ = X | R | A;
      *q++ = X | R;
      *q++ = X | A;
      *q++ = X | A;
      *q++ = 0;
      *q++ = X | A;
      *q++ = X | R;
      *q++ = X | R;
      *q++ = X;

      q = right;
      *q++ = X | R | A;		/* 0,-1 */
      *q++ = X | A;		/* 1,-1 */
      *q++ = A;			/* 0,0 */
      *q++ = A;			/* 1,0 */
      *q++ = X | R;		/* 0,1 */
      *q++ = X;			/* 1,1 */

      q = above;
      *q++ = X | R | A;		/* -1,0 */
      *q++ = R;			/*  0,0 */
      *q++ = X | A;		/*  1,0 */
      *q++ = X | R;		/* -1,1 */
      *q++ = R;			/*  0,1 */
      *q++ = X;			/*  1,1 */
      inited++;
    }

  rpolar = apolar = 0;
  if (p==NULL) printf("NULL pointer\n");
  if (p->where & PS_AT)
    {
      q = at;
      i1 = -1;
      j1 = -1;
    }
  else if (p->where & PS_RIGHT)
    {
      q = right;
      i1 = 0;
      j1 = -1;
      rpolar=polar;
    }
  else if (p->where & PS_ABOVE)
    {
      q = above;
      i1 = -1;
      j1 = 0;
      apolar=polar;
    }

  iz = p->iy;
  ir = p->ix;
  i2 = j2 = 1;
  ni = i2 - i1 + 1;
  for (j=j1, ls=ilist; j<=j2; j++)
    {
      for (i=i1; i<=i2; i++, q++)
	{
	  ii = i;
	  jj = j;
	  if (rpolar)
	    {
	      if (iz==0) continue;		/* no cases of this */
	      if (ir==mr-1 && i==1) ii = -ir + 0;
	      if (ii!=i) goto ADDIT;
	    }
	  else if (apolar)
	    {
	      if (ir==0 && i<0)	    ii = -ir + mr-1;
	      if (ir==mr-1 && i==1) ii = -ir + 0;
	      if (ii!=i) goto ADDIT;
	    }

	  if (ir == 0 && i < 0)
	    continue;
	  if (iz == 0 && j < 0)
	    continue;
	  if (ir == (mr-1) && (i==1 || (i==0 && *q==R)))
	    continue;
	  if (iz == (mz-1) && (j==1 || (j==0 && *q==A)))
	    continue;

	  if (!*q) continue;		/* skip at[0,0] */

	ADDIT:
	  ls->ix = (byte) (ir + ii);
	  ls->iy = (byte) (iz + jj);
	  ls->where = *q;
	  ls++;
	}
    }
  return (ls);
}

/*-----------------------------------------------------------------------------
|	has_psi
|	* see if psi occurs in ir,iz...ir+1,iz or ir,iz...ir,iz+1
|	* if yes, load into q (array list[], struct POINTS), increment q
|	* save_linear, save_extremum generate quad coefficients and
|	  increment q as a counter only!
-----------------------------------------------------------------------------*/
POINTS_ *has_psi(float psi, int ir, int iz, float x, float y, POINTS_ * q)
{
  float delt, delta, psi0, xa,ya,ra,ta;
  float x1, y1, xcalc, ycalc, *px_neighbor, *py_neighbor;
  char ctype[100], fmt1[200], fmt2[200];
  double dbl_delt;
  FloatSc ysc;

  int irz, verysmall, tellMe, tell_debug3;
  QCOEFF kx, ky;
  POINTS_ *q0, *q2, *qi, *p, *psame;
  BIN_ *b;

  tellMe = 0;

  irz = getp00(ir,iz);		/* load p00,p01,p10 etc pointers into buf[] */

  tell_debug3 = (debug3==1  && ir>0 && iz>0 && ir<(mr-1) && iz<(mz-1));
  if (tell_debug3) {
    fprintf(file3, "%10.4g,%10.4g ", x, y);
    fprintf(file3, "%5d %5d: ", ir, iz);
    //fprintf(file3, "%f %f %f %f--", *pm0, *p00, *p10, *p20);
    //fprintf(file3, "%f %f %f %f  ", *p0m, *p00, *p01, *p02);
    fprintf(file3, "%.8f %.8f %.8f %.8f--", *pm0, *p00, *p10, *p20);
    fprintf(file3, "%.8f %.8f %.8f %.8f  ", *p0m, *p00, *p01, *p02);
  }

  b = bin + irz;
  kx.flag = ky.flag = 0;
  q0 = q;
  psi0 = *p00;
  delta = psi - psi0;

  ysc = yscale;
  if( yscale<0.0 ) ysc=-yscale;
  *ctype = '\0';

  strcpy(fmt1, "%s: ir=%d iz=%d\n"
	  " x(ir,iz)= %g y(ir,iz)=%g\n" 
	  " x(ir+1,iz)= %g y(ir+1,iz)=%g\n"
	  " x(ir+delt,iz) =%g y(ir+delt,iz)=%g");

  strcpy(fmt2, "%s: ir=%d iz=%d\n" 
	  " x(ir,iz)= %g y(ir,iz)=%g \n" 
	  " x(ir,iz+1)= %g y(ir,iz+1)=%g \n"
	  " x(ir,iz+delt) =%g y(ir,iz+delt)=%g");

  /*------ Case 1: *p00 == psi, get q->x, q->y directly */

  if (*p00 == psi)
    {
      q->ix = (byte) ir;
      q->iy = (byte) iz;
      q->where = PS_AT;
      q->used = 0;

      if ((ftype != 5 ) && (ftype != 6)) {
        q->x = (ptFlt) ((float) ir * rscale + roffset);
        q->y = (ptFlt) ((float) iz * zscale + zoffset);
      }

      else {
          xcalc = (x - x_min)*xscale + wx0;
	  ycalc = (y_max - y)*ysc    + wy0;

	  if (xcalc != xcalc || ycalc != ycalc) {	/* nan */
	    overflow++;
	    return q;
	  }

	  if (fabs(xcalc) >= 32768. || fabs(ycalc) >= 32768.) {
	    overflow++;
	    return q;
	  }
          q->x = (ptFlt) xcalc;
	  q->y = (ptFlt) ycalc;
      } 

      q++;
      if (tell_debug3) fprintf(file3, "Eq\n");
      return (q);
    }

  /*------ Case 2: *p00 < psi, get kx.a ... ky.c */

  if (*p00 < psi)
    {
      if (polar && iz==0) ;
      else if (polar || ir<(mr-1))		/* r: look for *p10 > psi */
	{
	  if (*p10 == psi) ;
	  else if (*p10 > psi)
	    { q = save_linear(ir, mr, 1, &kx, q); strcat(ctype,"A-"); } /* right, asc, linear */
	  else if (ir == 0 || ir == mr - 2);
	  else if ((b - 1)->rsign <= 0 || (b + 1)->rsign >= 0);
	  else
	    { q = save_extremum(ir, mr, 1, psi, 1, &kx, qi=q);
	      if (q!=qi) strcat(ctype,"A>"); }
	}
      if (iz < (mz - 1))			/* z: look for *p01 > psi */
	{
	  if (*p01 == psi );
	  else if (*p01 > psi)
	    { q = save_linear(iz, mz, mr, &ky, q); strcat(ctype, "B-"); }
	  else if (iz == 0 || iz == mz - 2);
	  else if ((b - mr)->asign <= 0 || (b + mr)->asign >= 0);
	  else
	    { q = save_extremum(iz, mz, mr, psi, 1, &ky, qi=q); 
	      if (q!=qi) strcat(ctype,"B>"); }
	}
    }

  /*------ Case 3: *p00 > psi, get kx.a ... ky.c */

  else						/* ......... *p00 > psi */
    {
      if (polar && iz==0) ;
      else if (polar || ir<(mr-1))		/* r: look for *p10 < psi */
	{
	  if (*p10 == psi) ;
	  else if (*p10 < psi)
	    { q = save_linear(ir, mr, 1, &kx, q); strcat(ctype,"C-"); }
	  else if (ir == 0 || ir == mr - 2);
	  else if ((b - 1)->rsign >= 0 || (b + 1)->rsign <= 0);
	  else
	    { q = save_extremum(ir, mr, 1, psi, -1, &kx, qi=q);
	      if (q!=qi) strcat(ctype,"C>"); }
	}
      if (iz < (mz - 1))			/* z: look for *p01 < psi */
	{
	  if (*p01 == psi) ;
	  else if (*p01 < psi)
	    { q = save_linear(iz, mz, mr, &ky, q); strcat(ctype,"D-"); }
	  else if (iz == 0 || iz == mz - 2);
	  else if ((b - mr)->asign >= 0 || (b + mr)->asign <= 0);
	  else
	    { 
	      //if (ir==27) tell_save_extremum = 1;
	      q = save_extremum(ir, mr, mr, psi, -1, &ky, qi=q);
	      if (q!=qi) strcat(ctype,"D>"); 
	      if (q!=qi) printf("Quad D (r=%d, z=%d) of (%d, %d), psi=%.10f: p00=%.10f, p10= %.10f, p01= %.10f, asign %d, %d\n",
				ir, iz, mr, mz, psi, *p00, *p10, *p01, (b-mr)->asign, (b+mr)->asign);
	      tell_save_extremum = 0;
	    }
	}
    }

  /*------ Now get p->x, p->y for linear or extremum */
  /*       First get x1, y1, then convert to int.    */

  /*if (debug3 && (q!=q0)) fprintf(file3, "%c", ctype); */
  if (tell_debug3) fprintf(file3, "%s\n", ctype);

  q2 = q;			/* q is a count of how many to add */
  for (q=p=q0; q<q2; q++)	/* We have quad coeffs, now add to list */
    {				/* calculate delt, use ir+delt or iz+delt */
      p->ix = (byte) ir;
      p->iy = (byte) iz;
      p->used = 0;
      if (q->where == PS_RIGHT)
	{
	  dbl_delt = *p10 - psi0;
	  if (dbl_delt != 0.0) dbl_delt = (double)delta / dbl_delt;
	  delt = dbl_delt;
	  /*delt = (*p10 == psi0) ? FZ : delta / (*p10 - psi0);*/
	  if (kx.flag) delt = quadratic(&kx, psi, delt);
	}
      else
	{
	  //tellMe = (ir>=52 && ir<=54 && iz==10) ;
	  dbl_delt = *p01 - psi0;
	  if (dbl_delt != 0.0) dbl_delt = (double)delta / dbl_delt;
	  if (tellMe) {
	    printf ("\nHere %d,%d: %g/(%.10g-%.10g), %.10g\n",
		    ir, iz, delta, *p01, psi0, dbl_delt);
	    //debug_temp = 1;
	  }
	  delt = dbl_delt;
	  /*delt = (*p01 == psi0) ? FZ : delta / (*p01 - psi0);*/	/* CRASH HERE */
	  if (ky.flag) delt = quadratic(&ky, psi, delt);		/* and HERE */
	}

      /*if (fabs(delt)<.000001 || fabs(delt-(float)1.)<.000001)
	printf("%d,%d:%d has delt = %f\n", ir,iz,q->where,delt);*/

      verysmall = (delt <= .00001) || ((1.-delt) < .00001);
      p->where = q->where;

      /*------ Not Polar */

      if (!polar) {
	xa = (float)ir;
	ya = (float)iz;

	/*------ Not Polar, yes Multiblock */
          
	if ( (ftype==5) || (ftype == 6) ) {
	  if (p->where == PS_RIGHT) {			      /* -- worked */
	    px_neighbor = buf + xoff + iz*ystep + (ir+1)*xstep ;
	    py_neighbor = buf + yoff + iz*ystep + (ir+1)*xstep ;
	    x1 = x + delt*(*px_neighbor-x);
	    y1 = y + delt*(*py_neighbor-y);

#ifdef DEAD_CODE
	    printf(fmt1, "PS_RIGHT",ir,iz,x,y,*px_neighbor,*py_neighbor,x1,y1);
	    get_splr(&x1,&y1,ir,iz,delt);
	    /* printf (" After spline: x =%g y=%g \n", x1,y1); */
#endif
	  }

	  else if (p->where != PS_AT) {		      /* - worked!! */
	    px_neighbor = buf + xoff + (iz+1)*ystep + ir*xstep ;
	    py_neighbor = buf + yoff + (iz+1)*ystep + ir*xstep ;
	    x1 = x + delt*(*px_neighbor-x);
	    y1 = y + delt*(*py_neighbor-y);

#ifdef DEAD_CODE
	    sprintf(fmt2, "PS_AT",ir,iz,x,y,*px_neighbor,*py_neighbor,x1,y1);
	    get_splz(&x1,&y1,ir,iz,delt);
	    /* printf (" After spline: x =%g y=%g \n", x1,y1);*/
#endif	/* ... of DEAD_CODE */
	  }
	}

	/*------ Not Polar, not Multiblock */

	else { 
	  if (p->where == PS_RIGHT)   xa += delt;
	  else if (p->where != PS_AT) ya += delt;
	}

	/*------ Not Polar and we have x1, y1 so find p->x, p->y */

	if ((ftype == 5) || (ftype == 6)) {
	  xcalc = (x1-x_min) * xscale + wx0;
	  ycalc = (y_max-y1) * ysc    + wy0;

	  if (xcalc != xcalc || ycalc != ycalc) {	/* nan */
	    overflow++;
	    return q;
	  }

	  if (fabs(xcalc) >= 32768. || fabs(ycalc) >= 32768.) {
	    overflow++;
	    return q;
	  }
	  p->x = (ptFlt)xcalc;
	  p->y = (ptFlt)ycalc;

#ifdef DEAD_CODE
	  if (tell_debug3) {
	    fprintf(file3, "x calc = %g = %d using %lg, %g, %d\n", 
		    x1, (int)p->x, x_min, xscale, wx0 );
	    fprintf(file3, "y calc = %g = %d using %g, %lg, %d\n", 
		    y1, (int)p->y, y_max, ysc, wy0 );
	  }
#endif

	  /*p->x = (ptFlt)(xa * xscale + xoffset);
	    p->y = (ptFlt)((y_max-ya) * yscale + yoffset); */
	}

	else {
	  p->x = (ptFlt)(xa * rscale + roffset);
	  p->y = (ptFlt)(ya * zscale + zoffset);
	}

	if (verysmall && (psame = same_as(p)) )
	  { psame->where = PS_AT; p--; }
      }

      /*------ Yes Polar */	

      else {
	ta = (float)ir;
	ra = get_loopval(rlp, iz);
	if (p->where == PS_RIGHT) ta += delt;
	else if (p->where != PS_AT)
	  ra += delt * (get_loopval(rlp, iz+1) - ra);
	ta *= dtheta;
	xa = ra * cos(ta);
	ya = ra * sin(ta);
	p->x = (ptFlt)(xa * xscale + xoffset);
	p->y = (ptFlt)(ya * yscale + yoffset);
      }

      p++;
    }

  return (p);
}                                      /*... end of has_psi */

/*-----------------------------------------------------------------------------
|	same_as
-----------------------------------------------------------------------------*/
POINTS_ *same_as(POINTS_ *p)
{
  POINTS_ *q;
  int ix, iy;
  ix = p->ix;
  iy = p->iy;
  for(q=p-1; q>=list && q->iy>=iy-1; q--)
    {
      if (q->ix >= ix-1 && q->ix <= ix+1 && 
	  ((int)q->x==(int)p->x && (int)q->y==(int)p->y))
	return q;
    }
  return 0;
}

/*-----------------------------------------------------------------------------
|	getp00
-----------------------------------------------------------------------------*/
int getp00(int ir,int iz)
{
  int irz;
  long sign;
  irz = mr * iz + ir;
  
  p00 = buf + off0 + ir * xstep + iz * ystep;
  pm0 = p00 - xstep;
  p10 = p00 + xstep;
  p20 = p10 + xstep;

  p0m = p00 - ystep;
  p01 = p00 + ystep;
  p02 = p01 + ystep;
  if (!polar) return(irz);
			  /* "iz"=r, "ir"=theta */
  sign = (ir < mr/2) ? xstep : -xstep;
  if (iz==0)
    p0m = p01 + sign * mr/2;

  if (ir==0) pm0 = p00 + xstep * (mr -1);
  if (ir==mr-1) p10 = p00 - xstep * (mr-1);
  p20 = (ir!=mr-2) ? p10 + xstep : p00 - xstep * (mr-2);
  
  return(irz);
}

/*-----------------------------------------------------------------------------
|	save_linear
-----------------------------------------------------------------------------*/
POINTS_ *save_linear(int ir, int mr, int n, QCOEFF *kp, POINTS_ * q)
{
  int right;
  right = (n==1);
  if (quad_always && (ir > 0 && ir < mr - 2))
    {
      getquadcoeff(n, ir, mr, kp);
      kp->flag++;
    }
  else
    {
      kp->a = (float) 0;
      kp->b = right ? *p10 - *p00 : *p01 - *p00;
      kp->c = *p00;
    }
  q->where = right ? (byte) PS_RIGHT : (byte) PS_ABOVE;
  return (++q);
}

/*-----------------------------------------------------------------------------
|	save_extremum
|	how: 1 if *p00 < psi, -1 if *p00 > psi
|	we come here because slope at 0 has opposite sign from slope at 1
|	this means there's a peak or valley in 0,1
-----------------------------------------------------------------------------*/
POINTS_ *save_extremum(int ir, int mr, int n, float psi,
		       int how, QCOEFF *k, POINTS_ * q)
{
  int ok, psi_is_gt_p00;
  float x, psix;
  psi_is_gt_p00 = (how == 1);

  //if (ir==0 || ir==mr-2) return(q);
  //if (b1sign <= 0 || b2sign >= 0) return(q);

  //------ Find a,b,c such that psi = axx + bx + c in x=0,1

  getquadcoeff(n, ir, mr, k);

  //------ Is there a solution for THIS psi in 0,1?
  //       dpsi = 2ax + b = 0

  ok = 0;
  if (k->a > 0.) {
    x = -k->b / (2.0 * k->a);
    if (x>=0. && x<=1.) {
      psix = k->a * x * x + k->b * x + k->c;
      if (*p00<=psi && psix>=psi) ok = 1;
    }
  }

#ifdef DEAD_CODE
  k->g = k->c - k->b * k->b / (4. * k->a);	// contents of sqrt = -(k->g/a)
  k->flag++;

  //------ Hope for an easy test which says if there's a solution for x
  //       2a*psi + b = +,- sqrt(bb-4ac)
  //       (4aa*psi*psi + bb + 4a*psi*b) = bb - 4ac
  //       4a*psi(a*psi + b) = -4ac
  //       psi (a*psi + b) = c     ...DUH!! I just proved the quadratic equation

  ok = 0;
  if ((k->b*k->b - 4.0*k->a*k->c) > 0.) ok = 1;

  if (how == 1)							// *p00<psi, so solution is +
    ok = (k->b > FZ && (float) 2 * k->a < -k->b && k->g >= psi);
  else								// *p00>psi, so solution is -
    ok = (k->b < FZ && (float) 2 * k->a > -k->b && k->g <= psi);
#endif

  if (ok)
    {
      q->where = (n == 1) ? (byte) PS_RIGHT : (byte) PS_ABOVE;
      q++;

      if (tell_save_extremum) {
	printf("save_extremum for *p00=%.10g, mr=%d, n=%d\na = %.10g\nb = %.10g\nc = %.10g\ng = %.10g\n", 
	       *p00, mr,n, k->a, k->b, k->c, k->g);
	for (x=0; x<=1.1; x+=0.1) {
	  psix = k->a * x * x + k->b * x + k->c;
	  printf("%f. %.10f\n", x, psix);
	}
      }
    }
  return (q);
}

/*-----------------------------------------------------------------------------
|   getquadcoeff
-----------------------------------------------------------------------------*/
void getquadcoeff(int n, int ir, int mr, QCOEFF * k)
{
  float a1, a2, b1, b2, c1;
  static int cubic = 0;
  float af0,af1,af2,af3;
  
  af1 = *p00;
  if (n==1)
    {
      af2 = *p10;
      af3 = *p20;
      af0 = *pm0;
    }
  else
    {
      af2 = *p01;
      af3 = *p02;
      af0 = *p0m;
    }

  if (!cubic) goto QUAD;

  k->a = (af3 - 3 * af2 + 3 * af1 - af0) / 6.;
  k->b = (3 * af2 - 6 * af1 + 3 * af0) / 6.;
  k->c = (-af3 + 6 * af2 - 3 * af1 - 2 * af0) / 6.;
  return;
           
  //------ This calc is done by averaging two functions f = axx + bx + c,
  //       one centered on *p01 and one centered on *p02 ??
  //       I don't really follow this, but it was really clear when I first developed it

QUAD:
  if (tell_save_extremum) printf("getquadcoeff basis %.10g, %.10g, %.10g, %.10g\n", af0, af1, af2, af3);
#ifdef DEAD_CODE
  c1 = (af1 + af2) / 2.;		// avg value at 1.5
  b1 = (af2 - af0) / 4.;		// .5 * slope at 1.0
  b2 = (af3 - af1) / 4.;		// .5 * slope at 2.0
  a1 = (af2 - 2 * af1 + af0) / 4.;
  a2 = (af3 - 2 * af2 + af1) / 4.;	// ((af3-af2) - (af2-af1)) / 4

  k->c = c1 - b2 + a2;			// i.e. k->c = af1
  k->b = b1 + b2 - 2. * a2;		// i.e. k->b = ((af2-af0) + 2*(af2-af1)) / 4
  k->a = a1 + a2;			// i.e. k->a = ((af3-af1) - (af2-af0)) / 4
#endif

  //------ Solve f = a(x-0)(x-0) + b(x-0) + c, ie centered on x1
  //       and   f = a(x-1)(x-1) + b(x-1) + c, ie centered on x2

#ifdef DERIVATION
  c1 = af1;
  c2 = af2;

  b1 = (af2-af0) / 2;
  b2 = (af3-af1) / 2;

  a1 = ((af2-af1) - (af1-af0)) / 2;
  a2 = ((af3-af2) - (af2-af1)) / 2;

  2*f = a1*x*x + a2*(x-1)*(x-1) + b1*x + b2*(x-1) + c1 + c2;
  2*f = c1+c2 -b2 +a2 + x*x*(a1+a2) + x*(b1+b2-a2-a2);

  4*c = af1+af1 + af2+af2 + -af3+af1 + af3-af2 -af2+af1               = 4*af1;
  4*a = af2-af1 -af1+af0 +af3-af2 -af2+af1                            = (af3-af2) - (af1-af0);

  4*b = af2-af0 +af3-af1 -af3+af2 -af3+af2 +af2-af1 +af2-af1;
      = af2-af1 + af1-af0 +af3-af2 +af2-af1 -2*(af3-af2) +2*(af2-af1) = 4*(af2-af1) +af1-af0 -af3+af2;

  check x=0: 4*f = 4*af1;
  check x=1: 4*f = 4*af1 + af3-af2-af1+af0 + 4*(af2-af1) +af1-af0 -af3+af2 = 4*af2;
#endif

  k->c = af1;
  k->a = ((af3-af2) - (af1-af0)) / 4.0;
  k->b = (af2-af1) + ( af1-af0-af3+af2) / 4.0;
}

/*-----------------------------------------------------------------------------
|	quadratic
|	* find point in (0..1) where value = psi
|	* x0 = linear likely guess (psi-psi1)/(psi2-psi1)
-----------------------------------------------------------------------------*/
float quadratic(QCOEFF * k, float psi, float x0)
{
  float f, df;
  float x, xx, dxx, err, x1, x2;
  double z;
  int count, tellMe;
  char c;

  count = 0;
  tellMe = 0;

  k->c -= psi;
  for (x = x0;;)
    {
      count++;
      if (count>200) {
	z = k->b * k->b - 4.0 * k->a * k->c;
	z = sqrt(z)/(2.0*k->a);
	x1 = -k->b + (float)z;
	x2 = -k->b - (float)z;
	if (x1>=0.0 && x1<=1.0) x = x1;
	else if (x2>=0.0 && x2<=1.0) x = x2;
	else x = x0;
	if (tellMe) {
	  printf("quadratic alternate, psi=%.10g, x0=%.10g, a=%.10g %.10g %.10g\n",
			 psi, x0, k->a, k->b, k->c);
	  printf("x1=%.10g, x2=%.10g, Press any key..", x1, x2); 
	  c=getchar(); 
	}
	goto QUAD_DONE;
      }
      f = k->a * x * x + k->b * x + k->c;
      if (f == FZ) goto QUAD_DONE;
      df = 2 * k->a * x + k->b;
      xx = x;
      dxx = f / df;
      x = xx - dxx;
      err = (x - xx) / (x + xx);
      if (err < FZ)
	err = -err;
      if (err < (float) .0001) goto QUAD_DONE;
    }
 QUAD_DONE:
  return (x);
}

#ifdef DEAD_CODE_CONT
/*-----------------------------------------------------------------------------
|	get_cont_minmax
-----------------------------------------------------------------------------*/
void get_cont_minmax(int how, int i, float *xn, float *xx,
				     float *yn, float *yx)
{
  float rfac, zfac;
  struct CONT *p;
  p = cdata + i;
  
  if (how == 1 && !polar)
    {
      *xn = *yx = (float) 0;
      *xx = (float) (p->mr - 1);
      *yn = (float) (p->mz - 1);
    }
  else if (how==1 && polar)
    {
      *xn = *yx = -(float)(p->mz-1);
      *xx = *yn =  (float)(p->mz-1);
    }
  else
    {
      rfac = (p->rmax - p->rmin) / (float) (p->mr - 1);
      zfac = (p->zmax - p->zmin) / (float) (p->mz - 1);
      *xn = *xn * rfac + p->rmin;
      *xx = *xx * rfac + p->rmin;
      *yn = *yn * zfac + p->zmin;
      *yx = *yx * zfac + p->zmin;
    }
}
#endif

/*-----------------------------------------------------------------------------
|	new_ncurve
|	* (return value is unused)
-----------------------------------------------------------------------------*/
int new_ncurve(CURVE_SET *cp, char how)
{
  //printf("new_ncurve, type %c, code %c\n", cp->gtype, how);
  if (cp->gtype == 'G') return(0);
  if (cp->gtype=='V' && (how=='D' || how=='H')) {
    if (how=='H') cp->vDensity++;
    else if (cp->vDensity>1) cp->vDensity--;
    //printf("new_ncurve, density %d\n", cp->vDensity);
  }

  else {
    if (how >= 'a' && how <= 'z') how += ('A'-'a');
    //printf("new_ncurve, %d --> ", cp->ncurve);
    if (how == 'D') {
      if (2*cp->ncurve > 5000) {
	printf("Maximum # curves (5000) would be exceeded.\n");
	return (1);
      }
      cp->ncurve *= 2;
    }
    else cp->ncurve /= 2;
    if (cp->ncurve < 1) cp->ncurve = 1;
    //printf("%d, ncurve=%d\n", cp->ncurve, ncurve);
  }
  redrawflag = 1;
  return (1);
}
/************************************************************/
/*-----------------------------------------------------------------------------
|	new_nvect
-----------------------------------------------------------------------------*/
int new_nvect(CURVE_SET *cp, char how)
{
  if (cp->gtype != 'V') return(0);
  if (how == 'I') cp->fskip >>= 1;
  else	if ( how == 'i' ) cp->fskip <<=1;
  else  if ( how == 'J') cp->lstep >>=1;
  else  if ( how == 'j') cp->lstep <<=1;
  if (cp->fskip < 1) cp->fskip = 1;
  if (cp->lstep < 1) cp->lstep = 1;  
  redrawflag = 1;
  return (1);
}




/*-----------------------------------------------------------------------------
|	contour_values
-----------------------------------------------------------------------------*/
void contour_values(CURVE_SET *cp)
{
   if (cp->gtype != 'G')
      show_values = redrawflag = 1;
}

/*=============================================================================
**                  DEBUG FUNCTIONS
**===========================================================================*/
/*-----------------------------------------------------------------------------
|   savegrid
|       write grid.dat = contains values on the grid.  if debug2 set.
-----------------------------------------------------------------------------*/
void savegrid(int ir, int iz)
{
  static FILE *f2;
  static int report_float=0;
  float xpsi;
  int ipsi,i,jz;
  static char text[] = "     ";
  static char format[] = "%5d";	/* # chars in max, below */
  static float mul = (float) 1000.;
  static float max = (float) 10000;
  /*static int irmin=0, irmax = 33;*/	/* limit region */
  static int irmin=0, irmax = 200;
  static int use_dots=0;
  float pmax;

  if (ir == mr - 1 && iz == mz - 1)
    {
      fclose(f2);
      debug2 = 0;
      return;
    }

  if (ir == 0 && iz == 0)
    {
      f2 = fopen("grid.dat", "wt");
      pmax = current_cp->iy.max;			/* (assumes values are all positive) */
      if (pmax > (float)999) mul=(float)1;
      else if (pmax > (float)99) mul=(float)10;
      else if (pmax > (float)9)  mul=(float)100;
      else mul=(float)1000;
      fprintf(f2, "GRID.DAT, mr=%d, mz=%d, multiplier=%d\n", mr, mz, (int)mul);
      fprintf(f2, "display psi for all ir, iz\n\n");
      fprintf(f2, "      ");				/* leading spaces for header line */
      for (i=0; i<=mz; i+=10) {
	if (report_float) fprintf(f2,"%9d", i);
	else fprintf(f2,"%5d", i);
	if (i<mz) fprintf(f2,"                                             ");
      }
      fprintf(f2, "\n");
    }

  if (ir == 0)
    fprintf(f2, "\n%4d. ", iz);
  /*fprintf(f2, "\n");*/

  if (ir >= irmin && ir < irmax)	/* (Can't fit ALL on the page) */
    {
      jz = mz-1-iz;
      if (!polar) i = jz * mr + ir;
      else if (jz==0) i=0;
      else i = 1 + (jz-1)*mr + ir;
      /*xpsi = *(buf + +off0+i) * mul;*/		/* mul ==> # decimal places keep */
      xpsi = *(buf +off0+i) * mul;		/* mul ==> # decimal places keep */
      if (xpsi > FZ)
	xpsi += (float) .5;
      else
	xpsi -= (float) .5;

      /*------ Print data in buf[] as float, in actual iz */

      if (report_float) {
	i = (iz-1)*mr + ir;
	fprintf(f2, "%f ", *(buf+off0+i));
      }

      /*------ Print dots and spaces */

      else if (use_dots==1) {
	if (xpsi<1700 || xpsi>1900)	/* a particular debug case */
	  fprintf(f2, " ");

	else {
	  ipsi = (int) (xpsi);
	  fprintf(f2, ".");
	}
      }

      /*------ Print data in buf[] * multiplier */

      else {
	if (xpsi >= max || xpsi <= -max)
	  fprintf(f2, text);

	else {
	  ipsi = (int) (xpsi);
	  fprintf(f2, format, ipsi);
	}
      }
    }
}

/*-----------------------------------------------------------------------------
|	savelevel -- from debug1, writes to disk all points at current level
|	how=1: "index"="nlist"=p2-list, output all of list
|	how=0: "index"=offset into list to calculate a p2, print "Done", then F/T
|	how=2: "index"=an integer to print if included in format
-----------------------------------------------------------------------------*/
void savelevel(int how, int index, char *format)
{
  POINTS *p, *p2;
  static FILE *file;
  char text[5], fname[20];

  if (how == 1)
    {
      //printf("savelevel 1, level %s\n", format);
      p2 = list + index;
      strcpy(fname,"debug");
      strcat(fname, format);
      strcat(fname, ".dat");
      file = fopen(fname, "wt");
      //file = fopen("debug.dat", "wt");
      fprintf(file, "Summary of array list[] (struct POINTS)\n");
      fprintf(file, "From drawcontour() calls savelevel(), after has_psi() loop\n");
      fprintf(file, "Level %s\n", format);
      fprintf(file, "The array is loaded in has_psi()\n\n");
      fprintf(file, "'where' values <VRAO> from p->where, & with 1,2,4,F8\n\n");
      fprintf(file, "      ir  iz where       x    y\n");
      fprintf(file, "      --  -- -----       -    -\n");

      /*fprintf(file, "Format is index: ix,iy,<where> ==> x,y\n");*/
      for (p = list; p < p2; p++)
	{
	  strcpy(text, "    ");
	  if (p->where & 1)  *(text) = 'V';
	  if (p->where & 2)  *(text + 1) = 'R';
	  if (p->where & 4)  *(text + 2) = 'A';
	  if (p->where & 0xf8) *(text + 3) = 'O';
	  fprintf(file, "%3d: %3d,%3d<%s> ==> %4d,%4d\n",
		  p - list, p->ix, p->iy, text, (int) p->x, (int) p->y);
	}
    }

  else if (how == 0)
    {
      p2 = list + index;
      fprintf(file, "\nDone scanning list[]\n");
      for (p = list; p < p2; p++)
	{
	  strcpy(text, "  ");
	  if (p->used & PS_FROM) *text = 'F';
	  if (p->used & PS_TO)   *(text + 1) = 'T';
	  fprintf(file, "%d:<%s> \n", p - list, text);
	}
      fclose(file);
    }

  else
    {
      fprintf(file, format, index);
      fclose(file);
      file = fopen("debug.dat", "at");
    }
}

#ifdef DEAD_CODE_CONT
/*-----------------------------------------------------------------------------
|   psicount
-----------------------------------------------------------------------------*/
void psicount(float *buf, long bufsize, float psi0, float dpsi,
	      int ncurve, float *count)
{
  int i;
  float *p, *buf2;
  long nvalues;

  for (i = 0; i <= ncurve; i++)
    count[i] = FZ;

  nvalues = bufsize / 4L;
  buf2 = buf + nvalues;

  for (p = buf; p < buf2; p++)
    {
      i = (int) ((*p - psi0) / dpsi);
      count[i]++;
    }
}
#endif

/* Spline fitted value for x,y - r axes  */

void  get_splr(float *x,float *y,
               int ir, int iz, float delta)
{

  int i;
  Float t,ht;

  float *xbuf0,*xbuf,*ybuf0,*ybuf;
  double *xspl,*yspl, *xvec,*yvec;

  xbuf0 = buf + xoff;
  ybuf0 = buf + yoff;

  xspl = (Float *)malloc(mr*sizeof(Float));
  yspl = (Float *)malloc(mr*sizeof(Float));

  xvec = (Float *)malloc(mr*sizeof(Float));
  yvec = (Float *)malloc(mr*sizeof(Float));

   for ( i=0,xbuf=xbuf0+iz*mr,ybuf=ybuf0+iz*mr;
         i<mr; i++)
   {
     *(xvec+i)=*(xbuf+i);
     *(yvec+i)=*(ybuf+i);
   }  
  
      splfit(rspl,xvec,xspl,mr-1);
      ht=(Float) (ir+delta)*hr;
      spleval(rspl,xvec,xspl,&t,&t,&t,&t,ht,mr-1,mr-1,1,0);
      *x=(float) t;

      splfit(rspl,yvec,yspl,mr-1);
      spleval(rspl,yvec,yspl,&t,&t,&t,&t,ht,mr-1,mr-1,1,0);
      *y=(float) t;

   free( xvec); free(yvec); free(xspl); free(yspl);
      
  }      
         

/* Spline fitted value for x,y  - for z    */

void  get_splz(float *x,float *y,
               int ir, int iz, float delta)
{

  int i;
  Float t,ht;

  float *xbuf0,*xbuf,*ybuf0,*ybuf;
  Float *xspl,*yspl, *xvec,*yvec;

  xbuf0 = buf + xoff;
  ybuf0 = buf + yoff;

  xspl = (Float *)malloc(mz*sizeof(Float));
  yspl = (Float *)malloc(mz*sizeof(Float));

  xvec = (Float *)malloc(mz*sizeof(Float));
  yvec = (Float *)malloc(mz*sizeof(Float));

   for ( i=0,xbuf=xbuf0+ir,ybuf=ybuf0+ir;
         i<mz; i++)
   {
     *(xvec+i)=*(xbuf+i*mr);
     *(yvec+i)=*(ybuf+i*mr);
   }  
  
      splfit(zspl,xvec,xspl,mz-1);
      ht=(Float) (iz+delta)*hz;
      spleval(zspl,xvec,xspl,&t,&t,&t,&t,ht,mz-1,mz-1,1,0);
      *x=(float) t;

      splfit(zspl,yvec,yspl,mz-1);
      spleval(zspl,yvec,yspl,&t,&t,&t,&t,ht,mz-1,mz-1,1,0);
      *y=t;

   free( xvec); free(yvec); free(xspl); free(yspl);
      
    }      

/*-----------------------------------------------------------------------------
|	redraw_mb.  draw contours for multiblock
|	xmin..ymax are limits according to current zoom, WORLD coordinates
-----------------------------------------------------------------------------*/
void redraw_mb(CURVE_SET *cp, float xmin, float xmax, float ymin, float ymax)
{
  BLOCK_XY *xy;
  BLOCK_F  *fl, *fl2, *ftmp;

  float ac_xmin,ac_xmax,ac_ymin,ac_ymax;

  CURVE_SET ac_cp;
  CVAR *xp, *yp, *zp;
  XTextProperty newName_TP;
  char newName[200];

  int k,iblock;
  int zpp;
  int ix,iy,iz;
  long ixy, nxy;
  int fstep;
  int itime,iqty;
  int iqty2;
  int imr,imz;
  int first_bl,last_bl;

  float zmin,zmax;
  float fmin,fmax;
  float x, y, z;

  int mod;
  char text[80];

  float l_scale;
  float tt1,tt2;
  int num_points;
  float lmax, lmax_bl;

  FILE *fd;
  if (debug_mb) { fd = fopen("debugmb.dat", "a+t"); fprintf(fd, "**redraw_mb\n"); }

  x_min = xmin; x_max = xmax;
  y_min = ymax; y_max = ymin;
  //printf("Redraw_mb, xmin..ymax = %g %g %g %g\n", xmin, xmax, ymin, ymax);
  //printf("Redraw_mb, width, height = %g %g\n", xmax - xmin, ymax - ymin);

  k=0;
  num_points =0;
  need_distmin_inited = 1;

  itime = cp->itime;
  iqty  = cp->iqty;
  flags = cp->flags;
  if (multi_topology==2 && read_entire_topology==0)	// itime is really...
    itime = 0;						// ...offset in buf to timestep

  if (cp->f_block) {
    first_bl = cp->f_block-1;
    last_bl  = cp->l_block;
  }
  else {
    first_bl = 0;
    last_bl = nnode;
  }

  fstep = 1;
  xp = &cp->ix;
  yp = &cp->iy;
  zp = &cp->iz;

  if (debug_M2_detail) printf("redraw_mb, view %d: blocks %d to %d at x_off=%lx\n",
			      cp->view_number, first_bl, last_bl, xp->off0);

  if (debug_m)
    printf("redraw_mb, time=%ld, qty=%ld, gtype=%c, blocks %d to %d\n",
	   itime, iqty, cp->gtype, first_bl, last_bl);

  fl = flist + itime *(nnode+1)*nqty + iqty*(nnode+1)+1 + first_bl;

  if( cp->gtype == 'V')			// drawxx.in plot line starts with two iqty's
    {
      iqty2 = cp->iqty2;
      fl2 = flist + itime *(nnode+1)*nqty + iqty2*(nnode+1)+1 + first_bl;
      if (show_values) {
	sprintf(text,"          VECTOR PLOT information\n");
	xprintf(text);
      }
    }
  ftmp = flist + itime *(nnode+1)*nqty + iqty*(nnode+1);
  
  xy = xylist + 1 + first_bl;

  fmin = zp->min = ftmp->fmin;
  fmax = zp->max = ftmp->fmax;

  //---- Add extrema to window title

  sprintf(newName,"%s,  extrema=(%10.3e,%10.3e)", cp->title, fmin, fmax);
  newName_TP.value = newName;
  newName_TP.encoding = XA_STRING;
  newName_TP.format = 8;
  newName_TP.nitems = strlen(newName);
  XSetWMName(mydisplay,cp->window, &newName_TP);

  /*---- if  vector plot option   - found maximum of length */
  /*     cp->gtype is normally 'M' */

  if ((cp->gtype == 'V' ) && (!(cp->flags & SINGLE)))
    {
      //printf("Vector - not single\n");
      iqty2 = cp->iqty2;
      fl2 = flist + itime *(nnode+1)*nqty + iqty2*(nnode+1)+1 + first_bl;
      lmax =0.0;
      num_points = 0;
      if(cp->ncurve == 0)cp->ncurve=default_ncurve;

      for (iblock = first_bl; iblock < last_bl; 
	   iblock++, xy+=1, fl+=fstep,fl2+=fstep)
	{
	  if (x_min > xy->xmax ) continue;
	  if (y_min > xy->ymax)  continue;
	  current_block_xy = xy;
	  current_iblock = iblock;

	  imr = xy->mr;
	  imz = (iblock < nnode_r)?xy->mz : 1;
	  lmax_bl = 0.0;

	  if ( cp->fskip >= imr) cp->fskip = 1;
	  if ( cp->lstep >= imz) cp->lstep = 1;
	  k = cp->fskip;
	  zpp = cp->lstep;
	 
          draw_v( imr, imz, xy->offset, xy->offset + imr*imz,
                  fl->offset, fl2->offset,
		  &lmax_bl, &k, &zpp, l_scale, 0, cp->vDensity );
	  
          num_points+=k;
          lmax = (lmax_bl > lmax )? lmax_bl:lmax;
	  if (show_values)
	    {
	      sprintf(text," Block number = %d \n number of points =  %d\n",
		      iblock,k);
	      xprintf (text);
	      sprintf(text," Max length of vectors =%g\n",
		      lmax_bl);
	      xprintf (text);
	    }
	}
    
      fl = flist + itime *(nnode+1)*nqty + iqty*(nnode+1)+1 + first_bl;
      fl2 = flist + itime *(nnode+1)*nqty + iqty2*(nnode+1)+1 + first_bl;
      ftmp = flist + itime *(nnode+1)*nqty + iqty*(nnode+1);
  
      xy = xylist+1+first_bl;
      zp->min = ftmp->fmin;
      zp->max = ftmp->fmax;

      tt1= (x_max-x_min)*xscale;
      tt2= (ymax -ymin)*yscale;
      l_scale =sqrt((tt1*tt1 +tt2*tt2)/ (num_points))*cp->ncurve/16;
      cp->force = l_scale;
      cp->forcemax = lmax;

    }		/* .. end of 'V' (vector) with lmax */

  if (show_values && cp->gtype == 'V' && !(cp->flags & SINGLE))
    {
      sprintf(text,"  Number of points =  %d\n", num_points);
      xprintf (text);
      sprintf(text," Max length of vectors =%g\n  length scale = %g\n",
	      lmax,l_scale);
      xprintf (text);
      show_values =0;
    }

  /*------ Scan blocks */

  for (iblock = first_bl; iblock < last_bl; iblock++, xy+=1, fl+=fstep,
	 fl2+=fstep)
    {
      if (show_values) 
	printf("**Block %d of %d: %g vs %g, %g vs %g\n", iblock, last_bl,
	       xmin, xy->xmax, ymax, xy->ymax);

      /*------ Test for skip block, don't see scenario for these conditions */

      if (xmin > xy->xmax ) continue;
      if (ymax > xy->ymax)  continue;
      
      imr = xy->mr;
      imz = (iblock < nnode_r)?xy->mz :1;
      //if (debug4) xprintf("iblock %d, imr=%d, imz=%d\n", iblock, imr, imz);
      if (debug_mb) fprintf(fd, "iblock %d, imr=%d, imz=%d\n", iblock, imr, imz);

      //printf("redraw_mb, block %d: %g,%g to %g,%g\n", 
      //iblock, xp->min, yp->min, xp->max, yp->max);
      current_block_xy = xy;
      current_iblock = iblock;

      xp->min = xy->xmin;
      xp->max = xy->xmax;
      xp->off0= xy->offset;

      yp->min = xy->ymin;
      yp->max = xy->ymax;
      yp->off0= xy->offset + imr*imz;

      zp->min = fmin;
      zp->max = fmax;
      zp->off0 = fl->offset;
      zp->dfaml = 1;
      zp->dstep = imr;

      loop[2].count = imr;
      loop[1].count = xy->mz;

      cp->i_block = iblock;

      /*---- Draw Dots - e.g. via 'xdraw m -v' */

      if (cp->multiDots == 1)
	{
	  if (debug_mb) fprintf(fd, "\nMultiDots, block %d, imr=%d, imz=%d\n", iblock, imr, imz);
	  nxy = imr * imz;
	  for (ixy=0; ixy<nxy; ixy++) {
	    x = *(buf + xp->off0 + ixy);
	    y = *(buf + yp->off0 + ixy);
	    z = *(buf + zp->off0 + ixy);
	    if (debug_mb) fprintf(fd, "%4d. %15.5g %15.5g %15.5g\n", ixy, x, y, z);
	  }
	}

      /*---- not Vector: call regular contour */

      if (cp->gtype !='V')
	{
	  if (debug_mb) xprintf("redraw_mb, call redraw1 for block %d, imr=%d, imz=%d\n"
			      , iblock, imr, imz);
	  //if (debug_M2_detail) printf("redraw_mb, block %d\n", iblock);
	  redraw1(cp,xmin,xmax,ymin,ymax);

	  if (finding_contour_value < 0) return;	/* set neg when dist=0 found */

	  if (finding_contour_value == 0) {
	    if ((flags & MARKERS)&&(iblock <nnode_r))
	      drawgrid(xmin,xmax,ymin,ymax);

	    if ((flags & MARKERS)&&(iblock >= nnode_r))
	      drawgrid_t(xmin,xmax,ymin,ymax,iblock,
			 xy->mr,xy->mz,xy->offset);
	  }
	}

      /*---- yes Vector: draw Vector */

      else
	{
	  if (flags & SINGLE)
	    {
	      //printf("Vector - yes SINGLE\n");
	      lmax = cp->forcemax;
	      l_scale = cp->force;
	      /*   draw a white vector grig */
	      if ( cp->i_single == imr )
		{
		  cp->i_single = 0;
		  cp->j_single +=1;
		  if (cp->j_single == imz )
		    cp->j_single =0;
		};

	      ix = imr*cp->j_single+cp->i_single;

	      /* test if points is in properly area */
	      /* zmin, zmax now are free            */

	      zmin = *(buf + xy->offset+ix);		/* x coordinate */ 
	      zmax = *(buf + xy->offset+ix+imr*imz);	/* y coordinate */ 

	      if ((zmin>=x_min )&&(zmin <= x_max)
		  && (zmax >= y_min) && ( zmax <=y_max))
		{
		  if (show_values)
		    {
		      sprintf(text," Vector coordinates are i=%d, j=%d\n",
			      cp->i_single,cp->j_single);
		      xprintf (text);
		    }

		  draw_v(ix, imz, xy->offset, xy->offset + imr*imz,
			 fl->offset, fl2->offset,
			 &lmax, &k, &zpp, l_scale, 2, cp->vDensity);
		}

	      else
		if(show_values)
		  {
		    sprintf(text," Vector (i=%d, j=%d) is out of range\n",
			    cp->i_single,cp->j_single);
		    xprintf (text);
		    show_values =0;
		  }
	    }

	  else				/* ... not SINGLE */
	    {
	      //printf("Vector, not SINGLE\n");
	      mod =1;
	      if ( flags & LABELF ) mod = 3;
	      if ( flags & MARKERS) mod = 4;
	      k = cp->fskip;
	      zpp = cp->lstep;
	      draw_v( imr, imz, xy->offset, xy->offset + imr*imz,
		      fl->offset, fl2->offset,
		      &lmax, &k, &zpp, l_scale, mod, cp->vDensity);
	    }

	  if( ps_modeon) ps_stroke();
	}				/*... end of yes Vector */
    }
  if (debug_mb) { fprintf(fd, "**done redraw_mb\n");  fclose(fd); }
}

/*-----------------------------------------------------------------------------
|	drawgrid
-----------------------------------------------------------------------------*/
void drawgrid(float xmin, float xmax, float ymin, float ymax)
{
  int ir, iz, ix, iy;
  char marker;

  BIN_ *b;
  BIN_ *b0;
  float *xbuf0,*xbuf,*ybuf0,*ybuf, fx, fy;
  FloatSc ysc;					/*NN*/

  //printf("drawgrid, marker is forced as #\n");
  marker = '#';		/* # for box, anything else for point */

  /*printf("Contour %d: %g\n", ii, psi);*/
  b0 = bin + kz0 * mr;
  xbuf0 = buf + xoff;
  ybuf0 = buf + yoff;
  xbuf = xbuf0;
  ybuf = ybuf0;
  if(!ps_modeon)
    XSetForeground(mydisplay, redraw_gc, WhitePixel(mydisplay, myscreen));
  else
    ps_color(FZ,FZ,FZ);

  if( yscale<0.0 ) ysc=-yscale;

  for (iz=kz0; iz<kzn; iz++,b0+= mr,xbuf0+=mr,ybuf0+=mr)
    {
      for (ir = kr0; ir < krn; ir++)
	{
	  xbuf = xbuf0 + ir;
	  ybuf = ybuf0 + ir;
	   
	  if ( *xbuf < x_min || *xbuf > x_max ||
	       *ybuf < y_min || *ybuf > y_max )
	    continue;
	   
	  b = b0 + ir;

	  if (!ps_modeon || ftype==6) {
	    fx = (*xbuf-x_min)*xscale+wx0;
	    fy = (y_max-*ybuf)*ysc +wy0;
	    ix = (int)((*xbuf-x_min)*xscale+wx0);
	    iy = (int)((y_max-*ybuf)*ysc +wy0);
	    if (marker != '#')		// currently forced as #
	      XDrawPoint(mydisplay, redraw_win, redraw_gc, ix, iy);
	    else
	      drawmarker(ix, iy, fx, fy);
	    }
	}
    }
  if (ps_modeon) ps_stroke();
}

/*-----------------------------------------------------------------------------
|	drawgrid triangle
-----------------------------------------------------------------------------*/
void drawgrid_t(float xmin, float xmax, float ymin, float ymax,int iblock,
		int mvert,int mcells,int xy_off)
{
  int  itr;

  int *p;
  int iv1,iv2,iv3;
  int ix1,ix2,ix3,iy1,iy2,iy3;
/*NN*/
  float *xbuf,*cellbuf,*ybuf;
  FloatSc ysc;

  xbuf = buf + xy_off;
  ybuf = xbuf + mvert;
  cellbuf = ybuf+ mvert;
  p    = (int *) cellbuf;
  XSetForeground(mydisplay, redraw_gc, WhitePixel(mydisplay, myscreen));

  if( yscale<0.0 ) ysc=-yscale;
  for (itr=1; itr <= mcells ; itr++,cellbuf+=1)
    {

      iv1 = *(int *)cellbuf;
      iv2 = *(int *)(cellbuf+ mcells);
      iv3 = *(int *)(cellbuf+ 2*mcells);

      ix1 = (int)((*(xbuf+iv1)-x_min)*xscale+wx0);
      ix2 = (int)((*(xbuf+iv2)-x_min)*xscale+wx0);
      ix3 = (int)((*(xbuf+iv3)-x_min)*xscale+wx0);
      iy1 = (int)((y_max - *(ybuf+iv1))*ysc+wy0);
      iy2 = (int)((y_max - *(ybuf+iv2))*ysc+wy0);
      iy3 = (int)((y_max - *(ybuf+iv3))*ysc+wy0);
      if (!ps_modeon)
	{
      
	  /* Draw grid lines if you want 
	     XDrawLine(mydisplay, redraw_win, redraw_gc,
	     ix1,iy1,ix2,iy2
	     );
	     XDrawLine(mydisplay, redraw_win, redraw_gc,
	     ix2,iy2,ix3,iy3
	     );
	     XDrawLine(mydisplay, redraw_win, redraw_gc,
	     ix3,iy3,ix1,iy1
	     );
	  */
		    
    
	  XDrawPoint(mydisplay, redraw_win, redraw_gc,
		     ix1,iy1
		     );
	  XDrawPoint(mydisplay, redraw_win, redraw_gc,
		     ix2,iy2
		     );
	  XDrawPoint(mydisplay, redraw_win, redraw_gc,
		     ix3,iy3
		     );
	} 
    }
}


/*-----------------------------------------------------------------------------
|	drawcontour_t
-----------------------------------------------------------------------------*/
void drawcontour_t(float psi, int ii,int xy_off,int f_off,
		   int mvert,int mcells)
{
  int  itr;
  int iv1,iv2,iv3;
  int u1,u2,u3;
  int ix1,ix2,ix3,iy1,iy2,iy3;

  float *xbuf,*cellbuf,*ybuf;
  float *fbuf;
  float cx1,cx2,cx3,cy1,cy2,cy3;
  FloatSc ysc;

  xbuf = buf + xy_off;
  ybuf = xbuf + mvert;
  fbuf = buf +  f_off;
  cellbuf = ybuf+ mvert;

  if( yscale<0.0 ) ysc=-yscale;

  for (itr=1; itr <= mcells ; itr++,cellbuf+=1)
    {
      iv1 = *(int *)cellbuf;
      iv2 = *(int *)(cellbuf+ mcells);
      iv3 = *(int *)(cellbuf+ 2*mcells);
      /*           if ( *xbuf < x_min || *xbuf > x_max ||
                *ybuf < y_min || *ybuf > y_max )
		*/
      u1=FindPoint(*(xbuf+iv1),*(ybuf+iv1),*(fbuf+iv1),
		    *(xbuf+iv2),*(ybuf+iv2),*(fbuf+iv2),
		    psi,
		    &cx1,&cy1);
	
      u2=FindPoint(*(xbuf+iv1),*(ybuf+iv1),*(fbuf+iv1),
		   *(xbuf+iv3),*(ybuf+iv3),*(fbuf+iv3),
		   psi,
		   &cx2,&cy2);

      u3=FindPoint(*(xbuf+iv3),*(ybuf+iv3),*(fbuf+iv3),
		   *(xbuf+iv2),*(ybuf+iv2),*(fbuf+iv2),
		   psi,
		   &cx3,&cy3);
      
	if ( u1 && u2 )
	  {
	    ix1 = (int)((cx1-x_min)*xscale+wx0);
	    iy1 = (int)((y_max - cy1)*ysc+wy0);
	    ix2 = (int)((cx2-x_min)*xscale+wx0);
	    iy2 = (int)((y_max - cy2)*ysc+wy0);

	    if ( !ps_modeon)     
	      XDrawLine(mydisplay, redraw_win, redraw_gc,
			ix1,iy1,ix2,iy2
			);
	    else
	      ps_line(ix1,iy1,ix2,iy2);

	  };
	if ( u1 && u3)
	  {
	    ix1 = (int)((cx1-x_min)*xscale+wx0);
	    iy1 = (int)((y_max - cy1)*ysc+wy0);
	    ix2 = (int)((cx3-x_min)*xscale+wx0);
	    iy2 = (int)((y_max - cy3)*ysc+wy0);
	   

	    if ( !ps_modeon)     
	      XDrawLine(mydisplay, redraw_win, redraw_gc,
			ix1,iy1,ix2,iy2
			);
	    else
	      ps_line(ix1,iy1,ix2,iy2);

	    /*
	 XDrawLine(mydisplay, redraw_win, redraw_gc,
		    ix1,iy1,ix2,iy2
		    );
		    */
	  };

	if ( u2 && u3 )
	  { 
	    ix1 = (int)((cx2-x_min)*xscale+wx0);
	    iy1 = (int)((y_max - cy2)*ysc+wy0);
	    ix2 = (int)((cx3-x_min)*xscale+wx0);
	    iy2 = (int)((y_max - cy3)*ysc+wy0);


	    if ( !ps_modeon)     
	      XDrawLine(mydisplay, redraw_win, redraw_gc,
			ix1,iy1,ix2,iy2
			);
	    else
	      ps_line(ix1,iy1,ix2,iy2);

	    /*

	 XDrawLine(mydisplay, redraw_win, redraw_gc,
		    ix1,iy1,ix2,iy2
		    );
		    */
	  };

    }
	return;
}
/*--------------------------------------------------------------------*/
/* FindPoint  finds a coordinates (x,y) for z which is between z1,z2  */
/*           =0 - no z between                                        */
/*           =1 - z is between z1,z2                                  */
/*           =2 - z is at z1=z2                                       */
/*--------------------------------------------------------------------*/       
int FindPoint( float x1, float y1, float z1,
	       float x2, float y2, float z2,
	       float z,
	       float *xp, float *yp   )
{
  float coef;

  if ( z1 == z )
    {
      *xp=x1; *yp=y1;
      return 2;
    }
  else if ( z2 == z )
    {
      *xp=x2; *yp=y2;
      return 2;
    }
  else if ( (z1 < z) && ( z < z2))
    {
      coef = (z - z1)/(z2 -z1 );
      *xp = (x2 -x1)*coef +x1;
      *yp = (y2 -y1)*coef +y1;
      return 1;
    }
  else if ( (z1 >  z) && ( z > z2))
    {
      coef = (z - z1)/(z2 -z1 );
      *xp = (x2 -x1)*coef +x1;
      *yp = (y2 -y1)*coef +y1;
      return 1;
    }
 else return 0;
}

/*-----------------------------------------------------------------------------
|	draw_vector for rectangle blocks
!	* 'mod' specifies color
-----------------------------------------------------------------------------*/
void draw_v(int mr, int mz, int x_off, int y_off, int q1_off, int q2_off,
	    float *lmax, int *npoints, int *mpoints,
	    float l_scale, int mod, int density)
{
  int ir, iz;

  float *xbuf0,*xbuf,*ybuf0,*ybuf;
  float *v10,*v20,*v1,*v2;
  FloatSc ysc;
  float l;
  float vx,vy;
  int ix1,iy1;
  XColor rColor;
  unsigned long col;
  int ijstep;
  int nij;
  char text[80];
  int rcount, zcount;
  XRectangle clip_area;

  /*printf("draw_v, density %d\n", density);*/

  xbuf0 = buf + x_off;
  ybuf0 = buf + y_off;
  xbuf = xbuf0;
  ybuf = ybuf0;
  v10  = buf + q1_off;
  v20  = buf + q2_off;
  ijstep = *npoints;
  *npoints =0;

  if (mod)
    {
      if ( mod ==5 ) /* white color for vector grid */
	XSetForeground(mydisplay, redraw_gc, WhitePixel(mydisplay, myscreen));

      else {
	XParseColor( mydisplay,cmap, "yellow",&rColor);
	XAllocColor( mydisplay,cmap,&rColor);
	col = rColor.pixel;
	XSetForeground (mydisplay, redraw_gc, col);
      }
    }

  if( yscale<0.0 ) ysc=-yscale;

  if( mod == 2) /* single vector plot */
    {
      xbuf = xbuf0 + mr;
      ybuf = ybuf0 + mr;
      v1    = v10  + mr;
      v2    = v20  + mr;

      if (show_values)
	{
	  sprintf(text," x=%g, y=%g \n", *xbuf,*ybuf);
	  xprintf (text);
	  sprintf(text," Vx=%g, Vy=%g \n", *v1, *v2);
	  xprintf (text);
	  show_values =0;
	}

      ix1 =  (int)((*xbuf-x_min)*xscale+wx0);
      iy1 =  (int)((y_max-*ybuf)*ysc +wy0);

      /*printf( " x=%f, y=%f, vx =%f, vy=%f \n", *xbuf,*ybuf,*v1, *v2 );*/
      vector( ix1,iy1, *v1, *v2, *lmax, l_scale,0);
    }

  else
    {
      nij = mr*(*mpoints); 
      for (iz=zcount=0; iz<mz; iz+=(*mpoints),v10+= nij,v20+=nij,xbuf0+=nij,
	     ybuf0+=nij, zcount++)
	{
	  for (ir=rcount=0; ir < mr; ir+=ijstep, rcount++)
	    {
	      xbuf = xbuf0 + ir;
	      ybuf = ybuf0 + ir;
	      v1    = v10  +ir;
	      v2    = v20  +ir;
	   
	      if ( *xbuf < x_min || *xbuf > x_max ||
		   *ybuf < y_min || *ybuf > y_max ) continue;
	   
	      vx = *v1;
	      vy = *v2;
	      l = vx*vx + vy*vy;
	      if ( !mod) /* mod = 0  then lmax calc */
		{
		  *npoints=(*npoints)+1;
		  *lmax = ( l > (*lmax))? l:(*lmax);
		  continue;
		}
	      else /* vect plot */
		{
		  ix1 =  (int)((*xbuf-x_min)*xscale+wx0);
		  iy1 =  (int)((y_max-*ybuf)*ysc +wy0);

		  if (density > 1) {
		    if (zcount%density != 0) continue;
		    if (rcount%density != 0) continue;
		  }
		  vector( ix1,iy1, vx,vy,*lmax,l_scale,mod-1);
		}
	    }
	}
    }
}

/*-------------------------------------------------------------------------*/
/* vector
/*-------------------------------------------------------------------------*/
void vector(    int ix1, int iy1, float vx, float vy, 
	   	float lmax, float l_scale, int mod
	    )
{
  float l;
  int ix2,iy2;
  int ix3, iy3,ix4,iy4;
  int ang_coef = 30;
  float sn,cs;
  float alpha, gam;
  float dl;
  float temp;
  float ar_length;
  /*printf("Vector plot \n");*/

    alpha = 3.1415926535897931/ang_coef;  
    sn = sin( alpha);
    cs = cos( alpha);
    ar_length = 4;
    /* if equal length then l =0.5 lmax */
    l = (mod <= 1)?(vx*vx +vy*vy):(0.5*lmax);
    temp = sqrt(l/(lmax))*l_scale;
     
    if(vx !=0 && vy !=0)
      {
	gam=atan((vy*vx>0)?(vy/vx):-(vy/vx));
	ix2 = ix1 + (int)(((vx>0)?1:-1)*temp *cos(gam));
	iy2 = iy1 + (int)(((vy>0)?-1:1)*temp *sin(gam));
      }
    else if (vy == 0)
      {
	iy2 = iy1;
	ix2 = ix1+(int)(((vx>0)?1:-1)*temp);
      }
    else if (vx ==0 )
      {
	ix2 =ix1;
	iy2 = iy1 + (int)((vy>0?-1:1) *temp);
      }
    else
      {
	ix2 =ix1;
	iy2 =iy1;
      }
    /* --------set scaled color---------------------   */
    if ( mod == 0 || mod == 2)
      {
	l = vx*vx + vy *vy;
	/*   Color scaled by length            */
	/*	temp = sqrt(l/(lmax))*l_scale; */
	temp = l/lmax;
	set_color_table((int) (temp *100 ), 100);
      }
     /*-----------------------------------------------   */
     if( !ps_modeon)
       XDrawLine(mydisplay, redraw_win, redraw_gc,
		 ix1,iy1,ix2,iy2
		 );
     else
       ps_line(ix1,iy1,ix2,iy2);

     /* Draw an arrow */
     /* Find (x3,y3),(x4,y4) */
     if ( ix2!=ix1)     
       {
	 ix3 = (int)(ix2 -(ix2-ix1)/ar_length);
	 iy3 = (int)((ix3 -ix1)*(iy2 - iy1)/(ix2-ix1))+iy1;
       }
     else
       {
	 ix3 = ix1;
       iy3 = (int)(iy2 -(iy2-iy1)/ar_length);
       }
     
     ix1 = (int)(cs*(ix3-ix2)-sn*(iy3-iy2) + ix2);
     iy1 = (int)(sn*(ix3-ix2)+cs*(iy3-iy2) + iy2);

     if ( !ps_modeon)
     XDrawLine(mydisplay, redraw_win, redraw_gc,
	       ix2,iy2,
	       ix1,iy1);
     else
       ps_line(ix1,iy1,ix2,iy2);
     
     ix1 = (int)(cs*(ix3-ix2)+sn*(iy3-iy2) + ix2);
     iy1 = (int)(sn*(ix2-ix3)+cs*(iy3-iy2) + iy2);
     
     if ( !ps_modeon)     
       XDrawLine(mydisplay, redraw_win, redraw_gc,
		 ix2,iy2,
		 ix1,iy1);
     else
       ps_line(ix1,iy1,ix2,iy2);
     
     

}/* vect plot */
