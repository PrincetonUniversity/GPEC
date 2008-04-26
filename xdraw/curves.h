/******************************************************************************
**  NAME      CURVES.H
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      structures LOOP and CURVE_SET xdraw.c, xinit.c
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/

typedef struct {		/*----- Loop structure */
    char ltype;			/* I,N,X,A,V,H */
    char use_sep;		/* if not '1', use nnode[] instead of count */
    int hcount;			/* # variables as header to loop */
    int ih0;			/* index of first header variable (.in uses) */
    int count;			/* # non-header elements in loop */
    long sep;			/* count to next index at this level */    
    int current;		/* temp: current (fixed) value */
    float *extra;
    char *labels;		/* individual family member names from .in */
    } LOOP;

typedef struct {		/*----- struct in CURVE_SET for var */
  int index;			/* index from draw.in; if 'i', add 0x1000 */
  int hloop;			/* it's a var in a header; index of loop */
  long off0, dstep, dfaml;
  float min,max;
  float avg,rms;
  char *label;
  int nitems;
  XTextItem  ml_text[10];
} CVAR;

typedef struct {
  int ncount;
  long nsep;
} NODE;

struct LABEL {
  int color, boxed;
  float x,y;
  char *text;
  struct LABEL *next;	/* We REALLY want (struct LABEL *), but no compile */
  int nitems;
  XTextItem ml_text[10];

  } ;

#define LBL_WHITE -1
#define LBL_ALIGNL -2
#define LBL_ALIGNC -3
#define LBL_ALIGNR -4

#define MAXLOOP 16

//------ Here's CURVE_SET
typedef struct {		/* Info about what's drawn in each window */
  char gtype;			/* N,G,g,C,S,D,B,K */
  char which;			/* '1' = 1st in the window, '2' for 2nd, 3rd etc. */
  Window window;
  int flags;			/* MARKERS,COLORS,ASPECT,FILL, etc */
  int use_log;
  int lstep;			/* index, offset in loop[], for step */
  int lfaml;			/* index, offset in loop[], for family */
  int fskip;
  int param[MAXLOOP];		/* >=0 ==>fixed value, -1=step, -2=family */
  CVAR ix,iy;
  CVAR iz;			/* for contours only, corresponds to lfaml */
  char *title;
  char *subtitle;
  struct LABEL *label;
  int mcount;			/* -1:markers all +; 0:cycle; n>0:n no line*/
  float force;			/* force one specific value for contour */
  float forcemin;		/* G: forcemin, forcemax used for ymin, ymax */
  float forcemax;
  float exclmin;
  float exclmax;
  float xcutoff;		/* N: lower limits in x and y for log */
  float ycutoff;
  int use_cutoff;

  int icurve[9];
  int ncurve, hzcurve;		/* c:#contour lines; G:1st,2nd family members */
  int multiDots;

  int f_block,l_block,i_block;	/* c: # of first and last block ; if =0 - all*/
  int iqty, iqty2;
  int itime, itime_rel, itime_abs;
  int i_single,j_single;
  int vDensity;

  /* Window description */
    int view_number; 
    int axis_width,curves_width; /* width of axis and curves */
  } CURVE_SET;

typedef struct {	/*----- struct for each variable */
  int index;		/* index assigned in draw.in; if 'i', add 0x1000 */
  char *name;
  int nitems;
  XTextItem ml_text[10];
} VARIABLE_ID;

typedef struct {
  char  dataType;	// T if topology start, F = func values for timestep
  int   nx, ny, nt;	// current topology's #coordinates, #timesteps
  float xmin, xmax;
  float ymin, ymax;
  int   iTime_rel;	// time relative to this topology (abs time is i_m2)
  long  iTop;		// offset in xx.bin of this topology
  long  iThis;		// offset in xx.bin of this timestep
} M2_ARRAY;

// array m2_array contains pointers to all topologies (header and timesteps)
// in the xx.bin file.  But array buf only contains one topology 
// (top header + assoc timesteps) or one timestep (top header + timestep).
// n_m2 is total # of timesteps
// i_m2 is pointer in m2_array, i_m2_prev is previous pointer

/*----------------------------------------------------------------*/
/* Structures for Multiblocks */
/*----------------------------------------------------------------*/
typedef struct {	/*----- struct for each XY block */
  int mr;
  int mz;
  int nt;
  int mr_prv;
  int mz_prv;
  int nt_prv;
  long offset;         /* in buf */
  float xmin;
  float xmax;
  float ymin;
  float ymax;
} BLOCK_XY;

typedef struct {	/*----- struct for each F block */
  float fmin;
  float fmax;
  long offset;
} BLOCK_F;

/*-----------------------------------------------------------------*/

/* CURVE_SET
 * type: N=unused, G=Graph, C=Contour,
 *       S=Solid, D=Dots on contour, B=Big
 */
#define MARKERS 0x1		/* flag bits */
#define COLORS	0x2		/* Combined graphs in window change color */
#define ASPECT	0x4		/* Graphs in window use data aspect ratio */
#define DOTS	0x8
#define POLAR   0x10
#define FORCED	0x20
#define FILL    0x40
#define ANIMATE 0x80
#define REVERSE 0x100
#define SINGLE  0x200
#define LABELF	0x400
#define SPLINE	0x800
#define E_VIS	0x1000		/* Extrema from visible family members only */
#define E_IGNORE 0x2000		/* Extrema: use / ignore -E option if any */
/*DD*/
#define FILES   0x4000		/* Extrema: use / ignore -E option if any */
/*DD*/

#define NNAME 50
#define LNAME 80
#define MAXWINDOW 9
#define NCHAR (LNAME*(NNAME+MAXWINDOW+1))	/* size of namebuf[] */
#define NCS 90					/* size of curveset[] */

#define LOOPBIT 0x1000
#define INDBITS 0x0FFF

/* Window attributes constant */
#define W_VIS 0x1
#define W_MENU 0x2

#define W_POZ 0x4
#define W_GEO 0x8
/*DD*/
/* Files */
#define LLNAME (LNAME*10)

typedef struct { 
  char name[LLNAME];
  int curve_begin;
  int curve_end;
} BINFILE; 

/*DD*/
#define DATA_LIMITS 1.0E-9	/* was 1.0e-10 */
