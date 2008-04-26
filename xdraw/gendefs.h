/******************************************************************************
**  NAME      GENDEFS.H
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      structures LOOP and CURVE_SET xdraw.c, xinit.c
**
**  Copyright (c) GlassWare 1994.  All rights reserved.
******************************************************************************/
#define Int int		// these were 'long' if not UNIX
#define Long int
typedef int byte;

#define FZ (float)0
#define FH (float).5

#define addfloat(bx,ix0) bx+ix0;
#define float_to_int(fx1) (int) ( (fx1>=FZ) ? (fx1+FH) : (fx1-FH) )

#define dprintf printf
#define NOX (mydisplay==NULL)

#if defined(MOTIF)
//#elif defined(MSWINDOWS)
#else
#define XCALLS
//#ifndef UNIX
//#define DOS
//#endif
#endif

#ifndef XPRINTF			/* dialog fns in setcolor.c */
extern void xprintf(const char *, ...);
#else
extern void xprintf(const char *, int, int, int, int, int, int, int, int);
#endif				/* end XPRINTF */

extern int xinput(char *, int *);

