/******************************************************************************
**  NAME        XTOOLS
**  AUTHOR        Sheryl M. Glasser
**
**  DESCRIPTION
**     Supporting functions for xwindows
**
**  Copyright (c) CounterPoint Graphics 1993.  All rights reserved.
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#if !defined(__MACH__)
#include <malloc.h>
#endif
#include <string.h>
#include <fcntl.h>
#include <math.h>

#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>

#ifdef MOTIF
#include <Xm/Xm.h>
#include <Xm/DrawingA.h>
#endif

#include "gendefs.h"
#include "curves.h"		/* (xcontour.h uses CURVE_SET) */
#include "xtools.h"
#include "xdraw.h"
#include "xcontour.h"
#include "xinit.h"
#include "xedit.h"
#include "setcolor.h"
#include "ps.h"

#ifdef MOTIF
Widget topLevel;
XtAppContext app;
extern Widget main_window(Widget toplevel, int nwindow, Widget *work_area,
			  Dimension width, Dimension height);
#endif
/*=============================================================================
**            DEFINITIONS AND VARIABLES
**===========================================================================*/

VIEW view[14];			// windows normally 3 x 3 or 4 x 3
extern Display *mydisplay;
extern int myscreen;

static Window parent;
static XSizeHints myhint;
unsigned long myforeground, mybackground;


Window root;
Font font;
Font font_greek, font_latin;

XFontStruct *font_struct;
XFontStruct *font_struct1;
extern int font_height;
extern int font_height1;
Colormap cmap;

extern int ps_only;
extern int background;

int nwindow = 0;
static char *maintitle;
static int xsep = 20, ysep = 45;
unsigned int xscreen=0, yscreen=0;		/* set by XGetGeometry() */
int nRow, nCol, base_dx, base_dy, maxWindow;
static int iWin=0;

static int argc;
static char **argv;

extern char datafile[];
extern CURVE_SET curveset[];

#define FPT1 12
#define FPT2 18

/*=============================================================================
**                    DISPLAY, WINDOW, GC INITIALIZATION
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	give_command_args -- save argc, argv for Motif
-----------------------------------------------------------------------------*/
void give_command_args(int argcc, char **argvv)
{
  argc = argcc;
  argv = argvv;
}

/*-----------------------------------------------------------------------------
|   opendisplay. initializes display_in, etc.
-----------------------------------------------------------------------------*/
opendisplay(char *title_in)
{
  int x0, y0;
  unsigned int bwidth, depth;
  char finalname[100];

  /*------- copy input data */
  maintitle = title_in;

  if(ps_only) {
    mydisplay = NULL;
    printf("Printing into .ps files. \n"); 
    goto NO_XTERM;
  }

  /*------- initialization */

#if defined(XCALLS)
  mydisplay = XOpenDisplay("");
#elif defined(MOTIF)
  XtSetLanguageProc(NULL, NULL, NULL);
  topLevel = XtVaAppInitialize(&app, "XMenu", NULL,0,
			       &argc, argv, NULL, NULL);
  mydisplay = XtDisplay(topLevel);
#endif

  if (NOX) {
    printf("WARNING -- XTerm display not found.\n");
    goto NO_XTERM;
  }

  myscreen = DefaultScreen(mydisplay);
  cmap = DefaultColormap(mydisplay, myscreen);

  /*------- set up font */

  if (!getgoodfont(finalname, FPT1, FPT2,	/* MSWINDOWS: dummy */
         "times-medium-r-normal")) return 0;
  font_struct = XLoadQueryFont(mydisplay, finalname);
  font = (*font_struct).fid;
  font_latin = font;
  if (!getgoodfont(finalname, FPT1, FPT2,	/* MSWINDOWS: dummy */
         "symbol-medium-r-normal")) return 0;
  font_struct1 = XLoadQueryFont(mydisplay, finalname);
  font_greek = (*font_struct1).fid;

  /*------- default pixel values */

  mybackground = !background ? BlackPixel(mydisplay, myscreen): 
    WhitePixel(mydisplay, myscreen);
  myforeground = !background? WhitePixel(mydisplay, myscreen) :
    BlackPixel(mydisplay, myscreen);

  //------ Geometry (xscreen, yscreen)

  root = parent = DefaultRootWindow(mydisplay);
  XGetGeometry(mydisplay, root, &root, &x0, &y0,
	       &xscreen, &yscreen, &bwidth, &depth);

  //yscreen = 800;		// test 4x3 screen on 1024 x 1280 screen
  nRow = nCol = 3;
  if ((float)yscreen / (float)xscreen < 0.7) nCol = 4;

NO_XTERM:
  //ncol = nCol_input;
  maxWindow = nRow * nCol;
  if (nCol == 1) {
    base_dx = xscreen * 2 / 3;
    base_dy = yscreen * 2 / 3;
  }
  else {
    base_dx = (xscreen + 1 - (nCol + 1) * xsep) / nCol;
    base_dy = (yscreen + 1 - (nRow + 1) * ysep) / nRow;
  }
  return (!NOX);
}

/*-----------------------------------------------------------------------------
|	get_screensize
-----------------------------------------------------------------------------*/
void get_screensize(unsigned int *dx, unsigned int *dy)
{
  *dx = xscreen; *dy = yscreen;
}

/*-----------------------------------------------------------------------------
|	closedisplay
-----------------------------------------------------------------------------*/
void closedisplay()
{
#ifndef MOTIF
  int i;
  for (i = nwindow - 1; i >= 0; i--)
    {
      XDestroyWindow(mydisplay, view[i].window);
      XFreeGC(mydisplay, view[i].gc);
    }
  XCloseDisplay(mydisplay);
#endif
}

/*-----------------------------------------------------------------------------
|   makewindow. creates window and gc
-----------------------------------------------------------------------------*/
Window makewindow(char *title, int has_menu)
{
  int argc = 0;		/* for XSetStandardProperties */
  int xwindow, ywindow;
  char **argv, *p;
  unsigned long fore, back;
  Window topwin, win;
  GC gc;
  VIEW *v;

  XSetWindowAttributes attrs;

#ifdef MOTIF
  Widget shell, area;
  Widget area1;
  char text[8];
  extern void redraw_m(Widget widget, 
		       XtPointer client_data, XtPointer call_data);
#endif
  extern GC dialog_gc;

  /*------- default program-specified window position and size */

  nextwindow(&xwindow, &ywindow);
  iWin++;
  myhint.width =  base_dx;
  myhint.height = base_dy;
  myhint.x = xwindow;	/* windows are positioned... */
  myhint.y = ywindow;	/* ..incrementally */
  topwin = parent;
  back = (has_menu == WIN_MNGR)? myforeground : mybackground;
  fore = (has_menu == WIN_MNGR)? mybackground : myforeground;
  myhint.flags = PPosition | PSize;

  /*-------- window creation */

  v = view + has_menu;
  p = title;
  if (!NOX)
    {
#if defined(XCALLS)
      argv = NULL;		/* (unused) */
      XWarpPointer(mydisplay, None, root, 0,0,0,0, myhint.x, myhint.y);
      win = v->window =
	XCreateSimpleWindow(mydisplay, topwin, myhint.x, myhint.y,
			    myhint.width, myhint.height, 5, fore, back);
      XSetStandardProperties(mydisplay, win,
			   p, p, None,	/* None,icon titles; icon pixmap */
			   argv, argc, &myhint);
#elif defined(MOTIF)
      shell =
	XtVaAppCreateShell(NULL, "XDraw", topLevelShellWidgetClass,
			   mydisplay, XmNx, myhint.x, XmNy, myhint.y,
			   XmNwidth, myhint.width, XmNheight, myhint.height,
			   XmNtitle, title, NULL);
      area1=main_window(shell,nwindow,&area,myhint.width,myhint.height);
      sprintf(text, "Win%d", nwindow);

      enableWidget(nwindow, shell, area1, area,text);

      win = v->window = XtWindowOfObject(area1);
      v->work_area =XtWindowOfObject(area);
      
      /*
      printf(" MAKEWINDOW\n");
      printf(" Window (view)= %d, mainwindow =%x, work_area = %x\n",
	nwindow, win, v->work_area);
      test_window ("Mainwindow", win, area1);
      test_window ("Work_area", v->work_area, area);
      printf(" End of test MAKEWINDOW \n");
      */
      /* 
     v->work_area--;
     */
      win = v->window =v->work_area;
#endif
    }
  nwindow++;
  if (NOX)
    return ((Window) 0);

  /*---------- gc creation and initialization */

  gc = v->gc = XCreateGC(mydisplay, win, (long) 0, 0);
  if (has_menu==WIN_MNGR) dialog_gc = gc;
  XSetBackground(mydisplay, gc, back);	/* (MSWINDOWS these unused) */
  XSetForeground(mydisplay, gc, fore);
  XSetFont(mydisplay, gc, font);

  /*---------- input event selection */

  set_winmask(win);

#ifndef MOTIF
  /*---------- generate expose event */
  XMapRaised(mydisplay, win);
#endif

  attrs.bit_gravity = ForgetGravity;
  XChangeWindowAttributes (mydisplay, win,
			   CWBitGravity, &attrs);


  return (win);
  /*return (v->work_area);*/
}
/*-----------------------------------------------------------------------------
|	nextwindow *** ERROR here
-----------------------------------------------------------------------------*/
void nextwindow(int *xw, int *yw)
{
  int i,n;
  n = nRow * nCol;
  i = (iWin >= n) ? n-1 : iWin;
  *xw = (i % nCol) * (base_dx + xsep);
  *yw = (i / nCol) * (base_dy + ysep);
}

/*-----------------------------------------------------------------------------
|   getview, getwindow
-----------------------------------------------------------------------------*/
int getview(Window w)
{
  VIEW *v;
  int i;
  for (i = 0, v = view; i < nwindow && v->window != w; i++, v++);
  return (i);
}

Window getwindow(int i)
{
  return((view+i)->window);
}

unsigned long white()
{
  return (myforeground);
}

/*-----------------------------------------------------------------------------
|   get_properties -- called from redraw
-----------------------------------------------------------------------------*/
void get_properties(Window w, unsigned int *pw, unsigned int *ph,
		    unsigned int *pd)
{
  Window root;
  int x0, y0;
  unsigned int bwidth, status;

  if (NOX)
    {
      *pw = myhint.width;
      *ph = myhint.height;
     }
  else
    status = XGetGeometry(mydisplay, w, &root,&x0, &y0, pw, ph, &bwidth, pd);
}

/*-----------------------------------------------------------------------------
|	textwidth
-----------------------------------------------------------------------------*/
int textwidth(char *string, int n)
{
  if (string==0) return 0;
  if (font_struct==0) return strlen(string);
  return (XTextWidth(font_struct, string, n));
}
















