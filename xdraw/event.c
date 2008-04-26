/******************************************************************************
**  NAME	EVENT.C
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
#include <fcntl.h>
#include <math.h>

#include <unistd.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#define zprintf if (zdebug) printf

//extern void xwin_to_integer(float, float, int*, int*);
extern void trigger_next_window();
int tiny_zoom(int inc, float fxavg, float *pfx1, float *pfx2, char which);

extern int dialogwindow;
extern Window dialog_win;
int input_enabled = 0;
int time_stride = 1;
int fam_stride = 1;

extern int debug_m;
extern int debug_scale;
static int debug_scale2=0;	/* print result of tiny_zoom */
static int debug_scale3=0;	/* print values for next special zoom */
extern float *buf;

#ifdef MOTIF
#include <Xm/Xm.h>

#include <Xm/MainW.h>
#include <Xm/DrawingA.h>
#include <Xm/Label.h>
#include <Xm/RowColumn.h>
#include <Xm/CascadeBG.h>
#include <Xm/CascadeB.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/FileSB.h>
#include <Xm/SeparatoG.h>
#include <Xm/RowColumn.h>
#include <Xm/Command.h>
#include <Xm/Text.h>

#endif

#include "gendefs.h"
#include "curves.h"		/* (xcontour.h uses CURVE_SET) */
#include "xtools.h"
#include "xdraw.h"
#include "xcontour.h"
#include "xinit.h"
#include "xedit.h"
#include "setcolor.h"

extern void postscript(char *);	/* near heap problem, try elim .h's*/
extern float get_world(int, int, XRectangle, float, float, float, float,
		       float*, float*);
extern int nearest_M_value(float, float, CURVE_SET *, char *, double *);

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

static XEvent myevent;
static KeySym mykey;
static Window redraw_win, event_win=-1;
static GC redraw_gc, event_gc;

extern int nwindow;
extern VIEW view[];
extern int font_height;
extern unsigned long myforeground, mybackground;
extern CURVE_SET curveset[];
extern LOOP loop[];
extern int exitflag, redrawflag, titledrawn;
extern unsigned int xscreen, yscreen;
int xhair_type = XHAIR_XOR;
extern int ftype;
extern int ntime;
extern float data_limits;
int special_debug = 0;		// set via '!'
int report_zoom = 0;		// set via '!'

extern int drawing_all_windows;

/*=============================================================================
**                      VARIABLES FOR ZOOM
**===========================================================================*/
int zoom_on = 0, coord_on = 0;
float zoomfac = (float).05, xfromcoord = 0, yfromcoord=0;
static int zoom_count = 0, coord_count = 0;
int zoomi=0, coordi=0;
static int rubber = 1;
Window zoomwin=0;
static GC zoomgc;
static float fxbox1, fybox1, fdxbox, fdybox, zoomdx, zoomdy;
static int xcurs, ycurs;
static int cursor_defined = 0;
static Cursor curs1, curs2;
static int special_zoom_flag = 0;
void special_zoom(int option);

#define CURS1 XC_draft_large
#define CURS2 XC_top_left_arrow
/*#define CURS1 XC_top_left_corner*/
/*#define CURS2 XC_bottom_right_corner*/

#ifdef USE_MENU
#define PRESSMASK (ButtonPressMask | ButtonReleaseMask | KeyPressMask)
#define WHEREMASK (LeaveWindowMask | EnterWindowMask)

#else			/* ...of USE_MENU */
#define PRESSMASK (ButtonPressMask | KeyPressMask)
#define WHEREMASK LeaveWindowMask
#endif			/* ...of USE_MENU */

#define MOTIONMASK (PointerMotionMask | PointerMotionHintMask)
#define EVENTMASK3 (PRESSMASK | WHEREMASK | ExposureMask) 
/*#define EVENTMASK3 (PRESSMASK | WHEREMASK | ExposureMask|StructureNotifyMask)*/
#define EVENTMASK4 (EVENTMASK3 | MOTIONMASK)

static int tell_event=0;

/*=============================================================================
**		FUNCTIONS FOR MENU VERSION
**===========================================================================*/
#ifndef USE_MENU
void set_menu_visibility(int i) { i++; }
void menu_setwindows(int i) { i++; }

#else			/* .....of #ifndef USE_MENU, so USE_MENU is defined */
/*-----------------------------------------------------------------------------
|	init_menuXlibrary
-----------------------------------------------------------------------------*/
void init_menuXlibrary()
{
  unsigned long back;
  back = glx_index(get_menubkgd());
  init_glx(mydisplay, root, myforeground, back,
	   font, font_struct, &myevent, xscreen, yscreen);
}
#endif			/* ......of USE_MENU undefined, defined */


/*-----------------------------------------------------------------------------
|	init_menus, clear_if_needed: non-Motif
-----------------------------------------------------------------------------*/
#ifndef MOTIF
#ifndef USE_MENU
  void init_menus() {}
#endif

void clear_if_needed(Window w, GC gc, unsigned int dx, unsigned int dy)
		    {;}


/*-----------------------------------------------------------------------------
|	init_menus, clear_if_needed: Motif
-----------------------------------------------------------------------------*/
#else			/* .....MOTIF */
void init_menus()
{
  int x, y;
  nextwindow(&x, &y);
  execute_motif(x, y, nwindow, mydisplay);
}

void clear_if_needed(Window w, GC gc, unsigned int dx, unsigned int dy)
{
  xrectangle(dx,dy,BlackPixel(mydisplay, myscreen),
	     WhitePixel(mydisplay, myscreen));
}

void test_window(char *s, Window win, Widget wid)
{
  unsigned int dx, dy, b, d, nchildren;
  int x0, y0;
  Window root, parent, *children;
  Position x,y;
  Dimension w, h;
  int i;

  XGetGeometry(mydisplay, win, &root, &x0, &y0, &dx, &dy, &b, &d);
  printf ("\n%s: Window %lx, \n its root %lx \n Its coordinates x,y = %d %d,\n size %d %d border %d\n",
	  s, win, root, x0, y0, dx, dy,b);
  if (wid != 0)
    {
      XtVaGetValues(wid, XmNx, &x, XmNy, &y, XmNwidth, &w, XmNheight, &h, NULL);
      printf ("%s: Widget %lx, Window %lx has x,y = %d %d, size %d %d\n",
	      s, wid, win, x, y, w, h);
    }
  XQueryTree(mydisplay, win, &root, &parent, &children, &nchildren);
  printf("%s: window=%lx, root=%lx, parent=%lx, nchildren=%d\n",
	 s, win, root, parent, nchildren);
  XGetGeometry(mydisplay, parent, &root, &x0, &y0, &dx, &dy, &b, &d);
  printf ("%s: Window %lx has x,y = %d %d, size %d %d\n",
	  s, parent, x0, y0, dx, dy);
  for(i=0; i<nchildren;i++)
    {
      XGetGeometry(mydisplay, children[i], &root, &x0, &y0, &dx, &dy, &b, &d);
      printf ("%s: Children[%d] %lx has root %lx \n x,y = %d %d,\n size %d %d\n border %d\n",
	  s, i,children[i],root, x0, y0, dx, dy,b);
      
    }

  XFree(children);
}

/*-----------------------------------------------------------------------------
|	set_expose: for calls to get_expose_info
-----------------------------------------------------------------------------*/
Window set_expose(int i)
{
  myevent.xexpose.display = mydisplay;
  /*
  myevent.xexpose.window = view[i].window;
  return view[i].window;
  */
  myevent.xexpose.window = view[i].work_area;
  view[i].window = view[i].work_area;
  return view[i].work_area;
}
#endif			/* ......of MOTIF undefined, defined */


/*=============================================================================
**                        EVENTS
**===========================================================================*/

/*-----------------------------------------------------------------------------
|	get_expose_info
-----------------------------------------------------------------------------*/
void get_expose_info(Display **d, Window *w, GC *gc, int *fonthp)
{
  int iwindow;

  if (NOX) return;
  mydisplay = *d = myevent.xexpose.display;
  redraw_win = *w = myevent.xexpose.window;
  iwindow = getview(redraw_win);	/* get which window is redrawing */
  if (iwindow < nwindow)
    redraw_gc = *gc = view[iwindow].gc;
  *fonthp = font_height;
}

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
void set_winmask(Window win)  
{
#ifndef MOTIF
  XSelectInput(mydisplay, win, EVENTMASK3);
#endif
}
/*-----------------------------------------------------------------------------
|	readkey, parsekey
-----------------------------------------------------------------------------*/
int parsekey()
{
  int i, ic;
  char text[20];
  text[0] = text[1] = 0;
  i = XLookupString(&myevent.xkey, text, 10, &mykey, 0);
  if (i < 1) ic=0;
  /*else if (i==0) ic = *(int *)&mykey; (extra line from xtools.c--ok? */
  else ic = *(int *)text;
  return(ic);
}

int readkey()
{
  zprintf("In readkey()\n");
  for(;;)
    {
      XNextEvent(mydisplay, &myevent);	/* read the next event */
      if (myevent.type == KeyPress) break;
    }
  return(parsekey());
}

/*-----------------------------------------------------------------------------
|   xmessage -- display or erase a message in the current window during redraw
-----------------------------------------------------------------------------*/
void xmessage(int x, int y, char *s)
{
  if (NOX)
    return;
  set_xor(redraw_gc, 1);
  XSetForeground(mydisplay, redraw_gc, myforeground);
  XDrawString(mydisplay, redraw_win, redraw_gc, x, y, s, strlen(s));
  set_xor(redraw_gc, 0);
}

void xrectangle(unsigned int dx, unsigned int dy, int bkgd, int fore)
{
  XSetForeground(mydisplay, redraw_gc, BlackPixel(mydisplay, myscreen));
  XFillRectangle(mydisplay, redraw_win, redraw_gc, 0, 0, dx, dy);
  XSetForeground(mydisplay, redraw_gc, WhitePixel(mydisplay, myscreen));
}

/*-----------------------------------------------------------------------------
|	print_event
-----------------------------------------------------------------------------*/
void print_event(int etype, int iwin, long event_win)
  {
    static int inited=0;
    XWindowAttributes att;
    if (!inited) { printf("Display=%lx\n",mydisplay); inited=1; }
    printf("------- ");
    if      (etype==MotionNotify)  printf("Motion ");
    else if (etype==Expose)        printf("Expose ");
    else if (etype==ButtonPress)   printf("Press  ");
    else if (etype==ButtonRelease) printf("Release");
    else if (etype==EnterNotify)   printf("Enter  ");
    else if (etype==LeaveNotify)   printf("Leave  ");
    else if (etype==KeyPress)      printf("Key");
    else if (etype==ConfigureNotify) printf("Configure  ");
    else                           printf("Unknown");

    XGetWindowAttributes(mydisplay, event_win, &att);
    printf(" in window %d, ID %lx, masks %lx,%lx\n",
	   iwin, event_win, att.your_event_mask, att.all_event_masks);
    fflush(stdout);
  }

/*-----------------------------------------------------------------------------
|	keystrokes
-----------------------------------------------------------------------------*/
static void keystrokes(void);
static void keystrokes()
{
  xprintf("***** Keystroke Summary *********************************\n");
  xprintf("|  q ......quit\n");
  xprintf("|  p ......print window (create .ps file)\n");
  xprintf("|  t ......display descriptive title of window\n");
  xprintf("|  k ......this message\n");
  xprintf("|  \n");
  xprintf("|  m ......toggle markers: on/off\n");
  xprintf("|  M ......toggle marker size: 2..8\n");
  xprintf("|  a ......toggle aspect ratio: preserved/auto-scale\n");
  xprintf("|  e ......toggle Extrema, y from:  all curves/visible curves/-E values\n");
  xprintf("|  f ......toggle fill: on/off\n");
  xprintf("|  l ......toggle automatic labels on curves (random positions)\n");
  xprintf("|  1 ......draw one curve of family (first,next,...)\n");
  xprintf("|  2,3..9..draw a second, third.., ninth curve of family (first,next,...)\n");
  xprintf("|  0 ......draw all curves of family\n");
  xprintf("|  # ......prompt for one curve of family, or time stride for type M\n");
  xprintf("|  ^.......prompt for stepsize for 1,2\n");
  xprintf("|  >,<.....draw next, previous step of outer loop (time for MG)\n");
  xprintf("|  . , ....< or > for ALL windows\n");
  xprintf("|  @ ......edit mode (partially implemented)\n");
  xprintf("|  \n");
  xprintf("|  z ......zoom (use crosshairs)\n");
  xprintf("|  o ......original, unzoomed size\n");
  xprintf("|  +,- ....zoom up or down slightly\n");
  xprintf("|  (bksp)..restore to previous zoom level (only works once)\n");
  xprintf("|  \n");
  xprintf("|  c ......get-coordinates mode on/off (use crosshairs)\n");
  xprintf("|  s ......get-slope mode on/off (use crosshairs)\n");
  xprintf("|  r ......get-ratio mode on/off (use crosshairs)\n");
  xprintf("|  g ......get-gradient mode on/off (use crosshairs)\n");
  xprintf("|  \n");
  xprintf("|  d,h ....double, halve number of lines in contour plot\n");
  xprintf("|  D,H ....increase, decrease density of vectors in vector plot\n");
  xprintf("|  b ......draw one block of MG contour lines \n");
  xprintf("|  B ......draw all blocks of MG contour lines \n");
  xprintf("|  v ......print value of contour lines\n");
  xprintf("|  (spc) ..print a blank line\n");
  xprintf("********************************************************\n");
}

#ifndef USE_MENU

#ifndef MOTIF
/*-----------------------------------------------------------------------------
|	abort_zoom
-----------------------------------------------------------------------------*/
 int abort_zoom(void);
 int abort_zoom()
{
  int retval;
  retval = 1;
  if      (zoom_on  && event_win != zoomwin) zoom(zoomwin, ZABORT);
  else if (coord_on && event_win != zoomwin) coord(coord_on);
  else retval = 0;

  return retval;
}

/*----------------------------------------------------------------------------
|	event -- for keystroke version
|	* wait for an event, process it, return
-----------------------------------------------------------------------------*/
void event()
{
  char c;
  int iwin, ic, dummy, n, *result, old_input_type, jwin, i;
  int xw, yw, xr, yr, x, use_xsync;
  unsigned int keys_buttons;
  Window rw, cw, oldw;
  static char editing=0;
  static char inited=0;
  static int tell_event = 0;
  static int expose_event_in_process = 0;
  CURVE_SET *cp;
  LOOP *qf;
  char text[200];

  if (NOX) exit(0);
  /*
  printf(" Number of events = %d\n", dummy=XQLength(mydisplay));
  */
  fflush(stdout);
  
  XNextEvent(mydisplay, &myevent);	/* read the next event */
  oldw = event_win;
  event_win = myevent.xany.window;

  iwin = getview(event_win);

  if (tell_event) print_event(myevent.type, getview(event_win), event_win);
  if (expose_event_in_process && myevent.type != KeyPress) return;

  switch (myevent.type)
    {

    case Expose:		/* repaint window on expose events */
      //if (myevent.xexpose.count == 0)
      while (XCheckTypedWindowEvent(mydisplay,event_win,Expose,&myevent));
	{
	  expose_event_in_process = 1;
	  XClearArea(mydisplay, event_win, 0, 0, 0, 0, False);

	  if (event_win == dialog_win) redraw_dialog();
	  else redraw();

	  if (!inited && (dialogwindow==0 || event_win==dialog_win))
	      { xprintf("Press <k> for Keystroke summary.\n"); inited=1; }
	  expose_event_in_process = 0;

	  /*---- Problem - spurious redraws of current expose were occurring */
	  /*     e.g. 'd' for lots of contours */
	  /*     originally used XSync to clear event buffer - failed for eg uncover 9 windows */
	  /*     WARN.  Current logic may not break out of loop */

	  use_xsync = 0;
	  if (use_xsync) {
	    //if (iwin == (nwindow-1)) drawing_all_windows = 0;
	    //if (!drawing_all_windows) XSync(mydisplay, True);
	    // else for (;;) {
	    for (;;) {
	      XPeekEvent(mydisplay, &myevent);	/* if event buffer is empty, waits! */
	      event_win = myevent.xany.window;
	      jwin = getview(event_win);
	      if (myevent.type == Expose && jwin == iwin)
		XNextEvent(mydisplay, &myevent);
	      else break;
	    }
	  }
	}
      break;
            
    case ConfigureNotify:
      while (XCheckTypedWindowEvent(mydisplay,event_win,ConfigureNotify,&myevent));
      XClearArea(mydisplay, event_win, 0, 0, 0, 0, False);

      if (event_win == dialog_win) redraw_dialog();
      else redraw();
      break;
                 
    case MappingNotify:			/* process keyboard mapping changes */
      XRefreshKeyboardMapping(&myevent.xmapping);
      break;

    case ButtonPress:			/* process mouse-button presses */
      if (abort_zoom()) ;		/* click in another win: disable */
      else if (zoom_on)
	zoom(zoomwin, ZGET);
      else if (coord_on)
	coord(8);
      break;

    case LeaveNotify:
      if (zoom_on) zoom(zoomwin, ZIGNORE);	/* Gives message for debug */
      break;

    case MotionNotify:
      if ((zoom_on || coord_on) && event_win == zoomwin)
      {
        XQueryPointer(mydisplay, event_win,
		    &rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);
        crosshair(xcurs, ycurs);
        xcurs = xw;
        ycurs = yw;
        crosshair(xw, yw);
      }
      break;

    case KeyPress:		/* process keyboard input */
      ic = parsekey();
      c = *(char *)&ic;

      //if (expose_event_in_process && c != 'q') return;

      n = get_graph(event_win, 0);
      cp = curveset + n;
      iwin = getview(event_win);
      redrawflag = 0;
      
      abort_zoom();		/* key in another win: disable */

      /*------ Call xinput for new char (if \n, manage input) */
      /*       first call to xinput sets input_enabled = 1 */


      if (input_enabled) {			/* retrieve data entered via xinput */
	*text = c; *(text+1) = 0;		/* make keyPress char into a string */
	old_input_type = input_enabled;

	if (xinput(text, &ic) == 1) {
	  //printf("Input Type %d\n", old_input_type);

	  if (old_input_type == 99) {	/* !: debug stuff */
	    if (ic == 1) { 
	      tell_event = 1 - tell_event; 
	      printf("tell_event toggled to %d\n", tell_event); 
	    }
	    if (ic == 2) special_debug = 1 - special_debug;
	    if (ic == 3) report_zoom = 1 - report_zoom;
	    if (ic == 4) special_zoom(1);
	    if (ic == 5) set_extra_psi(1, 0.0);
	    if (ic == 6) tell_buffer();
	    if (ic == 7) tell_array();
	  }

	  else if (old_input_type == 5) {	/* #: stride for type M */
	    time_stride = (ic > ntime)? 1:ic;
	    //printf("Input 5: %d\n", ic);
	  }

	  //else if (ftype == 6)			/* #: time stride for type M */
	  //time_stride= ( ic > ntime )? 1: ic; 

	  else if (old_input_type == 2) {
	    qf = loop + cp->lfaml;		/* ^: time stride */
	    if (ic == 0) ic=1;
	    if (abs(ic) >= qf->count) 
	      xprintf("Invalid entry.  Value must be in plus-or-minus %d\n", qf->count-1);
	    else fam_stride = ic;
	  }

	  else if (old_input_type == 4) {	/* #: index# for type G eg */
	    cp->flags |= SINGLE;
	    for(i=0; i<9; i++) cp->icurve[i] = -1;
	    cp->icurve[0] = ic - 1;
	    toggle_single(cp, '1');
	  }

	}
      }

      /*------ Enable xinput to enter a value (leaves focus in display window) */

      else if (c == '#') {
	if ( ftype != 6 ) {	/* type G: draw one specified curve */
	  printf("Enter desired index: ");
	  input_enabled = 4;
	}

	else {
	  printf(" Enter desired time stride: ");
	  input_enabled = 5;
	}
      }

      else if (c == '^') {
	if (ftype == 0) {
	  qf = loop + cp->lfaml;
	  printf("Enter stride for '1'..'9' (# curves=%d): ", qf->count);
	  input_enabled = 2;
	}
      }

      else if (c == '!') {
	printf("Enter 1=tell event, 2=special debug, 3=report zoom, 4=special zoom\n");
	printf("      5=force extra psi, 6=tell buffer, 7=tell array: ");
	input_enabled = 99;
      }

      else if (c==' ') xprintf("\n");

      /*------ Quit, Print, Keystrokes, Title */

      else if (c == 'q') exitflag = 1;
      else if (c == 'p' && iwin < nwindow)
	{
	  parse_title(cp, cp->title, text);
	  if (cp->subtitle)
	    {
	      strcat(text, "; ");
	      parse_subtitle(cp, text+strlen(text));
	    }
	  postscript(text);
	}

      else if (c=='k') keystrokes();
      else if (c=='t') print_title(iwin, cp);

      /*--------- Zoom commands */

      else if (zoom_on  && event_win==zoomwin)
      {
        if (c == 'z') zoom(zoomwin, ZABORT);
      }
      else if (coord_on && event_win==zoomwin &&
	       (c!='c' && c!='s' && c!='r' && c!='g'));

      else if (c == 'z') zoom(event_win, ZENABLE);
      else if (c == 'o') zoom(event_win, ZORIGINAL);
      else if (c == 'c') coord(1);
      else if (c == 's') coord(2);
      else if (c == 'r') coord(3);
      else if (c == 'g') coord(4);	/* gradient for generalized contour plots */

      else if (c =='+' || c =='-') addzoom(event_win, c=='+');
      else if (c == '\010' )       addzoom(event_win, 2);	/* backspace */

      /*----------- Display commands */

      else if (c == 'l') toggle_flabel(cp);
      else if (c == 'd' || c=='h') new_ncurve(cp, c);
      else if (c == 'D' || c=='H') new_ncurve(cp, c);
      else if (c == 'I' || c=='i' ||
	       c=='J' || c== 'j') new_nvect(cp, c);
      else if (c == 'v') contour_values(cp);
      else if (c == 'm') toggle_markers(cp, &dummy);
      else if (c == 'M') toggle_markersize(cp);
      else if (c == 'a') toggle_aspect(cp, &dummy);
      else if (c == 'f') toggle_fill(cp);
      else if (c == 'e') toggle_extrema_use(cp);
      else if (c == 'b' || c == 'B') toggle_block(cp,c);

      /*------ Advancing commands */

      //else if (c=='1' || c=='2' || c=='0') toggle_single(cp,c);
      else if (c>='0' && c<='9') toggle_single(cp,c);
      else if (c=='>' || c=='<') toggle_single(cp,c);
      else if (c=='.' || c==',') toggle_single(cp,c);

      /*--------- Editor commands */

      else if (c == '@')
	{
	  if (editcurve(iwin,cp,NULL))
	    XClearArea(mydisplay, event_win, 0, 0, 0, 0, True);
	}

      if (redrawflag)
	XClearArea(mydisplay, event_win, 0, 0, 0, 0, True);
      break;
    }				/* switch (myevent.type) */
}				/* while (done==0) */

#else			/* ...of not MOTIF */
/*-----------------------------------------------------------------------------
|	event -- for MOTIF version
|	* All action in execute_motif(), called from init_menus()
-----------------------------------------------------------------------------*/
void event() {  exitflag = 1; }

#endif			/* ...of not MOTIF */

#else			/* ...USE_MENU is defined */

/*-----------------------------------------------------------------------------
|	
-----------------------------------------------------------------------------*/
int get_eventmask()
{
return EVENTMASK3;
}

Window get_op_window(Window win)
{
  return(win_is_menu(event_win) ? win : event_win);
}

/*-----------------------------------------------------------------------------
|	enable_window
|	* Preparation for Zoom, Coord, Slope
|	* If cursor is in menu, need to warp to selected window.
|	  This will first generate a LeaveNotify, EnterNotify event
|	* THEN enable appropriate function
-----------------------------------------------------------------------------*/
Window enable_window(Window win)
{
  Window root;
  int x0, y0, xw, yw;
  unsigned int dx, dy, bwidth, depth;
  if (win != event_win)
  {
    XGetGeometry(mydisplay, win, &root,
	  &x0, &y0, &dx, &dy, &bwidth, &depth);
    xw = dx / 2;
    yw = dy / 2;
    XWarpPointer(mydisplay, None, win, 0,0,0,0, xw,yw);
    do {
      XNextEvent(mydisplay, &myevent);
    }
    while (myevent.type != EnterNotify);
    event_win = myevent.xany.window;
  }
  return(event_win);
}

/*-----------------------------------------------------------------------------
|	event - for menu version
|		call when enable expose for each window, then call repeatedly
|	Verify: QueryPointer returns same as .xbutton.x etc.
|		value of .button
|		only gives events if in a valid window
|		ButtonPress only when it's a change
-----------------------------------------------------------------------------*/
#define redraw_exposed() redraw()	// had different defs for UNIX, not UNIX

void event()
{
  char c;
  int iwin, i, ic, etype, enter;
  Window oldw;
  short int state;
  int button;

  if (NOX) exit(0);
  zprintf("waiting.."); fflush(stdout);
  XNextEvent(mydisplay, &myevent);	/* read the next event */
  oldw = event_win;
  event_win = myevent.xany.window;
  etype = myevent.type;
  iwin = getview(event_win);		/* index to window (NOT curveset) */
  //print_event(etype,iwin,event_win);	// WARN! could be menu window!

  if (etype==Expose && myevent.xexpose.count==0)
    {
      if (win_is_menu(event_win)) redrawmenu(NULL, 0);
      else redraw_exposed();	/* getproperties() sets redraw_win */
    }

  else if (etype==MappingNotify)	/* process keyboard mapping changes */
    XRefreshKeyboardMapping(&myevent.xmapping);

  else if (etype==ButtonPress || etype==ButtonRelease)
    {
      button = myevent.xbutton.button;
      if (button==Button1 || button!=Button2)
	button = LEFTMOUSE;
      else
	button = RIGHTMOUSE;
      state = (etype==ButtonPress) ? DOWN : UP;
      if (menu_event(button, state)) ;
      else if (state != DOWN) ;
      else if (zoom_on)
        {
	  if (button==LEFTMOUSE && event_win==zoomwin)	/* WARN! can be Abort */
	    zoom(zoomwin, ZGET);
	  else zoom(zoomwin, ZDISABLE);
	}
      else if (button==LEFTMOUSE) set_selected_iwin(iwin);
    }

  else if (etype==EnterNotify || etype==LeaveNotify)
    {
      /*menu_event(INPUTCHANGE, iwin);*/
      enter = (etype==EnterNotify);
      menu_ievent(enter ? event_win : 0L);
    }

  else if (etype==MotionNotify && zoom_on && event_win==zoomwin)
    {
      crosshair(xcurs, ycurs);
      xcurs = myevent.xbutton.x;
      ycurs = myevent.xbutton.y;
      crosshair(xcurs, ycurs);
    }

  else if (etype==KeyPress)		/* process keyboard input */
    {
      ic = parsekey();
      c = *(char *)&ic;
      menu_event(KEYBD, (short int)c);
    }
}

/*-----------------------------------------------------------------------------
|	clear_redraw
-----------------------------------------------------------------------------*/
void clear_redraw(int iwin)
{
  XClearArea(mydisplay, view[iwin].window, 0, 0, 0, 0, True);
}

#endif			/* ...of USE_MENU */

/*=============================================================================
**				ZOOM AND COORD
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	zoom
|	* 1=ZENABLE, 2=ZGET=process button hit while enabled, 3=ZABORT
|	  0=ZDISABLE, -1=ZORIGINAL, -2=ZIGNORE=Leave event during zoom
|	  4=ZPREVIOUS
-----------------------------------------------------------------------------*/
void zoom(Window nextwin, int enable)
{
  /*Window nextwin;*/
  int i, xw, yw, xr, yr;
  unsigned keys_buttons;
  Window rw, cw;
  VIEW *v;

  zprintf("Zoom(%d): zoomwin=%lx, myevent.xmotion.win=%lx\n",
	  enable, zoomwin, nextwin /*myevent.xmotion.window*/);
  if (enable==ZIGNORE) return;		/* ignore LeaveEvent */

  /*nextwin = myevent.xkey.window;*/
  i = getview(nextwin);
  if (i == nwindow) return;
  v = view + i;
  zoomi = i;

  XQueryPointer(mydisplay, myevent.xmotion.window,
		&rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);

  /*printf("Zoom %d\n", enable);*/

  if (enable==ZDISABLE || enable==ZORIGINAL || enable==ZABORT)
    {
      if (zoom_on && zoomwin != 0)
	{
	  XUndefineCursor(mydisplay, zoomwin);
	  XSelectInput(mydisplay, zoomwin, EVENTMASK3);
	}
      else if (zoom_on && zoomwin == 0)
	zprintf("Wierd! zoom_on TRUE, zoomwin==0!\n");
      if (enable == ZORIGINAL)
      {
        v->clipped = 0;
	v->zoom_limit =0;
        XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);
      }
      else if (enable==ZABORT)
      {
	/*newclip(-1, xw, yw);*/
        if (zoom_count > 0)
	  XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);
	else crosshair(xcurs, ycurs);
      }
      zoom_on = zoom_count = 0;
    }

  else if (enable == ZGET)	/* Process button hit while enabled */
    {
      if (nextwin != zoomwin) return;
      ztest(&xw, &yw);	/* WARN! */
      xcurs = xw;
      ycurs = yw;
      if (special_zoom_flag) special_zoom(2);

      if (report_zoom) printf("zoom: xz[n] = %d; yz[n] = %d; n++;\n", 
			      xcurs, ycurs);	/* for values for special_zoom */

      if (zoom_count++ == 0)
	{
	  XDefineCursor(mydisplay, zoomwin, curs2);
	  if (xhair_type==XHAIR_CURSOR) xhair_type = XHAIR_1ST_LINE;
	  newclip(0, xcurs, ycurs);
	  crosshair(-1,0);
	  crosshair(xw, yw);	/* "erase" leaves crosshair on move */
	  if (xhair_type==XHAIR_1ST_LINE) xhair_type = XHAIR_CURSOR;
	}
      else
	{
	  zoom(zoomwin,ZDISABLE);	/* (no need to erase crosshair) */
	  newclip(1, xcurs, ycurs);
	  XClearArea(mydisplay, zoomwin, 0, 0, 0, 0, True);
	}
    }

  else				/* Enable crosshair */
    {
      //if( v->zoom_limit) { printf("BEEP%c",'\007'); return; }
      if (!cursor_defined) define_cursor(nextwin);
      zoom(zoomwin,ZDISABLE);

      zoomwin = nextwin;
      zoomi = getview(zoomwin);
      zoomgc = view[zoomi].gc;

      XDefineCursor(mydisplay, zoomwin, curs1);
      XSelectInput(mydisplay, zoomwin, EVENTMASK4);
      crosshair(-1, 0);	/* Initialize */
      crosshair(xcurs = xw, ycurs = yw);
      zoom_on = 1;
    }
}

/*-----------------------------------------------------------------------------
|	special_zoom
|	* find problem zooming into crit2
-----------------------------------------------------------------------------*/
void special_zoom(int option)
{
  static int xz[20], yz[20];
  static int i, n;
  if (option == 1) {
    n = i = 0;
    special_zoom_flag = 1;
    /*------ Begin zoom for 'slice' */
    xz[n] = 127; yz[n] = 128; n++;	// after zoom, '-' 3 times
    xz[n] = 139; yz[n] = 146; n++;	// 'crash in ticks'
    xz[n] = 108; yz[n] = 103; n++;
    xz[n] = 125; yz[n] = 114; n++;
    xz[n] = 113; yz[n] = 95; n++;
    xz[n] = 136; yz[n] = 112; n++;
    xz[n] = 131; yz[n] = 174; n++;
    xz[n] = 171; yz[n] = 201; n++;
    xz[n] = 143; yz[n] = 68; n++;
    xz[n] = 171; yz[n] = 95; n++;
    xz[n] = 169; yz[n] = 122; n++;
    xz[n] = 198; yz[n] = 152; n++;
    /*------ End 'crit' zoom */
  }

  else {
    xcurs = xz[i]; ycurs = yz[i]; i++;
    if (i == n) {
      special_zoom_flag = 0;
      printf("\nDONE SPECIAL ZOOM%c\n", '\007');
    }
  }
}

/*-----------------------------------------------------------------------------
|	define_cursor
-----------------------------------------------------------------------------*/
void define_cursor(Window nextwin)
{
   static Pixmap crosspix=0;
#define XX 0x80
   static unsigned char cross[] = {
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0xff,0xff,0xff,0xff,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0, 0,XX,0,0,
      0,XX,0,0, 0,XX,0,0, 0,XX,0,0
    };

   XColor fore, bkgd;
   unsigned bestdx, bestdy;

   if (xhair_type==XHAIR_CURSOR)
   {
      XQueryBestCursor(mydisplay, nextwin, xscreen, yscreen,
		       &bestdx, &bestdy);
      /*printf("Best cursor size %d %d\n", bestdx, bestdy);*/
      if (bestdx > 32) bestdx = 32;
      crosspix = XCreatePixmapFromBitmapData(mydisplay, nextwin, (char *)cross,
					     31, 31, 1, 0, 1);
      bkgd.pixel = mybackground;
      fore.pixel = myforeground;
      XQueryColor(mydisplay, cmap, &bkgd);
      XQueryColor(mydisplay, cmap, &fore);
      curs1 = curs2 = XCreatePixmapCursor(mydisplay, crosspix, crosspix,
					  &fore, &bkgd, 15, 15);
    }
   else
   {
      curs1 = XCreateFontCursor(mydisplay, CURS1);
      curs2 = XCreateFontCursor(mydisplay, CURS2);
   }
   cursor_defined = 1;
   /*   nextwin;*/
}

/*-----------------------------------------------------------------------------
|	tell_zoom -- prepare to call functions from motif
-----------------------------------------------------------------------------*/
void tell_zoom(Window zwin, int which, int flag)
{
   zoomi = which;
   zoomwin = view[which].window;
   zoomgc  = view[which].gc;
   coord_on = flag;
   /*myevent.xany.window = zwin;*/	/* should be same as zoomwin */
   /* zwin;*/
   myevent.xany.window = zoomwin;	/* WARN! zwin=0 from xmotif */
}

void tell_key(XEvent *p)
{
  myevent.xany.window = p->xany.window;
  myevent.xkey = p->xkey;
}

/*-----------------------------------------------------------------------------
|   newclip - called by zoom() to get current clip  & toggle_aspect 
|   * input xcurs,ycurs is relative to window, i.e. 0..width
|   * iflag: 0=1st point, 1=2nd point, -1=original
|   * gx, gy = fraction (0.0 ... 1.0) within inner box
|   * see DEAD_CODE in bak/src0704 for original method	
|   * v->f[i] is fractional position in box of cursor *from previous zoom*
-----------------------------------------------------------------------------*/
int newclip(int iflag, int xcurs, int ycurs)
{
  Window root;
  int x0, y0;
  VIEW *v;
  unsigned int bwidth, depth, zdx, zdy;
  double x2,y2,x1,y1,dx,dy, dx21, dy21;
  static double gx0, gy0;
  double gx, gy;
  float xmin,xmax,ymin,ymax, t;
  byte reportMe;
  int iResultx, iResulty;
  char text[80];
  float fx1, fy1, fx2, fy2, fx1_new, fy1_new, fx2_new, fy2_new;
  long lx1, ly1, lx2, ly2;
  double avgx, avgy, diffx, diffy;
  int samex, samey;

  v = view + zoomi;
  if (iflag == -1) { v->clipped = 0; return 0; }

  XGetGeometry(mydisplay, zoomwin, &root,	/* Geometry can change! */
	       &x0, &y0, &zdx, &zdy, &bwidth, &depth);
  zoomdx = (float)zdx;
  zoomdy = (float)zdy;
  if (v->zoom_limit) return 0;

  /*----- get fractional position of cursor relative to unzoomed -----*/

  gx = (double) (xcurs - fxbox1 * zoomdx) / (fdxbox * zoomdx);
  gy = (double) (ycurs - fybox1 * zoomdy) / (fdybox * zoomdy);

  if (debug_scale) {
    if (iflag==0) printf("\n");
    if (v->clipped && iflag==0) 
      printf("DS previous fraction %g,%g and %g,%g\n",
	     v->f[0],v->f[3],v->f[2],v->f[1]);
    printf("DS Cursor fractionally at %g, %g in box", gx,gy);
  }

  if (v->clipped)	/* gx, gy is fraction relative to unzoomed window */
    {			/* so using f[] is like one-step zoom from original */
      gx = v->f[0] + gx * (v->f[2] - v->f[0]);
      gy = v->f[3] + gy * (v->f[1] - v->f[3]);
      if (debug_scale) printf(" --> %g, %g", gx, gy);
    }

  if (debug_scale) printf("\n");

  /*------ First point: get gx0, gy0, data_limits_x, data_limits_y */
  /*       the data_limits are an appropriate 'epsilon' for dx, dy */
  /*  	   NOTE.  see also 'nan detected' in redraw0 */

  if (iflag == 0) {
    gx0 = gx;
    gy0 = gy;
    return 1;
  }

  /*------ Second point.  Check new zoom boundaries */

  samex = samey = reportMe = 0;
  fx1 = v->f[0];  fy1 = v->f[1];
  fx2 = v->f[2];  fy2 = v->f[3];

  /*------ Check x limits */

  if (v->zoom_limit & 1) samex++;

  else {
    if (gx0 < gx) { fx1_new = gx0; fx2_new = gx; }
    else          { fx1_new = gx;  fx2_new = gx0; }
    avgx  = fabs(fx1_new + fx2_new) / 2.0;
    diffx = fabs(fx2_new - fx1_new) / 2.0;

    if (diffx < (1.0e-7 * avgx)) {
      tiny_zoom(1, (float)avgx, &fx1_new, &fx2_new, 'X');
      samex++;
      strcat(text, "X");
      reportMe++;
    }
  }

  /*------ Check y limits */

  if (v->zoom_limit & 2) samey++;

  else {
    if (gy0 > gy) { fy1_new = gy0; fy2_new = gy; }
    else          { fy1_new = gy;  fy2_new = gy0; }
    avgy  = fabs(fy1_new + fy2_new) / 2.0;
    diffy = fabs(fy2_new - fy1_new) / 2.0;

    if (diffy < (1.0e-7 * avgy)) {
      tiny_zoom(1, (float)avgy, &fy2_new, &fy1_new, 'Y');
      samey++;
      if (v->zoom_limit & 1 || samex) strcat(text, " and ");
      strcat(text, "Y");
      reportMe++;
    }
  }

  /*------ Report zoom limits */

  strcpy(text, "Zoom limit reached in ");
  if (samex || (v->zoom_limit & 1)) { strcat(text, "X"); samex = 1; }
  if (samey || (v->zoom_limit & 2)) {
    if (samex) strcat(text, " and ");
    strcat(text, "Y");
  }

  if (reportMe) printf("%s%c\n", text, '\007');

  if (!(v->zoom_limit & 1) && !samex) {
    v->f1[0] = v->f[0]; v->f[0] = fx1_new;
    v->f1[2] = v->f[2]; v->f[2] = fx2_new;
    v->clipped = 1;
    redrawflag=1;
  }

  if (!(v->zoom_limit & 2) && !samey) {
    v->f1[1] = v->f[1]; v->f[1] = fy1_new;
    v->f1[3] = v->f[3]; v->f[3] = fy2_new;
    v->clipped = 1;
    redrawflag=1;
  }

  if (debug_scale) {
    lx1 = *(long *)&v->f[0];
    lx2 = *(long *)&v->f[2];
    ly1 = *(long *)&v->f[3];
    ly2 = *(long *)&v->f[1];
    printf("\nDS Zoom, fnew  = %.10g, %.10g to %.10g, %.10g\n",
	   v->f[0], v->f[3], v->f[2], v->f[1]);
    printf  ("DS Zoom, lfnew = %lx, %lx to %lx, %lx\n",
	   lx1, ly2, lx2, ly1);
  }

#ifdef DEAD_CODE

  dx = v->xmax - v->xmin;		/* current zoom (float to double, World) */
  dy = v->ymax - v->ymin;

  x2 = v->xmin + gx * dx;		/* coord of 2nd point (as double) */
  y2 = v->ymin + gy * dy;

  x1 = v->xmin + gx0 * dx;		/* coord of 1st point (as double) */
  y1 = v->ymin + gy0 * dy;

  if (debug_scale)
    printf("DS Zoom is %.10lg, %.10lg to %.10lg %.10lg\n", x1, y1, x2, y2);

  long lx1, ly1, lx2, ly2;
  float fx, fy;
  fx = x2; lx2 = *(long *)&x2;
  fy = y2; lx1 = *(long *)&y2;

  dx21 = fabs(x2 - x1);
  dy21 = fabs(y2 - y1);

  if (fabs(dx21) > 1.0e-11 && fabs(dy21) > 1.0e-11)
    {
      /* printf (" Newclip: DATA_LIMITS \n");*/
      v->f1[0] = v->f[0];      v->f[0] = gx0;
      v->f1[1] = v->f[1];      v->f[1] = gy0;
      v->f1[2] = v->f[2];      v->f[2] = gx;
      v->f1[3] = v->f[3];      v->f[3] = gy;

      if (v->f[0] > v->f[2]) fswap(0, 2);
      if (v->f[1] < v->f[3]) fswap(1, 3);
      printf("f[1]=%g, f[3]=%g\n", v->f[1], v->f[3]);
      v->clipped = 1;
    }

  else
    {
      printf("Zoomed limits are reached %c\n", '\007');
      printf("%.10lg: %.10lg, %.10lg\n", dx21, x2, x1);
      printf("%.10lg: %.10lg, %.10lg\n", dy21, y2, y1);
    }

  /*printf ("f[0],f[1],f[2],f[3]=\n %g %g %g %g \n",
    v->f[0],v->f[1],v->f[2],v->f[3]);*/
#endif
  return 1;
}

/*-----------------------------------------------------------------------------
|   givebox -- called by redraw() in preparation for possible zoom
|   * zoomdx etc are static, depend on window, not on current zoom
-----------------------------------------------------------------------------*/
void givebox(int x0, int y0, int w, int h, int ww, int wh)
{
  zoomdx = (float)ww;
  zoomdy = (float)wh;
  fxbox1 = (float) x0 / zoomdx;
  fybox1 = (float) y0 / zoomdy;
  fdxbox = (float) w / zoomdx;
  fdybox = (float) h / zoomdy;
}

/*-----------------------------------------------------------------------------
|	set_aspected_limits
|	* IN:  xmin1, xmax1 = autoscale data limits, xmin, xmax = aspected
|	* which: 0=x, 1=y
-----------------------------------------------------------------------------*/
void set_aspected_limits(Window win, int which, float xmin1, float xmax1,
			 float xmin, float xmax)
{
  VIEW *v;
  float min, max, temp;
  v = view + getview(win);
  if (! (v->clipped & 2))
    {
      min = (xmin - xmin1) / (xmax1 - xmin1);
      max = (xmax - xmin1) / (xmax1 - xmin1);
      if (which==1) { temp = min; min = max; max = temp; }
      v->f[0+which] = min;
      v->f[2+which] = max;
      if (which==1) v->clipped |= 2;
#ifdef DEAD_CODE_ASP
      printf("Aspect%s:\
 %f ... %f --> %f ... %f ==> Fraction = %f ... %f\n",
	     (which==0)?"x":"y", xmin1,xmax1,xmin,xmax,
	     v->f[0+which], v->f[2+which]);
      if (which==0) printf("f[] before: ");
      else printf("f[] after: ");
      printf("%f ... %f\n", v->f[1], v->f[3]);
#endif
    }
}

/*-----------------------------------------------------------------------------
|   get_zoomedclip - called by redraw() to get current clip
|   e.g. xmin --> xmin + fx1 * (xmax-xmin)
-----------------------------------------------------------------------------*/
int get_zoomedclip(Window w, float *fx1, float *fy1, float *fx2, float *fy2)
{
  VIEW *v;
  int i;

  i = getview(w);
  v = view + i;
  if (i == nwindow || !v->clipped)
    return (0);

  if (fx1 != NULL)
    {
      *fx1 = v->f[0];
      *fx2 = v->f[2];
      *fy1 = v->f[1];
      *fy2 = v->f[3];
    }
  return (v->clipped);
}

/*-----------------------------------------------------------------------------
|	tiny_zoom
|	* use IEEE format of float: 23 bits value, 8 bits exponent, 1 bit sign
-----------------------------------------------------------------------------*/
int tiny_zoom(int inc, float fxavg, float *pfx1, float *pfx2, char which)
{
  float fx1, fx2, f_one, fx1_new, fx2_new;
  double ffx1, ffx2, d_one;
  long lx1, lx2, lx1_new, lx2_new, himask, lomask, simask, l_one;

  simask = 0x80000000;
  himask = 0x7ff80000;
  lomask = 0x0007ffff;		/* 23 bits */
  fx1 = *pfx1; lx1 = *(long *)pfx1; ffx1 = (double)fx1;
  fx2 = *pfx2; lx2 = *(long *)pfx2; ffx2 = (double)fx2;

  /*------ Create a "one" which will add 1 to the 'long' version of x */

  l_one = lx1 & himask;
  l_one = l_one >> 23;
  l_one = l_one - 23;
  l_one = l_one << 23;
  f_one = *(float *)&l_one;
  d_one = (double)f_one;

  //if (inc) { ffx1 += d_one;  ffx2 -= d_one; }
  //else     { ffx1 -= d_one;  ffx2 += d_one; }

  if (inc) { ffx1 = fxavg - d_one;  ffx2 = fxavg + d_one; }
  else     { ffx1 -= d_one;  ffx2 += d_one; }
  //else     { ffx1 = fxavg - d_one;  ffx2 = fxavg + d_one; }

  fx1_new = (float)ffx1; lx1_new = *(long *)&fx1_new;
  fx2_new = (float)ffx2; lx2_new = *(long *)&fx2_new;

  if (debug_scale || debug_scale2) {
    printf("\n'One' = %lx = %g, inc=%d\n", l_one, f_one, inc);
    printf("f%c1 = %.10g = %lx, new = %.10g = %lx\n", which, fx1, lx1, fx1_new, lx1_new);
    printf("f%c2 = %.10g = %lx, new = %.10g = %lx\n", which, fx2, lx2, fx2_new, lx2_new);
  }

  *pfx1 = fx1_new;
  *pfx2 = fx2_new;
  return 0;
}

/*-----------------------------------------------------------------------------
|	addzoom
|	* inc = 1 for zoom in - bigger image, smaller surrounding box,
|	        0 for zoom out, 2 for previous
|	* f[]'s are the fractional positions of the previous zoom
|	  coordinates relative to unzoomed data
-----------------------------------------------------------------------------*/
void addzoom(Window w, int inc)
{
  double df, avgx, avgy, diffx, diffy;
  float fx1, fy1, fx2, fy2;
  float fx1_new, fy1_new, fx2_new, fy2_new;
  long lx1, lx2, ly1, ly2;
  long lx1_new, lx2_new, ly1_new, ly2_new;
  char text[80];
  int iResultx, iResulty;
  byte reportMe;
  VIEW *v;

  df = inc ? -zoomfac : zoomfac;	/* e.g. plus-or-minus 0.05 */
  df = (float) 1 + df / (float)2;	/* e.g. 1.025 or 0.975 */
  v = view + getview(w);

  if (!v->clipped)
    {
       v->f[0] = v->f[3] = (float)0;
       v->f[2] = v->f[1] = (float)1;
       v->clipped = 1;
    }

  if (v->zoom_limit == 3 && inc == 1 )
    return;

  if (v->clipped && inc == 2)		/* from key '\010' */
    {
      v->f[0] = v->f1[0];
      v->f[1] = v->f1[1];
      v->f[2] = v->f1[2];
      v->f[3] = v->f1[3];

      v->zoom_limit = 0;
      redrawflag=1;
      return;
    }

  fx1 = v->f[0];  fy1 = v->f[1];
  fx2 = v->f[2];  fy2 = v->f[3];

  avgx  = (fx1 + fx2) / 2.0;
  avgy  = (fy1 + fy2) / 2.0;

  diffx = (fx2 - fx1) / 2.0;		/* both positive numbers */
  diffy = (fy1 - fy2) / 2.0;

  /*------ Calculate new fractional distances */

  reportMe = 0;
  iResultx = iResulty = 0;
  strcpy(text, "Zoom limit reached in ");

  if (inc == 1 && v->zoom_limit & 1) {
    fx1_new = fx1; fx2_new = fx2;
  }

  else {
    if (debug_scale) printf("\nDS Addzoom, X avg %.10g, diff %.10g, delta %.10g\n", 
			    (float)avgx, (float)diffx, (float)(df * diffx));
    fx1_new = avgx - df * diffx;
    fx2_new = avgx + df * diffx;
    if (fx1_new == fx1 || fx2_new == fx2 || fx1_new == fx2_new) {
      iResultx = tiny_zoom(inc,(float)avgx, &fx1_new, &fx2_new, 'X');
      if (inc) {
	v->zoom_limit |= 1;
	strcat(text, "X");
	reportMe++;
      }
    }
    else v->zoom_limit &= 2;
  }

  if (inc == 1 && v->zoom_limit & 2) {
    fy1_new = fy1; fy2_new = fy2;
  }
  else {
    if (debug_scale) printf("\nDS Addzoom, Y avg %.10g, diff %.10g, delta %.10g\n", 
			    (float)avgy, (float)diffy, (float)(df * diffy));
    fy1_new = avgy + df * diffy;
    fy2_new = avgy - df * diffy;
    if (fy1_new == fy1 || fy2_new == fy2 || fy1_new == fy2_new) {
      iResulty = tiny_zoom(inc, (float)avgy, &fy2_new, &fy1_new, 'Y');
      if (inc) {
	v->zoom_limit |= 2;
	if (v->zoom_limit & 1) strcat(text, " and ");
	strcat(text, "Y");
	reportMe++;
      }
    }
    else v->zoom_limit &= 1;
  }

  /*------ Report zoom limits */

  if (reportMe) printf("%s%c\n", text, '\007');

  if (fx1_new!=fx1 || fx2_new!=fx2 || fy1_new!=fy1 || fy2_new!=fy2) {
    v->f[0] = fx1_new;    v->f[2] = fx2_new;
    v->f[1] = fy1_new;    v->f[3] = fy2_new;
    redrawflag=1;
  }

  if (debug_scale) {
    lx1 = *(long *)&fx1_new;
    lx2 = *(long *)&fx2_new;
    ly1 = *(long *)&fy1_new;
    ly2 = *(long *)&fy2_new;
    printf("\nDS Addzoom, new  = %.10g, %.10g to %.10g, %.10g\n",
	   fx1_new, fy2_new, fx2_new, fy1_new);
    printf("DS Addzoom, %lx, %lx to %lx, %lx\n",
	   lx1, ly2, lx2, ly1);

    fx1 = (float)(df*diffx); lx1 = *(long *)&fx1;
    fy1 = (float)(df*diffy); ly1 = *(long *)&fy1;
    printf("DS Addzoom, diff = %.10g, %.10g = %lx, %lx\n",
	   fx1, fy1, lx1, ly1);
  }
}

/*-----------------------------------------------------------------------------
|	show_coord
|	* coord_on: 1=coordinate 2=slope, 3=ratio,
|	  4=gradient (generalized contour plots)
-----------------------------------------------------------------------------*/
void show_coord(int xcurs, int ycurs, int count, CURVE_SET *cp)
{
  float xmin, ymin, xmax, ymax, dx,dy, x,y, delta;
  float psixy, gradx, grady, min_dist, psi_seg;
  double psi0, psix, psiy, xtell, ytell;
  static float x1,y1;
  XRectangle bigbox, clipbox;
  static char *point[] = { "point on ", "y-value for " };
  static char *slope[] = { "slope", "ratio" };
  static char *which[] = { "1st ", "2nd "};
  char text[500], vtext[200], gtext[200], *pt;
  int i, j, n, iError, iview;
  XColor xcolor;
  
  /*------ Calculate actual coordinate x,y from xcurs,ycurs */
  
  iview = getview(zoomwin);
  get_box_info(zoomwin, iview, &bigbox, &clipbox, &xmin, &ymin, &xmax, &ymax);
  get_world(xcurs, ycurs, clipbox, xmin, ymin, xmax, ymax, &x, &y);

  /*------ Create 1st part of text */

  i = coord_on - 2;		/* 0: slope, 1: ratio */
  j = (count==1) ? 0 : 1;
  *text = 0;
  if (coord_on==1 || coord_on == 4)
    sprintf(text, "Coordinate = ");
  else if (coord_on==2 || coord_on==3)
    sprintf(text, "%s%s%s = ", which[j], point[i], slope[i]);

  pt = text + strlen(text);
  xtell = x; ytell = y;

  if (coord_on == 1) {
    if (cp->use_log & 1) xtell = (float)pow(10.0, (double)x);
    if (cp->use_log & 2) ytell = (float)pow(10.0, (double)y);
  }

  if (coord_on == 1 || coord_on==2 || coord_on==4) {
    //printf("x=%g --> %g, y=%g --> %g\n", x, xtell, y, ytell);
    //sprintf(pt, "%.10g, %.10g", xtell, ytell);
    sprintf(pt, "%g, %g", xtell, ytell);
  }
  else if (coord_on==3) sprintf(pt, "%.4g", y);
  else sprintf(pt, "Coord type %d", coord_on);

  /*------ If contour, find value at point */

  *vtext = '\0';
  if ((coord_on == 1 || coord_on == 4) && cp->gtype=='M') {
    /*printf("\ncalling nearest_m_value for coord at %s\n", pt);*/
    crosshair(xcurs, ycurs);			/* turn crosshair off */
    pt = text + strlen(text);			/* append value to text */

    iError = 
      nearest_M_value(x, y, cp, vtext, &psi0);		/* find nearest contour value */

    if (coord_on == 4) {
      delta = 0.00001;
      x1 = (1.0 + delta) * x;
      y1 = (1.0 + delta) * y;

      iError = 
	nearest_M_value(x1, y, cp, gtext, &psix);		/* find nearest contour value */
      iError =
	nearest_M_value(x, y1, cp, gtext, &psiy);		/* find nearest contour value */

      psix -= psi0; dx = (x1 - x);
      psiy -= psi0; dy = (y1 - y);

      gradx = psix / (double)dx;
      grady = psiy / (double)dy;
      sprintf(gtext, "     Grad  = %g, %g", gradx, grady);
      //sprintf(gtext + strlen(gtext), " (from %lg/%g, %lg/%g)", psix, dx, psiy, dy);
    }

    xcolor.pixel = getcolor(65535, 65535, 65535);
    XSetForeground(mydisplay, zoomgc, xcolor.pixel);	/* got changed during redraw1 */
    crosshair(xcurs, ycurs);				/* turn crosshair back on */
  }

  /*------ Print message */

  xfromcoord = x;
  yfromcoord = y;
  strcat(text, "\n");
  xprintf(text);

  /*------ Print 2nd part of message */

  if (coord_on == 1 || coord_on==4) {	/* 'c' or 'g' mode */
    if (cp->gtype=='M') {
      xprintf("%s\n", vtext);
      if (coord_on == 4) xprintf("%s\n", gtext);
    }
  }

  else if (count==2)			/* 's' or 'r' mode */
  {
    n = strlen(which[j]) + strlen(point[i]);
    pt = text;
    //for (j=0,*pt++='\n'; j<n; j++,*pt++=' ');
    for (j=0; j<n; j++,*pt++=' ');
    sprintf(pt, "%s = ", slope[i]);
    pt += strlen(pt);
    if (coord_on==2)
    {
      dx = x - x1;
      dy = y - y1;
      if (dx == FZ) sprintf(pt, "infinity");
      else sprintf(pt, "%g", dy / dx);
    }

    else if (coord_on==3)
    {
      if (y==FZ) sprintf(pt, "infinity");
      else sprintf(pt, "%g", (float)y1 / (float)y);
    }

    strcat(text, "\n");
    printf(text);

#ifndef MOTIF
    if (coord_on==2 || coord_on==3)	/* disable 's' or 'r', redraw win */
      coord(coord_on);
#endif
  }
  x1 = x;
  y1 = y;
}

/*-----------------------------------------------------------------------------
|	coord
|	* 1=enable coord, 2=enable slope, 3=enable ratio,
|	  4=enable gradient, 8=process hit
-----------------------------------------------------------------------------*/
void coord(int enable)
{
  Window nextwin;
  static int cursor_defined = 0;
  int i, old_coord_on;
  int xw, yw, xr, yr;
  unsigned keys_buttons;
  Window rw, cw;
  VIEW *v;
  CURVE_SET *cp;

  nextwin = myevent.xkey.window;
  i = getview(nextwin);
  if (i == nwindow) return;
  v = view + i;
  cp = curveset+i;

  if (enable != 8 && coord_on)		/* it's a toggle, on-->off */
    {
      XUndefineCursor(mydisplay, zoomwin);
      XSelectInput(mydisplay, zoomwin, EVENTMASK3);
      old_coord_on = coord_on;
      coord_on = coord_count = 0;
      /*XClearArea(mydisplay, nextwin, 0, 0, 0, 0, True);*/
      crosshair(xcurs, ycurs);
      if (old_coord_on == enable) return;
    }

  XQueryPointer(mydisplay, myevent.xmotion.window,
		&rw, &cw, &xr, &yr, &xw, &yw, &keys_buttons);
  xcurs = xw;
  ycurs = yw;

  if (enable == 8)		/* Process button hit while enabled */
    {
      if (nextwin != zoomwin) return;
      if (coord_count >= 2) coord_count = 0;
      coord_count++;
      show_coord(xcurs,ycurs, coord_count, cp);
    }

  else				/* Enable! (off-->on) */
  {
    zoomwin = nextwin;
    coordi = getview(zoomwin);
    zoomgc = view[coordi].gc;
    coord_on = enable;
    if (!cursor_defined) define_cursor(nextwin);
    XDefineCursor(mydisplay, zoomwin, curs1);
    XSelectInput(mydisplay, zoomwin, EVENTMASK4);
    crosshair(-1,0);	/* Initialize */
    crosshair(xcurs, ycurs);
  }
}

/*-----------------------------------------------------------------------------
|	crosshair
|	* normally called in pairs: erase the old, draw the new
|	* crosshair types: 1=normal, all others for cases where no xor
|	*   0=save under, draw line; 3=user cursor
-----------------------------------------------------------------------------*/
void crosshair(int x, int y)
{
  int x1, y1, x2, y2;
  unsigned int dx, dy;
  static int state = -1;
  static XImage *xih=0, *xiv=0;
  static XImage *xiH=0, *xiV=0;

  if (!rubber) return;
  if (x==-1) { state = 0; return; }

  x1 = (int)(fxbox1 * zoomdx);
  y1 = (int)(fybox1 * zoomdy);
  x2 = (int)((fxbox1 + fdxbox) * zoomdx);
  y2 = (int)((fybox1 + fdybox) * zoomdy);
  dx = x2-x1+1;
  dy = y2-y1+1;
  if (x < x1 || x > x2 || y < y1 || y > y2) return;

  if (xhair_type==XHAIR_CURSOR) ;
  else if (xhair_type==XHAIR_XOR || xhair_type==XHAIR_1ST_LINE)
    {
      if (xhair_type==XHAIR_XOR) set_xor(zoomgc, 1);
      XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
      XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
    }

  else
    {
  /* This method involved saving H & V lines as images, then
   * using XPutImage.  Problem if resize window: didn't resize lines
   */
      if (xiH == 0)
	{
	  xih = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiv = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
	  XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
	  xiH = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiV = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  XPutImage(mydisplay, zoomwin, zoomgc, xih, 0,0,x1,y, dx,1);
	  XPutImage(mydisplay, zoomwin, zoomgc, xiv, 0,0,x,y1, 1,dy);
	}

#ifdef TRY_XORP
  /* This image was TRY_SETP with the hope that XPutImage() responds
   zoom* to set_xor().  It didn't!
   */
      set_xor(zoomgc, 1);
      XPutImage(mydisplay, zoomwin, zoomgc, xiH, 0,0,x1,y, dx,1);
      XPutImage(mydisplay, zoomwin, zoomgc, xiV, 0,0,x,y1, 1,dy);
#endif
/*
Note about bad response time on quetzal.  Can have "state=0" here and
nothing on screen for a while.  Then quick flash on an intermediate crosshair,
then the final one steady.  Means it takes a while for the X-function
events to be implemented.
*/
      if (state==1)
	{
	  XPutImage(mydisplay, zoomwin, zoomgc, xih, 0,0,x1,y, dx,1);
	  XPutImage(mydisplay, zoomwin, zoomgc, xiv, 0,0,x,y1, 1,dy);
	  state=0;
	}
      else
	{
	  xih = XGetImage(mydisplay, zoomwin, x1,y, dx,1, AllPlanes,XYPixmap);
	  xiv = XGetImage(mydisplay, zoomwin, x, y1,1, dy,AllPlanes,XYPixmap);
	  if (xhair_type==XHAIR_LINE)
	    {
	      XDrawLine(mydisplay, zoomwin, zoomgc, x1, y, x2, y);
	      XDrawLine(mydisplay, zoomwin, zoomgc, x, y1, x, y2);
	    }
	  else
	    {
	      XPutImage(mydisplay, zoomwin, zoomgc, xiH, 0,0,x1,y, dx,1);
	      XPutImage(mydisplay, zoomwin, zoomgc, xiV, 0,0,x,y1, 1,dy);
	    }
	  state=1;
	}
    }

  set_xor(zoomgc, 0);
}
