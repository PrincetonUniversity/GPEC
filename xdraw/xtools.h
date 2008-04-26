/*-----------------------------------------------------------------------------
|	xtools.h
-----------------------------------------------------------------------------*/
extern void give_command_args(int, char**);
extern int opendisplay(char *title);
extern void closedisplay(void);
extern void get_screensize(unsigned int *dx, unsigned int *dy);
extern void init_menuXlibrary(void);
extern int get_eventmask(void);
extern int getview(Window);
extern Window makewindow(char *title, int has_menu);
extern void nextwindow(int *xw, int *yw);
extern Window getwindow(int);
extern Window get_op_window(Window win);
extern void get_expose_info(Display **d,Window *w,GC *g, int *font_heightp);
extern void get_properties(Window w, unsigned int *wp, unsigned int *hp,
			   unsigned int *dp);
extern Window enable_window(Window win);

extern void clear_if_needed(Window, GC, unsigned int, unsigned int);
extern unsigned long white(void);
extern int textwidth(char *s, int n);
extern void event(void);
extern int readkey(void);
extern int parsekey(void);
extern void clear_redraw(int);
extern void xmessage(int x, int y, char *s);
extern int  get_zoomedclip(Window w, float *x1, float *y1, float *x2, float *y2);
extern void tell_zoom(Window, int, int);
extern void tell_key(XEvent *);
extern void addzoom(Window w, int inc);
#ifndef MOTIF
extern void zoom(Window w, int enable);
#endif
extern void coord(int enable);
extern void define_cursor(Window);
extern int   newclip(int i, int xcurs,int ycurs);
extern void set_aspected_limits(Window, int, float, float, float, float);
extern void givebox(int wx0, int wy0, int w, int h, int width, int height);
extern void crosshair(int x,int y);
extern Window set_expose(int i);
extern void set_winmask(Window);
extern void xrectangle(unsigned int dx, unsigned int dy, int bkgd, int fore);

#define ZGET	   2
#define ZENABLE    1
#define ZDISABLE   0
#define ZORIGINAL -1
#define ZIGNORE   -2
#define ZABORT	   3
#define ZPREVIOUS  4

#define XHAIR_XOR    0
#define XHAIR_LINE   1
#define XHAIR_CURSOR 2
#define XHAIR_1ST_LINE 8

typedef struct			/* In xtools.c and event.c */
  {
    Window window;
    Window work_area;
    GC gc;
    float f[4];			/* x1,y1,x2,y2, as read from cursor */
    float f1[4];
    float old_f[4], old_f1[4];
    int clipped;

    /* Window attributes */
    int window_flags;		/* Visibility,MENU,Status,SIZE,POZITION*/
    int x,y;			/* Position*/
    int width,height;		/* Geometry*/
    float xmin,xmax,ymin,ymax;  /* xmin,xmax,ymin,ymax currently used */
    float old_xmin, old_xmax;	/* save for restore from nan eg */
    float old_ymin, old_ymax;
    byte zoom_limit;
  }
VIEW;
#define MENU    0x1
#define STATUS  0x2
#define VISIBILITY 0x4

#define WIN_MNGR MAXWINDOW+1








