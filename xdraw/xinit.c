/******************************************************************************
**  NAME      XINIT.C
**  AUTHOR    Sheryl M. Glasser
**
**  DESCRIPTION
**      Read in data from .in and .bin, load structures
**
**  Copyright (c) GlassWare 1993.  All rights reserved.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <fcntl.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <X11/Xlib.h>
#define O_BINARY 0
#define PATHSEP '/'
#define STATENAME "state.xdraw"
#define LATIN 1
#define GREEK 2
#define NO_END_STR 1
#define L_TBL 20 /* number of symbols in table */

struct tlang 
{
  int length;
  char *token;
  int code;
};

typedef struct  tlang TLANG;

static TLANG table[L_TBL]={
{6,"\\alpha",'\141'},                    /* 1 */
{5,"\\beta",'\142'},
{6,"\\gamma",'\147'},
{3,"\\pi",'p'},
{4,"\\psi",'y'},
{4,"\\rho",'r'},
{6,"\\theta",'q'},
{4,"\\phi",'f'},
{5,"\\phi1",'j'},
{6,"\\sigma",'s'},
{3,"\\xi",'x'},
{4,"\\chi",'c'},
{6,"\\omega",'w'},
{6,"\\delta",'d'},
{3,"\\mu",'m'},
{6,"\\Delta",'D'},
{6,"\\Gamma",'G'},
{6,"\\Omega",'W'},
{6,"\\Theta",'Q'},
{7,"\\Lambda",'L'}         /*20 - */
};

#ifndef O_RDONLY
#define O_RDONLY 0			/* for penning */
#endif

#include "gendefs.h"
#include "curves.h"
#include "xinit.h"
#include "ps.h"
#include "xcontour.h"
#include "xtools.h"
#include "menuwin.h"
#include "xdraw.h"
#include "xedit.h"		/* for init_label */

static void freadf(FILE *, char *, char *);
static void readstate(void);
static FILE *get_inname(int argc, char *argv[]);
static int parsevar(char *);
static void testvar(CURVE_SET *);
static void read_asc_header(FILE *frd);
static void use_asc_header(void);
static int var_exists(int index);
static void set_indexed_value(char *p, char *fname);
static char *set_subtitle(char *, CURVE_SET *, char, char *);

#define QUOTE '\"'

extern float *buf;
static unsigned Long bufsize;
static char namebuf[NCHAR];
static char newtitle[100];
static char newlabel[100];
static char basepath[100];
static char xpolar[20]="X",ypolar[20]="Y";

extern int ftype;	/* 0=graph, 1=contour, 2=self-determined, 3,4 */
extern int multi_topology, n_m2, i_m2;
extern int background;  /* 0=background = black */
extern int gtype;	/* graph in window: gGCQLBD etc */
static int version;	/* 0=Alan's draw.in, 1=Sheryl's */
extern int default_ncol, default_ncurve;
extern int pscolor, pstitle, splinetype, showcomment, figures;
extern int labeldy, labelsep;
extern int ps_only, axis_width,curves_width;
extern float psfontfac;
extern int xhair_type;
extern int dialogwindow;
extern Font font_greek, font_latin;
extern int nwindow;
extern VIEW view[];
extern int tellMe;
extern int force_invert_bytes;

static char datapath[100];
char bin_file_name[LLNAME];
char prefix[20];
int zdebug = 0;			/* Only while zoom error - see xtools */
int ztestoption=0;
extern int debug_index, ldebug, verifying;

extern char title[];

extern LOOP loop[];		/* Loop structure of the data */
extern int nloop;
static int iloop, ih0=0;
static int counting_outer=0, i_outer=0;
int param[MAXLOOP];
int outerloop_added = 0;
int outerloop_count[] = { 0,0,0,0,0,0,0,0,0,0 };

static int use_hcount=0, stringbytes=0;
static char *dummytitle = " ";
extern NODE *nodelist;
extern int nnode, inode, ivar, ncount_equal;
extern int debug_mb;

#define NBLOCK 20

#define MAXVAR 128
VARIABLE_ID varid[MAXVAR];	/* Variables used for plotting */
int nvar=0;

extern CURVE_SET curveset[];	/* Describes the plots */
extern CURVE_SET *cp, *cp2;
extern int ncset;

/*DD*/
#define MAX_FILES 20
static int files_count=0;
static BINFILE binfile[MAX_FILES]; 
/*DD*/

#define is_white(c) (c==' ' || c=='\t')
#define is_eol(c)   (c=='\0' || c=='\n')
#define is_crlf(c)  (c=='\n' || c=='\r' || c==27)

#ifdef USEMAIN
/*-----------------------------------------------------------------------------
|	main (for debugging)
-----------------------------------------------------------------------------*/
void main(int argc, char *argv[])
{
  init(argc,argv);
}
#endif

/*-----------------------------------------------------------------------------
|	ztest
|	* testing = 0 via codeview --> use xtest[], ytest[] as zoom
|	* ztestoption: 1=print zoom coords, 0=no option,
|	  2 ==> array of xtest[], ytest[] (xw = char *)
-----------------------------------------------------------------------------*/
void ztest(int *xw, int *yw)
{
   static int testing=-1;	/* normal value -1 */
   static int nprint=0;
   char *p,*q;
   int i;

   static int xtest[12], ytest[12];
   if (ztestoption==1)
   {
      xprintf("%d %d   ", *xw, *yw);
      if (nprint == 1) xprintf("\n");
      nprint = nprint ? 0 : 1;
   }
   else if (ztestoption==2)
   {
      for(i=0, p=(char *)xw; q=skip_to_arg(p,2*i+2,0); i++)
        sscanf(q, "%d %d", xtest+i, ytest+i);
      *(xtest+i) = -1;
      testing = ztestoption = 0;
   }
   else if (testing>=0)
   {
      if (xtest[testing]==-1) testing=0;
      if (testing==0) xprintf("Warning! autozoom in effect\n");
      *xw = xtest[testing];
      *yw = ytest[testing];
      testing++;
   }
}

/*=============================================================================
**                  INITIALIZATION
**===========================================================================*/
/*-----------------------------------------------------------------------------
|	parsevar
-----------------------------------------------------------------------------*/
static int parsevar(char *p)
{
  int mask, i;
  mask = 0;
  if ( ftype != 6)
    {
      if (*p=='i' || *p=='I') { mask = LOOPBIT; p++; }
      sscanf(p,"%d", &i);
      i |= mask;
      return(i);
    }

  mask = LOOPBIT;
  switch (*p)
    {
    case 't':
    case 'T':
      return(0);
    case 'x':
    case 'X':
      return(1|mask);
    case 'y':
    case 'Y':
      return(2|mask);
    }
  sscanf(p,"%d", &i);
  if( i>=0 && i< 100)
    return(i);

  return(-1);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
static int variable_ok(int index)
{
  int i, result;
  result = 1;
  i = index & ~LOOPBIT;
  if (index & LOOPBIT)
    {
      if (i < 0 || i >= nloop) result = 0;
    }
  else if (i < 0 || i > loop[ivar].count + loop[ivar].ih0) result = 0;
  return result;
}

static void testvar(CURVE_SET *cp)
{
  int ok1, ok2;
  ok1 = variable_ok(cp->ix.index);
  ok2 = variable_ok(cp->iy.index);
  if (ok1 && ok2) return;
  i_abort("Inappropriate variable found",
	  ok2 ? cp->ix.index : cp->iy.index, cp - curveset);
}

/*-----------------------------------------------------------------------------
|	parse_title, parse_subtitle
|	* returns parsed version of title or subtitle in text
|	* <i0> variables expanded, SINGLE expanded
-----------------------------------------------------------------------------*/
int parse_title(CURVE_SET *cp, char *title, char *text)
{
  char *p, *p2;
  int var, i;
  var = -1;
  p = strchr(title,'<');
  if (p) p2 = strchr(p,'>');
  if (p && p2) var = 0;
  if (var != -1)
  {
    *p2 = '\0';
    i = parsevar(p+1);
    *p2 = '>';
    if (!variable_ok(i) || ! (i & LOOPBIT)) var = -1;
    else var = i & ~LOOPBIT;
  }
  if (var == -1) {
    strcpy(text, title);
    if (ftype != 0) sprintf(text+strlen(text),",  extrema=(%10.3e,%10.3e)", cp->iz.min, cp->iz.max);
  }

  else
  {
    strncpy(text, title, p-title);
    sprintf(text+(p-title), "%d", cp->param[var]);
    strcat(text, p2+1);
  }
  return( (var==-1) ? 0 : 1);
}

int parse_subtitle(CURVE_SET *cp, char *string)
{
  if (!cp->subtitle) return 0;
  if ((cp->flags & SINGLE) && load_single_label(cp, string)) ;
  else parse_title(cp, cp->subtitle, string);
  return 1;
}

/*-----------------------------------------------------------------------------
|	axis_label
-----------------------------------------------------------------------------*/
void axis_label(char *p, CVAR *xp)
{
  char text[20], *q;
  int n;
  static int inited=0;
  if (!inited) { *newlabel='\0'; inited++; }
  for(q=text; !is_white(*p) && !is_eol(*p); ) *q++ = *p++;
  *q = '\0';
  for(q=newlabel; *q && strcmp(q, text); q += strlen(q) +1) ;
  if (!*q) { strcpy(q, text); *(q+strlen(q)+1) = '\0'; }
  xp->label = q;

  m_lang_str(xp->label,xp->ml_text,&n);
  xp->nitems =n;
}

/*-----------------------------------------------------------------------------
|    freadf
|	version=0: read a header line, then a value
|	version=1: read a label on current line to (:), then a value
|	both cases skip over blank lines
-----------------------------------------------------------------------------*/
static void freadf(FILE * frd, char *format, char *s)
{
  char text[80], *p;
  int i, nline;
  nline = version ? 1 : 2;	/* version 0 has to read a header */
  for (i = 0; i < nline;)	/* Read 1 or 2 lines, plus blank lines */
    {
      fgets(text, 80, frd);
      *strchr(text, '\n') = 0;
      if (*text != 0)
	i++;
    }
  if (version == 1)
    {
      if (!strcmp(format,"%l")) { strcpy(s,text); return; }
      p = strchr(text, ':');
      for (p++; isspace(*p); p++);
      if (!strcmp(format, "%s"))
	strcpy(s, p);
      else if (!strcmp(format, "%c"))
	*s = *p;
      else if (!strcmp(format, "%d"))
	sscanf(p, format, (int *) s);
      else if (!strcmp(format, "%f"))
	sscanf(p, format, (float *) s);
    }
  else
    {
      if (!strcmp(format, "%s"))
	strcpy(s, text);
      else
	sscanf(text, format, s);
    }
}

static int blankline(char *s)
{
  char *p;
  for(p=s; *p==' ' || *p=='\t'; p++);
  return (*p < ' ');
}

/*-----------------------------------------------------------------------------
|	readstate
-----------------------------------------------------------------------------*/
static void readstate()
{
  FILE *frd;
  char ddelim,*p,fname[LNAME];
  char *inipath, ininame[LNAME+100];
  int i, verbose;
  //#define code(s) !strcmp(fname, s)

  verbose = 0;
  ddelim = PATHSEP;
  *datapath=0;				/* init variables */
  strcpy(prefix, "draw");

  p = ininame; *p = 0;
  inipath = getenv("XDRAWPATH");
  if (verbose) printf("inipath XDRAWPATH=%s\n", inipath);     // (null)

  if(!inipath) inipath = (char *) getenv("HOME");  // /home/sheryl or /root
  if (verbose) printf("inipath HOME=%s\n", inipath);

  if (inipath)
    {
      strcat(ininame, inipath);
      p = ininame + strlen(ininame);
      if (*(p-1) != PATHSEP) { *p++ = PATHSEP; *p = 0; }
    }

  strcat(ininame, "xdraw.ini");
  if (verbose) printf("INI file is %s\n", ininame);

  frd = fopen(ininame, "rt");

  if (frd==NULL)
    {
      if (verbose) printf("Can't find inifile %s, try %s\n", ininame, STATENAME);
      strcpy(p, STATENAME);
      if (verbose) printf("ini file %s\n", STATENAME);
      if ((frd = fopen(STATENAME, "rt"))==NULL) return;
    }
  //xprintf("**readstate %s\n", ininame);

  /*-----xdraw.ini or state.xdraw exists */

  for (;;)
    {
      p = fgets(fname, LNAME, frd);		/* read a line */
      if (!p || !strncmp(fname,"END",3))
	break;
      if (*fname==';') continue;		/* comment line */
      if (!(p=strchr(fname,':'))) continue;	/* look for (:) */

      for (*p++=0; is_white(*p); p++);		/* skip white */
      if (is_eol(*p)) continue;
      //printf("%s\n", fname);

      if (!strcmp(fname,"path"))		   	/* "path" */
	{
	  strcpy(datapath, p);
	  if (verbose) printf("datapath=%s\n",datapath);
	  p = datapath + strlen(datapath) - 1;
	  if (*p != ddelim)
	    {
	      *p++ = ddelim; *p = 0;
	    }
	  if (verbose) printf("datapath=%s\n",datapath);
        }
      else if (!strcmp(fname,"where_am_i")) 
	{ if (verbose) printf("INI file used ('where_am_i') is in %s\n",p); }
      else if (!strcmp(fname,"prefix"))			/* "prefix" */
	sscanf(p, "%s", prefix);
      else if (!strcmp(fname,"crosshair"))		/* "crosshair" */
	sscanf(p, "%d", &xhair_type);
      else if (!strcmp(fname,"splinetype"))		/* "splinetype" */
	sscanf(p, "%d", &splinetype);
      else if (!strcmp(fname,"figures"))		/* # "figures" for labels */
	sscanf(p, "%d", &figures);
      else if (!strcmp(fname,"markerstyle"))		/* "markerstyle" */
	setmarkerstyle(p);
      else if (!strcmp(fname,"showcomment"))		/* "showcomment" */
	sscanf(p, "%d", &showcomment);
      else if (!strcmp(fname,"menu"))			/* "menu" */
        {
	  sscanf(p, "%d", &i);
	  set_menu_visibility(i);
        }
      else if (!strcmp(fname,"pscolor"))		/* "pscolor" */
	sscanf(p, "%d", &pscolor);
      else if (!strcmp(fname,"psonly"))			/* "psonly" */
	sscanf(p, "%d", &ps_only);
      else if (!strcmp(fname,"pstitle"))		/* "pstitle" */
	sscanf(p, "%d", &pstitle);	      
      else if (!strcmp(fname,"psfontsize"))		/* "psfontsize" */
        {
	  sscanf(p, "%f", &psfontfac);
	  if (psfontfac < (float).5) psfontfac = (float).5;
	  else if (psfontfac > (float)1.5) psfontfac = (float)1.5;
        }
      else if (!strcmp(fname,"axiswidth"))		/* "axiswidth" */
        {
	  sscanf(p, "%d", &axis_width);
	  if (axis_width > 4) axis_width = 4;
	  else if (axis_width < 0) axis_width = 1;
        }
      else if (!strcmp(fname,"curveswidth"))		/* "curveswidth" */
        {
	  sscanf(p, "%d", &curves_width);
	  if (curves_width > 8) curves_width = 8;
	  else if (curves_width < 0) curves_width = 1;
        }
      else if (!strcmp(fname,"pstest"))
        sscanf(p, "%d %d", &labeldy, &labelsep);

      else if (!strcmp(fname,"ncontour"))		/* "ncontour" */
	sscanf(p, "%d", &default_ncurve);
      else if (!strcmp(fname,"big"))			/* "big, was fillscreen" */
        { sscanf(p, "%d", &default_ncol);
	  default_ncol = default_ncol ? 2 : 3; }
      else if (!strcmp(fname,"xpolar"))			/* "xpolar" */
	sscanf(p, "%s", xpolar);
      else if (!strcmp(fname,"ypolar"))			/* "ypolar" */
	sscanf(p, "%s", ypolar);

      else if (!strcmp(fname,"debug"))			/* "debug" for limits */
	sscanf(p, "%d", &ldebug);
      else if (!strcmp(fname,"test_data"))	       /* "test_data" for data test */
	sscanf(p, "%d", &debug_mb);

      else if (!strcmp(fname,"zdebug"))			/* "zdebug" for zoom */
	sscanf(p, "%d", &zdebug);
      else if (!strcmp(fname,"debugindex"))		/* "debugindex" */
	sscanf(p, "%d", &debug_index);		/* >=0: index to savelevel*/
      else if (!strcmp(fname,"verify"))			/* "verify" */
	sscanf(p, "%d", &verifying);
      else if (!strcmp(fname,"autozoom"))		/* "autozoom" */
      {
	sscanf(p, "%d", &ztestoption);
	if (ztestoption > 1) ztest((int *)p, NULL);
      }
      else if (!strcmp(fname,"dialog"))
	{ sscanf(p, "%d", &dialogwindow); if (dialogwindow) dialogwindow=1; }
    }
  fclose(frd);
}

/*-----------------------------------------------------------------------------
|	get_inname
|	* parse command line args, get file suffix, filename
|	* read 1st line of draw*.in, get ftype
-----------------------------------------------------------------------------*/
static FILE *get_inname(int argc, char *argv[])
{
  int i,suffix;
  FILE *frd;
  char *p, filename[100], text[150];

  /*------ Command Line Args */

  *basepath = 0;

  for (i = 1, ftype = suffix = 0; i < argc; i++)
    {
      p = argv[i];
      if (*p++ == '-')				/* ------ Dash options in command line */
	{
	  if      (*p == 'c')  ftype = 1;
	  else if (*p == 'C')  pscolor = 1;
	  else if (*p == 'P')  ps_only = 1;
          else if (*p == 'B')  background =1;
	  else if (*p == 'p')  ps_suffix(p + 1);
	  else if (*p == 'P')  {ps_only = 1;}
	  else if (*p == 'd')  dialogwindow = 1;
  	  else if (*p == 'b')  { strcpy(basepath, p + 1); strcat(basepath,"/"); }
	}
      else suffix = i;
    }

  /*------ Open draw*.in */

  strcpy(filename, datapath);
  if (*prefix=='*' || *prefix=='?')
    {  if (!suffix) strcat(filename, "draw"); }
  else strcat(filename, prefix);
  if (suffix)
    strcat(filename, argv[suffix]);
  strcat(filename, ".in");			/* ------ Reading drawx.in */
  if (!(frd = fopen(filename, "r")))
    s_abort("Cannot open %s", filename,0,0);

  /*------ Read 1st line of draw.in */

  fgets(text, LNAME, frd);			// "filename" or "Type"
  version = (strncmp(text, "Type", 4)) ? 0 : 1;	// 1 if "Type"
  if (tellMe) xprintf ("**get_inname %s, version %d\n", filename, version);

  if (version == 1)
    {
      ftype = 0;
      multi_topology = 0;
      p = strchr(text, ':');			/* get type G,C,I etc */
      if (p != NULL)
	for (p++; *p == ' ' || *p == '\t'; p++);
      for(; *p && *p!='\n'; p++)
	{
	  if (*p>='a' && *p<='z') *p -= ('a'-'A');

	  if (ftype==6 && *p=='2') 		/* if find M, then 2, "multi_topology" */
	    multi_topology = 1;
	  else if (*p == 'C') ftype=1;		/* Contour data format */
	  else if (*p == 'G') ftype=0;		/* original Graphics format */
	  else if (*p == 'I') ftype=2;		/* multi-dim (Indexed) */
	  else if (*p == 'H') ftype=3;		/* 'G' with Header in .in */
	  else if (*p == 'L') ftype=4;		/* 'G' with 1 header rec */
          else if (*p == 'N') ftype=5;          /* 'N' Contour Gener */
          else if (*p == 'M') ftype=6;          /* 'M' Multi Block Cont. Gener */

	  else s_abort("Undefined char in 1st line %s\n",filename,0,0);
	}
    }
  return(frd);
}

/*-----------------------------------------------------------------------------
|   skip_to_arg -- skip to the n-th argument in line text (1st arg: n=1)
|	           if endprompt non-zero (eg ':') skip it first
|		   returns NULL if no more
-----------------------------------------------------------------------------*/
char *skip_to_arg(char *text, int n, char endprompt)
{
  register char *p;
  int i, open_quote;

  if (endprompt) p = strchr(text, endprompt);	/* skip past endprompt */
  else p = NULL;
  if (!p) p = text;
  else p++;

  for (i = 0; i < n; i++)
    {
      open_quote = 0;
      for (; is_white(*p); p++);		/* skip leading spaces */
      if (is_eol(*p)) return (NULL);
      if (i == n - 1) return (p);
      for (; !is_white(*p) && !is_eol(*p); p++)		/* skip over arg */
	{
	  if (*p=='"') for(p++; *p!='"' && !is_eol(*p); p++) ;
	}
    }
  return (NULL);
}

/*-----------------------------------------------------------------------------
|	readline
-----------------------------------------------------------------------------*/
static int readline(FILE *frd, char *fname, int size);
static int readline(FILE *frd, char *fname, int size)
{
  char *p, *q;
  int leading_done;
  for(*fname = 0; ;)
    {
      p = fname + strlen(fname);
      if (fgets(p, p-fname+size, frd)==0) return 0;
      q = strchr(fname, '\n');
      if (q) *q=0;
      for(q=p, leading_done=0; *q; q++)
	{
	  if (!leading_done && (*q==' ' || *q=='\t')) continue;
	  leading_done = 1;
	  *p++ = *q;
	}
      *p = 0;
      if (p==fname) return 0;
      if (*(p-1) != '\\') return 1;
      *(p-1) = 0;
    }
}

/*-----------------------------------------------------------------------------
|	arg_suffix
-----------------------------------------------------------------------------*/
void arg_suffix(char *p, int axis)
{
  char *q, *qw;
  for(qw=p+1; !is_white(*qw) && !is_eol(*qw); qw++) ;
  for(q=qw; q>=p; q--)
    if (*(q-1)>='0' && *(q-1)<='9') break;

  if (q < qw) for(; q<qw; q++) {
    if (*q == '.') cp->flags |= DOTS;
    if (*q == '&') cp->flags |= SPLINE;
    if (*q == 'L' && ftype==0) {
      cp->use_log |= axis;
    }
    //if (cp->flags & (DOTS|SPLINE)) *q =' ';
    *q = ' ';
  }
}

/*-----------------------------------------------------------------------------
|	read_one_or_two
-----------------------------------------------------------------------------*/
void read_one_or_two(char *p, float *arg1, float *arg2, int *use)
{
  char *q, *qs;
  int flag;
  flag = 0;
  q=strchr(p, ',');
  for(qs=p; *qs && *qs>' '; qs++) ;		/* qs=1st white */
  if (q && qs>q)				/* there's a comma, scan 1 or 2 values */
    {
      if (q > p  ) { sscanf(p,  "%f", arg1); flag |= 1; }
      if (qs> q+1) { sscanf(q+1,"%f", arg2); flag |= 2; }
    }
  if (use != NULL) *use = flag;
}

/*-----------------------------------------------------------------------------
|	init
-----------------------------------------------------------------------------*/
void init(int argc, char *argv[])
{

  char fname[LLNAME];
  char *nameptr, *p, *q, *qs;
  VARIABLE_ID *qv;
  CVAR *vp;
  int i,j;
  int more, next, already;
  int usecolor;
  char code, *p1, *p2, c;
  FILE *frd;
  int win_at=0;
  VIEW *v;

/*------ Read state.xdr or state.xdraw, open drawx.in */

  give_command_args(argc, argv);	/* save argc,argv for Motif */
  usecolor = 1;
  readstate();				/* Read state.xdr or state.xdraw */
  frd=get_inname(argc,argv);		/* open eg."draw1.in", read ftype */

  if (ftype==0)				// Read more header info if any
    read_asc_header(frd);		// e.g. outer loop stuff

  for(;;)				/* skip to "filename(s)" */
    {
      fgets(fname,LNAME,frd);
      if ((int)strlen(fname) > 1) break;
    }

/*------ read binary file names and open input binary files */
  more =0;  /*DD*/
  for (i=0,buf=NULL,bufsize=0;;i++)		/* read data file name(s) */
    {
      if (*basepath) strcpy(fname, basepath);	/* from -b */
      else strcpy(fname, datapath);
      //printf (" filename = %s  \n",fname);
      p = fname + strlen(fname);
      fgets(p,LNAME,frd);
      if (*p == ';') continue;			// ignore "commented-out" filename
      if ((int)strlen(p) <= 1) break;		// done if blank line
      *strchr(p, '\n') = 0;
      if (*p=='/')
	strcpy(fname,p);
      else if (*basepath && (q=strchr(p,'/')))
      {
	 for(q++; strchr(q,'/'); q=strchr(q,'/')+1);
	 strcpy(p, q);
      }

      if (i==0) ps_dataname(p);
      
      //printf("call binread for %s\n", fname);
      strcpy(bin_file_name, fname);
      buf = binread(fname, &bufsize);		/* read data (.bin file) */

      if (multi_topology==2) {			// read first block of multi topology
	if (n_m2==0 || i_m2<0) exit(1);		// (prev call merely initialized arrays) 
	buf = binread(fname, &bufsize);
      }

      if (*basepath && i>0) continue;		/* -b: only read 1 file */
      if (buf == NULL) {
        printf("Abort!\n");
	exit(1);
      }

      //------ Bookkeeping for multiple input data files

      //printf (" Nnode = %d  inode = %d \n",nnode,inode);
      binfile[files_count].curve_begin = more;/* number of curve beg in file*/
      binfile[files_count].curve_end = nnode-1;
      strcpy(binfile[files_count].name, fname);/*file name */
      more = nnode;
      files_count++;
    }
  if (outerloop_added) use_asc_header();
  loop_structure();

  /*------ Read rest of drawx.in */

  nameptr = namebuf;
  if (version == 0)				/* Read usecolor */
    freadf(frd, "%d", (char *) &usecolor);
  freadf(frd, "%s", title);			/* Read global title */
  for(p=title+strlen(title)+1,more=1;;more=0)
  {
    fgets(nameptr,LNAME,frd);
    if (blankline(nameptr)) break;
    if (more) { strcat(title, "\n"); p++; }
    strcpy(p, nameptr);
    p += strlen(p) +1;
  }
  *p = 0;

  /*------ Read variable names ------*/

  for(;;)				/* skip to "independent var names" */
    {
      fgets(fname,LNAME,frd);
      if ((int)strlen(fname) > 1) break;
    }
			/* Skip "variable names" */
  for (nvar=already=0,qv=varid; ;nvar++,qv++)
    {
      if (!already) fgets(nameptr, LNAME, frd);
      p1 = skip_to_arg(nameptr, 1, '\0');	/* 1st arg=0,1,2.. in seq */
      if (p1 == NULL) break;
      p2 = skip_to_arg(nameptr, 2, '\0');	/* 2nd arg=name, eg. x,y,theta*/
      if ( ftype == 6 )       /*  ERROR handling needed!! */
	{
	  i=parsevar(p1)&~LOOPBIT;
	  if( i<0 ) exit (1); /* ERROR information */
	  qv = varid+i;
	}
      qv->index = parsevar(p1);
      qv->name = p2;
      m_lang_str(qv->name,qv->ml_text,&j);
      qv->nitems =j;

      already = (p2=strchr(p2, '\\')) != NULL;
      if (!already)
        {
	  *(p = strchr(nameptr, '\n')) = 0;
	  nameptr = p + 1;
        }
      else
        {
	  *p2 = 0;
          nameptr = read_flabels(frd, p2+1, qv->index);
        }
    }

  /*---- ftype=6: read info about dependent variables --------- */

  if ( ftype == 6)
    {
      /*------ skip to "dependent var names" ------*/

      for(;;) {
	fgets(fname,LNAME,frd);
	if ((int)strlen(fname) > 1) break;
      }

      /*------ Skip "variable names", get them ------*/

      for (nvar=3,already=0,qv=varid+3; ;nvar++,qv++)
	{
	  if (!already) fgets(nameptr, LNAME, frd);
	  p1 = skip_to_arg(nameptr, 1, '\0');	/* 1st arg=0,1,2.. in seq */
	  if (p1 == NULL) break;
	  p2 = skip_to_arg(nameptr, 2, '\0');	/* 2nd arg=name, eg. x,y,theta*/
	  
	  /*  ERROR handling needed!! */
	  i=parsevar(p1)&~LOOPBIT;
	  i+=2;
	  if( i<0 ) exit (1); /* ERROR information ?????*/
	  
	  qv = varid+i;
	  qv->index = i;
	  qv->name = p2;
	  m_lang_str(qv->name,qv->ml_text,&j);
	  qv->nitems =j;
	  
	  already = (p2=strchr(p2, '\\')) != NULL;
	  if (!already)
	    {
	      *(p = strchr(nameptr, '\n')) = 0;
	      nameptr = p + 1;
	    }
	  else
	    {
	      *p2 = 0;
	      nameptr = read_flabels(frd, p2+1, qv->index);
	    }
	}
      
    }	// .. done with variable names for ftype == 6

  /*------ read info about each plot ------*/

  for(;;)		/* skip to "iqty option title or ix iy" */
    {
      fgets(fname,LNAME,frd);
      if ((int)strlen(fname) > 1) break;
    }

  for(ncset=more=0; ncset<NCS; ncset++)
    {
      cp=curveset+ncset;
      if (more==0)
	{
	  if (!readline(frd, fname, LLNAME)) break;
	  next = 2;
	}
      if (*fname==';')			/* skip comment line */
	{ ncset--; continue; }
      p = fname;

      if (*p>='0' && *p<='9') ;			// iqty must start with number or i or I
      else if (*p=='i' || *p=='I') ;		// abort if eg from old-style Type=C
      else
        s_abort("Error specifying x-variable for graph spec:\n%s\n",
		fname,ncset,0);

      /*---- Initialize cp */

      cp->title = cp->subtitle = NULL;
      cp->label = NULL;
      cp->window = 0;
      cp->flags = 0;
      cp->mcount = -1;
      cp->gtype = 'G';				/* if contour, need -c */
      cp->multiDots = 0;
      
      cp->which = (char)(more? '2' : '1');
      cp->lstep = cp->lfaml = -1;
      cp->window = 0;
      cp->fskip = 1;
      for(i=0; i<nloop; cp->param[i++]=0) ;
      cp->iz.index = -1;
      cp->forcemin = cp->forcemax = (float)0;

      cp->f_block = 0;
      cp->l_block = 0;
      cp->iqty     =0;
      cp->itime    =0;
      cp->itime_rel = cp->itime_abs = 0;
      cp->exclmax = cp->exclmin = 0.0;
      cp->xcutoff = cp->ycutoff = 0.0;
      cp->use_cutoff = 0;
      cp->view_number=0;    /* window: -1 - nowin, numb -nomb 0f view*/
      win_at = (STATUS | MENU | VISIBILITY);
      cp->axis_width = axis_width;
      cp->curves_width = curves_width;	
      cp->vDensity = 1;
      cp->use_log = 0;
      vp = &cp->ix;
      p=fname;

      /*------ Read ix (init stuff if M) */

      arg_suffix(p, 1);

      if(ftype == 6)
	{
	  cp->gtype = 'M';
	  vp->index = 1|LOOPBIT;
          vp = &cp->iy;
          vp->index = 2|LOOPBIT;
          next =1;
	}
      else
	{   
	  vp->index = parsevar(p);
	  vp->hloop = v_is_header(vp->index);
	}

      /*------ Read iy suffix (or iqty if M) */

      p = skip_to_arg(fname,next,0);
      arg_suffix(p, 2);

      /*------ Read iy (n/a for M) */

      if( ftype != 6 )
	{
	  vp = &cp->iy;
	  vp->index = parsevar(p);
	  vp->hloop = v_is_header(vp->index);
	  testvar(cp);
	}

      else		/* Type 6: test if M --> V */
	{
	  sscanf(p,"%d",&cp->iqty);
	  cp->iqty--;
          next = 2;
	  p = skip_to_arg(fname,next,0);
	  if (p!=NULL && *p!='-' && *p!=QUOTE) /* Vplot option if exists iqty2, e.g. "26 27"*/
	    {
	      cp->gtype = 'V';
	      sscanf(p,"%d",&cp->iqty2);
	      cp->iqty2--;
	      cp->fskip =1;
	      cp->lstep =1;
	    }
	  else
	    {
	      next--;
	    }
	}
	  
      for(j=next+1,more=0;;j++)			/* If more iy's,get... */
        {					/* ...more,next for next pass */
	  p = skip_to_arg(fname,j,0);
	  if (p==NULL || *p=='-' || *p==QUOTE) break;
	  if (*p=='i' || *p=='I');
	  else if (*p>='0' && *p<='9');

	  else break;
	  if (!more) { more=1; next=j; }
	}

      /*---- Scan parameters (dash options) for this cp */

      initparam(NULL);				/* param[i]=0 */
      if (p && *p=='-')				/*........ any dash options? */
	for(;;j++)				/* get next dash option */
        {
	  p = skip_to_arg(fname,j,0);
	  if (p==NULL || *p != '-') break;
	  p++;					/* p on 1st char after dash*/
	  code = *p++;				/* p on 1st char after code*/
	  /*printf("curve option %c\n", code);*/
	  if	  (code=='s') cp->lstep = parsevar(p) & INDBITS;
	  else if (code=='f') cp->lfaml = parsevar(p) & INDBITS;
	  else if (code=='j') sscanf(p,"%d", &cp->fskip);
	  else if (code=='a') cp->flags |= ASPECT;
	  else if (code=='C') cp->flags |= COLORS;
	  else if (code=='p') cp->flags |= POLAR;
	  else if (code=='M')
	    { if (isdigit(*p)) sscanf(p,"%d", &cp->mcount); else cp->mcount=0; }
	  else if (code=='x') axis_label(p, &cp->ix);
	  else if (code=='y') axis_label(p, &cp->iy);
	  else if (code=='z') cp->iz.index = parsevar(p) & INDBITS;
	  else if (code=='T')
	    nameptr = set_subtitle(p, cp, code, nameptr);
	  else if (code=='X') { sscanf(p,"%f",&cp->force); cp->flags|=FORCED; }
	  else if (code=='E') 		/* -E<min>,<max>, can have either... */
	    {				/* ...absent, or max < min */
	      read_one_or_two(p, &cp->forcemin, &cp->forcemax, NULL);
	    }
	  else if (code=='g')		/* -g<min>,<max>, can have either... */
	    {				/* ...absent, or max < min */
	      read_one_or_two(p, &cp->exclmin, &cp->exclmax, NULL);
	    }
	  else if (code=='N')		/* -N<min>,<max>, can have either... */
	    {				/* ...absent, or max < min */
	      read_one_or_two(p, &cp->xcutoff, &cp->ycutoff, &cp->use_cutoff);
	    }
	  else if (code=='e')   sscanf(p,"%d",&cp->hzcurve);
	  else if (code=='L')   sscanf(p,"%d",&cp->ncurve);

	  /* NN multiblock*/
	  else if (code == 'B') {
	    sscanf(p,"%d",&cp->f_block);
	    cp->l_block = cp->f_block;
	  }

	  else if (code == 'q')
	    { sscanf(p,"%d",&cp->iqty);}
	  else if (code == 't')
	    { sscanf(p,"%d",&cp->itime);}

	  else if (code=='F') cp->flags |= FILL;
	  else if (code=='A') { sscanf(p,"%d",&cp->ncurve);/* timestep here */
				cp->flags |= ANIMATE; }
	  else if (code=='#') { sscanf(p,"%d",&cp->ncurve);
				cp->flags |= SINGLE; }
	  else if (code=='r') cp->flags |= REVERSE;

	  else if (code=='c')		/* can be -cS (solid), -cD (dots) etc */
	    {
	      if (ftype == 6) cp->multiDots = 1;
	      else {
		cp->gtype = 'C';
		c = toupper(*p);
		if (c=='S' || c=='D' || c=='L' || c=='N') 
		  cp->gtype = *p++;
		if (cp->iz.index < 0) cp->iz.index = 0;
	      }
	    }

	  else if (code=='i') set_indexed_value(p, fname);
	  else if (code=='l')			/* labels on curves */
	    {
	      cp->flags |= LABELF;
	      if (*p > ' ') init_label(cp, p, 0);
	    }
	  else if (code=='b' && *p>' ')		/* boxed labels */
	    {
	      cp->flags |= LABELF;
	      init_label(cp, p, 1);
	    }

	  /* Window attributes settings */
          else if (code =='m') win_at^= MENU;
	  else if (code =='u') win_at^= STATUS;
          else if (code =='v') cp->view_number=-1; /* non visible*/
	}

      getparam(cp);			/* choose lstep, set cp->param[] */
     
      if ((ftype != 6) && !curvetest(cp))		/* test if all OK */
	d_abort("Illegal combination found in curve %d",cp-curveset,0);

      if (p /*&& cp->which=='1'*/)		/* Read title, if any */
	{
	  if (*p == QUOTE) p++;			/* skip quote char */
	  strcpy(nameptr, p);
	  p = nameptr;
	  nameptr = strchr(nameptr, 0) + 1;
	}
      else if (cp->title != NULL)
	p = cp->title;
      else if (title != NULL)			/* no explicit title, use.. */
	p = title;				/* .. general label */
      else p = dummytitle;
	       
      /*if (cp->which=='1')*/ cp->title = p;

      if (cp->flags & POLAR)
	{
	  if (cp->ix.label==NULL) cp->ix.label = xpolar;
	  if (cp->iy.label==NULL) cp->iy.label = ypolar;
	}

      get_limits(cp);
      //printf("After get_limits, extrema %g, %g\n", cp->iz.min, cp->iz.max);

      /* Window check & install*/
      if (cp->which =='2') 
        {
          cp->view_number=(cp-1)->view_number;
        }
      else if ((nwindow < MAXWINDOW) && (!cp->view_number))
        {
	  cp->view_number = nwindow;
	  win_at |= VISIBILITY;
          v=view+nwindow++;
	  v->window_flags = win_at;
	  v->clipped = 0;
	  v->zoom_limit  = 0;
	  /* !! Add SIZE + POZITION - later */
        }
    }
  cp2 = cp;
  ncset = cp2 - curveset;
  nwindow =0;
}

/*-----------------------------------------------------------------------------
|	set_subtitle
-----------------------------------------------------------------------------*/
static char *set_subtitle(char *p, CURVE_SET *cp, char code, char *nameptr)
{
  char *q, delim;
  int n;
  if (*p==' ' || *p=='\t') return p;
  if (*p=='"') delim = *p++;
  else delim = '\0';
  for(q=p++;; p++)
    {
      if (*p==delim || *p=='\n') break;
      else if (!delim && (*p==' ' || *p=='\t')) break;
    }
  n = p - q;
  strncpy(nameptr, q, n);
  *(nameptr + n) = '\0';
  if (code=='t') cp->subtitle = nameptr;
  else cp->title = nameptr;
  p = nameptr + n + 1;
  return p;
}

/*-----------------------------------------------------------------------------
|	set_indexed_value
-----------------------------------------------------------------------------*/
static void set_indexed_value(char *p, char *fname)
{
  char *q, delim, text[20];
  int index, maxind;

  q = strchr(p,':');
  if (q==NULL) q=strchr(p,'=');
  if (q==NULL) d_abort("-i option in curve %d must include a (:) or (=)",
		                    cp - curveset,0);
  delim = *q; *q = 0;
  sscanf(p, "%d", &index);
  sprintf(text, "%s%c", p-2, delim);
  *q = delim;
  maxind = (inode==-1) ? nloop-1 : inode-1;
  if (index<0 || index>maxind)
    {
      sprintf(fname,"%s option in curve %d uses bad index %d (max value %d)",
	      text,cp - curveset,index,maxind);
      s_abort("%s",fname,0,0);
    }
  sscanf(q+1, "%d", &param[index]);
  if (param[index] >= loop[index].count)
    {
      sprintf(fname,"Illegal value %s%d option in curve %d (max value %s%d)",
	      text, param[index],cp-curveset,text,loop[index].count-1);
      s_abort("%s",fname,0,0);
    }
}

/*-----------------------------------------------------------------------------
|	v_is_header
-----------------------------------------------------------------------------*/
int v_is_header(int index)
{
  int i,j;
  LOOP *lp;
  if (index & LOOPBIT) return(-1);
  for(i=0,lp=loop; i<nloop; i++,lp++)
    {
      j = index - (int)lp->ih0;
      if (j < 0 ||  j >= (int)lp->hcount) continue;
      return(i);
    }
  return(-1);
}

/*-----------------------------------------------------------------------------
|	use_asc_header
-----------------------------------------------------------------------------*/
static void use_asc_header()
{
  LOOP *lp;
  int i, n, total;
  for(i=total=1; i<=outerloop_added; i++)
  {
    n = loop[i].count = outerloop_count[i];
    total *= n;
  }
  loop[0].count /= n;
}

/*-----------------------------------------------------------------------------
|	read_asc_header
-----------------------------------------------------------------------------*/
static void read_asc_header(FILE *frd)
{
  char text[150], *p, *q;
  int index, value, i;
  force_invert_bytes = -1;

  for(;;)
    {
      fgets(text, 150, frd);
      p = strchr(text, ':');
      if (p==NULL) return;

      if (!strncmp(text, "InvertBytes", 11)) 
	sscanf(p+1, "%d", &force_invert_bytes);

      else if (!strncmp(text, "Outer", 5))
        {
          sscanf(p+1, "%d", &outerloop_added);
	  for(i=0; i<=outerloop_added; i++)
	    outerloop_count[i] = i ? 1 : 0;
        }

      else if (!strncmp(text, "Loop", 4) && outerloop_added)
        {
	  for(p++; p && (q=strchr(p,'='));)
	    {
	      *q = ' ';
	      sscanf(p, "%d %d", &index, &value);
	      p = strchr(q+1, ' ');
	      if (index <= outerloop_added && value > 0)
		outerloop_count[index] = value;
	    }
	}
    }
}

/*=============================================================================
**			Edit
**===========================================================================*/
#define outtext(s) printf(s); fflush(stdout);
void settextcolor(int i) { i++; }

#define TEXT0 2
#define TEXT1 8
#define TEXT2 14

/*-----------------------------------------------------------------------------
|	curvetest
-----------------------------------------------------------------------------*/
int curvetest(CURVE_SET *cp)
{
  int indexf, error, i;
  CVAR *xp, *yp, *zp;
  char text[100];

  settextcolor(TEXT2);

/*------ Error test */
  error = 0;
  xp = &cp->ix;
  yp = &cp->iy;
  zp = &cp->iz;
  if (cp->gtype != 'G')
    {
      if (!(xp->index & LOOPBIT) || !(yp->index & LOOPBIT))
	error |= 0x10;
      if (zp->index & LOOPBIT) error |= 20;
    }
  else
    {
      indexf = cp->lfaml | LOOPBIT;
      if (indexf==xp->index) error |= 1;
      if (indexf==yp->index) error |= 2;
      if ((xp->index & LOOPBIT) && (yp->index & LOOPBIT)) error |= 4;
    }

  for(i=1; i<0x40; i*=2)
    {
      printvar(0, NULL, text);			/* Init text */
      if (!(error & i)) continue;
      printvar(-1, "WARNING! ", text);		/* Append to text */
      if (i & 3)
	{
	  printvar(indexf, "%s is used for both family and ", text);
	  if (i==1) printvar(-1, "X.\n", text);
	  else if (i==2) printvar(-1, "Y.\n", text);
	}
      else if (i==4)
	{
	  printvar(yp->index, "Y and X (%s and ", text);
	  printvar(xp->index, "%s) are both index variables.\n", text);
	}
      else if (i==0x10)
	strcat(text,"Contour x and y must both be index variables.\n");
      else if (i==0x20)
	strcat(text,"Contour must have a non-index z specified.\n");
      else continue;
      outtext(text);
    }

  settextcolor(TEXT0);
  return(error==0);
}

/*-----------------------------------------------------------------------------
|	printvar -- create a section of newtitle[]
|	using label assoc. w/ 'index'.  format==NULL ==> initialize
-----------------------------------------------------------------------------*/
int printvar(int index, char *format, char *newtitle)
{
  int i;
  VARIABLE_ID *q;
  char s[10],name[30];
  char *p;

  if (format==NULL) { *newtitle=0; return(0); }
  p = newtitle + strlen(newtitle);
  if (index==-1) { strcpy(p,format); return(0); }

  if (index & LOOPBIT) sprintf(s,"i%d",index & INDBITS);
  else sprintf(s,"%d",index);

  for(i=0,q=varid; i<nvar && q->index!=index; i++,q++) ;
  if (i==nvar) sprintf(name,"<%s>",s);
  else strcpy(name, (varid+i)->name);

  sprintf(p, format, name);
  return(i);
}


int m_lang_str(char *str,XTextItem *pnt_lang,int *n_items)
{
  int lang=LATIN;
  char      *tmp, *t_beg;
/*  XTextItem *pnt_lang;*/
  int end_of_str =0;
  int i;
  int j=0;
  int k=0;/* number of items*/
  char c1,c2;
  int res;

/* count the number of lang items */
  tmp=str;
 
  k=(*str=='$')?0:1;
  c1='a'; /* any */
/*  many's $ are ignored */
  c2=*str;
  if( c2 =='$') 
    {
      lang = GREEK;
      while ( *tmp == '$' ) tmp++; /* $ filtering */
    }

  while( c2 !='\0' )
    {
      if((c2 =='$') && c1 !='\\' &&(*(tmp+1)!='$')
	 &&(*(tmp+1)!='\0'))
	 k++;
      c1=c2; tmp++;
      c2=*tmp;
    }
/*   if (!k) return 0;  */ /* Only LATIN coding */
/*-----------------------------------------------*/
/* ------------------------------------------------*/
/*  create a tXTEXTITEM mas                       */
/*-----------------------------------------------*/
/*  pnt_lang = (XTextItem *) malloc(k*sizeof(XTextItem));*/
  lang = LATIN;
  tmp = t_beg =str;
  i=0;
  c2 = *str;
  if( c2 =='$') 
    {
      lang = GREEK;
      while ( *str == '$' ) str++; /* $ filtering */
    }
  while (!end_of_str && (i!=k))
    {
      t_beg = tmp = str;
      pnt_lang[i].chars = str;
      pnt_lang[i].delta = (lang==LATIN)?LATIN:GREEK;
      pnt_lang[i].font = (lang==LATIN)?font_latin :font_greek;
      j=0;
      while( (*str !='$')&&(*str!='\0'))
	    {
	      if (*str != '\\')
		{
		  *tmp++ =*str++;
		}
	      else
		if (lang == GREEK)
		  {
		    str = to_greek(str,&res);
		    if(res)
		      *tmp++=res;
		    else
		      *tmp++ =*str++;
		  }
		else
		  *tmp++=*str++;
	      j++;
	    }
      pnt_lang[i].nchars =j;
      while(*str =='$')str++;
      lang=(lang==LATIN)?GREEK:LATIN; i++;
      if(*str == '\0') end_of_str = 1;
    }
  *n_items =k;
  return k;
}

/*-----------------------------------------------*/
char *to_greek(char *str,int *res)
{
  int i=0;
  
  *res =0;
  for ( i=0; i<L_TBL; i++)
    {
      if(!strncmp(str,table[i].token,table[i].length))
	{
	  *res =table[i].code;
	  str +=table[i].length;
	  break;
	}
    }
  return str;	  
}      
