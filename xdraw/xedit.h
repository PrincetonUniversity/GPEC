/******************************************************************************
**	NAME		XEDIT.H
**	AUTHOR		Sheryl M. Glasser
**
**	DESCRIPTION
**		
**
**	Copyright (c) Toptools SCF 1990.  All rights reserved.
******************************************************************************/
extern int editcurve(int iwin, CURVE_SET *cp, char *mname2);
extern void print_title(int iwin, CURVE_SET *p);
extern char *get_descriptive_title(int iwin, CURVE_SET *p);
extern void init_label(CURVE_SET *, char *, int);

struct LABEL *get_label(int k, struct LABEL *q, int *lx, int *ly,
	       char **text, int *nline, int *align,
	       FloatSc xscale, FloatSc xoffset, FloatSc yscale, FloatSc yoffset,
	       int *nitems, XTextItem **item);


extern void get_labelbox(int *, char **);
extern char *get_newline(char *s);








