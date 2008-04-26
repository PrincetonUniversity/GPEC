/*-----------------------------------------------------------------------------
|	xcontour.h
-----------------------------------------------------------------------------*/
extern  void redraw1(CURVE_SET *cp,
		     float xmin, float xmax, float ymin, float ymax );
extern void redraw_mb(CURVE_SET *cp,
                     float xmin, float xmax, float ymin, float ymax);
extern	int new_ncurve(CURVE_SET *cp, char how);
extern  void contour_values(CURVE_SET *cp);
extern	void get_contlim(float *xlim, float *ylim);

#ifdef XCONTOUR
extern  void drawcontour(float psi,int ii);
extern  FLAG *surround(FLAG *ilist, POINTS *p);
extern	POINTS_ *has_psi(float psi, int ir, int iz, float x, float y,
                        POINTS_ *q);
extern	POINTS_ *same_as(POINTS_ *p);
extern	POINTS_ *nextpoint(int k,POINTS_ *q,POINTS_ *qln);
extern	double anglecmp(double a1, double a2);
extern	POINTS_ *save_linear(int ir,int mr,int n,QCOEFF *kp,POINTS_ *q);
extern	POINTS_ *save_extremum(int ir,int mr,int n,float psi,int how,
	QCOEFF *k,POINTS_ *q);
extern	void getquadcoeff(int n,int ir,int mr, QCOEFF *k);
extern  int getp00(int,int);
extern	float quadratic(QCOEFF *k, float psi, float x0);
extern	void savelevel(int how, int index, char *format);
extern  void savegrid(int ir,int iz);

extern void drawgrid(float xmin, float xmax, float ymin, float ymax);
extern void drawgrid_t(float xmin, float xmax, float ymin, float ymax,
                       int iblock,int mvert,int mcells,int xy_off);
 
extern void psicount(float *buf, long bufsize,
		     float psi0,float dpsi,int ncurve,float *count);
extern void  get_splr(float *x,float *y,
               int ir, int iz, float delta);
extern void  get_splz(float *x,float *y,
               int ir, int iz, float delta);
extern int FindPoint( float x1, float y1, float z1,
	       float x2, float y2, float z2,
	       float z,
	       float *xp, float *yp   );
extern void drawcontour_t(float psi, int ii,int xy_off,int f_off,
		   int mvert,int mcells);
int new_nvect(CURVE_SET *cp, char how);
void draw_v(int mr,int mz,int x_off,int y_off,int q1_off,int q2_off,
	    float *lmax, int *npoints, int *zpoints,
	    float l_scale, int mod, int density);
void vector(int ix1, int iy1, float vx, float vy, 
	    float lmax, float l_scale, int mod);

#define ptFlt float
void test_this_segment(float psi, ptFlt x2, ptFlt y2, ptFlt x1, ptFlt y1,
		       int ir2, int iz2, byte where2,
		       int ir1, int iz1, byte where1);
void transfer_distmin();
void tell_contour_value(float *px, float *py, float *pz, int mr,
			float *psixy, float *psi_seg,  
			float *gradx, float *grady);

#endif
void set_extra_psi(int method, float psi);
void enable_contour_test(int);

typedef struct
{
  double distmin;
  double psi;
  int irmin1, izmin1, irmin2, izmin2;
  byte where_min1, where_min2;
  double xmin1, xmin2, ymin1, ymin2;
  double xNearest, yNearest;
  BLOCK_XY *block_xy;
} NEAR_SEG;











