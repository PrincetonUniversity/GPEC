/*-----------------------------------------------------------------------------
|	setcolor.h
-----------------------------------------------------------------------------*/
extern void set_xor(GC gc, int enable);
extern	void setcolor(int k, int ncurve, int how);
extern	void set_color_table(int k, int ncurve);
extern	float fancycolor(int k, int n,int i);
extern	void testpalette(int x, int y, int dx, int dy);
extern unsigned long getcolor(int r, int g, int b);
extern unsigned long  hlsrgb(float hue, float light, float sat);
static float value(float n1, float n2, float hue);
extern int getgoodfont (char *fonttitle,int lowpoint,int highpoint,
			char *searchfontname);
int SetShades(double hue, double chroma,long *pixels, int n);
int AllocShades(long *pixels, int n, int k);




