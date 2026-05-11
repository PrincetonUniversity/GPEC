#define Float double

void splfit(Float *xs, Float *fs, Float *fs1, int m);
void spleval(Float *xs, Float *fs, Float *fs1,
	     Float *f, Float *f1, Float *f2, Float *f3,
	     Float x, int m, int mmax, int nqty, int mode);
