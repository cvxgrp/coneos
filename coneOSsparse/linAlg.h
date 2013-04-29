// x = b*a
extern void setAsScaledArray(double *x, const double * a,const double b,int len);

// a*= b
extern void scaleArray(double * a,const double b,int len);

// x'*y
extern double innerProd(const double * x, const double * y, int len);

extern double calcNormDiff(double *a, double *b, int l); 


// ||v||_2^2
extern double calcNormSq(const double * v,int len);

// ||v||_2
extern double calcNorm(const double * v,int len);

// ||v||_inf
extern double calcNormInf(const double *v, int len);

// saxpy
extern void addScaledArray(double * a, const double * b, int n, const double sc);
