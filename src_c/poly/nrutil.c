#include "../ptv.h"


#ifndef NRUTIL_C
#define NRUTIL_C
#define NR_END 1
#define FREE_ARG char*

static double sqrarg, swap;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

#ifndef _INC_MATH
    #include <math.h>
#endif

#define NRANSI
#include "nrutil.h"

void gaussj(double **a, int n, double **b, int m)
{
    int *indxc,*indxr,*ipiv;
    int icol=0;
    int irow=0;
    int i,j,k,l,ll;
    double big,dum,pivinv,temp;

    indxc=ivector(1,n);
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if (ipiv[j] != 1)
                for (k=1;k<=n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    } else if (ipiv[k] > 1) 
                        nrerror("gaussj: Singular Matrix-1");
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
            for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a[icol][icol] == 0.0) 
            nrerror("gaussj: Singular Matrix-2");
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=1;l<=n;l++) a[icol][l] *= pivinv;
        for (l=1;l<=m;l++) b[icol][l] *= pivinv;
        for (ll=1;ll<=n;ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
            }
    }
    for (l=n;l>=1;l--) {
        if (indxr[l] != indxc[l])
            for (k=1;k<=n;k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 45%. */


void covsrt(double **covar, int ma, int ia[], int mfit)
{
    int i,j,k;
    double temp;

    for (i=mfit+1;i<=ma;i++)
        for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit;
    for (j=ma;j>=1;j--) {
        if (ia[j]) {
            for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
            for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
            k--;
        }
    }
}
#undef SWAP

/* (C) Copr. 1986-92 Numerical Recipes Software 45%. */


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *nrvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_vector(double *v, long nl)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
  nch =0; nrh = 0;
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


void lfit(double x[], double y[], double sig[], int ndat, double a[], int ia[],
	      int ma, double **covar, double *chisq, 
		  void (__cdecl *funcs)(double index, double *iodata, int iosize) )

{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	void gaussj(double **a, int n, double **b, int m);
	int i,j,k,l,m,mfit=0;
	double ym,wt,sum,sig2i,**beta,*afunc;

	beta=dmatrix(1,ma,1,1);
	afunc=dvector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) nrerror("lfit: no parameters to be fitted");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}
	for (i=1;i<=ndat;i++) {
		funcs(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/SQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	gaussj(covar,mfit,beta,1);
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=beta[++j][1];
	*chisq=0.0;
	for (i=1;i<=ndat;i++) {
		funcs(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,ia,mfit);
	free_dvector(afunc,1,ma);
	free_dmatrix(beta,1,ma,1,1);
}

/* (C) Copr. 1986-92 Numerical Recipes Software 45%. */

#endif // NRUTIL_C
