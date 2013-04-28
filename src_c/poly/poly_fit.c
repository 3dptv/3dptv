/***************************************************************************
Package: Calibration and mapping functions, based on polynomial fitting.  

Author:    Ad Holten

Address:  Department of Applied Physics
          Eindhoven University of Technology
          Den Dolech 2
          5612AZ Eindhoven, The Netherlands

Creation Date:  April 2010

****************************************************************************/

// most declarations can be found in polyfit.h

#include "../ptv.h"
#include <math.h> 
#include "nrutil.h"

FITGLOB fitGlob;

void   save_3Dfitparameters3(char* pathname, char* mode, char* title, FitInfo *pXfit, FitInfo *pYfit);
void   build_indexedname(char* pathname, char* basename, int index);
double poly_xmin, poly_xmax, poly_ymin, poly_ymax, poly_zmin, poly_zmax;
int    poly_order, poly_crossorder;


void FitFunction2D(double index, double *iodata, int iosize)
{
    int order      = fitGlob.order;
    int crossorder = fitGlob.crossorder;
    double x, y;

    x = fitGlob.in2D[(int)index-1].x;
    y = fitGlob.in2D[(int)index-1].y;

    // getIOdata2D() returns iodata[] starting at index 0. 
    // However, this function is called by the numerical recipies function lfit()
    // which expects iodata[] to start at index 1.
    // By incrementing iodata, getIOdata2D() will start filling iodata[1].

    iodata++;    // for the iodata returned, iodata[1] will be the first parameter
    getIOdata2D(iodata, x, y, order, crossorder);
}

void InvFitFunction3D(double index, double *iodata, int iosize)
{
    int order      = fitGlob.order;
    int crossorder = fitGlob.crossorder;
    double X, Y, Z;

    X = fitGlob.pnt3[(int)index-1].x;
    Y = fitGlob.pnt3[(int)index-1].y;
    Z = fitGlob.pnt3[(int)index-1].z;

    // getIOdata2D() returns iodata[] starting at index 0. 
    // However, this function is called by the numerical recipies function lfit()
    // which expects iodata[] to start at index 1.
    // By incrementing iodata, getIOdata2D() will start filling iodata[1].

    iodata++;
    getIOdata3D(iodata, X, Y, Z, order, crossorder);
}

double dotmult(double* v1, double *v2, int nelem)
{
    int i;
    double result = 0;
    for (i = 0; i < nelem; i++)
        result += v1[i] * v2[i];
    return result;
}

// ===== 2D and 3D fitting tools ====

double Map2D_X(double x, double y, FitInfo* pfit)
{
    return Mapper2D(x, y, pfit->order, pfit->crossorder, 
                      pfit->parms1, pfit->nparms1);
}

double Map2D_Y(double x, double y, FitInfo* pfit)
{
    return Mapper2D(x, y, pfit->order, pfit->crossorder, 
                      pfit->parms2, pfit->nparms2);
}

int getNumParms2D(int order, int crossorder)
{
    int ix, iy, k, cnt;

    cnt = 2*order + 1;
    for (k = 0; k <= crossorder; k++)
        for (ix = 0; ix < k; ix++)
            for (iy = 0; iy < k; iy++)
                if ((ix + iy) == k)
                    cnt++;
    return cnt;
}


int getNumParms3D(int order, int crossorder)
{
    int ix, iy, iz, k, cnt;

    cnt = 3*order + 1;
    for (k = 0; k <= crossorder; k++)
        for (ix = 0; ix < k; ix++)
            for (iy = 0; iy < k; iy++)
                for (iz = 0; iz < k; iz++)
                    if ((ix + iy + iz) == k)
                        cnt++;
    return cnt;
}

    //   meaning of ia[] for order=3, crossorder=3
    //   1    1
    //   2    x
    //   3    x^2
    //   4    x^3
    //   5    y
    //   6    y^2
    //   7    y^3
    //   8    z
    //   9    z^2
    //  10    z^3
    //  11    y * z
    //  12    x * z
    //  13    x * y
    //  14    y * z^2
    //  15    y^2 * z
    //  16    x * z^2
    //  17    x * y * z
    //  18    x * y^2
    //  19    x^2 * z
    //  20    x^2 * y

void Fit3D_to_2D_body(FitInfo* pfi, coord_3d *pnt3, coord_2d *pnt2, int ndata)
{
	int    i, j, nparams, loop;
	double *x, *yu, *yv, *sig, chisq, dx, dy;
	int    *ia;
	double **covar, result, *xdata, *ydata, *params1, *params2, msdist;
	int    order, crossorder;

	order = crossorder = 3;			// working with different weight values now

	fitGlob.order      = order;
	fitGlob.crossorder = crossorder;
   	fitGlob.pnt3       = pnt3;
	fitGlob.pnt2       = pnt2;

	nparams = getNumParms3D(order, crossorder);

	// The lfit() function of the numerical recipies library expect arrays to
	// start at index 1, (the element at index 0 is never used) 
	x       = (double*) malloc(sizeof(double) * (ndata+1));	// index list
	yu      = (double*) malloc(sizeof(double) * (ndata+1));	// x-values to fit
	yv      = (double*) malloc(sizeof(double) * (ndata+1));	// y-values to fit
	sig     = (double*) malloc(sizeof(double) * (ndata+1));	// error in input values
	ia      = (int*)    malloc(sizeof(int)    * (nparams+1));
	covar   = (double**) dmatrix(1,nparams,1,nparams);

	// The rest of the program is using parameter lists starting from element 0.
	// By decrementing the pointers to the allocated param. lists, lfit will fill
	// the lists from index 0.
	params1 = (double*) malloc(sizeof(double) * nparams) - 1;
	params2 = (double*) malloc(sizeof(double) * nparams) - 1;

	for (i=1; i <= ndata; i++) {
		x[i]  = i;								
		yu[i] = pnt2[i-1].x;
		yv[i] = pnt2[i-1].y;
		sig[i] = 1.0;		
	}
	for (i=1; i <= nparams; i++) {
		params1[i] = 0.0;
		params2[i] = 0.0;
		ia[i]   = 1;
		for (j = 1; j <= nparams; j++) covar[i][j] = 0;
	}

   	// --- Fitting u and v ---
	if (usingZplanes) {	// now it is safe to fit with the requested order and cross order
		lfit(x, yu, sig, ndata, params1, ia, nparams, covar, &chisq, InvFitFunction3D);
		lfit(x, yv, sig, ndata, params2, ia, nparams, covar, &chisq, InvFitFunction3D);
	}
	else {
		//             1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7
		int cases[] = {0,1,0,1,0,1,0,0,1,2,3,0,1,2,3,0,1};
		// for (loop=0; loop<sizeof(cases)/sizeof(int); loop++) {
		for (loop=0; loop<7; loop++) {
			switch (cases[loop]) {
				case 0:
					for (i=1; i <= nparams; i++)
						ia[i] = 0;
					ia[ 1] = 1;     // const
					ia[ 2] = 1;     // x
					ia[ 5] = 1;     // y
					ia[ 8] = 1;     // z
					ia[11] = 1;     // y*z
					ia[12] = 1;     // x*z
					ia[13] = 1;     // x*y
					break;
				case 1:
					for (i=1; i <= nparams; i++)
						ia[i] = 0;
					ia[3] = 1;		// x^2
					ia[6] = 1;		// y^2
					ia[9] = 1;		// z^2
					break;
				case 2:
					for (i=1; i <= nparams; i++)
						ia[i] = 0;
					ia[14] = 1;		// z^2 * y
					ia[15] = 1;		// y^2 * z
					ia[16] = 1;		// z^2 * x
					ia[17] = 1;		// x * y * z
					ia[18] = 1;		// y^2 * x
					ia[19] = 1;		// z^2 * z
					ia[20] = 1;		// x^2 * y
					break;
				case 3:
					for (i=1; i <= nparams; i++)
						ia[i] = 0;
					ia[4]  = 1;		// x^3
					ia[7]  = 1;		// y^3
					ia[10] = 1;		// z^3
					break;
			}
       		// --- Fitting u and v ---
			lfit(x, yu, sig, ndata, params1, ia, nparams, covar, &chisq, InvFitFunction3D);
			lfit(x, yv, sig, ndata, params2, ia, nparams, covar, &chisq, InvFitFunction3D);
		}
	}
	msdist = 0;
	result = 0;
	xdata = (double*) dvector(1,nparams);
	ydata = (double*) dvector(1,nparams);
	for(i = 0; i<ndata; i++)	{
		InvFitFunction3D(i+1, xdata, nparams);
		InvFitFunction3D(i+1, ydata, nparams);
		dx = pnt2[i].x - dotmult(params1+1, xdata+1, nparams); 
		dy = pnt2[i].y - dotmult(params2+1, ydata+1, nparams);
		msdist += dx*dx + dy*dy;
	}
	free_dvector(xdata,1,nparams);
	free_dvector(ydata,1,nparams);

	free(sig);
	free(ia);
	free(x);
	free(yu);
	free(yv);
	free_dmatrix(covar,1,nparams,1,nparams);

	pfi->order      = order;
	pfi->crossorder = crossorder;
	pfi->nparms1    = nparams;
	pfi->nparms2    = nparams;
	pfi->parms1     = params1 + 1;
	pfi->parms2     = params2 + 1;
	pfi->result     = sqrt(msdist/ndata);

	printf("result 3D->2D fit = %.2lf\n", pfi->result);
}


void Fit3D_to_2D(FitInfo* pfi, int degree[2], coord_3d *pnt3, coord_2d *pnt2, int ndata)
{
	int    i, j, nparams;
	double *x, *yu, *yv, *sig, chisq, dx, dy;
	int    *ia;
	double **covar, result, *xdata, *ydata, *params1, *params2, msdist;

	fitGlob.order      = degree[0];
	fitGlob.crossorder = degree[1];
   	fitGlob.pnt3       = pnt3;
	fitGlob.pnt2       = pnt2;

	nparams = getNumParms3D(degree[0], degree[1]);

	// The lfit() function of the numerical recipies library expect arrays to
	// start at index 1, (the element at index 0 is never used) 
	x       = (double*) malloc(sizeof(double) * (ndata+1));	// index list
	yu      = (double*) malloc(sizeof(double) * (ndata+1));	// x-values to fit
	yv      = (double*) malloc(sizeof(double) * (ndata+1));	// y-values to fit
	sig     = (double*) malloc(sizeof(double) * (ndata+1));	// error in input values
	ia      = (int*)    malloc(sizeof(int)    * (nparams+1));
	covar   = (double**) dmatrix(1,nparams,1,nparams);

	// The rest of the program is using parameter lists starting from element 0.
	// By decrementing the pointers to the allocated param. lists, lfit will fill
	// the lists from index 0.
	params1 = (double*) malloc(sizeof(double) * nparams) - 1;
	params2 = (double*) malloc(sizeof(double) * nparams) - 1;

	for (i=1; i <= ndata; i++) {
		x[i]  = i;								
		yu[i] = pnt2[i-1].x;
		yv[i] = pnt2[i-1].y;
		sig[i] = 1.0;		
	}
	for (i=1; i <= nparams; i++) {
		params1[i] = 0.0;
		params2[i] = 0.0;
		ia[i]   = 1;
		for (j = 1; j <= nparams; j++) covar[i][j] = 0;
	}
   	// --- Fitting u and v ---
	lfit(x, yu, sig, ndata, params1, ia, nparams, covar, &chisq, InvFitFunction3D);
	lfit(x, yv, sig, ndata, params2, ia, nparams, covar, &chisq, InvFitFunction3D);

	msdist = 0;
	result = 0;
	xdata = (double*) dvector(1,nparams);
	ydata = (double*) dvector(1,nparams);
	for(i = 0; i<ndata; i++)	{
		InvFitFunction3D(i+1, xdata, nparams);
		InvFitFunction3D(i+1, ydata, nparams);
		dx = pnt2[i].x - dotmult(params1+1, xdata+1, nparams); 
		dy = pnt2[i].y - dotmult(params2+1, ydata+1, nparams);
		msdist += dx*dx + dy*dy;
	}
	free_dvector(xdata,1,nparams);
	free_dvector(ydata,1,nparams);

	free(sig);
	free(ia);
	free(x);
	free(yu);
	free(yv);
	free_dmatrix(covar,1,nparams,1,nparams);

	pfi->order      = degree[0];
	pfi->crossorder = degree[1];
	pfi->nparms1    = nparams;
	pfi->nparms2    = nparams;
	pfi->parms1     = params1 + 1;
	pfi->parms2     = params2 + 1;
	pfi->result     = sqrt(msdist/ndata);

	printf("result 3D->2D fit = %.2lf\n", pfi->result);
}

void draw_fiterror(Tcl_Interp* interp, double x1, double y1, int ic, FitInfo fit, coord_3d fix)
{
	int ix2, iy2;
	double x2, y2;

	x2 = fitdata_Map3D_x(fix.x, fix.y, fix.z, &fit);
	y2 = fitdata_Map3D_y(fix.x, fix.y, fix.z, &fit);
	
	ix2 = (int)(x1 + 10*(x2-x1));
	iy2 = (int)(y1 + 10*(y2-y1));

	drawvector(interp, (int)x1, (int)y1, ix2, iy2, 1, ic, "red");
}

void draw_fitted_fixpoints2d(Tcl_Interp* interp, int iimg, FitInfo fit, sort_point* fix, int nfix)
{
	int ifix, ix, iy;
	// char *dotcol[] = {"yellow", "green", "cyan", "red"};
	// k = loop > 3 ? 3 : loop;
	for (ifix=0; ifix<nfix; ifix++) {
		ix = (int) fitdata_Map2D_x(fix[ifix].x, fix[ifix].y, &fit);
		iy = (int) fitdata_Map2D_y(fix[ifix].x, fix[ifix].y, &fit);
        drawcross(interp, ix, iy, cr_sz+2, iimg, "yellow");
        draw_pnr(interp, ix, iy, fix[ifix].pnr, iimg, "yellow");
	}
}

void find_grid_points(coord_3d *fixori, coord_2d *pixori, int *pcnt, double eps, 
					  FitInfo fit, coord_3d* fix, int nfix, target* pix, int npix)
{
	double xp, yp, dx, dy, dmin, d;
	int ifix, jpix, used, n, j, k, jmn;

	k = 0; 
	for (ifix=0; ifix<nfix; ifix++) {
		xp = fitdata_Map3D_x(fix[ifix].x, fix[ifix].y, fix[ifix].z, &fit);
		yp = fitdata_Map3D_y(fix[ifix].x, fix[ifix].y, fix[ifix].z, &fit);
		dmin = 1e10;
		for (jpix=0; jpix<npix; jpix++) {	// look for the nearest target
			dx = xp - pix[jpix].x;
			dy = yp - pix[jpix].y;
			d = sqrt(dx*dx + dy*dy);
			if (d < dmin && d < eps) {
				dmin   = d;
				fixori[k] = fix[ifix];
				pixori[k].x = pix[jpix].x;
				pixori[k].y = pix[jpix].y;
				pixori[k].pnr = jpix;
				jmn = jpix;
			}
		}
		if (dmin < 1e10)
			k++;
	}
	n = k;

	// targets connected to the same 3D-point are thrown out of the list
	j = k = 0;
	for (ifix=0; ifix<n; ifix++) {
		used = 0;
		for (jpix=0; jpix<n; jpix++)
			if (pixori[ifix].pnr == pixori[jpix].pnr) used++;
		if (used == 1) {		// point ok
			fixori[k  ] = fixori[ifix];
			pixori[k++] = pixori[ifix];
		}
		else j++;
	}
	printf("Found: %d valid points, skipped %d\n", k,j);
	*pcnt = k;
}

void find_grid_points2d(coord_2d *fixori, coord_2d *pixori, int *pcnt, double eps, 
					  FitInfo fit, coord_2d* fix, int nfix, target* pix, int npix)
{
	double xp, yp, dx, dy, dmin, d;
	int ifix, jpix, used, n, j, k;

	k = 0;
	for (ifix=0; ifix<nfix; ifix++) {
		xp = fitdata_Map2D_x(fix[ifix].x, fix[ifix].y, &fit);
		yp = fitdata_Map2D_y(fix[ifix].x, fix[ifix].y, &fit);
		dmin = 1e10;
		for (jpix=0; jpix<npix; jpix++) {	// look for the nearest target
			dx = xp - pix[jpix].x;
			dy = yp - pix[jpix].y;
			d = sqrt(dx*dx + dy*dy);
			if (d < dmin && d < eps) {
				dmin   = d;
				fixori[k]   = fix[ifix];
				pixori[k].x = pix[jpix].x;
				pixori[k].y = pix[jpix].y;
				pixori[k].pnr = jpix;
			}
		}
		if (dmin < 1e10)
			k++;
	}
	n = k;

	// targets connected to the same 3D-point are thrown out of the list
	j = k = 0;
	for (ifix=0; ifix<n; ifix++) {
		used = 0;
		for (jpix=0; jpix<n; jpix++)
			if (pixori[ifix].pnr == pixori[jpix].pnr) used++;
		if (used == 1) {		// point ok
			fixori[k  ] = fixori[ifix];
			pixori[k++] = pixori[ifix];
		}
		else j++;
	}
	printf("Found: %d valid points, skipped %d\n", k,j);
	*pcnt = k;
}

void GetFitPixels_Planes(FitInfo* fit, double* Zlev, int nz, Rect* prect, FitInfo* pfit3d2d)
{
	double  xmin,xmax,ymin,ymax,zmin,zmax;
	double  X,Y,Z,xp,yp, dx,dy,dz, x0,y0;
	int     i,ix,iy,iz, nx,ny;
	coord_2d Pplane[11*11], Ppix[11*11];
	int     degree[2];
	BOOL    ok = TRUE;

	// measuring volume in mm
    xmin = X_lay[0];    xmax = X_lay[1];
    ymin = Y_lay[0];    ymax = Y_lay[1];
    zmin = Zmin_lay[0]; zmax = Zmax_lay[1];

	nx = 11; ny = 11;
	dx = (xmax-xmin)/nx;
	dy = (ymax-ymin)/ny;
	dz = (zmax-zmin)/nz;

	x0 = xmin;
	y0 = ymin;
	for (iz=0; iz<nz; iz++) {
		Z = zmin + iz*dz;
		i = 0;
		for (ix=0; ix<nx; ix++)
			for (iy=0; iy<ny; iy++) {
				X = x0 + ix*dx;
				Y = y0 + iy*dy;
                xp = fitdata_Map3D_x(X,Y,Z, pfit3d2d);
                yp = fitdata_Map3D_y(X,Y,Z, pfit3d2d);

				// bounding rect of the pixel coordinates
				if (i == 0) {
					prect->left = prect->right = xp;
					prect->top = prect->bottom = yp;
				}
				else {
					if (prect->left   > xp) prect->left   = xp;
					if (prect->right  < xp) prect->right  = xp;
					if (prect->bottom > yp) prect->bottom = yp;
					if (prect->top    < yp) prect->top    = yp;
				}
				Pplane[i].x = X;
				Pplane[i].y = Y;
				Ppix[i].x   = xp;
				Ppix[i].y   = yp;
				i++;
			}
		Zlev[iz] = Z;
		degree[0] = 3;		// order = Poly.planesfit.order;
		degree[1] = 3;		// cross = Poly.planesfit.cross;
		ok = Fit2D_to_2D(&fit[iz], degree, Ppix, Pplane, nx*ny);
	}
}

void GetFit_Pixels_Lines(FitInfo* pxFitinfo, FitInfo* pyFitinfo, FitInfo* fit, double* Zlev, int nz, Rect* prect)
{
	double  x1,x2,y1,y2,dx,dy, x,y;
	double  X,Y,Z;
	int     iz,nx,ny,cnt;
	int     degree[2] = {3,3};
	coord_2d *xin, *xout, *yin, *yout;

	double sumx, sumy, sumz, sumzx, sumzy, sumzz;

	x1 = prect->left;	x2 = prect->right;
	y1 = prect->bottom;	y2 = prect->top;
	
	// normalize boundaries
	if (x1 == x2) x2 = x1+100;
	if (y1 == y2) y2 = y1+100;
	if (x1 > x2) x=x2, x2=x1, x1=x;
	if (y1 > y2) y=y2, y2=y1, y1=y;

	nx = ny = 11;
	dx = (x2-x1)/(nx-1);
	dy = (y2-y1)/(ny-1);

	xin   = (coord_2d*) malloc(nx*ny * sizeof(coord_2d));
	xout  = (coord_2d*) malloc(nx*ny * sizeof(coord_2d));
	yin   = (coord_2d*) malloc(nx*ny * sizeof(coord_2d));
	yout  = (coord_2d*) malloc(nx*ny * sizeof(coord_2d));

 	cnt = 0;
	for (x=x1; x<=x2; x+=dx)
		for (y=y1; y<=y2; y+=dy)
	{
		// Get the line parameters (linear regression)
		double xof, xrc, yof, yrc;
	    sumz = sumx = sumy = sumzz = sumzx = sumzy = 0;
		for(iz=0; iz < nz; iz++) {
			X = Map2D_X(x, y, &fit[iz]);
			Y = Map2D_Y(x, y, &fit[iz]);
			Z = Zlev[iz];
			sumx  += X;
			sumy  += Y;
			sumz  += Z;
			sumzx += X*Z;
			sumzy += Y*Z;
			sumzz += Z*Z;
		}
		xrc = (nz*sumzx - sumz*sumx)/(nz*sumzz - sumz*sumz);
		xof = (sumx - xrc*sumz)/nz;		
		yrc = (nz*sumzy - sumz*sumy)/(nz*sumzz - sumz*sumz);
		yof = (sumy - yrc*sumz)/nz;		

		xin[cnt].x  = x;    yin[cnt].x  = x; 
		xin[cnt].y  = y;    yin[cnt].y  = y; 
		xout[cnt].x = xof;  yout[cnt].x = yof; 
		xout[cnt].y = xrc;  yout[cnt].y = yrc; 
		cnt++;
	}
	Fit2D_to_2D(pxFitinfo, degree, xin, xout, cnt);
	Fit2D_to_2D(pyFitinfo, degree, yin, yout, cnt);

	free(xin);
	free(xout);
	free(yin);
	free(yout);
}

void freefitparms(FitInfo* pfi)
{
	free(pfi->parms1);
	free(pfi->parms2);
}

BOOL Fit2D_to_2D(FitInfo* pfi, int degree[2], coord_2d* Pin, coord_2d* Pout, int ndata)
{
	int    i, j, nparams;
	double *x, *yu, *yv, *sig, chisq, dx, dy;
	int    *ia;
	double **covar, msdist;
	double *xdata, *ydata;
	double *params1, *params2;

	fitGlob.order      = degree[0];
	fitGlob.crossorder = degree[1];
	fitGlob.in2D       = Pin;
	fitGlob.out2D      = Pout;
	fitGlob.ndata      = ndata;

	nparams = getNumParms2D(degree[0], degree[1]);

	// The lfit() function of the numerical recipies library expect arrays to
	// start at index 1, (the element at index 0 is never used) 
	x       = (double*) malloc(sizeof(double) * (ndata+1));	// index list
	yu      = (double*) malloc(sizeof(double) * (ndata+1));	// x-values to fit
	yv      = (double*) malloc(sizeof(double) * (ndata+1));	// y-values to fit
	sig     = (double*) malloc(sizeof(double) * (ndata+1));	// error in input values
	ia      = (int*)    malloc(sizeof(int)    * (nparams+1));
	covar   = (double**) dmatrix(1,nparams,1,nparams);

	// The rest of the program is using parameter lists starting from element 0.
	// By decrementing the pointers to the allocated parameter lists, 
	// lfit will fill the lists from index 0.
	params1 = (double*) malloc(sizeof(double) * nparams) - 1;
	params2 = (double*) malloc(sizeof(double) * nparams) - 1;

	for (i=1; i <= ndata; i++) {
		x[i]  = i;								
		yu[i] = Pout[i-1].x;
		yv[i] = Pout[i-1].y;
		sig[i] = 1.0;		
	}
	for (i=1; i <= nparams; i++) {
		params1[i] = 0.0;
		params2[i] = 0.0;
		ia[i]      = 1;
		for (j = 1; j <= nparams; j++) covar[i][j] = 0;
	}
	// --- Fitting x ---
	lfit(x, yu, sig, ndata, params1, ia, nparams, covar, &chisq, FitFunction2D);
	// --- Fitting y ---
	lfit(x, yv, sig, ndata, params2, ia, nparams, covar, &chisq, FitFunction2D);

	msdist = 0;
	xdata = (double*) dvector(1,nparams);
	ydata = (double*) dvector(1,nparams);
	for(i = 0; i<ndata; i++)	{
		FitFunction2D(i+1, xdata, nparams);
		FitFunction2D(i+1, ydata, nparams);
		dx = Pout[i].x - dotmult(params1+1, xdata+1, nparams);
		dy = Pout[i].y - dotmult(params2+1, ydata+1, nparams);
		msdist += dx*dx + dy*dy;
	}
	free_dvector(xdata,1,nparams);
	free_dvector(ydata,1,nparams);

	free(sig);
	free(ia);
	free(x);
	free(yu);
	free(yv);
	free_dmatrix(covar,1,nparams,1,nparams);

	pfi->order      = degree[0];
	pfi->crossorder = degree[1];
	pfi->nparms1    = nparams;
	pfi->nparms2    = nparams;
	pfi->parms1     = params1 + 1;
	pfi->parms2     = params2 + 1;
	pfi->result     = sqrt(msdist/ndata);

	printf("result fit 2D->2D = %lf\n", pfi->result); 
	return TRUE;
}


// ==== 2D-point to 3D-lines ===

void getIOdata2D(double* iodata, double x, double y, int order, int crossorder)
{
	int i, ix, iy, k, co;

	// iodata++;		// nr functions are using indices starting from 1
	iodata[     0] = 1;
	iodata[     1] = x;
	iodata[order+1] = y;

	for(i = 1; i<order; i++) {
		iodata[     i+1] = iodata[     i] * x;
		iodata[order+i+1] = iodata[order+i] * y;
	};
	co = 2*order+1;
 	for (k = 2; k <= crossorder; k++) {	// xy x^2y ... y^2x ...
#ifdef REVXYZ
		for (iy = 1; iy < k; iy++)
			for (ix = 1; ix < k; ix++)
#else
		for (ix = 1; ix < k; ix++)
			for (iy = 1; iy < k; iy++)
#endif
				if ((ix + iy) == k) {
					iodata[co] = pow(x,ix) * pow(y,iy);
					co++;
				}
	}
}

double Mapper2D(double x, double y, int order, int crossorder, double* parms, int nparms)
{
   int i;
   double* iodata, result;

   iodata = (double*) malloc(sizeof(double) * nparms);
   getIOdata2D(iodata, x, y, order, crossorder);

   result = 0;
   for (i = 0; i < nparms; i++)
      result += iodata[i] * parms[i];

   free(iodata);
   return result;
}

double Mapxy_X(double x, double y, POLYFIT* pfit)
{
	return Mapper2D(x, y, pfit->order, pfit->crossorder, 
		              pfit->Xparms1, pfit->nxparms1);
}

double Mapxy_Xrc(double x, double y, POLYFIT* pfit)
{
	return Mapper2D(x, y, pfit->order, pfit->crossorder, 
		              pfit->Xparms2, pfit->nxparms2);
}

double Mapxy_Y(double x, double y, POLYFIT* pfit)
{
	return Mapper2D(x, y, pfit->order, pfit->crossorder, 
		              pfit->Yparms1, pfit->nyparms1);
}

double Mapxy_Yrc(double x, double y, POLYFIT* pfit)
{
	return Mapper2D(x, y, pfit->order, pfit->crossorder, 
		              pfit->Yparms2, pfit->nyparms2);
}

void getIOdata3D(double* iodata, double x, double y, double z, int order, int crossorder)
{
	int i, ix, iy, iz, k, co;

	iodata[       0]  = 1;
	iodata[       1]  = x;
	iodata[  order+1] = y;
	iodata[2*order+1] = z;

	for(i = 1; i<order; i++) {
		iodata[       i+1] = iodata[       i] * x;
		iodata[  order+i+1] = iodata[  order+i] * y;
		iodata[2*order+i+1] = iodata[2*order+i] * z;
	};
	co = 3*order+1;
 	for (k = 2; k <= crossorder; k++) {
		for (ix = 0; ix < k; ix++)
			for (iy = 0; iy < k; iy++)
				for (iz = 0; iz < k; iz++)
					if ((ix + iy + iz) == k) {
						// printf("co=%d, ix=%d iy=%d iz=%d\n", co, ix,iy,iz);
						if (ix == 0)
							iodata[co] = pow(y,iy) * pow(z,iz);
						else if (iy == 0)
							iodata[co] = pow(x,ix) * pow(z,iz);
						else if (iz == 0)
							iodata[co] = pow(x,ix) * pow(y,iy);
						else
							iodata[co] = pow(x,ix) * pow(y,iy) * pow(z,iz);
						co++;
					}
	}
}

double Mapper3D(double X, double Y, double Z, int order, int crossorder, double* parms, int nparms)
{
	int i;
	double* iodata, result;

	iodata = malloc(sizeof(double) * nparms);
	getIOdata3D(iodata, X, Y, Z, order, crossorder);

	result = 0;
	for (i = 0; i < nparms; i++)
	  result += iodata[i] * parms[i];

	free(iodata);
	return result;
}

double Map3D_x(double X, double Y, double Z, POLYFIT* pfit)
{
	return Mapper3D(X, Y, Z, pfit->order, pfit->crossorder, 
		              pfit->Xparms1, pfit->nxparms1);
}

double Map3D_y(double X, double Y, double Z, POLYFIT* pfit)
{
	return Mapper3D(X, Y, Z, pfit->order, pfit->crossorder, 
		              pfit->Yparms1, pfit->nyparms1);
}

// Using FitInfo instead of POLYFIT
double fitdata_Map3D_x(double X, double Y, double Z, FitInfo* pfit)
{
	return Mapper3D(X, Y, Z, pfit->order, pfit->crossorder, 
						pfit->parms1, pfit->nparms1);
}

double fitdata_Map3D_y(double X, double Y, double Z, FitInfo* pfit)
{
	return Mapper3D(X, Y, Z, pfit->order, pfit->crossorder, 
						pfit->parms2, pfit->nparms2);
}

double fitdata_Map2D_x(double X, double Y, FitInfo* pfit)
{
	return Mapper2D(X, Y, pfit->order, pfit->crossorder, 
						pfit->parms1, pfit->nparms1);
}

double fitdata_Map2D_y(double X, double Y, FitInfo* pfit)
{
	return Mapper2D(X, Y, pfit->order, pfit->crossorder,
						pfit->parms2, pfit->nparms2);
}

void parsepolycalibfile(char* fname, POLYFIT *pfit3dpix, POLYFIT *pfitpixepi)
{
	FILE*    fp; 
	char     buf[256];
	int      i, order, crossorder, nXparms1, nXparms2, nYparms1, nYparms2;
	int      ok, value, part;
	POLYFIT* pfit;

	if (fopen_s(&fp, fname, "r") != 0) {
		printf("Can't open %s\n", fname);
		return;
	}
	for (part=1; part<=2; part++) {
		pfit = (part == 1) ? pfit3dpix : pfitpixepi;
		i = 0;
		while (fgets(buf, 256, fp) && i < 6) {
			if (sscanf_s(buf, "%d", &value) == 1) switch (++i)
			{
				// case 1: ncameras  = value; break;
				case 1: order      = value; break;
				case 2: crossorder = value; break;
				case 3: nXparms1   = value; break;
				case 4: nXparms2   = value; break;
				case 5: nYparms1   = value; break;
				case 6: nYparms2   = value; break;
			}
		}
		if (i != 6) {
			printf("Bad header in %s\n", fname);
			return;
		}
		// prepare the POLYFIT structures
		// *ppolyfit = fit = malloc(sizeof(POLYFIT) * ncameras);

		pfit->order      = order;
		pfit->crossorder = crossorder;
		pfit->ncameras  = n_img;			// necessary ???
		pfit->nxparms1  = nXparms1;
		pfit->nxparms2  = nXparms2;
		pfit->nyparms1  = nYparms1;
		pfit->nyparms2  = nYparms2;
		pfit->Xparms1   = malloc(sizeof(double) * nXparms1);
		pfit->Xparms2   = malloc(sizeof(double) * nXparms2);
		pfit->Yparms1   = malloc(sizeof(double) * nYparms1);
		pfit->Yparms2   = malloc(sizeof(double) * nYparms2);

		// read all the parameters1 and parameters2 for X and Y, for this camera.
		ok = 1;
		for (i = 0; ok && i < nXparms1; i++)
			ok = fscanf_s(fp, "%le", &pfit->Xparms1[i]) == 1;

		for (i = 0; ok && i < nXparms2; i++)
			ok = fscanf_s(fp, "%le", &pfit->Xparms2[i]) == 1;

		for (i = 0; ok && i < nYparms1; i++)
			ok = fscanf_s(fp, "%le", &pfit->Yparms1[i]) == 1;

		for (i = 0; ok && i < nYparms2; i++)
			ok = fscanf_s(fp, "%le", &pfit->Yparms2[i]) == 1;
	}
	fclose(fp);
	if (!ok)
		printf("Not all the parameters are found in the file: %s\n", fname);
}


// --- functions used by the shaking algorithm ---

void save_3Dfitparameters2(char* pathname, char* title, FitInfo *pXfit, FitInfo *pYfit)
{
	save_3Dfitparameters3(pathname, "w", title, pXfit, pYfit);
	printf("Saved result to: %s\n", pathname);
}

void save_3Dfitparameters(char* basename, char* title, int ncam, FitInfo *Xfit, FitInfo *Yfit)
{
	int ic;
	char pathname[256];

	printf("save_3Dfitparameters() ---> OLD FUNCTION, DON'T USE !!!!\n");
	for (ic=0; ic<ncam; ic++) {
		build_indexedname(pathname, basename, ic+1);
		if (Yfit)
			save_3Dfitparameters2(pathname, title, &Xfit[ic], &Yfit[ic]);
		else
			save_3Dfitparameters2(pathname, title, &Xfit[ic], NULL);
	}
}

void saveparms_3D_to_pix(char* pathname, char* mode, FitInfo *pfit)
{
	char* title = "# fit parameters: 3D -> pixels";
	if (*mode != 'w' && *mode != 'a') {
		printf("saveparms_3D_to_pix(), BAD mode: %s\n", mode);
		return;
	}
	save_3Dfitparameters3(pathname, mode, title, pfit, NULL);
	printf("%s result to: %s\n", *mode == 'w' ? "Saved" : "Added", pathname);
}

void saveparms_2D_to_rays(char* pathname, char* mode, FitInfo *pXfit, FitInfo *pYfit)
{
	char* title = "# fit parameters: 2D -> epilpolar lines";
	if (*mode != 'w' && *mode != 'a') {
		printf("saveparms_2D_to_rays(), BAD mode: %s\n", mode);
		return;
	}
	save_3Dfitparameters3(pathname, mode, title, pXfit, pYfit);
	if (*mode == 'w')
		printf("Saved parameters in: %s\n", pathname);
	else
		printf("Added parameters to: %s\n", pathname);
}

void saveparms_2D_to_2D(char* pathname, char* mode, FitInfo *pfit)
{
	char* title = "# fit parameters: 2D -> 2D";
	if (*mode != 'w' && *mode != 'a') {
		printf("saveparms_2D_to_2D(), BAD mode: %s\n", mode);
		return;
	}
	save_3Dfitparameters3(pathname, mode, title, pfit, NULL);
	if (*mode == 'w')
		printf("Saved parameters in: %s\n", pathname);
	else
		printf("Added parameters to: %s\n", pathname);
}

void add_3Dfitparameters2(char* pathname, char* title, FitInfo *pXfit, FitInfo *pYfit)
{
	// Append the data to pathname = fopen mode "a"
	save_3Dfitparameters3(pathname, "a", title, pXfit, pYfit);
	printf("Added the parameter values to: %s\n", pathname);
}

void add_2D_to_epipolar(char* pathname, FitInfo *pXfit, FitInfo *pYfit)
{
	// Append the data to pathname = fopen mode "a"
	char* title = "# fit parameters: 2D -> epilpolar lines";
	save_3Dfitparameters3(pathname, "a", title, pXfit, pYfit);
	printf("Added the parameter values to: %s\n", pathname);
}

void save_3Dfitparameters3(char* pathname, char* mode, char* title, FitInfo *pXfit, FitInfo *pYfit)
{
	int i;
	FILE *fp;

	if (fopen_s(&fp, pathname, mode) != 0) {
		printf("Data NOT SAVED\nCan't open file: %s\n", pathname);
		return;
	}
	fprintf(fp, title);
	fprintf(fp, "\n\n");
	// fprintf(fp, "%4d ncameras\n",        ncam);  
	fprintf(fp, "%4d order\n",           pXfit->order);  
	fprintf(fp, "%4d crossorder\n",      pXfit->crossorder);  
	fprintf(fp, "%4d nXparameters1\n",   pXfit->nparms1);  
	fprintf(fp, "%4d nXparameters2\n",   pYfit ? pXfit->nparms2 : 0);
	fprintf(fp, "%4d nYparameters1\n",   pYfit ? pXfit->nparms1 : pXfit->nparms2); 
	fprintf(fp, "%4d nYparameters2\n\n", pYfit ? pXfit->nparms2 : 0);

	for (i = 0; i < pXfit->nparms1; i++)
		fprintf(fp, "%.14le\n", pXfit->parms1[i]);
	for (i = 0; i < pXfit->nparms2; i++)
		fprintf(fp, "%.14le\n", pXfit->parms2[i]);
	for (i = 0; pYfit && i < pYfit->nparms1; i++)
		fprintf(fp, "%.14le\n", pYfit->parms1[i]);
	for (i = 0; pYfit && i < pYfit->nparms2; i++)
		fprintf(fp, "%.14le\n", pYfit->parms2[i]);
	fprintf(fp,"\n");
	fclose(fp);
}

void destroy_POLYFIT_list(POLYFIT *fitlist)
{
	int ic;
	int ncam = fitlist[0].ncameras;
	for (ic=0; ic<ncam; ic++) {
		free(fitlist[ic].Xparms1);
		free(fitlist[ic].Xparms2);
		free(fitlist[ic].Yparms1);
		free(fitlist[ic].Yparms2);
	}
}

void ray_tracing_poly(double x, double y, /*int pnr,*/ int icam, Line3D* pline)
{
	double xrc = Mapxy_Xrc(x, y, &fitpixepi[icam]);
	double yrc = Mapxy_Yrc(x, y, &fitpixepi[icam]);

	pline->xoff = Mapxy_X(x, y, &fitpixepi[icam]);
	pline->yoff = Mapxy_Y(x, y, &fitpixepi[icam]);
	pline->zoff = 0;
	pline->xrc  = xrc;
	pline->yrc  = yrc;
}

void find_candidate_poly_epiline(coord_2d crd[], target pix[], int num, 
								 double xa, double ya, double xb, double yb, 
								 double eps, target pt1, candidate cand[], int* count, const char** argv)

// target pt1				a target of the first camera
// Line3D line				3D-line belonging to this target
// target pix[]				all targets found by the second camera
// int num,					number of targets in pix[]
// Line3D candlines[]		3D-lines belonging to the targets in pix[] 
// coord_2d crd[]
// double xa, ya, xb, yb	limits of the search area, in pixels 
// double eps				epipolar band width
// candidate cand[]			containe for the candidates found
// int* count)				number of the candidates found


/*  binarized search in a x-sorted coord-set, exploits shape information  */
{
	int		j, dummy, j0, dj; //, n,nx,ny,sumg;
	double	m, b, d, temp, qn, qnx, qny, qsumg, corr, xmin, xmax, ymin, ymax;
	target	pt2;
	double	tol_band_width, particle_size;
	int		dumbbell = 0, subcount = 0;

	// Beat Mai 2010 for dumbbell
	if (atoi(argv[1])==3)
	    dumbbell=1;

	if (dumbbell==0) {
	    /////here is new Beat version of April 2010
	    if (pt1.nx>pt1.ny) particle_size = pt1.nx;
	    else               particle_size = pt1.ny;
	    tol_band_width = eps*0.5*(pix_x+pix_y)*particle_size;
	}
	else
	    tol_band_width = eps;
    
	if (tol_band_width < 0.06)
		tol_band_width = 0.06;


	// The tol_band_width is in mm on the ccd. However, all coordintaes are
	// in pixels, so the tol_band_width is converted to mm.
	tol_band_width /= 0.5*(pix_x+pix_y);

	/* define sensor format for search interrupt */
	//xmin = (-1) * pix_x * imx/2;  xmax = pix_x * imx/2;
	//ymin = (-1) * pix_y * imy/2;  ymax = pix_y * imy/2;
	/* define sensor format for search interrupt in pixels */
	xmin = 0; xmax = imx;
	ymin = 0; ymax = imy;

	for (j=0; j<4; j++) {
		cand[j].pnr = -999;
		cand[j].tol = -999;
		cand[j].corr = -999;
	}

	/* line equation: y = m*x + b */
	if (xa == xb)  xa += 1e-10;
	m = (yb-ya)/(xb-xa);  b = ya - m*xa;

	if (xa > xb) { temp = xa;  xa = xb;  xb = temp; }
	if (ya > yb) { temp = ya;  ya = yb;  yb = temp; }

	// n = pt1.n;	nx = pt1.nx;  ny = pt1.ny; sumg = pt1.sumg;
	if (xb>xmin && xa<xmax && yb>ymin && ya<ymax)	/* in sensor area */
	{
		 /* now inflate the sensor area by the tol_band_width */
		xa -= tol_band_width;;
		ya -= tol_band_width;;
		xb += tol_band_width;;
		yb += tol_band_width;;

		/* binarized search for start point of candidate search */
		for (j0=num/2, dj=num/4; dj>1; dj/=2) {
			if (crd[j0].x < xa) j0 += dj;
			else                j0 -= dj;
		}
		j0 -= 12;  if (j0 < 0) j0 = 0;   /* due to trunc */

		for (j=j0, *count=0; j<num; j++) {				/* candidate search */
			if (crd[j].x > xb) return;					/* finish search */

			if (crd[j].x>xa && crd[j].x<xb && crd[j].y>ya && crd[j].y<yb)
			{
				d = fabs ((crd[j].y - m*crd[j].x - b) / sqrt(m*m+1));

				// Beat: modified in April 2010 to allow for better treatment of 
				// different sized traced particles, in particular colloids and tracers
				// old : if (d < eps) {
				// new : // if (nx>ny) particle_size=nx;
				//		 // else	   particle_size=ny;
				//		 if (d < tol_band_width) {

				if (d < tol_band_width) {
					pt2 = pix[crd[j].pnr];

					if (pt1.n < pt2.n)       qn    = (double) pt1.n/pt2.n;
					else		       	     qn    = (double) pt2.n/pt1.n;
					if (pt1.nx < pt2.nx)     qnx   = (double) pt1.nx/pt2.nx;
					else		       	     qnx   = (double) pt2.nx/pt1.nx;
					if (pt1.ny < pt2.ny)     qny   = (double) pt1.ny/pt2.ny;
					else		       	     qny   = (double) pt2.ny/pt1.ny;
					if (pt1.sumg < pt2.sumg) qsumg = (double) pt1.sumg/pt2.sumg;
					else                     qsumg = (double) pt2.sumg/pt1.sumg;

					// empirical correlation coefficient from shape and 
					// brightness parameters 
					corr = (4*qsumg + 2*qn + qnx + qny);
              			
					/* create a tendency to prefer those matches with brighter targets */
					corr *= (double)(pt1.sumg + pt2.sumg);

					if (qn>=cn && qnx>=cnx && qny>=cny && qsumg>csumg) {
						if (*count < maxcand) {
							cand[*count].pnr   = j;
							cand[*count].tol   = d;
							cand[*count].corr  = corr;
							(*count)++;
						}
						else {
							dummy = (int)maxcand;
							printf("in find_candidate_plus: count > maxcand\n");
						}
					}
				}
			}
		}
	}
	else  *count = -1;		   /* out of sensor area */
}

void mid_point_poly(Line3D line1, Line3D line2, 
                    double* pdist, double* pX, double *pY, double *pZ)
{
    Waist w;
    Line3D_mindist(line1, line2, &w);
    *pdist = w.dist;
    *pX = w.x;
    *pY = w.y;
    *pZ = w.z;
}

double mindist(Line3D *line1, Line3D *line2)
{
    // line 1 : r1 + a1.s1
    // line 2 : r2 + a2.s2
    //
    // n = cross(s1,s2);
    // n = n / abs(n);
    // D = abs(dot(n,r1) - dot(n,r2));

    double rx,ry, ax,ay, sx,sy, bx,by, nx,ny,nz, N, d;

	rx = line1->xoff;
	ry = line1->yoff;
	// rz = 0;
	ax = line1->xrc;
	ay = line1->yrc;
	// az = 1;;

	sx = line2->xoff;
	sy = line2->yoff;
	// sz = 0;
	bx = line2->xrc;
	by = line2->yrc;
	// bz = 1;

    // n = cross(a,b)
	nx = ay - by;           // nx = ay * bz - az * by;
	ny = bx - ax;           // ny = az * bx - ax * bz;
	nz = ax*by - ay*bx;     // nz = ax * by - ay * bx;

    N  = sqrt(nx*nx + ny*ny + nz*nz);
    d  = nx*(rx - sx) + ny*(ry - sy);
    return fabs(d)/N;
}

void Line3D_mindist(Line3D line1, Line3D line2, Waist *pwaist)
{
	double r1x,r1y,a1x,a1y,r2x,r2y,a2x,a2y;
	double h1,h2,h3,h4,h5,h7,h8,h9,h10,h11,h12,h13;
	double p,q, Px,Py, Qx,Qy;

	// The target lines calculated by polynomials are going through 
	// the point (xoff, yoff, 0) with direction vector (xrc, yrc, 1)
	r1x = line1.xoff;
	r1y = line1.yoff;
	// r1z = 0;
	a1x = line1.xrc;
	a1y = line1.yrc;
	// a1z = 1;

	r2x = line2.xoff;
	r2y = line2.yoff;
	// r2z = 0;
	a2x = line2.xrc;
	a2y = line2.yrc;
	// a2z = 1;

	h1 = a1y - a2y;					// h1 = a1y * a2z - a1z * a2y;
	h2 = a2x - a1x;					// h2 = a1z * a2x - a1x * a2z;
	h3 = a1x * a2y - a1y * a2x;
 
	h4 = r1x - r2x;
	h5 = r1y - r2y;
	// h6 = r1z - r2z;

	h7  = h3*a2y - h2;			// h7  = h3*a2y - h2*a2z;
	h8  = h3*h5;						// h8  = h3*h5 - h2*h6;
	h12 = (h3*a1y - h2)/h7;	//	h12 = (h3*a1y - h2*a1z)/h7;
	h13 = h8/h7;

	h9  = h2*h4  - h1*h5;
	h10 = h2*a1x - h1*a1y;
	h11 = h2*a2x - h1*a2y;

	p = (h9 - h11*h13)/(h11*h12 - h10);
	q = h13 + p*h12;
	Px = r1x + p*a1x; Py = r1y + p*a1y; // Pz = r1z + p*a1z; 
	Qx = r2x + q*a2x; Qy = r2y + q*a2y; // Qz = r2z + q*a2z; 

	pwaist->dist = sqrt((Px-Qx)*(Px-Qx) + (Py-Qy)*(Py-Qy) + (p-q)*(p-q));
	pwaist->x    = (Px + Qx)/2;
	pwaist->y    = (Py + Qy)/2;
	pwaist->z    = (p + q)/2;
}

// original call: volumedimension(), in multimed.c
int volumedimension_poly(double* xmax, double* xmin, 
						 double* ymax, double* ymin, double* zmax, double* zmin)
{
	int	ip, iz, icam;
	double X,Y,Z, Zmin, Rmax=0,Zmax;
	critparameters cpar;

	/* find extrema of imaged object volume */
	if (! parse_polycriteria_par2("parameters/poly_criteria.par", &cpar))
		return FALSE;
	*xmin = *xmax = (cpar.xlay[0] + cpar.xlay[1])/2;
	*ymin = *ymax = (cpar.ylay[0] + cpar.ylay[1])/2;
	*zmin = Zmin = cpar.zminlay[0];
	*zmax = Zmax = cpar.zmaxlay[1];

	for (icam=0; icam<n_img; icam++) {
		coord_2d P[4] = {0, 0.0, 0.0, 
						 1, 0.0, (double)imy,
						 2, (double)imx, 0.0,
						 3, (double)imx, (double)imy};
		// intersect with image vertices rays
		// corner points 0(0,0), 1(imx,0), 2(0,imy) and 3(imx,imy)
		for (ip = 0; ip < 4; ip++) {
			Line3D line;
			// ray_tracing_poly(P[ip].x, P[ip].y, ip, icam, &line);
			ray_tracing_poly(P[ip].x, P[ip].y, icam, &line);

			for (iz = 0; iz < 2; iz++) {
				Z = iz == 0 ? Zmin : Zmax;
				X = line.xoff + (Z-line.zoff)*line.xrc;
				Y = line.yoff + (Z-line.zoff)*line.yrc;

				if (X < *xmin) *xmin=X;
				if (X > *xmax) *xmax=X;
				if (Y < *ymin) *ymin=Y;
				if (Y > *ymax) *ymax=Y;
			}
		}
	}
	return TRUE;
}

int InCubeEsp(Waist* w1, Waist* w2, double eps)
{
	return fabs(w1->x - w2->x) <= eps && 
				 fabs(w1->y - w2->y) <= eps &&
				 fabs(w1->z - w2->z) <= eps ;
}

int inVolume(double x, double y, double z)
{
    return (x >= X_lay[0] && x <= X_lay[1] &&
            y >= Y_lay[0] && y <= Y_lay[1] &&
            z >= Zmin_lay[0] && z <= Zmax_lay[1]);
}

int epi_mm_poly(double xp1, double yp1, int ic1, int ic2, 
				double* pxmin, double* pymin, double* pxmax, double* pymax)
/*
	double xp1, yp1;                        input coord 
	int    ic1, ic2;                        camera indices
	double *pxmin, *pymin, *pxmax, *pymax;  output search window
*/
{
	double x1,y1,z1,x2,y2,z2,t,xa,ya,za,xb,yb,zb;
	double Xo,Xr,Yo,Yr,Zo, xmin, xmax, ymin, ymax, zmin, zmax;
	Line3D line;  
	int    found1, found2;

    xmin = X_lay[0];    xmax = X_lay[1];
    ymin = Y_lay[0];    ymax = Y_lay[1];
    zmin = Zmin_lay[0]; zmax = Zmax_lay[1];

	ray_tracing_poly(xp1, yp1, ic1, &line);
	Xo = line.xoff;		// intersection of the epipolar line with the plane z = 0
	Yo = line.yoff;
	Zo = 0;
	Xr = line.xrc;
	Yr = line.yrc;

	// determine the intersection points of the epipolar line with the bottom and
	// top of the illuminated volume.
    z1 = zmin;
    x1 = Xo + (z1-Zo)*Xr;
    y1 = Yo + (z1-Zo)*Yr;
    z2 = zmax;
    x2 = Xo + (z2-Zo)*Xr;
    y2 = Yo + (z2-Zo)*Yr;
	if (x1 > x2) {		// choose the points so that x1<x2
        t = x2; x2 = x1; x1 = t;
        t = y2; y2 = y1; y1 = t;
        t = z2; z2 = z1; z1 = t;
	}
	
	// First intersection point
    found1 = 0;
	if (inVolume(x1,y1,z1)) {
        xa = x1;
        ya = y1;
        za = z1;
        found1 = 1;
	}
	if (!found1 && x1<xmin) {
        xa = xmin;
        za = Zo + (xa - Xo)/Xr;
        ya = Yo + (za - Zo)*Yr;
        if (inVolume(xa,ya,za))
            found1 = 1;
	}
	if (!found1 && y1<ymin) {
        ya = ymin;
        za = Zo + (ya - Yo)/Yr;
        xa = Xo + (za - Zo)*Xr;
        if (inVolume(xa,ya,za))
            found1 = 1;
	}
	if (!found1 && y1>ymax) {
        ya = ymax;
        za = Zo + (ya - Yo)/Yr;
        xa = Xo + (za - Zo)*Xr;
        if (inVolume(xa,ya,za))
            found1 = 1;
	}

	// second intersection point
    found2 = 0;
	if (inVolume(x2,y2,z2)) {
        xb = x2;
        yb = y2;
        zb = z2;
        found2 = 1;
	}
	if (!found2 && x2>xmax) {
        xb = xmax;
        zb = Zo + (xb - Xo)/Xr;
        yb = Yo + (zb - Zo)*Yr;
        if (inVolume(xb,yb,zb))
            found2 = 1;
	}
	if (!found2 && y2<ymin) {
        yb = ymin;
        zb = Zo + (yb - Yo)/Yr;
        xb = Xo + (zb - Zo)*Xr;
        if (inVolume(xb,yb,zb))
            found2 = 1;
	}
	if (!found2 && y2>ymax) {
        yb = ymax;
        zb = Zo + (yb - Yo)/Yr;
        xb = Xo + (zb - Zo)*Xr;
        if (inVolume(xb,yb,zb))
            found2 = 1;
	}
	if (found1 && found2) {
		*pxmin = Map3D_x(xa,ya,za, &fit3dpix[ic2]);
		*pymin = Map3D_y(xa,ya,za, &fit3dpix[ic2]);
		*pxmax = Map3D_x(xb,yb,zb, &fit3dpix[ic2]);
		*pymax = Map3D_y(xb,yb,zb, &fit3dpix[ic2]);
		return 1;	// ie, intersects the illuminated volume
	}
	return (0);
}

void det_lsq_3d_poly(double xp[], double yp[], int np, double *pX, double *pY, double *pZ)
{
    int     i,count_inner=0,n,m,flag[4];
    double  d_inner=0, dist, X_pos[6], Y_pos[6], Z_pos[6], XX,YY,ZZ, si0,sqX,sqY,sqZ;
    Line3D  line[4];

    for (i = 0; i < 4; i++) {
        flag[i] = 0;
        if(xp[i] > -999){
            flag[i] = 1;
            ray_tracing_poly(xp[i], yp[i], i, &line[i]);
        }       
    }
    count_inner = 0;
    for (n=0; n<n_img-1; n++){
        for (m=n+1; m<n_img; m++){
            if (flag[n]==1 && flag[m]==1) {
                mid_point_poly(line[n], line[m], &dist, &XX, &YY, &ZZ);
                d_inner += dist;
                X_pos[count_inner]=XX;
                Y_pos[count_inner]=YY;
                Z_pos[count_inner]=ZZ;
                count_inner++;
            }
        }
    }
    d_inner/=(double)count_inner;       
    XX=0.;YY=0.;ZZ=0.;
    for(i=0;i<count_inner;i++){
       XX+=X_pos[i]; 
       YY+=Y_pos[i];
       ZZ+=Z_pos[i];
    }
    XX/=(double)count_inner;YY/=(double)count_inner;ZZ/=(double)count_inner;
    //end of new det_lsq
    *pX = XX;
    *pY = YY;
    *pZ = ZZ;

    //statistics
    si0=0.;sqX=0.;sqY=0.;sqZ=0.;
    for(i=0;i<count_inner;i++){
       si0+=pow(X_pos[i]-XX,2.)+pow(Y_pos[i]-YY,2.)+pow(Z_pos[i]-ZZ,2.);
       sqX+=pow(X_pos[i]-XX,2.); 
       sqY+=pow(Y_pos[i]-YY,2.);
       sqZ+=pow(Z_pos[i]-ZZ,2.);           
    }
    si0/=(double)count_inner;sqX/=(double)count_inner;sqY/=(double)count_inner;sqZ/=(double)count_inner;
    
    mean_sigma0 += pow(si0,0.5);
    rmsX += pow(sqX,0.5);
    rmsY += pow(sqY,0.5);
    rmsZ += pow(sqZ,0.5);
    //end of statistics
}


void boundingrect(Rect* prect, coord_2d* points, int npoints, int flag)
{
    int i;
    double x1,x2,y1,y2;

    if (flag == 0) {
        x1 = x2 = points[0].x;
        y1 = y2 = points[0].y;
    }
    else {
        x1 = prect->left;
        x2 = prect->right;
        y1 = prect->bottom;
        y2 = prect->top;
    }

    for (i=1; i<npoints; i++) {
        if (x1 > points[i].x) x1 = points[i].x;
        if (x2 < points[i].x) x2 = points[i].x;
        if (y1 > points[i].y) y1 = points[i].y;
        if (y2 < points[i].y) y2 = points[i].y;
    }
    prect->left   = x1;
    prect->right  = x2;
    prect->bottom = y1;
    prect->top    = y2;
}

void CreateZplaneFitdata(coord_2d *pixcrd, coord_2d *planecrd, FitInfo *pfit, double z,
					  double x1, double x2, int nx, double y1, double y2, int ny)
{
	int i;
	double x, y, dx, dy;

	dx = (x2-x1)/(nx-1);
	dy = (y2-y1)/(ny-1);

	i = 0;
	for (x = x1; x <= x2; x += dx)
		for (y = y1; y <= y2; y += dy) {
			pixcrd[i].x   = fitdata_Map3D_x(x, y, z, pfit);
			pixcrd[i].y   = fitdata_Map3D_y(x, y, z, pfit);
			planecrd[i].x = x;
			planecrd[i].y = y;
			i++;
		}
}
