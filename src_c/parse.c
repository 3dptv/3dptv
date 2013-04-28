// parse,c
//
// functions for reading and saving data from the various text files
// ad holten, 04-2013

#include "ptv.h"


int read_strings(char* fname, int nparms, char val[][256])
{
	int  i, ne, c;
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Can't open file: %s\n", fname);
		return 0;
	}
    for (i=0, ne=0; i<nparms; i++) {
        val[i][0] = 0;    // reset the string val[i]
        if (fscanf (fp, "%255[^\n]\n", val[i]) == 1) ne++;
		
		// if \n has not been found, skip the rest of the line
		if (strlen(val[i]) >= 255)
			while ((c=getc(fp))>0 && c!='\n');
    }
    fclose(fp);
    if (ne < nparms) {
        printf("%s : missing %d parameters\n", fname, nparms-ne);
		return 0;
	}
	return 1;
}

int read_allstrings(char* fname, char val[][256], int maxcnt)
{
	int  i,c;
	FILE *fp = fopen(fname, "r");
	if (!fp) {
		printf("Can't open file: %s\n", fname);
		return 0;
	}
    for (i=0; i<maxcnt; i++)    // clears the strings buffer
		val[i][0] = 0;

	for (i=0; fscanf (fp, "%255[^\n]\n", val[i]) == 1; i++) {
		// if \n has not been found, skip the rest of the line
		if (strlen(val[i]) >= 255)
			while ((c=getc(fp))>0 && c!='\n');
	}
	fclose(fp);
	return i;
}

void get_processingmethod(int* pMapping, int* pUsingplanes)
{
    // Determines:
    // 1: the method of particle mapping, ETHZ or polynomials
    // 2: what kind of calibration method is used, a 3D-body or planes 

    char  s1[256] = "ETHZ", s2[256] = "body";
    FILE* fp = fopen("parameters/method.par","r");
    if (fp) {
        fscanf(fp, "%s%s", s1, s2);
        fclose(fp);
    }
    *pMapping     = toupper(s1[0]) == 'P' ? POLYN : ETHZ;  // ETHZ is default
    *pUsingplanes = tolower(s2[0]) == 'b' ? 0 : 1;         // 3D-body is default
}

BOOL read_3Dcoordinates(char *fixp_name, coord_3d **plist, int *pcount)
{
	FILE *fp;
	int i,k,n;
	coord_3d* fixp;

	fp = fopen_rp(fixp_name);
	if (!fp) return FALSE;
	
	// count number of points
	n = 0;
	while (fscanf(fp, "%d %*lf %*lf %*lf", &i) == 1)  n++;
	rewind(fp);

	// allocate and fill the 3D-coord. buffer
	fixp = (coord_3d*) malloc(n*sizeof(coord_3d));
	k = 0;
	while (fscanf (fp, "%d %lf %lf %lf", &fixp[k].pnr,
			&fixp[k].x, &fixp[k].y, &fixp[k].z) != EOF)  k++;

	fclose (fp);
	*plist  = fixp;
	*pcount = k;
	return TRUE;
}

BOOL read_2Dcoordinates(char *pixp_name, coord_2d **plist, int *pcount)
{
	int i,n;
	coord_2d *pt;

	FILE *fp = fopen_rp(pixp_name);
	if (!fp) return FALSE;
	
	// count number of points
	n = 0;
	while (fscanf(fp, "%d %*lf %*lf", &i) == 1)  n++;
	rewind(fp);

	// allocate and fill the 2D-coord. buffer
	pt = (coord_2d*) malloc(n*sizeof(coord_2d));
	i = 0;
	while (i<n && fscanf (fp, "%d%lf%lf", &pt[i].pnr, &pt[i].x, &pt[i].y) == 3)
		i++;

	fclose (fp);
	*plist  = pt;
	*pcount = n;
	return TRUE;
}

BOOL read_target_crds(char *pathname, coord_3d **pfixp, int *pcount)
{
    // Read the 3D-coordinates from a target file pathname.
    int i, n;
    coord_3d *pt;

	FILE *fp = fopen_rp(pathname);
    if (!fp) return FALSE;

	// count number of points
	n = 0;
	while (fscanf(fp, "%d %*lf %*lf %*lf", &i) == 1)  n++;
	rewind(fp);

	// allocate and fill the 3D-coord. buffer
	pt = (coord_3d*) malloc(n*sizeof(coord_3d));
	i = 0;
	while (i<n && fscanf (fp, "%d%lf%lf%lf", &pt[i].pnr, &pt[i].x, &pt[i].y, &pt[i].z) == 4)
		i++;
    
	fclose (fp);
	*pfixp  = pt;
    *pcount = n;
    return TRUE;
}


// line in distribution file:
//   col:1 2   3   4   5   6   7   8   9   10 11 12 13
//   index xp1 yp1 xp2 yp2 xp3 yp3 xp4 yp4 X  Y  Z  r (of d?)

void ReadDataFromDistribFile(char* pathname, int ic, int* pnelem, target* pixels, coord_3d* fixed)
{
    FILE* fp = fopen(pathname, "r");
    char buf[1024];
    int  i, cnt, nvalid;
    double *val, *p;

    // count lines
    for (cnt = 0; fgets(buf,1024,fp) != NULL; cnt++);
    rewind(fp);

    printf("cnt = %d\n", cnt);
    // read all the values
    p = val = (double*) malloc(cnt*13*sizeof(double));
    while (fscanf(fp, "%lf", p++) == 1);

    // count the valid points
    nvalid = 0;
    for (i=0; i<cnt; i++) {
        if (val[i*13 + 1 + 2*ic] >= 0)  // xp
            nvalid++;
    }
    printf("nvalid = %d\n", nvalid);
    if (nvalid > 10000) {
        printf("Too many points, aborted\n");
        return;
    }

    nvalid = 0;
    for (i=0; i<cnt; i++) {
        pixels[nvalid].x = val[i*13 + 1 + 2*ic];
        if (pixels[nvalid].x >= 0) {
            pixels[nvalid].y   = val[i*13 + 2 + 2*ic];
            pixels[nvalid].pnr = i;

            fixed[nvalid].x   = val[i*13 + 9];
            fixed[nvalid].y   = val[i*13 + 10];
            fixed[nvalid].z   = val[i*13 + 11];
            fixed[nvalid].pnr = i;
            nvalid++;
        }
    }
    *pnelem = nvalid;
    fclose(fp);
    free(val);
}


BOOL parse_ptv_par(char* fname)
{
    /* read from main parameter file */
    int ne, i;
	if (! (fpp = fopen_rp(fname)) ) return FALSE;
	
    ne = fscanf (fpp, "%d\n", &n_img);		// remind that fscanf could return -1, 0 and 1
    for (i=0; i<4; i++) {
        ne += fscanf (fpp, "%[^\n]\n", img_name[i]);
        ne += fscanf (fpp, "%[^\n]\n", img_cal[i]);
    }
    ne += fscanf(fpp, "%d",  &hp_flag);
    ne += fscanf(fpp, "%d",  &allCam_flag);
    ne += fscanf(fpp, "%d",  &tiff_flag);
    ne += fscanf(fpp, "%d",  &imx);
    ne += fscanf(fpp, "%d",  &imy);	
    ne += fscanf(fpp, "%lf", &pix_x);
    ne += fscanf(fpp, "%lf", &pix_y);
    ne += fscanf(fpp, "%d",  &chfield);
    ne += fscanf(fpp, "%lf", &mmp.n1);
    ne += fscanf(fpp, "%lf", &mmp.n2[0]);
    ne += fscanf(fpp, "%lf", &mmp.n3);
    ne += fscanf(fpp, "%lf", &mmp.d[0]);
    fclose(fpp);

    if (ne < 21)
        printf("%s : missing parameters, cnt < 21\n", fname);
    return ne == 21;
}

BOOL parse_polyptv_par(char* fname)
{
    /* read the main parameters needed for the polynomial mapping method */
    int ne, i;
	if (! (fpp = fopen_rp(fname)) ) return FALSE;

    ne = fscanf (fpp, "%d\n", &n_img);
    for (i=0; i<4; i++) {
        ne += fscanf (fpp, "%[^\n]\n", img_name[i]);
        ne += fscanf (fpp, "%[^\n]\n", img_cal[i]);
    }
    ne += fscanf(fpp, "%d",  &hp_flag);
    ne += fscanf(fpp, "%d",  &allCam_flag);
    ne += fscanf(fpp, "%d",  &tiff_flag);		// not needed
    ne += fscanf(fpp, "%d",  &imx);
    ne += fscanf(fpp, "%d",  &imy);
    ne += fscanf(fpp, "%lf", &pix_x);
    ne += fscanf(fpp, "%lf", &pix_y);
    fclose(fpp);

    if (ne < 16) {
        printf("%s : missing parameters, cnt < %d\n", fname, 16);
		return FALSE;
	}
    return TRUE;
}


BOOL parse_criteria_par(char* fname)
{
    /* read illuminated layer data */
    int ne;
    if (! (fpp = fopen_rp(fname)) ) return FALSE;

	ne  = fscanf (fpp, "%lf", &X_lay[0]);
	ne += fscanf (fpp, "%lf", &Zmin_lay[0]);
	ne += fscanf (fpp, "%lf", &Zmax_lay[0]);
	ne += fscanf (fpp, "%lf", &X_lay[1]);
	ne += fscanf (fpp, "%lf", &Zmin_lay[1]);
	ne += fscanf (fpp, "%lf", &Zmax_lay[1]);
	ne += fscanf (fpp, "%lf", &cnx);
	ne += fscanf (fpp, "%lf", &cny);
	ne += fscanf (fpp, "%lf", &cn);
	ne += fscanf (fpp, "%lf", &csumg);
	ne += fscanf (fpp, "%lf", &corrmin);
	ne += fscanf (fpp, "%lf", &eps0);

	//	if (map_method != ETHZ) {
	//		ne += fscanf (fpp, "%lf", &Y_lay[0]);
	//		ne += fscanf (fpp, "%lf", &Y_lay[1]);
	//		ne -= 2;	// so, if all went well, ne will be 12
	//	}
	fclose (fpp);

    if (ne < 12)
        printf("%s : missing parameters, cnt < 12\n", fname);

	mmp.nlay = 1;
    return ne == 12;
}

BOOL parse_polycriteria_par(char* fname)
{
    /* read illuminated volume data */
    int ne;
	if (! (fpp = fopen_rp(fname)) ) return FALSE;

	ne  = fscanf (fpp, "%lf", &X_lay[0]);	// remind that fscanf could return -1, 0 and 1
	ne += fscanf (fpp, "%lf", &Y_lay[0]);
	ne += fscanf (fpp, "%lf", &Zmin_lay[0]);
	ne += fscanf (fpp, "%lf", &X_lay[1]);
	ne += fscanf (fpp, "%lf", &Y_lay[1]);
	ne += fscanf (fpp, "%lf", &Zmax_lay[0]);
	Zmin_lay[1] = Zmin_lay[0];
	Zmax_lay[1] = Zmax_lay[0];

	ne += fscanf (fpp, "%lf", &cnx);
	ne += fscanf (fpp, "%lf", &cny);
	ne += fscanf (fpp, "%lf", &cn);
	ne += fscanf (fpp, "%lf", &csumg);
	ne += fscanf (fpp, "%lf", &corrmin);
	ne += fscanf (fpp, "%lf", &eps0);

	fclose (fpp);
    if (ne < 12)
        printf("%s : missing parameters, cnt < 12\n", fname);

	mmp.nlay = 1;
    return ne == 12;
}

BOOL parse_polycriteria_par2(char* fname, critparameters *cpar)
{
    /* read illuminated volume data */
    int ne;
	if (! (fpp = fopen_rp(fname)) ) return FALSE;

	//ne  = fscanf (fpp, "%lf", &cpar->xlay[0]);	// remind that fscanf could return -1, 0 and 1
	//ne += fscanf (fpp, "%lf", &cpar->ylay[0]);
	//ne += fscanf (fpp, "%lf", &cpar->zminlay[0]);
	//ne += fscanf (fpp, "%lf", &cpar->xlay[1]);
	//ne += fscanf (fpp, "%lf", &cpar->ylay[1]);
	//ne += fscanf (fpp, "%lf", &cpar->zmaxlay[0]);
	//cpar->zminlay[1] = cpar->zminlay[0];
	//cpar->zmaxlay[1] = cpar->zmaxlay[0];

	ne  = fscanf (fpp, "%lf", &cpar->xmin);	// remind that fscanf could return -1, 0 and 1
	ne += fscanf (fpp, "%lf", &cpar->ymin);
	ne += fscanf (fpp, "%lf", &cpar->zmin);
	ne += fscanf (fpp, "%lf", &cpar->xmax);
	ne += fscanf (fpp, "%lf", &cpar->ymax);
	ne += fscanf (fpp, "%lf", &cpar->zmax);


	ne += fscanf (fpp, "%lf", &cpar->cnx);
	ne += fscanf (fpp, "%lf", &cpar->cny);
	ne += fscanf (fpp, "%lf", &cpar->cn);
	ne += fscanf (fpp, "%lf", &cpar->csumg);
	ne += fscanf (fpp, "%lf", &cpar->corrmin);
	ne += fscanf (fpp, "%lf", &cpar->eps0);

	fclose (fpp);
    if (ne < 12)
        printf("%s : missing parameters, cnt < 12\n", fname);
	
	// traditional named parameters
	cpar->zminlay[1] = cpar->zminlay[0] = cpar->zmin;
	cpar->zmaxlay[1] = cpar->zmaxlay[0] = cpar->zmax;
	cpar->xlay[0] = cpar->xmin;
	cpar->xlay[1] = cpar->xmax;
	cpar->ylay[0] = cpar->ymin;
	cpar->ylay[1] = cpar->ymax;

	mmp.nlay = 1;
    return ne == 12;
}

BOOL parse_sequence_par(char* fname, char seq_name[][128], 
					    int *pseq_first, int *pseq_last)
{
    /* read sequence parameters */
    int i, ne;

    if (! (fpp = fopen_rp(fname)) ) return FALSE;

    ne = 0;
    for (i=0; i<4; i++)
        ne += fscanf(fpp, "%s\n", seq_name[i]);
    ne += fscanf (fpp,"%d\n", pseq_first);
    ne += fscanf (fpp,"%d\n", pseq_last);
    fclose (fpp);

    if (ne < 6)
        printf("parameters/sequence.par : missing parameters, cnt < 6\n");
    return ne == 6;
}

BOOL parse_body_calori_par(char* fname, char imgname[][256], char* fixpname, int deg3d2d[2], int deg2depi[2])
{
    int ncam, nplanes, i, n;
	char val[20][256];

	n = read_allstrings(fname, val, 20);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);

	if (nplanes<1 || n<2 + nplanes*(1+ncam) + 4) {
		printf("Not enough parameters in %s\n", fname);
		return FALSE;
	}
	n = 2;
	for (i=0; i<4; i++)
		strcpy(imgname[i], val[n++]);
	strcpy(fixpname, val[n++]);

	deg3d2d[0]  = atoi(val[n++]);
	deg3d2d[1]  = atoi(val[n++]);
	deg2depi[0] = atoi(val[n++]);
	deg2depi[1] = atoi(val[n++]);
	return TRUE;
}


BOOL parse_body_manori_par(char* fname, int list[][6], int* pncam, int* pnori)
{
    int ncam, nori, i, ic, n, nplanes;
	char val[100][256];

	n = read_allstrings(fname, val, 100);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);
	nori    = atoi(val[2]);

	if (ncam<1 || nori<2 || n < 3+ncam*nori) {
		printf("Not enough parameters in %s\n", fname);
		return FALSE;
	}
	n = 3;
	for (ic=0; ic<ncam; ic++) 
		for (i=0; i<nori; i++)
			list[ic][i] = atoi(val[n++]);
	*pncam = ncam;
	*pnori = nori;
	return TRUE;
}

BOOL parse_mult_calori_par(char* fname, int iplane, 
		int *pncam, int *pnplanes, char imgname[][256], 
		char* fixpname, int deg3d2d[2], int deg2depi[2])
{
    int ncam, nplanes, i, n;
	char val[100][256];

	n = read_allstrings(fname, val, 100);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);

	if (ncam<1 || nplanes<1 || n < 6+nplanes*(ncam+1)) {
		printf("Error in %s, not ebnough parameters\n", fname);
		return FALSE;
	}
	// iplane count from 0
	if (iplane<0 || iplane>=nplanes) {
		printf("plane index (=%d) out of range\n", iplane);
		return FALSE;
	}
	n = 2 + iplane*(ncam+1);
	strcpy(fixpname, val[n++]);
	for (i=0; i<4; i++)
		strcpy(imgname[i], val[n++]);

	n = 2 + nplanes*(1+ncam);
	deg3d2d[0]  = atoi(val[n++]);
	deg3d2d[1]  = atoi(val[n++]);
	deg2depi[0] = atoi(val[n++]);
	deg2depi[1] = atoi(val[n++]);

	*pncam    = ncam;
	*pnplanes = nplanes;
	return TRUE;
}

BOOL parse_multi_manori_par(char* fname, int iplane,  int list[][MaxOri], int *pncam, int *pnori)
{
    int ncam, nplanes, i, ic, n;
	char val[100][256];

	n = read_allstrings(fname, val, 100);
	ncam    = atoi(val[0]);
	nplanes = atoi(val[1]);
	nori    = atoi(val[2]);

	// rowlen : number of points per camera
	if (ncam<1 || nplanes<1 || n < 2+nplanes*ncam*nori) {
		printf("Error in %s, not ebnough parameters\n", fname);
		return FALSE;
	}
	// iplane count from 0
	if (iplane<0 || iplane>=nplanes) {
		printf("plane index (=%d) out of range\n", iplane);
		return FALSE;
	}
	n = 3 + iplane*ncam*nori;
	for (ic=0; ic<ncam; ic++) 
		for (i=0; i<4; i++)
			list[ic][i] = atoi(val[n++]);
	*pncam = ncam;
	*pnori = nori;
	return TRUE;
}


double get_sortradius(char *pathname, double defaultR)
{
	// Reads the radius used in the sorting functions.
	// If the file is not found, the default value is returned.
	double R = 0;
	FILE *fp = fopen("parameters/sortgrid.par", "r");
    if (fp) {
        fscanf (fp, "%lf", &R);
        fclose(fpp);
    }
    if (R == 0)
		R = defaultR;
	return R;
}

BOOL parse_track_par(char* fname, trackparameters *tpar)
{
	char val[20][256];

	if (! read_strings(fname, 9, val))	// function displays a message on error
		return FALSE;
	tpar->dvxmin = atof(val[0]);
	tpar->dvxmax = atof(val[1]);
	tpar->dvymin = atof(val[2]);
	tpar->dvymax = atof(val[3]);
	tpar->dvzmin = atof(val[4]);
	tpar->dvzmax = atof(val[5]);
	tpar->dangle = atof(val[6]);
	tpar->dacc   = atof(val[7]);
	tpar->add    = atoi(val[8]);
	return TRUE;
}

BOOL parse_shaking_par(char* fname, int *seq_first, int *seq_last, int *max_points, int *max_frames)
{
	char val[20][256];
	if (! read_strings(fname, 4, val))	// function displays a message on error
		return FALSE;
	*seq_first  = atoi(val[0]);
	*seq_last   = atoi(val[1]);
	*max_points = atoi(val[2]);
	*max_frames = atoi(val[3]);
	return TRUE;
}


// --- ptv_is file ---
// 730
//  -1  287     21.024     35.334     -1.597
//  -1   -2      8.221     43.613    -38.305
//  -1    0     12.323     39.568    -10.853
//  -1    6     -7.815    -33.691     -1.733


BOOL parse_ptv_is(char* pathname, ptv_is **plist, int *pcount)
{
    ptv_is *list, P;
    int i, n;
    BOOL ok;

    FILE *fp = fopen(pathname, "r");
    if (!fp) {
        printf("Can't open %s\n", pathname); 
		*plist  = NULL;
		*pcount = 0;
        return FALSE;
    }
    ok = fscanf(fp, "%d\n", &n) == 1;
    if (ok) {
        list = (ptv_is*) malloc(n * sizeof(ptv_is));
        for (i=0; ok && i<n; i++) {
            ok = fscanf(fp, "%d%d%lf%lf%lf", &P.prev, &P.next, &P.X, &P.Y, &P.Z) == 5;
            list[i] = P;
        }
    }
    fclose(fp);
    if (!ok) {
		free (list);
        printf("Error in reading %s\n", pathname);
		*pcount = 0;
		*plist  = NULL;
	}
	else {
		*pcount = n;
		*plist  = list;
	}
    return ok;
}



// --- rt_is file ---
// 254
//   1     0.604    -1.154    46.558  484  450  518  538
//   2    -0.271     0.835    70.853  404  389  455  459
//   3    -3.468     3.179    45.896  449  492  559  505

BOOL parse_rt_is(char* pathname, rt_is **plist, int *pcount)
{
    rt_is  P;
    int i, n;
    FILE* fp;
    BOOL ok;


    fp = fopen(pathname, "r");
    if (!fp) {
        printf("Can't open %s\n", pathname); 
		*plist  = NULL;
		*pcount = 0;
        return FALSE;
    }
    ok = fscanf(fp, "%d\n", &n) == 1;
    if (ok) {
        *plist = malloc(n * sizeof(rt_is));
        for (i=0; ok && i<n; i++) {
            ok = fscanf(fp, "%d%lf%lf%lf%d%d%d%d", &P.id, &P.X, &P.Y, &P.Z, 
                    &P.tnr[0], &P.tnr[1], &P.tnr[2], &P.tnr[3]) == 8;
            (*plist)[i] = P;
        }
    }
    fclose(fp);
    if (!ok) {
		free (plist);
        printf("Error in reading %s\n", pathname);
		*pcount = 0;
		*plist  = NULL;
	}
	else {
		*pcount = n;
	}
    return ok;
}

// --- the start of a target file ---
// 980
//   0  157.9141  122.0844    12     4     4   693    -1
//   1  428.3399  125.6797     9     3     3   612    -1
//   2  628.1592  128.5085    10     4     3   710    -1

BOOL parse_targetfile(char* pathname, target **tlist, int *pcount)
{
    target* T;
    int i, n;
    FILE* fp;
    BOOL ok;

	fp = fopen(pathname, "r");
    if (!fp) {
        printf("Can't open %s\n", pathname); 
		*tlist  = T = NULL;
		*pcount = 0;
        return FALSE;
    }
    ok = fscanf(fp, "%d\n", &n) == 1;
    if (ok)
        T = malloc(n * sizeof(target));
    for (i=0; ok && i<n; i++) {
        ok = fscanf(fp, "%d%lf%lf%d%d%d%d%d\n", 
      &T[i].pnr, &T[i].x, &T[i].y, &T[i].n, &T[i].nx, &T[i].ny, &T[i].sumg, &T[i].tnr) == 8;
    }
    fclose(fp);
	if (!ok) {
		free (T);
		printf("Error in reading %s\n", pathname);
		*tlist  = NULL;
		*pcount = 0;
	}
	else {
		*pcount = n;
		*tlist  = T;
	}
    return ok;
}

#include <direct.h>
#include "ERRNO.H"

BOOL createfolder(char* pathname)
{
	char path[256], *p;
	
	strcpy(path, pathname);
	for (p=path; *p; p++)
		if (*p == '\\') *p = '/';
	
	p = path;
	do {
		p = strchr(p,'/');
		if (p != NULL) *p = 0;
		if (_mkdir(path) != 0 && errno != EEXIST) {
			printf("Creating directory %s failed !!!\n", path);
			return FALSE;
		}
		if (p != NULL) *p++ = '/';
	} while (p != NULL);
	return TRUE;
}
