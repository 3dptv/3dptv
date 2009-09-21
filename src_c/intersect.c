/****************************************************************************

Routine:				intersect.c

Author/Copyright:		Hans-Gerd Maas

Address:				Institute of Geodesy and Photogrammetry
						ETH - Hoenggerberg
						CH - 8093 Zurich

Creation Date:			6.4.88
	
Description:			2 * (point + direction cosines) ==> intersection
	
Routines contained:		-

****************************************************************************/
#include "ptv.h"



/*-------------------------------------------------------------------------*/
/* 2 cameras */




void intersect_rt (X1,Y1,Z1,a1,b1,c1,X2,Y2,Z2,a2,b2,c2,X,Y,Z)

double  X1, Y1, Z1, a1, b1, c1, X2, Y2, Z2, a2, b2, c2, *X, *Y, *Z;

/* intersection, given two points with direction cosines */
/* only valid, if Z1 = Z2 = 0 , which is the case after ray tracing */

{
	if (fabs(b1-b2) > fabs(a1-a2))	*Z = (Y2-Y1) / ((b1/c1) - (b2/c2));
	else							*Z = (X2-X1) / ((a1/c1) - (a2/c2));
	
	*X = ( X1 + X2  +  *Z * (a1/c1 + a2/c2)) /2;
	*Y = ( Y1 + Y2  +  *Z * (b1/c1 + b2/c2)) /2;
}





void intersect (X1,Y1,Z1,a1,b1,c1,X2,Y2,Z2,a2,b2,c2,X,Y,Z)

double  X1, Y1, Z1, a1, b1, c1, X2, Y2, Z2, a2, b2, c2, *X, *Y, *Z;

/* intersection, given two points with direction cosines */

{
	if (fabs(b1-b2) > fabs(a1-a2))
		*Z = (Y2 - Y1 + Z1*(b1/c1) - Z2*(b2/c2)) / ((b1/c1) - (b2/c2));
	else
		*Z = (X2 - X1 + Z1*(a1/c1) - Z2*(a2/c2)) / ((a1/c1) - (a2/c2));
	
	*X = ( X1 + X2  +  *Z * (a1/c1 + a2/c2)) / 2;
	*Y = ( Y1 + Y2  +  *Z * (b1/c1 + b2/c2)) / 2;
}





/*-------------------------------------------------------------------------*/
/* 3 cameras */





void
intersect_rt_3 (X1,Y1,Z1,a1,b1,c1,X2,Y2,Z2,a2,b2,c2,X3,Y3,Z3,a3,b3,c3, X,Y,Z)

double	X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, X3,Y3,Z3,a3,b3,c3, *X,*Y,*Z;

/* intersection, given three points with direction cosines */

{
	int		i, n_max;
	double	base[6], max_base=0;

	/* find the maximum base component to decide, wether to take 
	   X or Y for determination of Z */
	   /* should be taken out as constant */
	base[0] = fabs (X1 + X2 - 2*X3);
	base[1] = fabs (X1 + X3 - 2*X2);
	base[2] = fabs (X2 + X3 - 2*X1);
	base[3] = fabs (Y1 + Y2 - 2*Y3);
	base[4] = fabs (Y1 + Y3 - 2*Y2);
	base[5] = fabs (Y2 + Y3 - 2*Y1);

	for (i=0; i<6; i++)	if (base[i] > max_base)
	{
		max_base = base[i];
		n_max = i;
	}

	switch (n_max)
	{
		case 0:	*Z = (2*X3 - X1 - X2) / ((a1/c1) + (a2/c2) - (2*a3/c3));
				break;
		case 1:	*Z = (2*X2 - X1 - X3) / ((a1/c1) + (a3/c3) - (2*a2/c2));
				break;
		case 2:	*Z = (2*X1 - X2 - X3) / ((a2/c2) + (a3/c3) - (2*a1/c1));
				break;
		case 3:	*Z = (2*Y3 - Y1 - Y2) / ((b1/c1) + (b2/c2) - (2*b3/c3));
				break;
		case 4:	*Z = (2*Y2 - Y1 - Y3) / ((b1/c1) + (b3/c3) - (2*b2/c2));
				break;
		case 5:	*Z = (2*Y1 - Y2 - Y3) / ((b2/c2) + (b3/c3) - (2*b1/c1));
				break;
	}

	*X = (X1+(*Z)*(a1/c1) + X2+(*Z)*(a2/c2) + X3+(*Z)*(a3/c3)) / 3;
	*Y = (Y1+(*Z)*(b1/c1) + Y2+(*Z)*(b2/c2) + Y3+(*Z)*(b3/c3)) / 3;
}





void
intersect_rt_3m (X1,Y1,Z1,a1,b1,c1,X2,Y2,Z2,a2,b2,c2,X3,Y3,Z3,a3,b3,c3, X,Y,Z)

double	X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, X3,Y3,Z3,a3,b3,c3, *X,*Y,*Z;

/* intersection, given three points with direction cosines */

{
	int		i, n_max1=0, n_max2;
	double	base[6], max_base1=0, max_base2=0, Za, Zb;

	/* find the maximum base component to decide, wether to take 
	   X or Y for determination of Z */
	   /* should be taken out as constant */
	base[0] = fabs (X1 + X2 - 2*X3);
	base[1] = fabs (X1 + X3 - 2*X2);
	base[2] = fabs (X2 + X3 - 2*X1);
	base[3] = fabs (Y1 + Y2 - 2*Y3);
	base[4] = fabs (Y1 + Y3 - 2*Y2);
	base[5] = fabs (Y2 + Y3 - 2*Y1);

	for (i=0; i<6; i++)
	{
		if (base[i] > max_base1)
		{
			max_base2 = max_base1;	max_base1 = base[i];
			n_max2 = n_max1;		n_max1 = i;
		}
		else if (base[i] > max_base2)
		{
			max_base2 = base[i];	n_max2 = i;
		}
	}

	/* compute *Z from the 2 longest base components,
	   what was empirically found to be the most accurate value */

	switch (n_max1)
	{
		case 0:	Za = (2*X3 - X1 - X2) / ((a1/c1) + (a2/c2) - (2*a3/c3));
				break;
		case 1:	Za = (2*X2 - X1 - X3) / ((a1/c1) + (a3/c3) - (2*a2/c2));
				break;
		case 2:	Za = (2*X1 - X2 - X3) / ((a2/c2) + (a3/c3) - (2*a1/c1));
				break;
		case 3:	Za = (2*Y3 - Y1 - Y2) / ((b1/c1) + (b2/c2) - (2*b3/c3));
				break;
		case 4:	Za = (2*Y2 - Y1 - Y3) / ((b1/c1) + (b3/c3) - (2*b2/c2));
				break;
		case 5:	Za = (2*Y1 - Y2 - Y3) / ((b2/c2) + (b3/c3) - (2*b1/c1));
				break;
	}
	switch (n_max2)
	{
		case 0:	Zb = (2*X3 - X1 - X2) / ((a1/c1) + (a2/c2) - (2*a3/c3));
				break;
		case 1:	Zb = (2*X2 - X1 - X3) / ((a1/c1) + (a3/c3) - (2*a2/c2));
				break;
		case 2:	Zb = (2*X1 - X2 - X3) / ((a2/c2) + (a3/c3) - (2*a1/c1));
				break;
		case 3:	Zb = (2*Y3 - Y1 - Y2) / ((b1/c1) + (b2/c2) - (2*b3/c3));
				break;
		case 4:	Zb = (2*Y2 - Y1 - Y3) / ((b1/c1) + (b3/c3) - (2*b2/c2));
				break;
		case 5:	Zb = (2*Y1 - Y2 - Y3) / ((b2/c2) + (b3/c3) - (2*b1/c1));
				break;
	}
	*Z = (Za + Zb) / 2;


	*X = (X1+(*Z)*(a1/c1) + X2+(*Z)*(a2/c2) + X3+(*Z)*(a3/c3)) / 3;
	*Y = (Y1+(*Z)*(b1/c1) + Y2+(*Z)*(b2/c2) + Y3+(*Z)*(b3/c3)) / 3;
}





void intersect_3 (X1,Y1,Z1,a1,b1,c1,X2,Y2,Z2,a2,b2,c2,X3,Y3,Z3,a3,b3,c3, X,Y,Z)

double	X1,Y1,Z1,a1,b1,c1, X2,Y2,Z2,a2,b2,c2, X3,Y3,Z3,a3,b3,c3, *X,*Y,*Z;

/* intersection, given three points with direction cosines */

{
	int		i, n_min;
	double	base[6], min_base=1e20;

	/* find the minimum base component to decide, wether to take 
	   X or Y for determination of Z */
	base[0] = fabs (a1-a2);  base[1] = fabs (a1-a3);  base[2] = fabs (a2-a3);
	base[3] = fabs (b1-b2);  base[4] = fabs (b1-b3);  base[5] = fabs (b2-b3);

	for (i=0; i<6; i++)	if (base[i] < min_base)
	{
		min_base = base[i];
		n_min = i;
	}
	
	if (n_min > 2)		/* min base in y-direction => take x */
		*Z = (X2 + X3 - 2*X1 + 2*Z1*(a1/c1) - Z2*(a2/c2) - Z3*(a3/c3))
		   / ((2*a1/c1) - (a2/c2) - (a3/c3));
	else				/* min base in y-direction => take x */
		*Z = (Y2 + Y3 - 2*Y1 + 2*Z1*(b1/c1) - Z2*(b2/c2) - Z3*(b3/c3))
		   / ((2*b1/c1) - (b2/c2) - (b3/c3));

	*X = (X1+(*Z-Z1)*(a1/c1) + X2+(*Z-Z2)*(a2/c2) + X3+(*Z-Z3)*(a3/c3)) / 3;
	*Y = (Y1+(*Z-Z1)*(b1/c1) + Y2+(*Z-Z2)*(b2/c2) + Y3+(*Z-Z3)*(b3/c3)) / 3;
}




