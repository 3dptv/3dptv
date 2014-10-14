/****************************************************************************

Routine:			ray_tracing

Author/Copyright:	Hans-Gerd Maas

Address:			Institute of Geodesy and Photogrammetry
					ETH - Hoenggerberg
					CH - 8093 Zurich

Creation Date:		21.4.88
	
Description:		traces one ray, given by image coordinates,
					exterior and interior orientation 
					through dufferent media
					(see Hoehle, Manual of photogrammetry)
	
Routines contained: 	   -

****************************************************************************/

#include "ptv.h"




void ray_tracing (x,y,Ex,I,mm,Xb2,Yb2,Zb2,a3,b3,c3)

double		x, y;
Exterior	Ex;
Interior	I;
mm_np		mm;
double		*Xb2, *Yb2, *Zb2, *a3, *b3, *c3;

/* ray-tracing, see HOEHLE and Manual of Photogrammetry */

{
	double	a1, b1, c1, a2, b2, c2, Xb1, Yb1, Zb1, d1, d2, cosi1, cosi2,
			vect1[3], vect2[3], factor, s2;

	s2 = sqrt (x*x + y*y + I.cc*I.cc);
	
	/* direction cosines in image coordinate system */
	vect1[0] = x/s2;  vect1[1] = y/s2;	  vect1[2] = -I.cc/s2;

	matmul (vect2, Ex.dm, vect1, 3,3,1);
	 
	/* direction cosines in space coordinate system , medium n1 */
	a1 = vect2[0];  b1 = vect2[1];  c1 = vect2[2];  
	
       	d1 = -(Ex.z0 - mm.d[0]) / c1;

	/* point on the horizontal plane between n1,n2 */
	Xb1 = Ex.x0 + d1*a1;  Yb1 = Ex.y0 + d1*b1;	Zb1 = Ex.z0 + d1*c1;
	
	cosi1 = c1;
	factor = cosi1 * mm.n1/mm.n2[0]
			 + sqrt (1 - (mm.n1*mm.n1)/(mm.n2[0]*mm.n2[0])
			 + (cosi1*cosi1)*(mm.n1*mm.n1)/(mm.n2[0]*mm.n2[0]));

	/* direction cosines in space coordinate system , medium n2 */
	a2 = a1 * mm.n1/mm.n2[0];
	b2 = b1 * mm.n1/mm.n2[0];
	c2 = c1 * mm.n1/mm.n2[0] - factor;
	
	d2 = -mm.d[0]/c2;

	/* point on the horizontal plane between n2,n3 */
	*Xb2 = Xb1 + d2*a2;  *Yb2 = Yb1 + d2*b2;  *Zb2 = Zb1 + d2*c2;
	
	cosi2 = c2;
	factor = cosi2 * mm.n2[0]/mm.n3 
			 + sqrt (1 - (mm.n2[0]*mm.n2[0])/(mm.n3*mm.n3)
			 + (cosi2*cosi2)*(mm.n2[0]*mm.n2[0])/(mm.n3*mm.n3));

	/* direction cosines in space coordinate system , medium mm.n3 */
	*a3 = a2 * mm.n2[0]/mm.n3;
	*b3 = b2 * mm.n2[0]/mm.n3;
	*c3 = c2 * mm.n2[0]/mm.n3 - factor;
}

