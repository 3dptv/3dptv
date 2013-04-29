/*********************************************************************

Author/Copyright:	 Jochen Willneff

Address:			 Institute of Geodesy and Photogrammetry
					 ETH - Hoenggerberg
					 CH - 8093 Zurich

Creation Date:		 September'97

Description:		display of image sequences

Routines contained: flow_demo_c

**********************************************************************/

#include "ptv.h"

int flow_demo_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
	int i, i_seq, nr, j,pft_version=3, num_points;
	char		   name[128];
	unsigned char  *imgf;
	Tk_PhotoHandle img_handle;
	Tk_PhotoImageBlock img_block;

	nr = atoi(argv[1]);

	fpp = fopen_rp ("parameters/sequence.par");		// replaced fopen_r, ad holten 12-2012
	if (!fpp) return TCL_OK;
	for (i=0; i<4; i++)
		fscanf (fpp, "%s\n", seq_name[i]);	   /* name of sequence */
	fscanf (fpp,"%d\n", &seq_first);
	fscanf (fpp,"%d\n", &seq_last);
	fclose (fpp);


	/* allocate memory */
	imgf = (unsigned char *) calloc (imgsize, 1);

	fpp = fopen ("parameters/pft_version.par", "r");
	if (fpp) {
		fscanf (fpp, "%d\n", &pft_version);
		pft_version = pft_version+3;
		fclose (fpp);
	}
	else {
		fpp = fopen ("parameters/pft_version.par", "w");
		fprintf(fpp,"%d\n", 0);
		fclose(fpp);
	}

	/* load and display images */
	for (i_seq=seq_first; i_seq<=seq_last; i_seq++) {
		compose_name_plus_nr (seq_name[nr], "", i_seq, name);
		fp1 = fopen_rp(name);				// replaced fopen_r, ad holten 12-2012
		if (!fp1) return TCL_OK;

		sprintf (buf, "display camera %d, image %d", nr+1, i_seq);
		Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
		Tcl_Eval(interp, ".text delete 2");
		Tcl_Eval(interp, ".text insert 2 $tbuf");

		read_image (interp, name, imgf);
		fclose (fp1);

		img_handle = Tk_FindPhoto( interp, "temp");
		Tk_PhotoGetImage (img_handle, &img_block);

		sprintf(buf, "newimage %d", nr+1, 0.5, 0.5, 1, 0);
		Tcl_Eval(interp, buf);

		if(pft_version==4) {
			sprintf (filename, "%s%s", name,"_targets");
			/* read targets of camera nr*/
			nt4[3][nr]=0;

			fp1 = fopen_rp(filename);			// replaced fopen, ad holten 12-2012
			if (!fp1) return TCL_OK;

			fscanf (fp1, "%d\n", &nt4[3][nr]);
			for (j=0; j<nt4[3][nr]; j++){
				fscanf (fp1, "%4d %lf %lf %d %d %d %d %d\n",
					&pix[nr][j].pnr,  &pix[nr][j].x,
					&pix[nr][j].y,    &pix[nr][j].n ,
					&pix[nr][j].nx ,  &pix[nr][j].ny,
					&pix[nr][j].sumg, &pix[nr][j].tnr);
			}
			fclose (fp1);
			num[nr] = nt4[3][nr];
			if (display) {
				for (j=0; j<num[nr]; j++)
					drawcross (interp, (int) pix[nr][j].x, (int) pix[nr][j].y,cr_sz, nr, "blue");
				printf ("drawing %d 2d ", num[nr]);
			}
			
			sprintf (filename, "res/rt_is.%d", i_seq);
			fp1 = fopen (filename, "r");
			if (fp1) {
				fscanf (fp1, "%d\n", &num_points);
				for (j=0; j<num_points; j++) {
					// ad holten, 12-2012, replaced the next lines with bugs
					//		if (n_img==4) {
					//			fscanf(fp1, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[j].x, &fix[j].y, &fix[j].z,
					//				&geo[0][j].pnr, &geo[1][j].pnr, &geo[2][j].pnr, &geo[3][j].pnr);
					//		}
					//		if (n_img==3) {
					//			fscanf(fp1, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[j].x, &fix[j].y, &fix[j].z,
					//				&geo[0][j].pnr, &geo[1][j].pnr, &geo[2][j].pnr);
					//		}
					//		if (n_img==2){ // Alex's patch. 24.09.09. Working on Wesleyan data of 2 cameras only
					//			fscanf(fp1, "%d %lf %lf %lf %d %d %d %d\n",
					//				&dumy, &fix[j].x, &fix[j].y, &fix[j].z,
					//				&geo[0][j].pnr, &geo[1][j].pnr);
					//		}
					// by:
					for (i=0; i<4; i++)	geo[i][j].pnr = -1;
					fscanf(fp1, "%*d %lf %lf %lf %d %d %d %d\n", 
						&fix[j].x, &fix[j].y, &fix[j].z, &geo[0][j].pnr, &geo[1][j].pnr, &geo[2][j].pnr, &geo[3][j].pnr);
					for (j=3; j>=n_img; j--) geo[i][j].pnr = -1;
				}
				fclose (fp1);
				if (display) {
					for (j=0; j<num_points; j++) {
						img_coord (fix[j].x,fix[j].y,fix[j].z, Ex[nr], I[nr], G[nr], ap[nr], mmp, &pix[nr][j].x,&pix[nr][j].y);
						metric_to_pixel (pix[nr][j].x,pix[nr][j].y, imx,imy, pix_x,pix_y, &pix[nr][j].x,&pix[nr][j].y, chfield);
						//if(geo[nr][j].pnr>-1){
						//	  drawcross (interp, (int) pix[nr][geo[nr][j].pnr].x, (int) pix[nr][geo[nr][j].pnr].y,cr_sz, nr, "yellow");
						//}
						drawcross (interp, (int) pix[nr][j].x, (int) pix[nr][j].y,cr_sz, nr, "red");
					}
					printf ("and %d corresponding 3d positions for frame %d\n", num_points,i_seq);
				}
			}	 
		}
		Tcl_Eval(interp, "update idletasks");
	}
	printf ("done\n\n");

	free (imgf);
	return TCL_OK;
}

