/*********************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	September'97

Description:	       	display of image sequences

Routines contained:    	flow_demo_c

**********************************************************************/
#include "ptv.h"

int flow_demo_c (ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, i_seq, nr, j,pft_version=3;
  char	        name[128];
  unsigned char	*imgf;
  Tk_PhotoHandle img_handle;
  Tk_PhotoImageBlock img_block;

  nr = atoi(argv[1]);

  fpp = fopen_r ("parameters/sequence.par");
  for (i=0; i<4; i++)
    fscanf (fpp, "%s\n", seq_name[i]);     /* name of sequence */
  fscanf (fpp,"%d\n", &seq_first);
  fscanf (fpp,"%d\n", &seq_last);
  fclose (fpp);


  /* allocate memory */
  imgf = (unsigned char *) calloc (imgsize, 1);

  fpp = fopen ("parameters/pft_version", "r");
  if (fpp){
      fscanf (fpp, "%d\n", &pft_version);
      fclose (fpp);
  }

  /* load and display images */
  for (i_seq=seq_first; i_seq<=seq_last; i_seq++){
      compose_name_plus_nr (seq_name[nr], "", i_seq, name);
      fp1 = fopen_r (name);	if (! fp1)	return TCL_OK;
      sprintf (buf, "display camera %d, image %d", nr+1, i_seq);
      Tcl_SetVar(interp, "tbuf", buf, TCL_GLOBAL_ONLY);
      Tcl_Eval(interp, ".text delete 2");
      Tcl_Eval(interp, ".text insert 2 $tbuf");

      read_image (interp, name, imgf);
      fclose (fp1);

      img_handle = Tk_FindPhoto( interp, "temp");
      Tk_PhotoGetImage (img_handle, &img_block);

      sprintf(buf, "newimage %d", nr+1);
      Tcl_Eval(interp, buf);

      
	  if(pft_version==4){
         sprintf (filename, "%s%s", name,"_targets");
	      /* read targets of camera nr*/
	     nt4[3][nr]=0;

	     fp1= fopen (filename, "r");
	     if (! fp1) printf("Can't open ascii file: %s\n", filename);

         fscanf (fp1, "%d\n", &nt4[3][nr]);
         for (j=0; j<nt4[3][nr]; j++){
	          fscanf (fp1, "%4d %lf %lf %d %d %d %d %d\n",
		      &pix[nr][j].pnr, &pix[nr][j].x,
		      &pix[nr][j].y, &pix[nr][j].n ,
		      &pix[nr][j].nx ,&pix[nr][j].ny,
		      &pix[nr][j].sumg, &pix[nr][j].tnr);
	     }
         fclose (fp1);

	     num[nr] = nt4[3][nr];

		 if (display){
	         for (j=0; j<num[nr]; j++){
	             drawcross (interp, (int) pix[nr][j].x, (int) pix[nr][j].y,cr_sz, nr, "blue");
	         }
		 }
	  }
      Tcl_Eval(interp, "update idletasks");
    }

  free (imgf);
  return TCL_OK;
}

