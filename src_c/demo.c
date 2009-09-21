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
  int i, i_seq, nr;
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

  /* load and display images */
  for (i_seq=seq_first; i_seq<=seq_last; i_seq++)
    {
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
      Tcl_Eval(interp, "update idletasks");
    }

  free (imgf);
  return TCL_OK;
}

