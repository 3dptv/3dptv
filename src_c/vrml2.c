/****************************************************************************

Author/Copyright:      	Jochen Willneff

Address:	      	Institute of Geodesy and Photogrammetry
		       	ETH - Hoenggerberg
		       	CH - 8093 Zurich

Creation Date:	       	end of 99 ...
	
Description:	       	creates VRML files from tracks and
                        detected particles
	
Routines contained:     vrmltracks_c, vrmldetections_c, vrmldettracks_c

****************************************************************************/
#include "ptv.h"

int vrmltracks_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, anz1, anz2, m, j, k;
  FILE *fp1, *fp2;
  char val[256];
  vector *line1, *line2;
  double color, ymin=0, ymax=0, cubes;
  double mx, my, mz, dx, dy, dz, du, dl, rotz, rotx;

  
  /* open file for line elements */
  fp2 = fopen ("tracks.wrl", "w");
  
  /* create header, background and coordsys for vrml-file */
  fprintf(fp2, "#VRML V2.0 utf8\n\n");
  fprintf(fp2, " Background { skyColor [ 1 1 1] }\n\n");
  
  /* create boundaries from object volume */
  volumedimension(&X_lay[1], &X_lay[0], &ymax, &ymin, &Zmax_lay[1], &Zmin_lay[0]);
  cubes=(Zmax_lay[1]-Zmin_lay[0])/200;
  
  /* create viewpoint */
  fprintf(fp2, "  DEF V1 Viewpoint {\n");
  fprintf(fp2, "   position   %7.3f %7.3f %7.3f\n",
	  (X_lay[0]+X_lay[1])/2,(ymax+ymin)/2,Zmax_lay[1]);      
  fprintf(fp2, "   orientation   1 0 0 0 \n"); 
  fprintf(fp2, "   description \"Start\" }\n\n");    
  
  /* create object volume */ 
  fprintf(fp2, "  Shape { geometry IndexedLineSet {\n");
  fprintf(fp2, "   coord Coordinate { point  [\n");
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0],ymin, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymin,Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymin, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax,Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f, ] }\n", X_lay[1], ymin, Zmax_lay[1]);
  fprintf(fp2, "   coordIndex [ \n");
  fprintf(fp2, "    0, 1, 2, 3, 0, -1,\n");
  fprintf(fp2, "    0, 4, -1,\n");
  fprintf(fp2, "    1, 5, -1,\n");
  fprintf(fp2, "    2, 6, -1,\n");
  fprintf(fp2, "    3, 7, -1,\n");
  fprintf(fp2, "    4, 5, 6, 7, 4, -1 ]\n }"); 
  fprintf(fp2, "    appearance Appearance {material Material { diffuseColor 0 0 1 } } }\n");

  
  fprintf(fp2, "\n\n# start trajectories\n\n");
  /* read trackfile from ptv and create vectorfield */
 
  for (i=seq_first; i<seq_last;i++)
    {
      if (i < 10)             sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/ptv_is.%2d",  i);
      else       sprintf (val, "res/ptv_is.%3d",  i);
  
      printf("Create VRML, read file: %s\n", val);     
      fp1 = fopen (val, "r");
      
      color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));

      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) calloc (anz1, sizeof (vector));
      for (j=0;j<anz1;j++) {
	fscanf (fp1, "%d\n", &line1[j].p);
	fscanf (fp1, "%d\n", &line1[j].n);
	fscanf (fp1, "%lf\n", &line1[j].x1);
	fscanf (fp1, "%lf\n", &line1[j].y1);
	fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      /* read next time step */     
      if (i+1 < 10)             sprintf (val, "res/ptv_is.%1d", i+1);
      else if (i+1 < 100)       sprintf (val, "res/ptv_is.%2d",  i+1);
      else       sprintf (val, "res/ptv_is.%3d",  i+1);
      
      fp1 = fopen (val, "r");     
      fscanf (fp1,"%d\n", &anz2);
      line2 = (vector *) malloc (anz2 * sizeof (vector));
      
      for (j=0;j<anz2;j++) {
	fscanf (fp1, "%d\n", &line2[j].p);
	fscanf (fp1, "%d\n", &line2[j].n);
	fscanf (fp1, "%lf\n", &line2[j].x1);
	fscanf (fp1, "%lf\n", &line2[j].y1);
	fscanf (fp1, "%lf\n", &line2[j].z1);
      }
      fclose (fp1);
      
      fprintf(fp2, "# time step %d\n", i);

	  fprintf(fp2, "  Shape { appearance DEF step%d Appearance {\n", i);
      fprintf(fp2, "   material Material {\n");
      fprintf(fp2, "    ambientIntensity 0.25\n");
      fprintf(fp2, "    diffuseColor 1 %.4f 0 } } }\n\n", color);
 
      for(j=0;j<anz1;j++) { 	
	m = line1[j].n;
	if (m >= 0) {	
	 
	  fprintf(fp2, "    Transform {translation %7.3f %7.3f %7.3f\n",
		  line1[j].x1, line1[j].y1, line1[j].z1);
	  fprintf(fp2, "    children [\n");
	  fprintf(fp2, "    Shape { geometry Box { size %3.2f %3.2f %3.2f }\n      appearance USE step%d }\n", 
		  cubes, cubes, cubes, i);
	  
	  fprintf(fp2, "    Shape { geometry Sphere { radius %3.2f }\n      appearance USE step%d }\n", cubes, i);
	  
	  fprintf(fp2, "    ] }\n");

	  /* create cubes of last time step */	 
	  if(i ==seq_last-1) {
	    
	  fprintf(fp2, "    Transform {translation %7.3f %7.3f %7.3f\n",
		  line2[m].x1, line2[m].y1, line2[m].z1);
	  fprintf(fp2, "    children [\n"); 
	  fprintf(fp2, "    Shape { geometry Box { size %3.2f %3.2f %3.2f }\n      appearance USE step%d }\n", 
		  cubes, cubes, cubes, i);
	   
	  fprintf(fp2, "    Shape { geometry Sphere { radius %3.2f }\n      appearance USE step%d }\n", cubes, i);
	  
	  fprintf(fp2, "    ] }\n");	  
	  }
	  
	  /* create last cubes of ending trajectories  */
	  
	  if(line2[m].n <= 0) {
	   
	  fprintf(fp2, "    Transform {translation %7.3f %7.3f %7.3f\n",
		  line2[m].x1, line2[m].y1, line2[m].z1);	  
	  fprintf(fp2, "    children [\n");
	  fprintf(fp2, "    Shape { geometry Box { size %3.2f %3.2f %3.2f }\n      appearance USE step%d }\n", 
		  cubes, cubes, cubes, i );
	  
	  fprintf(fp2, "    Shape { geometry Sphere { radius %3.2f }\n      appearance USE step%d }\n", cubes, i);
	  
	  fprintf(fp2, "    ] }\n");	  
	  }

	    mx=(line1[j].x1+line2[m].x1)/2;
	    my=(line1[j].y1+line2[m].y1)/2; 
	    mz=(line1[j].z1+line2[m].z1)/2; 
	    dx=line1[j].x1-line2[m].x1;
	    dy=line1[j].y1-line2[m].y1;
	    dz=line1[j].z1-line2[m].z1;
	    du=sqrt(dx*dx+dy*dy);
	    dl=sqrt(dx*dx+dy*dy+dz*dz);
	    
	    rotz=0;
	    if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
	    rotx=0;
	    if(du == 0.0) {rotx=M_PI/2;}

	    if(du != 0.0) {
	      if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
	      if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
	      if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
	      if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
	      if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
	      if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
	      if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
	      if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
	    }
	    
	    
	    fprintf(fp2, "    Transform { translation %7.3f %7.3f %7.3f\n",mx, my, mz);
	    fprintf(fp2, "    children [ Transform {rotation 0 0 1 %7.5f\n",rotz);
	    fprintf(fp2, "    children [ Transform {rotation 1 0 0 %7.5f \n",rotx);	
	    fprintf(fp2, "    children [ Shape { geometry Cylinder { radius %3.2f height %3.2f } \n      appearance USE step%d } ]\n",cubes/4, dl, i);
	    fprintf(fp2, "    children [ Shape { geometry Box {size %3.2f %3.2f %3.2f } \n      appearance USE step%d } ]\n\n",cubes/2, dl,cubes/2, i);
		fprintf(fp2, "    } ] } ] }\n");
		
	    /* end of cylinder */ 

	    
	  fprintf(fp2, "    Shape { geometry IndexedLineSet {\n");
	  fprintf(fp2, "     coord Coordinate { point [\n");
	  fprintf(fp2, "      %7.3f %7.3f %7.3f,\n",line1[j].x1, line1[j].y1, line1[j].z1);
	  fprintf(fp2, "      %7.3f %7.3f %7.3f, ] }\n", line2[m].x1, line2[m].y1, line2[m].z1);
	  fprintf(fp2, "     coordIndex [ 0, 1, -1] } \n      appearance USE step%d}\n\n", i);	  
	   
	}	
      }
  
      fprintf(fp2, "# end of time step %d\n\n", i);     
      strcpy(val, "");
      free(line1); free(line2);
    }  /* end of sequence loop */
  
  fprintf(fp2, "# trajectories finished\n");
  
  fclose(fp2); 
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 \"Tracks written to VRML-File: tracks.wrl\"");
  Tcl_Eval(interp, "update idletasks");  

  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");
  
  return TCL_OK;
}


int vrmldetections_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, anz1, j, k, dumy;
  FILE *fp1, *fp2;
  char val[256];
  vector *line1;
  double color, ymin=0, ymax=0, cubes;

  /* open file for line elements */
  fp2 = fopen ("detections.wrl", "w");

  /* create header, background and coordsys for vrml-file */
  fprintf(fp2, "#VRML V2.0 utf8\n\n");
  fprintf(fp2, " Background { skyColor [ 1 1 1] }\n\n");
  
  /* create boundaries from object volume */
  volumedimension(&X_lay[1], &X_lay[0], &ymax, &ymin, &Zmax_lay[1], &Zmin_lay[0]);
  cubes=(Zmax_lay[1]-Zmin_lay[0])/600;
  
  /* create viewpoint */
  fprintf(fp2, "  DEF V1 Viewpoint {\n");
  fprintf(fp2, "   position   %7.3f %7.3f %7.3f\n",
	  (X_lay[0]+X_lay[1])/2,(ymax+ymin)/2,Zmax_lay[1]);      
  fprintf(fp2, "   orientation   1 0 0 0 \n"); 
  fprintf(fp2, "   description \"Start\" }\n\n");    
  
  /* create object volume */ 
  fprintf(fp2, "  Shape { geometry IndexedLineSet {\n");
  fprintf(fp2, "   coord Coordinate { point  [\n");
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0],ymin, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymin,Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymin, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax,Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f, ] }\n", X_lay[1], ymin, Zmax_lay[1]);
  fprintf(fp2, "   coordIndex [ \n");
  fprintf(fp2, "    0, 1, 2, 3, 0, -1,\n");
  fprintf(fp2, "    0, 4, -1,\n");
  fprintf(fp2, "    1, 5, -1,\n");
  fprintf(fp2, "    2, 6, -1,\n");
  fprintf(fp2, "    3, 7, -1,\n");
  fprintf(fp2, "    4, 5, 6, 7, 4, -1 ]\n }"); 
  fprintf(fp2, "    appearance Appearance {material Material { diffuseColor 0 0 1 } } }\n\n\n"); 

  /* read trackfile from ptv and create particle positions */
  for (i=seq_first; i<=seq_last ;i++)
    { 
      if (i < 10)             sprintf (val, "res/rt_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/rt_is.%2d",  i);
      else       sprintf (val, "res/rt_is.%3d",  i);
      
      printf("Create VRML, read file: %s\n", val);         
      fp1 = fopen (val, "r");
      
      color = ((double)(i-seq_first))/((double)(seq_last+1-seq_first));
      
      fscanf (fp1,"%d\n", &anz1);
      line1 = (vector *) calloc (anz1, sizeof (vector));
      for (j=0;j<anz1;j++) {
	fscanf (fp1, "%d %lf %lf %lf %d %d %d %d\n",
		&line1[j].p, &line1[j].x1, &line1[j].y1,
		&line1[j].z1, &dumy, &dumy, &line1[j].type, &dumy);
      }
      
      fclose (fp1);
      
      fprintf(fp2, "  Shape { appearance DEF step%d Appearance {\n", i);
      fprintf(fp2, "   material Material {\n");
      fprintf(fp2, "    ambientIntensity 0.25\n");
      fprintf(fp2, "    diffuseColor 1 %.4f 0 } } }\n\n", color);
            
      for(j=0;j<anz1;j++) {


	  fprintf(fp2, "    Transform {translation %7.3f %7.3f %7.3f\n",
		  line1[j].x1, line1[j].y1, line1[j].z1);
	  fprintf(fp2, "    children [\n"); 
	  fprintf(fp2, "    Shape { geometry Box { size %3.2f %3.2f %3.2f } appearance USE step%d }\n", 
		  cubes, cubes, cubes, i);
	  /*
	  fprintf(fp2, "    Shape { geometry Sphere { radius %3.2f } appearance USE step%d }\n", cubes, i);
	  */
	  fprintf(fp2, "    ] }\n");	  
	  }

      fprintf(fp2, "\n# end of time step %d\n\n", i); 

      strcpy(val, "");
      free(line1);

    }  /* end of sequence loop */
  
  fprintf(fp2, "# detections finished\n");
  
  fclose(fp2); 
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 \"Detections written to VRML-File: detections.wrl\"");
  Tcl_Eval(interp, "update idletasks");

  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");
  
  return TCL_OK;
}


int vrmldettracks_c(ClientData clientData, Tcl_Interp* interp, int argc, const char** argv)
{
  int i, anz1, anz2, m, j, k;
  FILE *fp1, *fp2;
  char val[256];
  vector *line1, *line2;
  double color, ymin=0, ymax=0, cubes;
  double mx, my, mz, dx, dy, dz, du, dl, rotz, rotx;
  
  /* open file for line elements */
  fp2 = fopen ("dt.wrl", "w");

  /* create boundaries from object volume */
  volumedimension(&X_lay[1], &X_lay[0], &ymax, &ymin, &Zmax_lay[1], &Zmin_lay[0]);
  cubes=(Zmax_lay[1]-Zmin_lay[0])/500;
  cubes=(ymax-ymin)/800;
  
  /* create header, background and coordsys for vrml-file */
  fprintf(fp2, "#VRML V2.0 utf8\n\n");
  fprintf(fp2, " Background { skyColor [ 1 1 1] }\n\n");
  
  /* create boundaries from object volume */
  volumedimension(&X_lay[1], &X_lay[0], &ymax, &ymin, &Zmax_lay[1], &Zmin_lay[0]);
  cubes=(Zmax_lay[1]-Zmin_lay[0])/600;
  
  /* create viewpoint */
  fprintf(fp2, "  DEF V1 Viewpoint {\n");
  fprintf(fp2, "   position   %7.3f %7.3f %7.3f\n",
	  (X_lay[0]+X_lay[1])/2,(ymax+ymin)/2,Zmax_lay[1]);      
  fprintf(fp2, "   orientation   1 0 0 0 \n"); 
  fprintf(fp2, "   description \"Start\" }\n\n");    
  
  /* create object volume */ 
  fprintf(fp2, "  Shape { geometry IndexedLineSet {\n");
  fprintf(fp2, "   coord Coordinate { point  [\n");
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0],ymin, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax, Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymin,Zmin_lay[0]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymin, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[0], ymax, Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f,\n", X_lay[1], ymax,Zmax_lay[1]);
  fprintf(fp2, "    %5.2f %5.2f %5.2f, ] }\n", X_lay[1], ymin, Zmax_lay[1]);
  fprintf(fp2, "   coordIndex [ \n");
  fprintf(fp2, "    0, 1, 2, 3, 0, -1,\n");
  fprintf(fp2, "    0, 4, -1,\n");
  fprintf(fp2, "    1, 5, -1,\n");
  fprintf(fp2, "    2, 6, -1,\n");
  fprintf(fp2, "    3, 7, -1,\n");
  fprintf(fp2, "    4, 5, 6, 7, 4, -1 ]\n }"); 
  fprintf(fp2, "    appearance Appearance {material Material { diffuseColor 0 0 1 } } }\n");

  fprintf(fp2, "\n\n# start trajectories\n\n");
  /* read trackfile from ptv and create vectorfield */
 
  for (i=seq_first; i<=seq_last;i++)
    {
      if (i < 10)             sprintf (val, "res/ptv_is.%1d", i);
      else if (i < 100)       sprintf (val, "res/ptv_is.%2d",  i);
      else       sprintf (val, "res/ptv_is.%3d",  i);
  
      printf("Create VRML, read file: %s\n", val);     
      fp1 = fopen (val, "r");
      
      color = ((double)(i-seq_first))/((double)(seq_last-1-seq_first));

      fscanf (fp1,"%d\n", &anz1);
      
      line1 = (vector *) calloc (anz1, sizeof (vector));
      for (j=0;j<anz1;j++) {
	fscanf (fp1, "%d\n", &line1[j].p);
	fscanf (fp1, "%d\n", &line1[j].n);
	fscanf (fp1, "%lf\n", &line1[j].x1);
	fscanf (fp1, "%lf\n", &line1[j].y1);
	fscanf (fp1, "%lf\n", &line1[j].z1);
      }

      strcpy(val, "");     
      fclose (fp1);

      if (i<seq_last) {
	/* read next time step */     
	if (i+1 < 10)             sprintf (val, "res/ptv_is.%1d", i+1);
	else if (i+1 < 100)       sprintf (val, "res/ptv_is.%2d",  i+1);
	else       sprintf (val, "res/ptv_is.%3d",  i+1);
	
	fp1 = fopen (val, "r");	
	fscanf (fp1,"%d\n", &anz2);
	line2 = (vector *) malloc (anz2 * sizeof (vector));
	
	for (j=0;j<anz2;j++) {
	  fscanf (fp1, "%d\n", &line2[j].p);
	  fscanf (fp1, "%d\n", &line2[j].n);
	  fscanf (fp1, "%lf\n", &line2[j].x1);
	  fscanf (fp1, "%lf\n", &line2[j].y1);
	  fscanf (fp1, "%lf\n", &line2[j].z1);
	}
	fclose (fp1);
      }
      
      
      fprintf(fp2, "  Shape { appearance DEF step%d Appearance {\n", i);
      fprintf(fp2, "   material Material {\n");
      fprintf(fp2, "    ambientIntensity 0.25\n");
      fprintf(fp2, "    diffuseColor 1 %.4f 0 } } }\n\n", color);
       
      for(j=0;j<anz1;j++) /* if( line1[j].z1 > -22) */ {
	
	
	fprintf(fp2, "    Transform {translation %7.3f %7.3f %7.3f\n",
		line1[j].x1, line1[j].y1, line1[j].z1);
	fprintf(fp2, "    children [\n");
	/*
	fprintf(fp2, "    Shape { geometry Sphere { radius %3.2f } appearance USE step%d }\n", cubes, i );
	*/
	fprintf(fp2, "    Shape { geometry Box {size %3.2f %3.2f %3.2f }  appearance USE step%d }\n",cubes, cubes, cubes, i);
	fprintf(fp2, " ] }\n\n");

	if (i<seq_last) {
	  m = line1[j].n;
	  if (m >= 0) {	

	  fprintf(fp2, "    Shape { geometry IndexedLineSet {\n");
	  fprintf(fp2, "     coord Coordinate { point [\n");
	  fprintf(fp2, "      %7.3f %7.3f %7.3f,\n",line1[j].x1, line1[j].y1, line1[j].z1);
	  fprintf(fp2, "      %7.3f %7.3f %7.3f, ] }\n", line2[m].x1, line2[m].y1, line2[m].z1);
	  fprintf(fp2, "     coordIndex [ 0, 1, -1] } \n      appearance USE step%d }\n\n", i);	  
		   

	    /* cylinder/cube to mark link */
	  /*
	    mx=(line1[j].x1+line2[m].x1)/2;
	    my=(line1[j].y1+line2[m].y1)/2; 
	    mz=(line1[j].z1+line2[m].z1)/2; 
	    dx=line1[j].x1-line2[m].x1;
	    dy=line1[j].y1-line2[m].y1;
	    dz=line1[j].z1-line2[m].z1;
	    du=sqrt(dx*dx+dy*dy);
	    dl=sqrt(dx*dx+dy*dy+dz*dz);
	    
	    rotz=0;
	    if(dy == 0.0) {rotz=-M_PI/2;} else {rotz = -atan(dx/dy);}
	    
	    rotx=0;
	    if(du == 0.0) {rotx=M_PI/2;}
	    
	    if(du != 0.0) {
	      if(dx>=0.0 && dy>=0.0 && dz> 0.0) {rotx =  atan(dz/du);}
	      if(dx>=0.0 && dy< 0.0 && dz> 0.0) {rotx = -atan(dz/du);}
	      if(dx< 0.0 && dy> 0.0 && dz> 0.0) {rotx =  atan(dz/du);}
	      if(dx< 0.0 && dy<=0.0 && dz> 0.0) {rotx = -atan(dz/du);}
	      if(dx>=0.0 && dy>=0.0 && dz< 0.0) {rotx =  atan(dz/du);}
	      if(dx>=0.0 && dy< 0.0 && dz< 0.0) {rotx = -atan(dz/du);}
	      if(dx< 0.0 && dy> 0.0 && dz< 0.0) {rotx =  atan(dz/du);}
	      if(dx< 0.0 && dy<=0.0 && dz< 0.0) {rotx = -atan(dz/du);}
	    }

  	    
	    fprintf(fp2, "    Transform { translation %7.3f %7.3f %7.3f\n",mx, my, mz);
	    fprintf(fp2, "    children [ Transform {rotation 0 0 1 %7.5f\n",rotz);
	    fprintf(fp2, "    children [ Transform {rotation 1 0 0 %7.5f \n",rotx);	
	    fprintf(fp2, "    children [ Shape { geometry Cylinder { radius %3.2f height %3.2f } \n      appearance USE step%d } ]\n",cubes/4, dl, i);
		fprintf(fp2, "    children [ Shape { geometry Box {size %3.2f %3.2f %3.2f } \n      appearance USE step%d } ]\n\n",cubes/2, dl,cubes/2, i);		
		fprintf(fp2, "    } ] } ] }\n");
		*/
	    /* end of cylinder */ 
	  }
	}
      }
      fprintf(fp2, "\n\n");     
      fprintf(fp2, "# end of time step %d\n\n", i);     
      strcpy(val, "");
      free(line1); free(line2);
    }  /* end of sequence loop */
  
  fprintf(fp2, "# trajectories finished\n");
  
  fclose(fp2); 
  Tcl_Eval(interp, ".text delete 2");
  Tcl_Eval(interp, ".text insert 2 \"Tracks/Detections written to VRML-File: dt.wrl\"");
  Tcl_Eval(interp, "update idletasks");  

  sprintf(val, "...done");
  Tcl_SetVar(interp, "tbuf", val, TCL_GLOBAL_ONLY);
  Tcl_Eval(interp, ".text delete 3");
  Tcl_Eval(interp, ".text insert 3 $tbuf");
  
  return TCL_OK;
}
