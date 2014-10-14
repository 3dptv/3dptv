//////////////////////////////////////////////////////////////////////////////
//
// this code was written by Beat Luthi at IfU, ETH Zürich, Okt 2007
//
// is represents an attempt to have ONE clean non-GUI version of the postPorcessing codes
// that float around in various Borland versions
//
// luethi@ifu.baug.ethz.ch
//
// last update/change: August 2011 Marc Wolf
//
//////////////////////////////////////////////////////////////////////////////

/*

This software links 3D particle positions of consequtivee time steps. 

Copyright (C) 2006 Beat Luthi, Risø, Nat. Lab, Denmark

This program is free software; you can redistribute it and/or modify it under 
the terms of the GNU General Public License v2, as published by the Free 
Software Foundation, provided that the above copyright notices and this 
permission notice appear in all copies of the software and related documentation.
You may charge a fee for the physical act of transferring a copy, and you may at 
your option offer warranty protection in exchange for a fee.

You may not copy, modify, sublicense, or distribute the Program except as 
expressly provided under this License.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR 
A PARTICULAR PURPOSE.

*/


#include "stdafx.h"


TpointList pointList;

FILE * input;
FILE *fpp;
int n;
float e;

static void flushline(FILE * fp);
static void map_slices_to_cycles();
static void readPTVFile(int n, int index);
static void read_scanning_PTVFile(int n, int index); // added by Beat March 2013 for scanning
static void prepare_fast_search();
static void doCubicSplines(bool single,int number);
static void setAllMatrixesToZero(int size);
static void makeAT(int n, int m);
static void makeATA(int n, int m);
static void makeATY(int n, int m,int wh);
static bool solve(int n, int m);
static void writeXUAPFile(int t);
static void followTrajPoint(FILE *fpp, int t,int startPoint);
static void readXUAPFile(int n, bool firstTime);
static void readXUAGFile(int n, bool firstTime);


int main(int argc, char *argv[])
{
	char garb[10];
	char pa[256];
	char name[256];
	int c;
	int deltaFrames,numCycles;

	//begin of read in control parameters
	///////////////////////////////////////////////////////////////////////////////////
	if (argc == 1) {
		//if (NULL == (input = fopen("C:/input.inp","r"))){ 
		//if (NULL == (input = fopen("D:/PTV/version_March_2013_scanning/input_2905.txt","r"))){ 
		if (NULL == (input = fopen("D:\ScanPTV_publish\version_March_2013_scanning/input_41.txt","r"))){ 	
		    cout<< "\ndid not find *.inp file";
	    }
	    else{
		    cout<< "\nautomatically and succesfully opened *.inp file \n";
	    }
	}
	else{
	    if (NULL == (input = fopen(argv[1],"r"))){
		    cout<< "\ndid not find *.inp file";
	    }
	    else{
		    cout<< "\nsuccesfully opened *.inp file \n";
	    }
	}
	//what should be done?
	fscanf(input,"%i",&n); flushline(input); if(n==1){pointList.xuap=true;}			else{pointList.xuap=false;}
	fscanf(input,"%i",&n); flushline(input); if(n==1){pointList.traj_point=true;}		else{pointList.traj_point=false;}
	fscanf(input,"%i",&n); flushline(input); if(n==1){pointList.derivatives=true;}		else{pointList.derivatives=false;}
	fscanf(input,"%i",&n); flushline(input); if(n==1){pointList.pressure=true;}		else{pointList.pressure=false;}
	fscanf(input,"%i",&n); flushline(input); if(n==1){pointList.Hessian=true;}		else{pointList.Hessian=false;}


	//data
    fscanf(input,"%s",pa); flushline(input);sprintf (pointList.path,pa);
    fscanf(input,"%i",&n); flushline(input);pointList.firstSFile               = n; // the code will compute the cycle numbers,i.e. the firstFile and lastFile by itself
	fscanf(input,"%i",&n); flushline(input);pointList.lastSFile                = n; // the code will compute the cycle numbers,i.e. the firstFile and lastFile by itself
	fscanf(input,"%i",&n); flushline(input);pointList.numSlices                = n;

    //fact
	fscanf(input,"%f",&e); flushline(input);pointList.deltaT_between_slice     = e;
	fscanf(input,"%f",&e); flushline(input);pointList.deltaT                   = e; // between scans, or for non-scanning
	fscanf(input,"%f",&e); flushline(input);pointList.viscosity                = e;

	//controls xuap
	fscanf(input,"%i",&n); flushline(input);pointList.PL                       = n;
	fscanf(input,"%i",&n); flushline(input);pointList.minLeftRight             = n;
	fscanf(input,"%f",&e); flushline(input);pointList.tolMaxVel                = e;

	//controls traj_accc
	fscanf(input,"%f",&e); flushline(input);pointList.maxRadius                = e;
	fscanf(input,"%f",&e); flushline(input);pointList.weDiv                    = e; // weighting for 2Q+diva error
	fscanf(input,"%f",&e); flushline(input);pointList.weAcc                    = e; // weighting for acceleration error
	fscanf(input,"%f",&e); flushline(input);pointList.weVel                    = e; // weighting for divu error, addeded by Marc, 14.07.2011
	flushline(input);
	fscanf(input,"%i",&n); flushline(input);pointList.minTrajLength            = n;
	fscanf(input,"%i",&n); flushline(input);pointList.polyConst                = n;
	fscanf(input,"%f",&e); flushline(input);pointList.c1                       = e;
	fscanf(input,"%f",&e); flushline(input);pointList.c2                       = e;
	fscanf(input,"%i",&n); flushline(input);pointList.maxRank                  = n;
	fscanf(input,"%i",&n); flushline(input);pointList.numOfFrames              = n;
	fscanf(input,"%i",&n); flushline(input);pointList.max_grid_X               = n;
	fscanf(input,"%i",&n); flushline(input);pointList.max_grid_Y               = n;
	fscanf(input,"%i",&n); flushline(input);pointList.max_grid_Z               = n;
	fscanf(input,"%i",&n); flushline(input);pointList.max_grid_C               = n;
	flushline(input);
	fscanf(input,"%f",&e); flushline(input);pointList.xminChamber              = e; //added by Markus, 20.07.2009
	fscanf(input,"%f",&e); flushline(input);pointList.xmaxChamber              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.xminChannel              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.xmaxChannel              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.zminChamber              = e; //added by Markus, 20.07.2009
	fscanf(input,"%f",&e); flushline(input);pointList.zmaxChamber              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.zminChannel              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.zmaxChannel              = e;
	fscanf(input,"%f",&e); flushline(input);pointList.yChamberChannel          = e;
	
	//end of read in control parameters
	///////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////
	// begin of #slice > 1 treatment
	if(pointList.numSlices>1){
		deltaFrames=pointList.lastSFile-pointList.firstSFile+1;
		numCycles=int((double)deltaFrames/(double)pointList.numSlices);
		pointList.firstFile=pointList.firstSFile;
		pointList.lastFile=pointList.firstSFile+numCycles-1;
	}
	else{
		//business as usual, like it was before scanning.
		pointList.firstFile=pointList.firstSFile;
		pointList.lastFile=pointList.lastSFile;
	}

	// end of #slice >1 treatment, 
	///////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////
	if(pointList.xuap){
	   pointList.PLh=int((double)pointList.PL/2.);
       pointList.count=0;
       pointList.maxVel=0.;
       pointList.meanVel=0.;
       pointList.meanAcc=0.;

	   if(pointList.numSlices>1){
		   // map slices frame and point id's to cycle frame and point id's
		   map_slices_to_cycles();
	   }

       for (int i=pointList.firstFile;i<pointList.lastFile+1;i++){
		   if(i % ((int)((double)20/(double)pointList.numSlices)+1) == 0){
             cout << "processing file ........."<<i<<"\n";
             cout << "max Vel.................."<<pointList.maxVel<<"\n";
	         cout << "mean Vel................."<<pointList.meanVel<<"\n";
	         cout << "mean Acc................."<<pointList.meanAcc<<"\n\n";
	      }
	      for (int ii=-pointList.PLh;ii<pointList.PLh+1;ii++){
			  if(pointList.numSlices>1){
				  // read in scanned ptv_is files
				  read_scanning_PTVFile(i,ii);
			  }
			  else{ //business as usual, no scanning
                  readPTVFile(i,ii);
			  }
          }
	      doCubicSplines(false,0);
          writeXUAPFile(i);
	   }
	}
    ///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	if(pointList.traj_point){

       pointList.count=0;
	   pointList.count2=0;
       pointList.count3=0;
	   pointList.count4=0;
	   pointList.count5=0;
	   pointList.count6=0;
       pointList.meanDiss=0.;
	   pointList.meanUSq=0.;

	   cout << "\npreparing grid for neighbor search\n";
	   prepare_fast_search();

       for (int i=pointList.firstFile;i<pointList.lastFile+1;i++){
		  if((double)pointList.count3/(double)pointList.count>0){
             cout << "point per sphere.............."<<(double)pointList.count3/(double)pointList.count<<"\n";
		     cout << "% rel. diva < 0.1............."<<100.*(double)pointList.count4/(double)pointList.count<<"\n";
	         cout << "% rel. acc  < 0.2............."<<100.*(double)pointList.count5/(double)pointList.count<<"\n";
			 cout << "% rel. divu < 0.1............."<<100.*(double)pointList.count6/(double)pointList.count<<"\n";
	         cout << "r.m.s. u [m/s]................"<<pow(pointList.meanUSq,0.5)<<"\n";
		     cout << "mean dissipation [m^2/s^3]...."<<pointList.meanDiss<<"\n\n";
		  }
		  cout << "processing file .............."<<i<<"\n";
		   
          c=sprintf (name, pointList.path);
	      c+=sprintf (name+c, "/trajPoint.");
          c+=sprintf (name+c, "%1d", i); 
          fpp = fopen(name,"w");
          followTrajPoint(fpp,i,0);
          fclose (fpp);
      
       }
	}
	///////////////////////////////////////////////////////////////////////////////////
	scanf("Please hit a key  %s", garb);  // to stop console

	return 0;
}

void flushline(FILE * fp)
{
    while(fgetc(fp)!='\n' && !feof(fp));
}

void map_slices_to_cycles()
{
    FILE *fpp;
	int c;
    int numOfPoints;
	int cid,old_cid;
    char name[256];


	old_cid=-1;
    for (int i=pointList.firstSFile;i<pointList.lastSFile+1;i++){

       if(i % 100 == 0){
             cout << "mapping slice file ..."<<i<<"\n";
	   }

	   // determine cycle id
	   cid=(int)( (double)(i-pointList.firstSFile)/(double)pointList.numSlices );
	   if(cid>old_cid){
		   old_cid=cid;
		   pointList.numPoints_per_cycle[cid]=0; //so now it is initiated for cummulative cycle point_id
	   }
  
       c=sprintf (name, pointList.path);
	   c+=sprintf (name+c, "/ptv_is.");
       c+=sprintf (name+c, "%1d", i); 
       
       fpp = fopen(name,"r");
       fscanf (fpp, "%d\0", &numOfPoints);
       
       for (int j=0; j<numOfPoints; j++){
		   // pointList.map_slice_cycle[i-pointList.firstSFile+1][j][0]=cid;//cycle id
		   pointList.map_slice_cycle[i-pointList.firstSFile][j]=j+pointList.numPoints_per_cycle[cid];//cummulated cycle point_id
       }
       fclose (fpp);

	   pointList.numPoints_per_cycle[cid]+=numOfPoints;
    }
}

void readPTVFile(int n, int index)
{
    FILE *fpp;
	int c;
    int numOfPoints;
	int left,right;
    double x,y,z,rmsDist;

    char name[256];

    if(n+index>pointList.firstFile-1 && n+index<pointList.lastFile+1){
       c=sprintf (name, pointList.path);
	   c+=sprintf (name+c, "/ptv_is.");
       c+=sprintf (name+c, "%1d", n+index); 
       
       fpp = fopen(name,"r");
       fscanf (fpp, "%d\0", &numOfPoints);
       pointList.point[index+pointList.PLh][0][0]=numOfPoints;
       for (int i=1; i<numOfPoints+1; i++){ //these lines 218-231 müssen geändert werden
           fscanf (fpp, "%d\0", &left);
           fscanf (fpp, "%d\0", &right);
           fscanf (fpp, "%lf\0", &x);
           fscanf (fpp, "%lf\0", &y);
           fscanf (fpp, "%lf\0", &z);
           rmsDist=0.005;
		   pointList.point[index+pointList.PLh][i][0]=left+1;//;//
		   pointList.point[index+pointList.PLh][i][1]=right+1;//;//

           pointList.point[index+pointList.PLh][i][2]=x*0.001;//;//
           pointList.point[index+pointList.PLh][i][3]=y*0.001;//;//
           pointList.point[index+pointList.PLh][i][4]=z*0.001;//;//
           pointList.point[index+pointList.PLh][i][15]=rmsDist;
      }
       fclose (fpp);
    }
    else{
       pointList.point[index+pointList.PLh][0][0]=0;
    }
}

void read_scanning_PTVFile(int n, int index)
{
    FILE *fpp;
	int c;
    int numOfPoints;
    int fid_left,left,fid_right,right,cid,old_cid,id,cpid,left_pid,right_pid;
    double x,y,z,rmsDist;

    char name[256];

    if(n+index>pointList.firstFile-1 && n+index<pointList.lastFile+1){

        pointList.point[index+pointList.PLh][0][0]=pointList.numPoints_per_cycle[n+index-pointList.firstFile];
		
		//now loop through slices etc
		for (int s=0;s<pointList.numSlices;s++){

		   c=sprintf (name, pointList.path);
	       c+=sprintf (name+c, "/ptv_is.");
           //c+=sprintf (name+c, "%1d", n+index); <- as was before scanning
		   // map n+index to proper slice frame id
		   id = ((n+index)-pointList.firstFile)*pointList.numSlices+s+pointList.firstSFile;
           c+=sprintf (name+c, "%1d", id); // <--and now we have cid mapped to id
           fpp = fopen(name,"r");
           fscanf (fpp, "%d\0", &numOfPoints);
           for (int i=0; i<numOfPoints; i++){ //these lines 218-231 müssen geändert werden
               fscanf (fpp, "%d\0",  &fid_left);
			   fscanf (fpp, "%d\0",  &left);
               fscanf (fpp, "%d\0",  &fid_right);
			   fscanf (fpp, "%d\0",  &right);
               fscanf (fpp, "%lf\0", &x);
               fscanf (fpp, "%lf\0", &y);
               fscanf (fpp, "%lf\0", &z);
               rmsDist=0.005;

			   cid       = n+index;
			   cpid      = 1+pointList.map_slice_cycle[id-pointList.firstFile][i];
			   if(fid_left>-1){
			      left_pid  = pointList.map_slice_cycle[fid_left-pointList.firstSFile ][left ];
			   }
			   else{
                  left_pid  = left; 
			   }
			   if(fid_right>-1){
			      right_pid = pointList.map_slice_cycle[fid_right-pointList.firstSFile][right];
			   }
			   else{
                  right_pid = right;
			   }

			   pointList.point[index+pointList.PLh][cpid][0]=left_pid+1;
		       pointList.point[index+pointList.PLh][cpid][1]=right_pid+1;

               pointList.point[index+pointList.PLh][cpid][2]=x*0.001;
               pointList.point[index+pointList.PLh][cpid][3]=y*0.001;
               pointList.point[index+pointList.PLh][cpid][4]=z*0.001;
               pointList.point[index+pointList.PLh][cpid][15]=rmsDist;
			   pointList.point[index+pointList.PLh][cpid][16]=pointList.deltaT_between_slice*((double)s-(double)pointList.numSlices/2); //delta t relative to time of middle slice 
		   
		   }
           fclose (fpp);
		}
    }
    else{
       pointList.point[index+pointList.PLh][0][0]=0;
    }
}

void doCubicSplines(bool single,int number)
{

   //int pointList.PLh=int((double)pointList.PL/2.);
   int nP=pointList.point[pointList.PLh][0][0];
   int ind[200]; //Marc & Beat: 27.04.2011 changed size from 21 to 200
   //double tolerance=0.15;//StrToFloat(paramForm->toleranceEdit->Text);
   double velocity;

   double weight,time;

   int start,end;
   if(!single){
      start=1;
      end=nP; 
   }
   else{
      start=number;
      end=number;
   }

   for(int i=start;i<end+1;i++){ 
      pointList.point[pointList.PLh][i][14]=0; //can be cubic splined
      int maxIndex=pointList.PLh;
      int minIndex=pointList.PLh;
      int index=pointList.PLh;
      int badCounter=0;
      ind[index]=i; 
      bool ok=true;

	  while(index>0 && ok){
          if(pointList.point[index][ind[index]][0]>0 && pointList.point[index][0][0]>0){ 
             ind[index-1]=pointList.point[index][ind[index]][0]; 
             index--;
             minIndex=index;
          }
          else{
             ok=false;
          }
      }
      index=pointList.PLh;
      ind[index]=i;
      ok=true;
      while(index<2*pointList.PLh && ok){
          if(pointList.point[index][ind[index]][1]>0 && pointList.point[index][0][0]>0){ 
             ind[index+1]=pointList.point[index][ind[index]][1]; 
             index++;
             maxIndex=index;
          }
          else{
             ok=false;
          }
      }



      //first do for x and u, then do for a
      if(maxIndex-minIndex>2+badCounter && maxIndex>pointList.PLh-1+pointList.minLeftRight && minIndex<pointList.PLh+1-pointList.minLeftRight){ 
//    if(maxIndex-minIndex>2+badCounter && maxIndex>9+minLength && minIndex<11-minLength){ 
      
	  //if(maxIndex-minIndex>2+badCounter ){ //ok (minIndex<10 && maxIndex>10){//
          pointList.point[pointList.PLh][i][14]=1;
          //x-Component
          setAllMatrixesToZero(4);
          for(int t=minIndex-pointList.PLh;t<maxIndex-pointList.PLh+1;t++){ //anpassen so dass t nicht mehr genau Zeit, sondern 'nur' master loop index ist
              weight     = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][15];//da muss man anpassen, ind,2
              weight     = 1.-1./(1.+exp(-300.*(weight-0.015)));
			  weight   = 1.; //Beat March 2009

			  //this is the only chamge due to scanning, works also if numSlice=0
			  if(pointList.numSlices>1){
                  time=(double)t*pointList.deltaT + pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][16];
			  }
			  else{
			      time=(double)t*pointList.deltaT;
			  }
			  //end of only chamge due to scanning

              pointList.A[t+pointList.PLh][0] = 1.*weight;
              pointList.A[t+pointList.PLh][1] = time*weight; 
              pointList.A[t+pointList.PLh][2] = pow(time,2.)*weight; // t is integer from e.g. -11 to 11
              pointList.A[t+pointList.PLh][3] = pow(time,3.)*weight;
              pointList.y[0][t+pointList.PLh] = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][2]*weight;
			  
          }
          makeAT(pointList.PL,4);
          makeATA(pointList.PL,4);
          makeATY(pointList.PL,4,0);
          solve(pointList.PL,4);

          pointList.point[pointList.PLh][i][5]=pointList.X[0];//filtered x position
          pointList.point[pointList.PLh][i][8]=pointList.X[1];//filtered velocity, derivative from filtered x
          pointList.point[pointList.PLh][i][11]=2.*pointList.X[2];// filtered acc, derivative from velocity

          //y-Component
          setAllMatrixesToZero(4);
          for(int t=minIndex-pointList.PLh;t<maxIndex-pointList.PLh+1;t++){
              weight     = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][15];
              weight     = 1.-1./(1.+exp(-300.*(weight-0.015)));
			  weight   = 1.; //Beat March 2009

			  //this is the only chamge due to scanning, works also if numSlice=0
			  if(pointList.numSlices>1){
                  time=(double)t*pointList.deltaT+pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][16];
			  }
			  else{
			      time=(double)t*pointList.deltaT;
			  }
			  //end of only chamge due to scanning

              pointList.A[t+pointList.PLh][0] = 1.*weight;
              pointList.A[t+pointList.PLh][1] = time*weight;
              pointList.A[t+pointList.PLh][2] = pow(time,2.)*weight;
              pointList.A[t+pointList.PLh][3] = pow(time,3.)*weight;
              pointList.y[0][t+pointList.PLh] = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][3]*weight;
          }
          makeAT(pointList.PL,4);
          makeATA(pointList.PL,4);
          makeATY(pointList.PL,4,0);
          solve(pointList.PL,4);

          pointList.point[pointList.PLh][i][6]=pointList.X[0]; //pointList.point[pointList.PLh][ind[pointList.PLh]][3];//
          pointList.point[pointList.PLh][i][9]=pointList.X[1]; //(1./(2.*pointList.pointList.deltaT))*(pointList.point[11][ind[11]][3]-pointList.point[9][ind[9]][3]);//
          pointList.point[pointList.PLh][i][12]=2.*pointList.X[2]; //(1./(pointList.pointList.deltaT*pointList.pointList.deltaT))*(pointList.point[11][ind[11]][3]-2.*pointList.point[pointList.PLh][ind[pointList.PLh]][3]+pointList.point[9][ind[9]][3]);//
          //z-Component
          setAllMatrixesToZero(4);
          for(int t=minIndex-pointList.PLh;t<maxIndex-pointList.PLh+1;t++){
              weight     = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][15];
              weight     = 1.-1./(1.+exp(-300.*(weight-0.015)));
			  weight   = 1.; //Beat March 2009

			  //this is the only chamge due to scanning, works also if numSlice=0
			  if(pointList.numSlices>1){
                  time=(double)t*pointList.deltaT+pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][16];
			  }
			  else{
			      time=(double)t*pointList.deltaT;
			  }
			  //end of only chamge due to scanning

              pointList.A[t+pointList.PLh][0] = 1.*weight;
              pointList.A[t+pointList.PLh][1] = time*weight;
              pointList.A[t+pointList.PLh][2] = pow(time,2.)*weight;
              pointList.A[t+pointList.PLh][3] = pow(time,3.)*weight;
              pointList.y[0][t+pointList.PLh] = pointList.point[t+pointList.PLh][ind[t+pointList.PLh]][4]*weight;
          }
          makeAT(pointList.PL,4);
          makeATA(pointList.PL,4);
          makeATY(pointList.PL,4,0);
          solve(pointList.PL,4);
          
          pointList.point[pointList.PLh][i][7]=pointList.X[0]; //pointList.point[pointList.PLh][ind[pointList.PLh]][4];//
          pointList.point[pointList.PLh][i][10]=pointList.X[1];//(1./(2.*pointList.pointList.deltaT))*(pointList.point[11][ind[11]][4]-pointList.point[9][ind[9]][4]);//
          pointList.point[pointList.PLh][i][13]=2.*pointList.X[2]; //(1./(pointList.pointList.deltaT*pointList.pointList.deltaT))*(pointList.point[11][ind[11]][4]-2.*pointList.point[10][ind[10]][4]+pointList.point[9][ind[9]][4]);//
          //max break!
          velocity=pow(pow(pointList.point[10][i][8],2.)+pow(pointList.point[10][i][9],2.)+pow(pointList.point[10][i][10],2.),0.5);
          if(velocity>pointList.tolMaxVel){
             pointList.point[pointList.PLh][i][14]=0;
          }
      }
   }
}

void setAllMatrixesToZero(int size)
{

    for(int i=0;i<500;i++){
       if(i<size){
          pointList.X[i]=0.;
          pointList.ATY[i]=0.;
          pointList.BTY[i]=0.;
          pointList.CTY[i]=0.;
       }
       pointList.Y[i]=0.;
       
       pointList.YuB[i]=0.;
       pointList.YvB[i]=0.;
       pointList.YwB[i]=0.;
      
       pointList.Yaz[i]=0.;
       pointList.Yay[i]=0.;
       pointList.Yax[i]=0.;

       for(int j=0;j<size;j++){
          pointList.A[i][j]=0.;
          pointList.AT[j][i]=0.;
          if(i<size){pointList.ATA[i][j]=0.;}
          pointList.B[i][j]=0.;
          pointList.BT[j][i]=0.;
          if(i<size){pointList.BTB[i][j]=0.;}
          pointList.C[i][j]=0.;
          pointList.CT[j][i]=0.;
          if(i<size){pointList.CTC[i][j]=0.;}
       }
    }
    
}
void makeAT(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
           pointList.AT[i][j]=pointList.A[j][i];
        }
     }
}
void makeATA(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
           pointList.ATA[i][j]=0.;
           for(int k=0;k<n;k++){
              pointList.ATA[i][j]=pointList.ATA[i][j]+pointList.AT[i][k]*pointList.A[k][j];
           }
        }
     }
}
void makeATY(int n, int m,int wh)
{
     for(int i=0;i<m;i++){
           pointList.ATY[i]=0.;
           for(int k=0;k<n;k++){
               pointList.ATY[i]=pointList.ATY[i]+pointList.AT[i][k]*pointList.y[wh][k];
           }
     }
}
bool solve(int n, int m)
{
    double faktor;
    bool ok=true;

    for(int i=1;i<m;i++){
       for(int j=i;j<m;j++){
          if(fabs(pointList.ATA[j][i-1])>0.){
             faktor=pointList.ATA[i-1][i-1]/pointList.ATA[j][i-1];
             for(int k=0;k<m;k++){
                pointList.ATA[j][k]=pointList.ATA[i-1][k]-faktor*pointList.ATA[j][k];
             }
             pointList.ATY[j]=pointList.ATY[i-1]-faktor*pointList.ATY[j];
          }
       }
    }
    for(int i=m-1;i>-1;i--){
       for(int j=i+1;j<m;j++){
          pointList.ATY[i]=pointList.ATY[i]-pointList.ATA[i][j]*pointList.X[j];
       }
       if(fabs(pointList.ATA[i][i])>0.){
          pointList.X[i]=pointList.ATY[i]/pointList.ATA[i][i];
       }
       else{
          ok=false;
       }
    }
    return ok;
}

////////////////////von Beat June 2011
void makeBT(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
           pointList.BT[i][j]=pointList.B[j][i];
        }
     }
}
void makeBTB(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
           pointList.BTB[i][j]=0.;
           for(int k=0;k<n;k++){
              pointList.BTB[i][j]=pointList.BTB[i][j]+pointList.BT[i][k]*pointList.B[k][j];
           }
        }
     }
}

void makeBTY(int n, int m,int wh)
{
     for(int i=0;i<m;i++){
           pointList.BTY[i]=0.;
           for(int k=0;k<n;k++){
			   switch (wh){
				   case 1:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YuB[k];
				       break;
				   case 2:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YvB[k];
				       break;
				   case 3:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YwB[k];
				       break;
			   }
           }
     }
}

void makeBTYa(int n, int m,int wh)
{
     for(int i=0;i<m;i++){
           pointList.BTY[i]=0.;
           for(int k=0;k<n;k++){
			   switch (wh){
				   case 1:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YaxB[k];
				       break;
				   case 2:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YayB[k];
				       break;
				   case 3:
                       pointList.BTY[i]=pointList.BTY[i]+pointList.BT[i][k]*pointList.YazB[k];
				       break;
			   }
           }
     }
}
bool solveB(int n, int m)
{
    double faktor;
    bool ok=true;

    for(int i=1;i<m;i++){
       for(int j=i;j<m;j++){
          if(fabs(pointList.BTB[j][i-1])>0.){
             faktor=pointList.BTB[i-1][i-1]/pointList.BTB[j][i-1];
             for(int k=0;k<m;k++){
                pointList.BTB[j][k]=pointList.BTB[i-1][k]-faktor*pointList.BTB[j][k];
             }
             pointList.BTY[j]=pointList.BTY[i-1]-faktor*pointList.BTY[j];
          }
       }
    }
    for(int i=m-1;i>-1;i--){
       for(int j=i+1;j<m;j++){
          pointList.BTY[i]=pointList.BTY[i]-pointList.BTB[i][j]*pointList.X[j];
       }
       if(fabs(pointList.BTB[i][i])>0.){
          pointList.X[i]=pointList.BTY[i]/pointList.BTB[i][i];
       }
       else{
          ok=false;
       }
    }
	return ok;
}

void makeCT(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
           pointList.CT[i][j]=pointList.C[j][i];
        }
     }
}
void makeCTC(int n, int m)
{
     for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
           pointList.CTC[i][j]=0.;
           for(int k=0;k<n;k++){
              pointList.CTC[i][j]=pointList.CTC[i][j]+pointList.CT[i][k]*pointList.C[k][j];
           }
        }
     }
}
void makeCTY(int n, int m,int wh)
{
     for(int i=0;i<m;i++){
           pointList.CTY[i]=0.;
           for(int k=0;k<n;k++){
               pointList.CTY[i]=pointList.CTY[i]+pointList.CT[i][k]*pointList.yC[wh][k];
           }
     }
}
bool solveC(int n, int m)
{
    double faktor;
    bool ok=true;

    for(int i=1;i<m;i++){
       for(int j=i;j<m;j++){
          if(fabs(pointList.CTC[j][i-1])>0.){
             faktor=pointList.CTC[i-1][i-1]/pointList.CTC[j][i-1];
             for(int k=0;k<m;k++){
                pointList.CTC[j][k]=pointList.CTC[i-1][k]-faktor*pointList.CTC[j][k];
             }
             pointList.CTY[j]=pointList.CTY[i-1]-faktor*pointList.CTY[j];
          }
       }
    }
    for(int i=m-1;i>-1;i--){
       for(int j=i+1;j<m;j++){
          pointList.CTY[i]=pointList.CTY[i]-pointList.CTC[i][j]*pointList.X[j];
       }
       if(fabs(pointList.CTC[i][i])>0.){
          pointList.X[i]=pointList.CTY[i]/pointList.CTC[i][i];
       }
       else{
          ok=false;
       }
    }
    return ok;
}
void writeXUAPFile(int t)
{

    FILE *fpp;
    char name[256];
    int c;

    c=sprintf (name, pointList.path);
    c+=sprintf (name+c, "/xuap.");
    c+=sprintf (name+c, "%1d", t);

    fpp = fopen(name,"w");
    
    for(int i=1;i<pointList.point[pointList.PLh][0][0];i++){
       if(pointList.point[pointList.PLh][i][14]>0){
           pointList.count++;
           double vel=pow( pow(pointList.point[pointList.PLh][i][8],2.)
                          +pow(pointList.point[pointList.PLh][i][9],2.)
                          +pow(pointList.point[pointList.PLh][i][10],2.),0.5);
           double acc=pow( pow(pointList.point[pointList.PLh][i][11],2.)
                          +pow(pointList.point[pointList.PLh][i][12],2.)
                          +pow(pointList.point[pointList.PLh][i][13],2.),0.5);
           pointList.meanVel=(pointList.meanVel*(double)(pointList.count-1)+vel)/(double)pointList.count;
           pointList.meanAcc=(pointList.meanAcc*(double)(pointList.count-1)+acc)/(double)pointList.count;
           if(vel>pointList.maxVel){
              pointList.maxVel=vel;
           }
        }
        for(int j=0;j<14;j++){
			if(j<5 || pointList.point[pointList.PLh][i][14]>0){
                fprintf(fpp, "%lf\t", pointList.point[pointList.PLh][i][j]);
			}
			else{
                fprintf(fpp, "%lf\t", 0.);
			}
        }
        fprintf(fpp, "%lf\n", pointList.point[pointList.PLh][i][14]);
    }
    fclose (fpp);
}

void followTrajPoint(FILE *fpp, int t,int startPoint)
{
     
	 int ind_X,ind_Y,ind_Z,ind_C,ind_count;
	 short ind_list[1000];
	
	 int pCounterA,pCounterB,pCounterC,numInTraj;
     int startT, startP;
     double dist,dx,dy,dz;
     double centerX,centerY,centerZ;
     double Liu[5],Liv[5],Liw[5],Liax[4],Liay[4],Liaz[4];
     double ux,uy,uz,vx,vy,vz,wx,wy,wz;
     double dix,diy,diz,absDi,Dx,Dy,Dz,lx,ly,lz,cx,cy,cz,refx,refy,refz;
     double w1,w2,w3,s11,s12,s13,s22,s23,s33,ww1,ww2,ww3,wwsij;
	 double s111,s222,s333,s112,s113,s221,s223,s331,s332,s123;
     double sijsjkski,wsq,twosijsij,R,Q,div,ref,diss,USq;
     int time;
	 double u[3];
	 double a[3];
	 double ref_diva,diva,reldiva,quality,reldivu;
    
     double minDistA[500];
     int minDistAIndex[500];
     double minDistB[500];
     int minDistBIndex[500];
     double minDistC[500];
     int minDistCIndex[500];
     
     double um,up,vm,vp,wm,wp;
     bool okc,contin;

     int rank;

     int start;
     int end;
     int minCounter;
	 int counter_f;

	 double avU[3];
	 double avA[3];

     bool ok;
     startT=t;

     if(t==pointList.firstFile){
		 cout << "\nreading initial xuap, may take some time....\n\n";
         readXUAPFile(t,true);
     }
     else{
        readXUAPFile(t,false);
     }
     

     start=1;
     end=(int)(pointList.point[2][0][0]+0.5); // 1. field tells number of rows
     

     int n;
     for(int nn=start;nn<end;nn++){
         time=2; // is set to 2 to calculate local acc
		 /*if(end>5000){
			 if(nn % 2000 == 0){
			     cout << "processing point ............."<<nn<<"\n";
			 }
		 }*/
//if((double)pointList.count3/(double)pointList.count>0){
//    cout << "point per sphere.............."<<(double)pointList.count3/(double)pointList.count<<"\n";
//}

         if(pointList.point[2][nn][11]>0. && !(pointList.occ[t-pointList.firstFile][nn]) ){
            startP=nn;
            ok=true;
            numInTraj=0;
            pointList.noDeriv=0;
            n=nn;
            while(ok){
				pointList.occ[t+time-2-pointList.firstFile][n]=true;
               //interpolieren und rausschreiben mit t,n (Zeit und Startpunkt)
               //%Da soll jetzt duidxj linear interpoliert werden
               //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               //%die nächsten Punkte zu Punkt x,y,z, finden
               pointList.count++;
               setAllMatrixesToZero(4);
			   contin=true;

			      if(pointList.derivatives){ 
				
                              
               centerX=pointList.point[time][n][2]; // filtered x,y,z
               centerY=pointList.point[time][n][3];
               centerZ=pointList.point[time][n][4];


               for(int i=0;i<500;i++){
                  minDistA[i]=1000; // initialization for the search radius, max. 1000 mm
                  minDistB[i]=1000;
                  minDistC[i]=1000;
				  minDistAIndex[i]=0;
				  minDistBIndex[i]=0;
				  minDistCIndex[i]=0;
               }
               /// Beat March 2013, these 3 for loops are VERY SLOW for large number of particles, i.e. for scanning
               //AAAAAAAAAAAA time step t-1
               // the following bit has been replaced by fast_search 
               /*for(int i=1;i<pointList.point[time-1][0][0]+1;i++){
			      dist=pow(pow(pointList.point[time-1][i][2]-centerX,2.)+pow(pointList.point[time-1][i][3]-centerY,2.)+pow(pointList.point[time-1][i][4]-centerZ,2.),0.5); // distance between measurement point and all points in the xuap file
				  if(dist<minDistA[pointList.maxRank] && pointList.point[time-1][i][11]>0.){ // check if cubic spline successful
                     rank=pointList.maxRank; // whole paragraph: sorting the points according to their distance to measurement point: k=0 closest point, k=500 last point
                     for(int k=pointList.maxRank;k>-1;k--){
                        if(dist<minDistA[k]){
                           rank=k;							
                        }
                     }
                     for(int l=pointList.maxRank;l>rank;l--){
                        minDistA[l]=minDistA[l-1];
                        minDistAIndex[l]=minDistAIndex[l-1];
                     }
                     minDistA[rank]=dist;
                     minDistAIndex[rank]=i;
                  }
			   }
			   for(int i=0;i<500;i++){
                  minDistA[i]=1000; // initialization for the search radius, max. 1000 mm
				  minDistAIndex[i]=0;
               }*/
			   //end of replaced stuff

			   //Beat March 2013 find i not as loop thorugh everything, but only through 3 x 3 surounding grid cells
			   ind_X=(int)((double)(centerX-pointList.minX)/pointList.dh_X);
			   ind_Y=(int)((double)(centerY-pointList.minY)/pointList.dh_Y);
			   ind_Z=(int)((double)(centerZ-pointList.minZ)/pointList.dh_Z);
			   if(ind_X<0 ){ind_X=0 ;}
			   if(ind_Y<0 ){ind_Y=0 ;}
			   if(ind_Z<0 ){ind_Z=0 ;}
			   if(ind_X>pointList.max_grid_X-1){ind_X=pointList.max_grid_X-1;}
			   if(ind_Y>pointList.max_grid_Y-1){ind_Y=pointList.max_grid_Y-1;}
			   if(ind_Z>pointList.max_grid_Z-1){ind_Z=pointList.max_grid_Z-1;}

			   ind_count=0;
			   for (int ii=-1;ii<2;ii++){
				   for (int jj=-1;jj<2;jj++){
					   for (int kk=-1;kk<2;kk++){
						   if(ind_X+ii>=0 && ind_Y+jj>=0 && ind_Z+kk>=0 && ind_X+ii<pointList.max_grid_X && ind_Y+jj<pointList.max_grid_Y && ind_Z+kk<pointList.max_grid_Z){
			                  ind_C=pointList.fast_search[time-1][ind_X+ii][ind_Y+jj][ind_Z+kk][pointList.max_grid_C];
						      if(ind_C>0){
							     for(int ind=0;ind<ind_C;ind++){
								     ind_list[ind_count]=pointList.fast_search[time-1][ind_X+ii][ind_Y+jj][ind_Z+kk][ind];
								     ind_count++;
					     	     }
							  }
						   }
					   }
				   }
			   }
			   for(int ind=0;ind<ind_count;ind++){						   
                  int i=ind_list[ind];
				  //end of i replacement

			      dist=pow(pow(pointList.point[time-1][i][2]-centerX,2.)+pow(pointList.point[time-1][i][3]-centerY,2.)+pow(pointList.point[time-1][i][4]-centerZ,2.),0.5); // distance between measurement point and all points in the xuap file
				  if(dist<minDistA[pointList.maxRank] && pointList.point[time-1][i][11]>0.){ // check if cubic spline successful
                     rank=pointList.maxRank; // whole paragraph: sorting the points according to their distance to measurement point: k=0 closest point, k=500 last point
                     for(int k=pointList.maxRank;k>-1;k--){
                        if(dist<minDistA[k]){
                           rank=k;							
                        }
                     }
                     for(int l=pointList.maxRank;l>rank;l--){
                        minDistA[l]=minDistA[l-1];
                        minDistAIndex[l]=minDistAIndex[l-1];
                     }
                     minDistA[rank]=dist;
                     minDistAIndex[rank]=i;
                  }
               }
               //BBBBBBBBBBBBBBB time step t
               // has been replaced by fast_search for(int i=1;i<pointList.point[time][0][0]+1;i++){
                  
               //Beat March 2013 find i not as loop thorugh everything, but only through 3 x 3 surounding grid cells
			   ind_count=0;
			   for (int ii=-1;ii<2;ii++){
				   for (int jj=-1;jj<2;jj++){
					   for (int kk=-1;kk<2;kk++){
						   if(ind_X+ii>=0 && ind_Y+jj>=0 && ind_Z+kk>=0 && ind_X+ii<pointList.max_grid_X && ind_Y+jj<pointList.max_grid_Y && ind_Z+kk<pointList.max_grid_Z){
			                  ind_C=pointList.fast_search[time][ind_X+ii][ind_Y+jj][ind_Z+kk][pointList.max_grid_C];
						      if(ind_C>0){
							     for(int ind=0;ind<ind_C;ind++){
								     ind_list[ind_count]=pointList.fast_search[time][ind_X+ii][ind_Y+jj][ind_Z+kk][ind];
								     ind_count++;
					     	     }
							  }
						   }
					   }
				   }
			   }
			   for(int ind=0;ind<ind_count;ind++){						   
                  int i=ind_list[ind];
				  //end of i replacement
			   
			      dist=pow(pow(pointList.point[time][i][2]-centerX,2.)+pow(pointList.point[time][i][3]-centerY,2.)+pow(pointList.point[time][i][4]-centerZ,2.),0.5);
                  if(dist<minDistB[pointList.maxRank] && pointList.point[time][i][11]>0.){
                     rank=pointList.maxRank;
                     for(int k=pointList.maxRank;k>-1;k--){
                        if(dist<minDistB[k]){
                           rank=k;
                        }
                     }
                     for(int l=pointList.maxRank;l>rank;l--){
                        minDistB[l]=minDistB[l-1];
                        minDistBIndex[l]=minDistBIndex[l-1];
                     }
                     minDistB[rank]=dist;
                     minDistBIndex[rank]=i;
                  }
               }
               //CCCCCCCCCCCCCCCCCCCCCCCCCC time step t+1
               // has been replaced by fast_search for(int i=1;i<pointList.point[time+1][0][0]+1;i++){

			   //Beat March 2013 find i not as loop thorugh everything, but only through 3 x 3 surounding grid cells
			   ind_count=0;
			   for (int ii=-1;ii<2;ii++){
				   for (int jj=-1;jj<2;jj++){
					   for (int kk=-1;kk<2;kk++){
						   if(ind_X+ii>=0 && ind_Y+jj>=0 && ind_Z+kk>=0 && ind_X+ii<pointList.max_grid_X && ind_Y+jj<pointList.max_grid_Y && ind_Z+kk<pointList.max_grid_Z){
			                  ind_C=pointList.fast_search[time+1][ind_X+ii][ind_Y+jj][ind_Z+kk][pointList.max_grid_C];
						      if(ind_C>0){
							     for(int ind=0;ind<ind_C;ind++){
								     ind_list[ind_count]=pointList.fast_search[time+1][ind_X+ii][ind_Y+jj][ind_Z+kk][ind];
								     ind_count++;
					     	     }
							  }
						   }
					   }
				   }
			   }
			   for(int ind=0;ind<ind_count;ind++){						   
                  int i=ind_list[ind];
				  //end of i replacement 

                  dist=pow(pow(pointList.point[time+1][i][2]-centerX,2.)+pow(pointList.point[time+1][i][3]-centerY,2.)+pow(pointList.point[time+1][i][4]-centerZ,2.),0.5);
                  if(dist<minDistC[pointList.maxRank] && pointList.point[time+1][i][11]>0.){
                     rank=pointList.maxRank;
                     for(int k=pointList.maxRank;k>-1;k--){
                        if(dist<minDistC[k]){
                           rank=k;
                        }
                     }
                     for(int l=pointList.maxRank;l>rank;l--){
                        minDistC[l]=minDistC[l-1];
                        minDistCIndex[l]=minDistCIndex[l-1];
                     }
                     minDistC[rank]=dist;
                     minDistCIndex[rank]=i;
                  }
               }
               
               pCounterA=0;
               pCounterB=0;
               pCounterC=0;
               
               int i;
                // Abfüllen der überbestimmten linearen Gleichungssyteme (3 Zeitschritte)
               for(int pointInd=0;(pointInd<pointList.maxRank) && (minDistA[pointInd]<pointList.maxRadius);pointInd++){
                  i=minDistAIndex[pointInd]; // for-loop through all points with radius smaller than maxRadius
                  if(pointList.point[time-1][i][11]>0.){ // Zeitschritt i-1
                     dist=pow(pow(pointList.point[time-1][i][2]-centerX,2.)+pow(pointList.point[time-1][i][3]-centerY,2.)+pow(pointList.point[time-1][i][4]-centerZ,2.),0.5); //distance (2. time)
                     dx=pointList.point[time-1][i][2]-centerX; //deviation from the meas. point
                     dy=pointList.point[time-1][i][3]-centerY;
                     dz=pointList.point[time-1][i][4]-centerZ;
                     pointList.A[pCounterA][0]=1.; // auffüllen A Matrix linearer Ansatz
                     pointList.A[pCounterA][1]=dx;
                     pointList.A[pCounterA][2]=dy;
                     pointList.A[pCounterA][3]=dz;
                        
                     pointList.y[1][pCounterA]=pointList.point[time-1][i][5]; // drei Komponenten Geschwindigkeit
                     pointList.y[2][pCounterA]=pointList.point[time-1][i][6];
                     pointList.y[3][pCounterA]=pointList.point[time-1][i][7];
                     pCounterA++;
                  }
               }
               for(int pointInd=0;(pointInd<pointList.maxRank) && (minDistB[pointInd]<pointList.maxRadius);pointInd++){
                  i=minDistBIndex[pointInd];
                  if(pointList.point[time][i][11]>0.){ // Zeitschritt i
                     dist=pow(pow(pointList.point[time][i][2]-centerX,2.)+pow(pointList.point[time][i][3]-centerY,2.)+pow(pointList.point[time][i][4]-centerZ,2.),0.5);
                     dx=pointList.point[time][i][2]-centerX;
                     dy=pointList.point[time][i][3]-centerY;
                     dz=pointList.point[time][i][4]-centerZ;
                     pointList.B[pCounterB][0]=1.;
                     pointList.B[pCounterB][1]=dx;
                     pointList.B[pCounterB][2]=dy;
                     pointList.B[pCounterB][3]=dz;
					 /*pointList.B[pCounterB][4]=dx*dx;
					 pointList.B[pCounterB][5]=dy*dy;
					 pointList.B[pCounterB][6]=dz*dz;
					 pointList.B[pCounterB][7]=dx*dy;
					 pointList.B[pCounterB][8]=dx*dz;
					 pointList.B[pCounterB][9]=dy*dz;*/
                     
                     pointList.YuB[pCounterB]=pointList.point[time][i][5];
                     pointList.YvB[pCounterB]=pointList.point[time][i][6];
                     pointList.YwB[pCounterB]=pointList.point[time][i][7];

					 pointList.YaxB[pCounterB]=pointList.point[time][i][8]; // acceleration
                     pointList.YayB[pCounterB]=pointList.point[time][i][9];
                     pointList.YazB[pCounterB]=pointList.point[time][i][10];
                     pCounterB++;
                  }
               }
               for(int pointInd=0;(pointInd<pointList.maxRank) && (minDistC[pointInd]<pointList.maxRadius);pointInd++){
                  i=minDistCIndex[pointInd];
                  if(pointList.point[time+1][i][11]>0.){
                     dist=pow(pow(pointList.point[time+1][i][2]-centerX,2.)+pow(pointList.point[time+1][i][3]-centerY,2.)+pow(pointList.point[time+1][i][4]-centerZ,2.),0.5);
                     dx=pointList.point[time+1][i][2]-centerX;
                     dy=pointList.point[time+1][i][3]-centerY;
                     dz=pointList.point[time+1][i][4]-centerZ;
                     pointList.C[pCounterC][0]=1.;
                     pointList.C[pCounterC][1]=dx;
                     pointList.C[pCounterC][2]=dy;
                     pointList.C[pCounterC][3]=dz;
                     
                     pointList.yC[1][pCounterC]=pointList.point[time+1][i][5];
                     pointList.yC[2][pCounterC]=pointList.point[time+1][i][6];
                     pointList.yC[3][pCounterC]=pointList.point[time+1][i][7];
                     pCounterC++;
                  }
               }
               //end loop through maxRank

               pointList.count3=pointList.count3+pCounterB;
               minCounter=3;
               
               if(pCounterA>minCounter && pCounterB>minCounter && pCounterC>minCounter){ // %jetzt wird endlich Punkt1 interpoliert
                  //%correct x,y,z with center of interpolation!
                                   
                  makeBT(pCounterB,4);  // Gegenwart: räumliche Ableitungen für Geschw. und Beschl.
                  makeBTB(pCounterB,4);
                  makeBTY(pCounterB,4,1);
                  solveB(pCounterB,4);
                  Liu[0]=pointList.point[time][n][5];
                  Liu[1]=pointList.X[1];
                  Liu[2]=pointList.X[2];
                  Liu[3]=pointList.X[3];                

                  makeBT(pCounterB,4);
                  makeBTB(pCounterB,4);
                  makeBTY(pCounterB,4,2);
                  solveB(pCounterB,4);
                  Liv[0]=pointList.point[time][n][6];
                  Liv[1]=pointList.X[1];
                  Liv[2]=pointList.X[2];
                  Liv[3]=pointList.X[3];                 

                  makeBT(pCounterB,4);
                  makeBTB(pCounterB,4);
                  makeBTY(pCounterB,4,3);
                  solveB(pCounterB,4);
                  Liw[0]=pointList.point[time][n][7]; // W-Komponente direkt von der Trajektorie
                  Liw[1]=pointList.X[1]; // Ableitungen, Konstanten C1, C2 und C3
                  Liw[2]=pointList.X[2];
                  Liw[3]=pointList.X[3]; 

				  double trace=Liu[1]+Liv[2]+Liw[3];

				  makeBT(pCounterB,4);
                  makeBTB(pCounterB,4);
                  makeBTYa(pCounterB,4,1);
                  solveB(pCounterB,4);
                  Liax[0]=pointList.point[time][n][8];
                  Liax[1]=pointList.X[1];
                  Liax[2]=pointList.X[2];
                  Liax[3]=pointList.X[3];                

                  makeBT(pCounterB,4);
                  makeBTB(pCounterB,4);
                  makeBTYa(pCounterB,4,2);
                  solveB(pCounterB,4);
                  Liay[0]=pointList.point[time][n][9];
                  Liay[1]=pointList.X[1];
                  Liay[2]=pointList.X[2];
                  Liay[3]=pointList.X[3];                 

                  makeBT(pCounterB,4);
                  makeBTB(pCounterB,4);
                  makeBTYa(pCounterB,4,3);
                  solveB(pCounterB,4);
                  Liaz[0]=pointList.point[time][n][10];
                  Liaz[1]=pointList.X[1];
                  Liaz[2]=pointList.X[2];
                  Liaz[3]=pointList.X[3];                
				  //end of old linear Ansatz stuff

                  //this is for du/dt                        
                  makeAT(pCounterA,4);
                  makeATA(pCounterA,4);
                  makeATY(pCounterA,4,1);
                  okc=solve(pCounterA,4);
				  if(!okc){contin=false;}
                  um=pointList.X[0]; // Interpolierter Geschwindigkeit auf centerX, um = u minus

                  makeCT(pCounterC,4);
                  makeCTC(pCounterC,4);
                  makeCTY(pCounterC,4,1);
                  okc=solveC(pCounterC,4);
				  if(!okc){contin=false;}
                  up=pointList.X[0]; // up = u plus
                  
				  Liu[4]=1./(pointList.deltaT)*(0.5*up-0.5*um); // central difference

				  //this is for dv/dt
                  makeAT(pCounterA,4);
                  makeATA(pCounterA,4);
                  makeATY(pCounterA,4,2);
                  okc=solve(pCounterA,4);
				  if(!okc){contin=false;}
                  vm=pointList.X[0];

                  makeCT(pCounterC,4);
                  makeCTC(pCounterC,4);
                  makeCTY(pCounterC,4,2);
                  okc=solveC(pCounterC,4);
				  if(!okc){contin=false;}
                  vp=pointList.X[0];
                  
				  Liv[4]=1./(pointList.deltaT)*(0.5*vp-0.5*vm);
                               
                  //this is for dw/dt
				  makeAT(pCounterA,4);
                  makeATA(pCounterA,4);
                  makeATY(pCounterA,4,3);
                  okc=solve(pCounterA,4);
				  if(!okc){contin=false;}
                  wm=pointList.X[0];
                                       
                  makeCT(pCounterC,4);
                  makeCTC(pCounterC,4);
                  makeCTY(pCounterC,4,3);
                  okc=solveC(pCounterC,4);
				  if(!okc){contin=false;}
                  wp=pointList.X[0];
                                          
                  Liw[4]=1./(pointList.deltaT)*(0.5*wp-0.5*wm);

				  ////////////////////

                  

                  //%omega,strain,div,ref
                  w1=Liw[2]-Liv[3];
                  w2=Liu[3]-Liw[1];
                  w3=Liv[1]-Liu[2];
                  s11=Liu[1];
                  s22=Liv[2];
                  s33=Liw[3];
                  s12=0.5*(Liu[2]+Liv[1]);
                  s13=0.5*(Liu[3]+Liw[1]);
                  s23=0.5*(Liv[3]+Liw[2]);
                  div=fabs(trace);
                  ref=fabs(s11)+fabs(s22)+fabs(s33);
                                    
                  //acceleration quality: Vorbereitung für polynomial fits
				  Dx=Liax[0];lx=Liu[4];cx=Liu[0]*Liu[1]+Liv[0]*Liu[2]+Liw[0]*Liu[3];
				  Dy=Liay[0];ly=Liv[4];cy=Liu[0]*Liv[1]+Liv[0]*Liv[2]+Liw[0]*Liv[3];
				  Dz=Liaz[0];lz=Liw[4];cz=Liu[0]*Liw[1]+Liv[0]*Liw[2]+Liw[0]*Liw[3];

                  dix=fabs(Dx-lx-cx);
                  diy=fabs(Dy-ly-cy);
                  diz=fabs(Dz-lz-cz);
				  refx=fabs(Dx)+fabs(lx)+fabs(cx);
				  refy=fabs(Dy)+fabs(ly)+fabs(cy);
				  refz=fabs(Dz)+fabs(lz)+fabs(cz);
                  absDi=(1./3.)*(dix/refx+diy/refy+diz/refz);

				  wsq=w1*w1+w2*w2+w3*w3;
                  twosijsij=2.*(s11*s11+s22*s22+s33*s33
                           +2.*(s12*s12+s13*s13+s23*s23));

                  Q=(1./4.)*(wsq-twosijsij);
				  diva=Liax[1]+Liay[2]+Liaz[3];
				  ref_diva=fabs((1./4.)*wsq)+fabs((1./4.)*twosijsij)+fabs(Liax[1])+fabs(Liay[2])+fabs(Liaz[3]); // Beat und Marc Juni 2011: (1./4.)*
					
				  // Vorbereitung für polynomial fit
                  pointList.traj[numInTraj][ 0]=pointList.point[time][n][2];
                  pointList.traj[numInTraj][ 1]=pointList.point[time][n][3];
                  pointList.traj[numInTraj][ 2]=pointList.point[time][n][4];
                  pointList.traj[numInTraj][ 3]=Liu[0];
                  pointList.traj[numInTraj][ 4]=Liv[0];
                  pointList.traj[numInTraj][ 5]=Liw[0];
                  pointList.traj[numInTraj][ 6]=Liax[0];
                  pointList.traj[numInTraj][ 7]=Liay[0];
                  pointList.traj[numInTraj][ 8]=Liaz[0];
                  pointList.traj[numInTraj][ 9]=w1;
                  pointList.traj[numInTraj][10]=w2;
                  pointList.traj[numInTraj][11]=w3;
                  pointList.traj[numInTraj][12]=s11;
                  pointList.traj[numInTraj][13]=s12;
                  pointList.traj[numInTraj][14]=s13;
                  pointList.traj[numInTraj][15]=s22;
                  pointList.traj[numInTraj][16]=s23;
                  pointList.traj[numInTraj][17]=s33;
                  pointList.traj[numInTraj][18]=Liu[4];//dudt
                  pointList.traj[numInTraj][19]=Liv[4];//dvdt
                  pointList.traj[numInTraj][20]=Liw[4];//dwdt

                  pointList.traj[numInTraj][21]=Liax[1];
                  pointList.traj[numInTraj][22]=Liax[2];
				  pointList.traj[numInTraj][23]=Liax[3];
				  pointList.traj[numInTraj][24]=Liay[1];
                  pointList.traj[numInTraj][25]=Liay[2];
				  pointList.traj[numInTraj][26]=Liay[3];
				  pointList.traj[numInTraj][27]=Liaz[1];
                  pointList.traj[numInTraj][28]=Liaz[2];
				  pointList.traj[numInTraj][29]=Liaz[3];
					
				  // Gewichtungen
                  if(ref_diva>0){
					  pointList.traj[numInTraj][30]=pointList.weDiv*fabs(2*Q+diva)/ref_diva+pointList.weAcc*absDi+pointList.weVel*div/ref;
                  }
                  else{
                      pointList.traj[numInTraj][30]=0.95; // Gewichtung für Punkte bei denen docubicspline und/oder die räumliche Interpolation nicht geklappt hat 
                  }
                  if(pointList.traj[numInTraj][30]>0.95){
                      pointList.traj[numInTraj][30]=0.95;
                  }
                  pointList.traj[numInTraj][31]=n;
                               
               }// end of if pCOunter>3 solve...

			   // Catch falls räumliche Interpolation versagt: es werden dennoch die Informationen verwendet, die da sind, 
			   // d.h. alle Punkte der Trajektorie werden verwendet
               if(!(pCounterA>minCounter && pCounterB>minCounter && pCounterC>minCounter ) || !(contin)){
                  pointList.traj[numInTraj][ 0]=pointList.point[time][n][2];
                  pointList.traj[numInTraj][ 1]=pointList.point[time][n][3];
                  pointList.traj[numInTraj][ 2]=pointList.point[time][n][4];
                  pointList.traj[numInTraj][ 3]=pointList.point[time][n][5];
                  pointList.traj[numInTraj][ 4]=pointList.point[time][n][6];
                  pointList.traj[numInTraj][ 5]=pointList.point[time][n][7];
                  pointList.traj[numInTraj][ 6]=pointList.point[time][n][8];
                  pointList.traj[numInTraj][ 7]=pointList.point[time][n][9];
                  pointList.traj[numInTraj][ 8]=pointList.point[time][n][10];
                  pointList.traj[numInTraj][ 9]=0.;
                  pointList.traj[numInTraj][10]=0.;
                  pointList.traj[numInTraj][11]=0.;
                  pointList.traj[numInTraj][12]=0.;
                  pointList.traj[numInTraj][13]=0.;
                  pointList.traj[numInTraj][14]=0.;
                  pointList.traj[numInTraj][15]=0.;
                  pointList.traj[numInTraj][16]=0.;
                  pointList.traj[numInTraj][17]=0.;
                  pointList.traj[numInTraj][18]=0.;//dudt
                  pointList.traj[numInTraj][19]=0.;//dvdt
                  pointList.traj[numInTraj][20]=0.;//dwdt

                  pointList.traj[numInTraj][21]=0;
                  pointList.traj[numInTraj][22]=0;
                  pointList.traj[numInTraj][23]=0;
                  pointList.traj[numInTraj][24]=0;
                  pointList.traj[numInTraj][25]=0;
                  pointList.traj[numInTraj][26]=0;
				  pointList.traj[numInTraj][27]=0;
                  pointList.traj[numInTraj][28]=0;
                  pointList.traj[numInTraj][29]=0;

                  pointList.traj[numInTraj][30]=1.;   //Wichtig
                  pointList.traj[numInTraj][31]=(double)n;
                  pointList.noDeriv++;

               }
			      }//end of derivatives

				  // gefilterte Infromationen werden direkt auf die Trajektorie geschrieben, gefiltert von cubicspline
				  else{
                      pointList.traj[numInTraj][ 0]=pointList.point[time][n][2];
                      pointList.traj[numInTraj][ 1]=pointList.point[time][n][3];
                      pointList.traj[numInTraj][ 2]=pointList.point[time][n][4];
                      pointList.traj[numInTraj][ 3]=pointList.point[time][n][5];
                      pointList.traj[numInTraj][ 4]=pointList.point[time][n][6];
                      pointList.traj[numInTraj][ 5]=pointList.point[time][n][7];
                      pointList.traj[numInTraj][ 6]=pointList.point[time][n][8];
                      pointList.traj[numInTraj][ 7]=pointList.point[time][n][9];
                      pointList.traj[numInTraj][ 8]=pointList.point[time][n][10];
					  pointList.traj[numInTraj][ 9]=(double)n;
				  }

               numInTraj++; //da sollte doch irgendein check sein, wenn numFrames überschritten wird, i.e. max trajlength erreicht ist.

               //schauen ob's einen nächsten gibt
			   if(pointList.point[time][n][1]>0 && time<pointList.lastFile-pointList.firstFile && time<pointList.numOfFrames-2){ //Beat March 2013, perhaps now the trajectories don't have these strange links
                   n=pointList.point[time][n][1];
                   time++; //da sollte doch irgendein check sein, wenn numFrames überschritten wird, i.e. max trajlength erreicht ist.
                   if( pointList.point[time][n][11]<1. ){  
                       ok=false;
                   }
               }
               else{
                  ok=false;
               }
            }//end while ok

			     if(pointList.derivatives){
                 
// Wenn eine Trajektorie fertig ist, dann wird weighted polynomial fit durchgeführt

            if(numInTraj-pointList.noDeriv>pointList.minTrajLength-1){   //Wichtig
               /////polynom business////////////////////////////////////////
               double su=0.;
               double x4[500],x5[500],x6[500];
               double x7[500],x8[500],x9[500],x10[500],x11[500];
               double x12[500],x13[500],x14[500];
               double x15[500],x16[500],x17[500];
			   double x18[500],x19[500],x20[500];
			   double x21[500],x22[500],x23[500],x24[500];

               double xp[500],yp[500],zp[500],up[500],vp[500],wp[500];
               double axp[500],ayp[500],azp[500];
               double w1p[500],w2p[500],w3p[500];
               double s11p[500],s12p[500],s13p[500],s22p[500],s23p[500],s33p[500];
               double utp[500],vtp[500],wtp[500];
               double daxdxp[500],daxdyp[500],daxdzp[500];
			   double daydxp[500],daydyp[500],daydzp[500];
			   double dazdxp[500],dazdyp[500],dazdzp[500];

               setAllMatrixesToZero(4);

			   // Bestimmung der Gewichtung (Beat's Diss eq. 2.29)
			   // Ordnung
               for(int ii=0;ii<numInTraj;ii++){
                  su=su+1-pointList.traj[ii][30];//reldiv(ii)
               }
               int order=(int)(su/pointList.polyConst+3.5);
               if(numInTraj<5){
                  order=2;
               }
               if(numInTraj<2){
                  order=1;
               }
                       
			   // Gewichtung
               for(int ii=0;ii<numInTraj;ii++){
                  pointList.we[ii]=1.-1./(1.+exp(-pointList.c1*(pointList.traj[ii][30]-pointList.c2)));//reldiv(ii) (Beat's Diss eq. 2.30)
               }

			   // Abfüllen der Matrix
               for(int ii=0;ii<numInTraj;ii++){
                  for(int ij=0;ij<order;ij++){
                     pointList.A[ii][ij]=pointList.we[ii]*pow((double)ii*pointList.deltaT+0.000000001,(double)(ij));
                  }
                 
                  pointList.y[4] [ii]=pointList.we[ii]*pointList.traj[ii][ 9];//w1(i)
                  pointList.y[5] [ii]=pointList.we[ii]*pointList.traj[ii][10];//w2(i)
                  pointList.y[6] [ii]=pointList.we[ii]*pointList.traj[ii][11];//w3(i)
                  pointList.y[7] [ii]=pointList.we[ii]*pointList.traj[ii][12];//s11(i)
                  pointList.y[8] [ii]=pointList.we[ii]*pointList.traj[ii][13];//s12(i)
                  pointList.y[9] [ii]=pointList.we[ii]*pointList.traj[ii][14];//s13(i)
                  pointList.y[10][ii]=pointList.we[ii]*pointList.traj[ii][15];//s22(i)
                  pointList.y[11][ii]=pointList.we[ii]*pointList.traj[ii][16];//s23(i)
                  pointList.y[12][ii]=pointList.we[ii]*pointList.traj[ii][17];//s33(i)

                  pointList.y[13][ii]=pointList.we[ii]*pointList.traj[ii][18];//dudt(i)
                  pointList.y[14][ii]=pointList.we[ii]*pointList.traj[ii][19];//dvdt(i)
                  pointList.y[15][ii]=pointList.we[ii]*pointList.traj[ii][20];//dwdt(i)

                  pointList.y[16][ii]=pointList.we[ii]*pointList.traj[ii][21];//daxdx(i)
                  pointList.y[17][ii]=pointList.we[ii]*pointList.traj[ii][22];//daxdy(i)
                  pointList.y[18][ii]=pointList.we[ii]*pointList.traj[ii][23];//daxdz(i)
				  pointList.y[19][ii]=pointList.we[ii]*pointList.traj[ii][24];//daydx(i)
                  pointList.y[20][ii]=pointList.we[ii]*pointList.traj[ii][25];//daydy(i)
                  pointList.y[21][ii]=pointList.we[ii]*pointList.traj[ii][26];//daydz(i)
				  pointList.y[22][ii]=pointList.we[ii]*pointList.traj[ii][27];//dazdx(i)
                  pointList.y[23][ii]=pointList.we[ii]*pointList.traj[ii][28];//dazdy(i)
                  pointList.y[24][ii]=pointList.we[ii]*pointList.traj[ii][29];//dazdz(i)
               }
               
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,4);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x4[ii]=pointList.X[ii];//w1
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,5);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x5[ii]=pointList.X[ii];//w2
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,6);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x6[ii]=pointList.X[ii];//w3
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,7);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x7[ii]=pointList.X[ii];//s11
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,8);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x8[ii]=pointList.X[ii];//s12
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,9);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x9[ii]=pointList.X[ii];//s13
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,10);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x10[ii]=pointList.X[ii];//s22
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,11);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x11[ii]=pointList.X[ii];//s23
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,12);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x12[ii]=pointList.X[ii];//s33
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,13);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x13[ii]=pointList.X[ii];//dudt
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,14);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x14[ii]=pointList.X[ii];//dvdt
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,15);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x15[ii]=pointList.X[ii];//dwdt
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,16);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x16[ii]=pointList.X[ii];//daxdx
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,17);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x17[ii]=pointList.X[ii];//daxdy
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,18);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x18[ii]=pointList.X[ii];//daxdz
               }
			   makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,19);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x19[ii]=pointList.X[ii];//daydx
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,20);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x20[ii]=pointList.X[ii];//daydy
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,21);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x21[ii]=pointList.X[ii];//daydz
               }
			   makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,22);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x22[ii]=pointList.X[ii];//dazdx
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,23);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x23[ii]=pointList.X[ii];//dazdy
               }
               makeAT(numInTraj,order);
               makeATA(numInTraj,order);
               makeATY(numInTraj,order,24);
               solve(numInTraj,order);
               for(int ii=0;ii<order;ii++){
                  x24[ii]=pointList.X[ii];//dazdz
               }

			   // Initialisierung
               for(int ii=0;ii<numInTraj;ii++){
                  w1p[ii]=0;
                  w2p[ii]=0;
                  w3p[ii]=0;
                  s11p[ii]=0;
                  s12p[ii]=0;
                  s13p[ii]=0;
                  s22p[ii]=0;
                  s23p[ii]=0;
                  s33p[ii]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][12]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][13]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][14]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][15]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][16]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][17]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][18]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][19]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][20]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][21]=0;//grdients?
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][22]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][23]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][24]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][25]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][26]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][27]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][28]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][29]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][30]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][31]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][32]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][33]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][34]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][35]=0;
				  pointList.point[ii+2][(int)pointList.traj[ii][31]][36]=0;
                  utp[ii]=0;
                  vtp[ii]=0;
                  wtp[ii]=0;
                  daxdxp[ii]=0;
                  daxdyp[ii]=0;
                  daxdzp[ii]=0;
				  daydxp[ii]=0;
                  daydyp[ii]=0;
                  daydzp[ii]=0;
				  dazdxp[ii]=0;
                  dazdyp[ii]=0;
                  dazdzp[ii]=0;

                  xp[ii] =pointList.traj[ii][0];
                  yp[ii] =pointList.traj[ii][1];
                  zp[ii] =pointList.traj[ii][2];
                  up[ii] =pointList.traj[ii][3];
                  vp[ii] =pointList.traj[ii][4];
                  wp[ii] =pointList.traj[ii][5];
                  axp[ii]=pointList.traj[ii][6];
                  ayp[ii]=pointList.traj[ii][7];
                  azp[ii]=pointList.traj[ii][8];
				  // Polynom in der Schlaufe
                  for(int ij=0;ij<order;ij++){
                     w1p[ii]= w1p[ii]+ x4[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));// change to get non-filtered data: pointList.traj[ii][ 9];//
                     w2p[ii]= w2p[ii]+ x5[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     w3p[ii]= w3p[ii]+ x6[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s11p[ii]=s11p[ii]+ x7[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s12p[ii]=s12p[ii]+ x8[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s13p[ii]=s13p[ii]+ x9[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s22p[ii]=s22p[ii]+ x10[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s23p[ii]=s23p[ii]+x11[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     s33p[ii]=s33p[ii]+x12[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     utp[ii]=utp[ii]+x13[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     vtp[ii]=vtp[ii]+x14[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     wtp[ii]=wtp[ii]+x15[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     daxdxp[ii]=daxdxp[ii]+x16[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     daxdyp[ii]=daxdyp[ii]+x17[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
					 daxdzp[ii]=daxdzp[ii]+x18[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
					 daydxp[ii]=daydxp[ii]+x19[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     daydyp[ii]=daydyp[ii]+x20[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
					 daydzp[ii]=daydzp[ii]+x21[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
					 dazdxp[ii]=dazdxp[ii]+x22[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                     dazdyp[ii]=dazdyp[ii]+x23[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
					 dazdzp[ii]=dazdzp[ii]+x24[ij]*pow((double)ii*pointList.deltaT+0.00001,(double)(ij));
                  }              
                  
               }// end for loop through traj

			   // Statistik für den Zuschauer
			   for(int ii=0;ii<numInTraj;ii++){                   
					 USq=up[ii]*up[ii]+vp[ii]*vp[ii]+wp[ii]*wp[ii];

                     ww1=w1p[ii]*s11p[ii]+w2p[ii]*s12p[ii]+w3p[ii]*s13p[ii];
                     ww2=w1p[ii]*s12p[ii]+w2p[ii]*s22p[ii]+w3p[ii]*s23p[ii];
                     ww3=w1p[ii]*s13p[ii]+w2p[ii]*s23p[ii]+w3p[ii]*s33p[ii];
                     wwsij=w1p[ii]*ww1+w2p[ii]*ww2+w3p[ii]*ww3;

                     s111=s11p[ii]*s11p[ii]*s11p[ii];
                     s222=s22p[ii]*s22p[ii]*s22p[ii];
                     s333=s33p[ii]*s33p[ii]*s33p[ii];
                     s112=s11p[ii]*s12p[ii]*s12p[ii];
                     s113=s11p[ii]*s13p[ii]*s13p[ii];
                     s221=s22p[ii]*s12p[ii]*s12p[ii];
                     s223=s22p[ii]*s23p[ii]*s23p[ii];
                     s331=s33p[ii]*s13p[ii]*s13p[ii];
                     s332=s33p[ii]*s23p[ii]*s23p[ii];
                     s123=s12p[ii]*s23p[ii]*s13p[ii];
                     sijsjkski=s111+s222+s333+3.*(s112+s113+s221+s223+s331+s332)+6.*s123; // mal überprüfen...

                     wsq=w1p[ii]*w1p[ii]+w2p[ii]*w2p[ii]+w3p[ii]*w3p[ii];
                     twosijsij=2.*(s11p[ii]*s11p[ii]+s22p[ii]*s22p[ii]+s33p[ii]*s33p[ii]
                                 +2.*(s12p[ii]*s12p[ii]+s13p[ii]*s13p[ii]+s23p[ii]*s23p[ii]));
			         diss=pointList.viscosity*twosijsij;

                     //rel diva quality
					 Q=(1./4.)*(wsq-twosijsij);
                     R=-(1./3.)*(sijsjkski+(3./4.)*wwsij);
					 diva=daxdxp[ii]+daydyp[ii]+dazdzp[ii];
				     ref_diva=2.*fabs((1./4.)*wsq)+2.*fabs((1./4.)*twosijsij)+fabs(daxdxp[ii])+fabs(daydyp[ii])+fabs(dazdzp[ii]);
                     if(ref_diva>0){
                         reldiva=fabs(2*Q+diva)/ref_diva; 
					 }
					 else{
                         reldiva=1.; 
					 }

                     //acceleration quality
					 ux=s11p[ii];
                     uy=s12p[ii]-0.5*w3p[ii];
                     uz=s13p[ii]+0.5*w2p[ii];
                     vx=s12p[ii]+0.5*w3p[ii];
                     vy=s22p[ii];
                     vz=s23p[ii]-0.5*w1p[ii];
                     wx=s13p[ii]-0.5*w2p[ii];
                     wy=s23p[ii]+0.5*w1p[ii];
                     wz=s33p[ii];
				     Dx=axp[ii];lx=utp[ii];cx=up[ii]*ux+vp[ii]*uy+wp[ii]*uz;
				     Dy=ayp[ii];ly=vtp[ii];cy=up[ii]*vx+vp[ii]*vy+wp[ii]*vz;
				     Dz=azp[ii];lz=wtp[ii];cz=up[ii]*wx+vp[ii]*wy+wp[ii]*wz;

                     dix=fabs(Dx-lx-cx);
                     diy=fabs(Dy-ly-cy);
                     diz=fabs(Dz-lz-cz);
				     refx=fabs(Dx)+fabs(lx)+fabs(cx);
				     refy=fabs(Dy)+fabs(ly)+fabs(cy);
				     refz=fabs(Dz)+fabs(lz)+fabs(cz);
					 if(refx>0 && refy>0 && refz>0){
                        absDi=(1./3.)*(dix/refx+diy/refy+diz/refz);
			         }
			         else{
                        absDi=1.;
			         }
                        

					 // rel divu quality
					 div=fabs(ux+vy+wz);
					 ref=fabs(ux)+fabs(vy)+fabs(wz);
					 if(ref>0){
                        reldivu=div/ref;
			         }
			         else{
                        reldivu=1.;
			         }

					 //totQuality
					 quality=pointList.weDiv*reldiva+pointList.weAcc*absDi+pointList.weVel*reldivu;
					    
					    //prepare for xuag files
                        pointList.point[ii+2][(int)pointList.traj[ii][31]][12]=ux;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][13]=uy;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][14]=uz;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][15]=vx;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][16]=vy;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][17]=vz;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][18]=wx;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][19]=wy;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][20]=wz;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][21]=1;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][22]=utp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][23]=vtp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][24]=wtp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][25]=daxdxp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][26]=daxdyp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][27]=daxdxp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][28]=daydxp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][29]=daydyp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][30]=daydzp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][31]=dazdxp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][32]=dazdyp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][33]=dazdzp[ii];
						pointList.point[ii+2][(int)pointList.traj[ii][31]][34]=quality;
						pointList.point[ii+2][(int)pointList.traj[ii][31]][35]=(double)numInTraj;	// 
						pointList.point[ii+2][(int)pointList.traj[ii][31]][36]=(double)ii;			// added by Beat and Marc, 31.08.2011
					    //end of prepare xuag files
                     
					 if(pointList.weDiv*fabs(2*Q+diva)/ref_diva+pointList.weAcc*absDi+pointList.weVel*reldivu){
                        pointList.count2++;
                        pointList.meanDiss=(pointList.meanDiss*(double)(pointList.count2-1)+diss)/(double)pointList.count2;
                        pointList.meanUSq=(pointList.meanUSq*(double)(pointList.count2-1)+USq)/(double)pointList.count2;                        
					 }
					 if(reldiva<0.1){
                        pointList.count4++;
					 }
					 if(absDi<0.2){
                        pointList.count5++;
					 }
					 if(reldivu<0.1){
                        pointList.count6++;
					 }
					 
                     fprintf(fpp, "%lf\t", xp[ii]);//1
                     fprintf(fpp, "%lf\t", yp[ii]);//2
                     fprintf(fpp, "%lf\t", zp[ii]);//3
                     fprintf(fpp, "%lf\t", up[ii]);//4
                     fprintf(fpp, "%lf\t", vp[ii]);//5
                     fprintf(fpp, "%lf\t", wp[ii]);//6
                     fprintf(fpp, "%lf\t", axp[ii]);//7
                     fprintf(fpp, "%lf\t", ayp[ii]);//8
                     fprintf(fpp, "%lf\t", azp[ii]);//9
                     fprintf(fpp, "%lf\t", w1p[ii]);//10
                     fprintf(fpp, "%lf\t", w2p[ii]);//11
                     fprintf(fpp, "%lf\t", w3p[ii]);//12
                     fprintf(fpp, "%lf\t", s11p[ii]);//13
                     fprintf(fpp, "%lf\t", s12p[ii]);//14
                     fprintf(fpp, "%lf\t", s13p[ii]);//15
                     fprintf(fpp, "%lf\t", s22p[ii]);//16
                     fprintf(fpp, "%lf\t", s23p[ii]);//17
                     fprintf(fpp, "%lf\t", s33p[ii]);//18
					 fprintf(fpp, "%lf\t", utp[ii]);//19
                     fprintf(fpp, "%lf\t", vtp[ii]);//20
                     fprintf(fpp, "%lf\t", wtp[ii]);//21
                     fprintf(fpp, "%lf\t", daxdxp[ii]);//22
                     fprintf(fpp, "%lf\t", daxdyp[ii]);//23
                     fprintf(fpp, "%lf\t", daxdzp[ii]);//24
					 fprintf(fpp, "%lf\t", daydxp[ii]);//25
                     fprintf(fpp, "%lf\t", daydyp[ii]);//26
                     fprintf(fpp, "%lf\t", daydzp[ii]);//27
					 fprintf(fpp, "%lf\t", dazdxp[ii]);//28
                     fprintf(fpp, "%lf\t", dazdyp[ii]);//29
                     fprintf(fpp, "%lf\t", dazdzp[ii]);//30
                     fprintf(fpp, "%lf\t", quality);//31 0=good, 1=bad
                     fprintf(fpp, "%lf\t", (double)(ii));//32 age along trajectory
					 fprintf(fpp, "%lf\n",pointList.traj[ii][9]);//33 reference to index in rt_is, ptv_is, xuap files
                     
                  }// end for
                  ////end of polynom business
               
            } //end if of polynom buisness
			    }//end of derivatives
				else{
                  double xp[500],yp[500],zp[500],up[500],vp[500],wp[500];
                  double axp[500],ayp[500],azp[500];
                  for(int ii=0;ii<numInTraj;ii++){ 
					  xp[ii] =pointList.traj[ii][0];
                      yp[ii] =pointList.traj[ii][1];
                      zp[ii] =pointList.traj[ii][2];
                      up[ii] =pointList.traj[ii][3];
                      vp[ii] =pointList.traj[ii][4];
                      wp[ii] =pointList.traj[ii][5];
                      axp[ii]=pointList.traj[ii][6];
                      ayp[ii]=pointList.traj[ii][7];
                      azp[ii]=pointList.traj[ii][8];
                     fprintf(fpp, "%lf\t", xp[ii]);//1
                     fprintf(fpp, "%lf\t", yp[ii]);//2
                     fprintf(fpp, "%lf\t", zp[ii]);//3
                     fprintf(fpp, "%lf\t", up[ii]);//4
                     fprintf(fpp, "%lf\t", vp[ii]);//5
                     fprintf(fpp, "%lf\t", wp[ii]);//6
                     fprintf(fpp, "%lf\t", axp[ii]);//7
                     fprintf(fpp, "%lf\t", ayp[ii]);//8
                     fprintf(fpp, "%lf\t", azp[ii]);//9
					 fprintf(fpp, "%lf\t", (double)(ii));//32 age along trajectory
					 fprintf(fpp, "%lf\n",pointList.traj[ii][9]);///n, reference to index in rt_is, ptv_is, xuap files
				  }
				}
            
         } // end if not occ und central
         
     }// end haupt for schlaufe
     
}

void readXUAPFile(int n, bool firstTime)
{
    int numOfPoints;
    double left,right,x,y,z,u,v,w,ax,ay,az,dummy,cubic;

    FILE *fpp;
    char name[256];
    int c;

	
    FILE *fpp_xuag;
    char name_xuag[256];
    int c_xuag;
    c_xuag=sprintf (name_xuag, pointList.path);
    c_xuag+=sprintf (name_xuag+c_xuag, "/xuag.");
    c_xuag+=sprintf (name_xuag+c_xuag, "%1d", n-1);

    int ind_X,ind_Y,ind_Z,ind_C,dummy_count;

	dummy_count=0;
    
    
       for(int i=0;i<pointList.numOfFrames;i++){

          
		   

          if(n-2+i>pointList.firstFile-1 && n-2+i<pointList.lastFile+1){
             if(i<pointList.numOfFrames-1 && !(firstTime)){
                 //write xuag
				 if(i==2){
				    fpp_xuag = fopen(name_xuag,"w");
				    for(int j=1;j<pointList.point[2][0][0]+1;j++){           
                       for(int arg=0;arg<36;arg++){
                          fprintf(fpp_xuag, "%lf\t", pointList.point[2][j][arg]);
                       }
                       fprintf(fpp_xuag, "%lf\n", pointList.point[2][j][36]);                        	
                    }
				    fclose (fpp_xuag);
				 }
			     //end write xuag
				 //rotate point entries
				 for(int j=0;j<pointList.point[i+1][0][0]+1;j++){           
					for(int k=0;k<37;k++){
                        pointList.point[i][j][k]=pointList.point[i+1][j][k];
                    }
                 }
				 //Beat March 2013 also rotate fast_search grid entries
				 for (int ii=0;ii<pointList.max_grid_X;ii++){
			         for (int jj=0;jj<pointList.max_grid_Y;jj++){
				         for (int kk=0;kk<pointList.max_grid_Z;kk++){
							 for (int cc=0;cc<pointList.max_grid_C+1;cc++){
					              pointList.fast_search[i][ii][jj][kk][cc]=pointList.fast_search[i+1][ii][jj][kk][cc];
							 }
				         }
			         }
		         }
             }
             else{
				 cout << "reading xuap ................."<<n-2+i<<"\n";

                //Beat March 2013, init fast_search counter
		        for (int ii=0;ii<pointList.max_grid_X;ii++){
			        for (int jj=0;jj<pointList.max_grid_Y;jj++){
				        for (int kk=0;kk<pointList.max_grid_Z;kk++){
					        pointList.fast_search[i][ii][jj][kk][pointList.max_grid_C]=0;
				        }
			        }
		        }

                numOfPoints=0;
                c=sprintf (name, pointList.path);
                c+=sprintf (name+c, "/xuap.");
                c+=sprintf (name+c, "%1d", n-2+i);
                fpp = fopen(name,"r");
                while(!feof(fpp)){
                   numOfPoints++;
                   fscanf (fpp, "%lf\0", &left);
                   fscanf (fpp, "%lf\0", &right);
                   fscanf (fpp, "%lf\0", &dummy); //measured x
                   fscanf (fpp, "%lf\0", &dummy); //measured y
                   fscanf (fpp, "%lf\0", &dummy); //measured z
                   fscanf (fpp, "%lf\0", &x); //cubic spline x
                   fscanf (fpp, "%lf\0", &y); //cubic spline y
                   fscanf (fpp, "%lf\0", &z); //cubic spline z
                   fscanf (fpp, "%lf\0", &u);
                   fscanf (fpp, "%lf\0", &v);
                   fscanf (fpp, "%lf\0", &w);
                   fscanf (fpp, "%lf\0", &ax);
                   fscanf (fpp, "%lf\0", &ay);
                   fscanf (fpp, "%lf\0", &az);
                   fscanf (fpp, "%lf\0", &cubic);
                   pointList.point[i][numOfPoints][0]=left;
                   pointList.point[i][numOfPoints][1]=right;
                   pointList.point[i][numOfPoints][2]=x;
                   pointList.point[i][numOfPoints][3]=y;
                   pointList.point[i][numOfPoints][4]=z;
                   pointList.point[i][numOfPoints][5]=u;
                   pointList.point[i][numOfPoints][6]=v;
                   pointList.point[i][numOfPoints][7]=w;
                   pointList.point[i][numOfPoints][8]=ax;
                   pointList.point[i][numOfPoints][9]=ay;
                   pointList.point[i][numOfPoints][10]=az;
                   pointList.point[i][numOfPoints][11]=cubic;

				   //Beat March 2013 prepare fast search
				   ind_X=(int)((double)(x-pointList.minX)/pointList.dh_X);
				   ind_Y=(int)((double)(y-pointList.minY)/pointList.dh_Y);
				   ind_Z=(int)((double)(z-pointList.minZ)/pointList.dh_Z);
				   if(ind_X<0 ){ind_X=0 ;}
				   if(ind_Y<0 ){ind_Y=0 ;}
				   if(ind_Z<0 ){ind_Z=0 ;}
				   if(ind_X>pointList.max_grid_X-1){ind_X=pointList.max_grid_X-1;}
				   if(ind_Y>pointList.max_grid_Y-1){ind_Y=pointList.max_grid_Y-1;}
				   if(ind_Z>pointList.max_grid_Z-1){ind_Z=pointList.max_grid_Z-1;}
				   if(cubic==1){
				      ind_C=pointList.fast_search[i][ind_X][ind_Y][ind_Z][pointList.max_grid_C];
					  pointList.fast_search[i][ind_X][ind_Y][ind_Z][ind_C]=numOfPoints;
					  pointList.fast_search[i][ind_X][ind_Y][ind_Z][pointList.max_grid_C]++;
			          if(pointList.fast_search[i][ind_X][ind_Y][ind_Z][pointList.max_grid_C]>=pointList.max_grid_C){
					      cout << "\nups, max in bin is reached, adpat in stdafx.h, read_scanning_PTVFile, folloTrajPoint\n";
					      pointList.fast_search[i][ind_X][ind_Y][ind_Z][pointList.max_grid_C]=pointList.max_grid_C;
						  dummy_count++;
				      }
				   }

                }
                fclose (fpp);
                pointList.point[i][0][0]=numOfPoints;
             }
          }
          else{
             pointList.point[i][0][0]=0;
          }
       }
}

void prepare_fast_search()
{
    double left,right,x,y,z,u,v,w,ax,ay,az,dummy,cubic;

    FILE *fpp;
    char name[256];
    int c;

	pointList.minX=1e6;
	pointList.minY=1e6;
	pointList.minZ=1e6;
	pointList.maxX=-1e6;
	pointList.maxY=-1e6;
	pointList.maxZ=-1e6;

    c=sprintf (name, pointList.path);
    c+=sprintf (name+c, "/xuap.");
	c+=sprintf (name+c, "%1d", pointList.firstFile+pointList.minLeftRight+2);
    fpp = fopen(name,"r");
    while(!feof(fpp)){
         fscanf (fpp, "%lf\0", &left);
         fscanf (fpp, "%lf\0", &right);
         fscanf (fpp, "%lf\0", &dummy); //measured x
         fscanf (fpp, "%lf\0", &dummy); //measured y
         fscanf (fpp, "%lf\0", &dummy); //measured z
         fscanf (fpp, "%lf\0", &x); //cubic spline x
         fscanf (fpp, "%lf\0", &y); //cubic spline y
         fscanf (fpp, "%lf\0", &z); //cubic spline z
         fscanf (fpp, "%lf\0", &u);
         fscanf (fpp, "%lf\0", &v);
         fscanf (fpp, "%lf\0", &w);
         fscanf (fpp, "%lf\0", &ax);
         fscanf (fpp, "%lf\0", &ay);
         fscanf (fpp, "%lf\0", &az);
         fscanf (fpp, "%lf\0", &cubic);
		 if(x<pointList.minX && cubic==1){pointList.minX=x;}
		 if(y<pointList.minY && cubic==1){pointList.minY=y;}
		 if(z<pointList.minZ && cubic==1){pointList.minZ=z;}
		 if(x>pointList.maxX && cubic==1){pointList.maxX=x;}
		 if(y>pointList.maxY && cubic==1){pointList.maxY=y;}
		 if(z>pointList.maxZ && cubic==1){pointList.maxZ=z;}
    }
    fclose (fpp);

	cout << "\nmin x ......................."<<pointList.minX<<"\n";//min max stuff
	cout << "max x ......................."<<pointList.maxX<<"\n";//min max stuff
	cout << "min y ......................."<<pointList.minY<<"\n";//min max stuff
	cout << "max y ......................."<<pointList.maxX<<"\n";//min max stuff
	cout << "min z ......................."<<pointList.minZ<<"\n";//min max stuff
	cout << "max z ......................."<<pointList.maxZ<<"\n";//min max stuff

	pointList.num_X=(int)((double)(pointList.maxX-pointList.minX)/pointList.maxRadius);
	pointList.num_Y=(int)((double)(pointList.maxY-pointList.minY)/pointList.maxRadius);
	pointList.num_Z=(int)((double)(pointList.maxZ-pointList.minZ)/pointList.maxRadius);

	cout << "\nopt. # grid cells in x-dir..."<<pointList.num_X<<"\n";//min max stuff 
	cout << "opt. # grid cells in y-dir..."<<pointList.num_Y<<"\n";//min max stuff
	cout << "opt. # grid cells in z-dir..."<<pointList.num_Z<<"\n";//min max stuff

	if(pointList.num_X>pointList.max_grid_X){pointList.num_X=pointList.max_grid_X;}
	if(pointList.num_Y>pointList.max_grid_Y){pointList.num_Y=pointList.max_grid_Y;}
	if(pointList.num_Z>pointList.max_grid_Z){pointList.num_Z=pointList.max_grid_Z;}

	cout << "\n# grid cells in x-dir........"<<pointList.num_X<<"\n";//min max stuff 
	cout << "# grid cells in y-dir........"<<pointList.num_Y<<"\n";//min max stuff
	cout << "# grid cells in z-dir........"<<pointList.num_Z<<"\n";//min max stuff

	pointList.dh_X=(double)(pointList.maxX-pointList.minX)/(double)pointList.num_X;
	pointList.dh_Y=(double)(pointList.maxY-pointList.minY)/(double)pointList.num_Y;
	pointList.dh_Z=(double)(pointList.maxZ-pointList.minZ)/(double)pointList.num_Z;

	
	
}

void readXUAGFile(int n, bool firstTime)
{
    int numOfPoints;
    double left,right,x,y,z,u,v,w,ax,ay,az,dummy,cubic;
	double ux,uy,uz,vx,vy,vz,wx,wy,wz,grad,ut,vt,wt,axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,quality;

    FILE *fpp;
    char name[256];
    int c;

       for(int i=0;i<pointList.numOfFrames;i++){
          if(n-2+i>pointList.firstFile-1 && n-2+i<pointList.lastFile){
             if(i<pointList.numOfFrames-1 && !(firstTime)){

				 for(int j=0;j<pointList.point[i+1][0][0]+1;j++){           
					for(int k=0;k<35;k++){
                        pointList.point[i][j][k]=pointList.point[i+1][j][k];
                    }
                 }
             }
             else{
                numOfPoints=0;
                c=sprintf (name, pointList.path);
                c+=sprintf (name+c, "/xuag.");
                c+=sprintf (name+c, "%1d", n-2+i);
                fpp = fopen(name,"r");
                while(!feof(fpp)){
                   numOfPoints++;
                   fscanf (fpp, "%lf\0", &left);
                   fscanf (fpp, "%lf\0", &right);
                   fscanf (fpp, "%lf\0", &x); //cubic spline x
                   fscanf (fpp, "%lf\0", &y); //cubic spline y
                   fscanf (fpp, "%lf\0", &z); //cubic spline z
                   fscanf (fpp, "%lf\0", &u);
                   fscanf (fpp, "%lf\0", &v);
                   fscanf (fpp, "%lf\0", &w);
                   fscanf (fpp, "%lf\0", &ax);
                   fscanf (fpp, "%lf\0", &ay);
                   fscanf (fpp, "%lf\0", &az);
                   fscanf (fpp, "%lf\0", &cubic);
				   fscanf (fpp, "%lf\0", &ux);
                   fscanf (fpp, "%lf\0", &uy);
                   fscanf (fpp, "%lf\0", &uz);
				   fscanf (fpp, "%lf\0", &vx);
                   fscanf (fpp, "%lf\0", &vy);
                   fscanf (fpp, "%lf\0", &vz);
				   fscanf (fpp, "%lf\0", &wx);
                   fscanf (fpp, "%lf\0", &wy);
                   fscanf (fpp, "%lf\0", &wz);
				   fscanf (fpp, "%lf\0", &grad);
				   fscanf (fpp, "%lf\0", &ut);
				   fscanf (fpp, "%lf\0", &vt);
				   fscanf (fpp, "%lf\0", &wt);
				   fscanf (fpp, "%lf\0", &axx);
				   fscanf (fpp, "%lf\0", &axy);
				   fscanf (fpp, "%lf\0", &axz);
				   fscanf (fpp, "%lf\0", &ayx);
				   fscanf (fpp, "%lf\0", &ayy);
				   fscanf (fpp, "%lf\0", &ayz);
				   fscanf (fpp, "%lf\0", &azx);
				   fscanf (fpp, "%lf\0", &azy);
				   fscanf (fpp, "%lf\0", &azz);
				   fscanf (fpp, "%lf\0", &quality);
                   pointList.point[i][numOfPoints][0]=left;//1
                   pointList.point[i][numOfPoints][1]=right;//2
                   pointList.point[i][numOfPoints][2]=x;//3
                   pointList.point[i][numOfPoints][3]=y;//4
                   pointList.point[i][numOfPoints][4]=z;//5
                   pointList.point[i][numOfPoints][5]=u;//6
                   pointList.point[i][numOfPoints][6]=v;//7
                   pointList.point[i][numOfPoints][7]=w;//8
                   pointList.point[i][numOfPoints][8]=ax;//9
                   pointList.point[i][numOfPoints][9]=ay;//10
                   pointList.point[i][numOfPoints][10]=az;//11
                   pointList.point[i][numOfPoints][11]=cubic;//12
				   pointList.point[i][numOfPoints][12]=ux;//13
                   pointList.point[i][numOfPoints][13]=uy;//14
                   pointList.point[i][numOfPoints][14]=uz;//15
				   pointList.point[i][numOfPoints][15]=vx;//16
                   pointList.point[i][numOfPoints][16]=vy;//17
                   pointList.point[i][numOfPoints][17]=vz;//18
				   pointList.point[i][numOfPoints][18]=wx;//19
                   pointList.point[i][numOfPoints][19]=wy;//20
                   pointList.point[i][numOfPoints][20]=wz;//21
				   pointList.point[i][numOfPoints][21]=grad;//22
				   pointList.point[i][numOfPoints][22]=ut;//23
                   pointList.point[i][numOfPoints][23]=vt;//24
                   pointList.point[i][numOfPoints][24]=wt;//25
				   pointList.point[i][numOfPoints][25]=axx;//26
				   pointList.point[i][numOfPoints][26]=axy;//27
				   pointList.point[i][numOfPoints][27]=axz;//28
				   pointList.point[i][numOfPoints][28]=ayx;//29
				   pointList.point[i][numOfPoints][29]=ayy;//30
				   pointList.point[i][numOfPoints][30]=ayz;//31
				   pointList.point[i][numOfPoints][31]=azx;//32
				   pointList.point[i][numOfPoints][32]=azy;//33
				   pointList.point[i][numOfPoints][33]=azz;//34
				   pointList.point[i][numOfPoints][34]=0;//35
				   pointList.point[i][numOfPoints][35]=0;//36
				   pointList.point[i][numOfPoints][36]=0;//37
				   pointList.point[i][numOfPoints][37]=quality;//38;
                }
                fclose (fpp);
                pointList.point[i][0][0]=numOfPoints;
             }
          }
          else{
             pointList.point[i][0][0]=0;
          }
       }
}
