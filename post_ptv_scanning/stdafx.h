// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#pragma once


#include <iostream>
#include <tchar.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <malloc.h>
#include <fstream>
#include <string.h>
#include <float.h>
#include <time.h>
using namespace std;

struct TpointList
{

  char experiment[256];

  int firstSFile;
  int lastSFile;
  int firstFile;
  int lastFile;
  int numSlices; //added by Beat March 2013
  short map_slice_cycle[40000][2000];//added by Beat March 2013 <----------das ist recht kritisch, kann es nicht weiter erhöhen 
  int numPoints_per_cycle[10000];//added by Beat March 2013
  int PL,minLeftRight,PLh;
  int count,count2,count3,count4,count5,count6;

  double maxVel;
  double meanVel;
  double meanAcc;
  double meanDiss;
  double meanUSq;
  double meanUisq;
  double meanDudxsq;
  double Re;
  
  bool xuap;
  bool traj_point;
  bool derivatives;
  bool pressure;
  bool Hessian;

  char path[256];

  double deltaT;
  double deltaT_between_slice;
  double tolMaxVel;
  int polyConst;
  int maxRank;
  double maxRadius;
  double weDiv;
  double weAcc;
  double weVel;
  double viscosity;

  double point[200][25000][37]; ////// Beat March 2013, da liegt die einzige 'Schwäche': Max traj length = 166,.... ohne xuag könnten es 500 sein
  bool occ[10000][25000]; // occ[10000][25000];
  short fast_search[200][18][18][9][26];//Beat March 2013, should roughly match the aspect ratio of the observation domain
  int max_grid_X,max_grid_Y,max_grid_Z,max_grid_C;//Beat March 2013, should roughly match the aspect ratio of the observation domain
  double maxX,minX,maxY,minY,maxZ,minZ,dh_X,dh_Y,dh_Z;
  int num_X,num_Y,num_Z;

  int numOfFrames;
  bool changed;

  double u,ux,uy,uz;
  double v,vx,vy,vz;
  double w,wx,wy,wz; 
  double ax,axx,axy,axz;
  double ay,ayx,ayy,ayz;
  double az,azx,azy,azz; 

  double meaU[2];
  double meaV[2];
  double meaW[2];
  double meaAx[2];
  double meaAy[2];
  double meaAz[2];

  int maxRowIndex;

  double A  [500][100];
  double AT [100][500];
  double ATA[100][100];
  double B  [500][100];
  double BT [100][500];
  double BTB[100][100];
  double C[500][100];
  double CT [100][500];
  double CTC[100][100];
	
  double Y  [500];
  double y[50][500];
  double yC[5][500];

  double X  [500];
  double ATY[500];
  double BTY[500];
  double CTY[500];
  
  double Yu[500];
  double Yv[500];
  double Yw[500];
    
  double YuB[500];
  double YvB[500];
  double YwB[500];

  double YaxB[500];
  double YayB[500];
  double YazB[500];
    
  double Yaz[500];
  double Yay[500];
  double Yax[500];

  double YuxB[500];
  double YuyB[500];
  double YuzB[500];
  double YvxB[500];
  double YvyB[500];
  double YvzB[500];
  double YwxB[500];
  double YwyB[500];
  double YwzB[500];

  double YpxB[500];
  double YpyB[500];
  double YpzB[500];

  double Aij[3][3];
  double Aaij[3][3];
  double uxx[9];
  double pij[3][3];

  double pointPerRadius[500][2]; //0.1mm resolution up to 10mm.
  int dis[500];
  int disA[500];
  int disB[500];
  int disC[500];
  double traj[500][65];
  double we[500];
  int minTrajLength;
  int noDeriv;
  double c1;
  double c2;

  double xminChamber; //added by Markus, 20.07.2009
  double xmaxChamber;
  double yminChannel;
  double ymaxChannel;
  double zminChamber;
  double zmaxChamber;
  double zminChannel;
  double zmaxChannel;
  double yChamberChannel;

  double xminChannel;
  double xmaxChannel;


  FILE *fpp;
};



// TODO: reference additional headers your program requires here
