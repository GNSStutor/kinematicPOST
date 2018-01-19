////////////////////////////////////////////////////////////////////////////////
//
// 基準局、ユーザ位置から、衛星の仰角・方位角を計算
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include "global_extern.h"

const double PI = 3.1415926535898;//

void calc_direction(int rcvn, int iter)
{
/////////////////////////////////////////////////////////////////////
//	
/////////////////////////////////////////////////////////////////////
	int i,k;
	double xpos1[PRN],ypos1[PRN],zpos1[PRN];
	double xpos2[PRN],ypos2[PRN],zpos2[PRN];
	double lat,lon;
	double kk = 0.5;
	double SVx2[RCVN][PRN]={0},SVy2[RCVN][PRN]={0},SVz2[RCVN][PRN]={0};

	double a[6]={0};

	//仰角、方位角を算出する基準となる位置を設定
	if(rcvn==1){//基準局
		a[0]=POSrcvlat[rcvn];
		a[1]=POSrcvlon[rcvn];
		a[2]=Ref_pos[0][rcvn];
		a[3]=Ref_pos[1][rcvn];
		a[4]=Ref_pos[2][rcvn];
	}
	else{//ユーザ側
		a[0]=User_lat;
		a[1]=User_lon;
		a[2]=User_pos[0];
		a[3]=User_pos[1];
		a[4]=User_pos[2];
	}

	for(i=0;i<SATn[rcvn];i++){//number of satellite
		k = SVn[rcvn][i];//satellite number

		//user latitude and longitude
		lat=(90-a[0])*2.0*PI/360.0;
		lon=a[1]*2.0*PI/360.0;
				
///////////////////////////////////////////////////
//　自身の位置（基準面、北、東、真上）に対して衛星を回転させる
		xpos1[k]=cos(lon)*(SVx[rcvn][k]-a[2])
			+sin(lon)*(SVy[rcvn][k]-a[3]);
		ypos1[k]=-sin(lon)*(SVx[rcvn][k]-a[2])
			+cos(lon)*(SVy[rcvn][k]-a[3]);
		zpos1[k]=SVz[rcvn][k]-a[4];
				
		//
		xpos2[k]=cos(lat)*xpos1[k]-sin(lat)*zpos1[k];
		ypos2[k]=ypos1[k];
		zpos2[k]=sin(lat)*xpos1[k]+cos(lat)*zpos1[k];
				
		//
		SVx2[rcvn][k]=cos(kk*PI)*xpos2[k]+sin(kk*PI)*ypos2[k];
		SVy2[rcvn][k]=-sin(kk*PI)*xpos2[k]+cos(kk*PI)*ypos2[k];
		SVz2[rcvn][k]=zpos2[k];
///////////////////////////////////////////////////
	}					

////////////////////////////////////////////////////////////////	
//	上記結果を利用して、仰角・方位角を求める
////////////////////////////////////////////////////////////////

	double root;
	for(k=0;k<=130;k++){
		Elevation[rcvn][k]=0;
		Azimuth[rcvn][k]=0;
	}

	for(i=0;i<SATn[rcvn];i++){
		k=SVn[rcvn][i];

		root=sqrt(SVx2[rcvn][k]*SVx2[rcvn][k]+
			SVy2[rcvn][k]*SVy2[rcvn][k]);		
		Elevation[rcvn][k] = atan(SVz2[rcvn][k]/root);
		Elevation[rcvn][k] = 360.0*Elevation[rcvn][k]/2.0/PI;
				
		//0<azi<90
		if(SVx2[rcvn][k]>=0 && SVy2[rcvn][k]>=0){
			Azimuth[rcvn][k] = atan(SVx2[rcvn][k]/SVy2[rcvn][k]);
			Azimuth[rcvn][k] = 180.0*Azimuth[rcvn][k]/PI;
		}
		//90<azi<180
		if(SVx2[rcvn][k]>=0 && SVy2[rcvn][k]<=0){
			Azimuth[rcvn][k] = PI+atan(SVx2[rcvn][k]/SVy2[rcvn][k]);
			Azimuth[rcvn][k] = 180.0*Azimuth[rcvn][k]/PI;
		}
		//180<azi<270
		if(SVx2[rcvn][k]<=0 && SVy2[rcvn][k]<=0){
			Azimuth[rcvn][k] = PI+atan(SVx2[rcvn][k]/SVy2[rcvn][k]);
			Azimuth[rcvn][k] = 180.0*Azimuth[rcvn][k]/PI;
		}
		//270<azi<360
		if(SVx2[rcvn][k]<=0 && SVy2[rcvn][k]>=0){
			Azimuth[rcvn][k] = 2.0*PI+atan(SVx2[rcvn][k]/SVy2[rcvn][k]);
			Azimuth[rcvn][k] = 180.0*Azimuth[rcvn][k]/PI;
		}			
	}


	//測位衛星毎で利用判断を行う　利用しないものは仰角を強制的に0
	for(i=0;i<SATn[rcvn];i++){	
		prn=SVn[rcvn][i];
		if(Cflag==0){if(prn>=71&&prn<101)
				Elevation[rcvn][prn]=0.0;
		}if(Rflag==0){if(prn>=101)
				Elevation[rcvn][prn]=0.0;
		}if(Gflag==0){if(prn<=32)
				Elevation[rcvn][prn]=0.0;
		}if(Jflag==0){if(prn>=33 && prn<=39)
				Elevation[rcvn][prn]=0.0;
		}if(Eflag==0){if(prn>=41&&prn<=70)
				Elevation[rcvn][prn]=0.0;
		}
	}

}