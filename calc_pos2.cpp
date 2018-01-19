//////////////////////////////////////////////////////////////////////////
//
// 単独測位計算（1回目の計算で自身の時計誤差を修正したあと再計算）
//
//////////////////////////////////////////////////////////////////////////


///#include <iostream>
#include <stdio.h>
#include <math.h>
#include "global_extern.h"

void trans_coordinates(int);
void least_square(int,double[],int,int);
int minver(double[],double);
#define PI 3.1415926535897932384626433832795

const double cs     = 299792458.0;//光速
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//円周率
const double f      = -4.442807633e-10;//
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;

//測位演算部分のトップ
void calc_pos2(int rcvn, int iter_main, int pos)
{
	int i,j;
	int iter = 0;
	double init[7];
	double ref[3]={0};

	init[0] = 0;//X（ECEF）の初期値
	init[1] = 0;//Y（ECEF）の初期値
	init[2] = 0;//Z（ECEF）の初期値
//	init[0] = -3960000.0;//X（ECEF）の初期値
//	init[1] = 3350000.0;//Y（ECEF）の初期値
//	init[2] = 3698000.0;//Z（ECEF）の初期値
	init[3] = 0.0;//受信機クロック誤差
	init[4] = 0.0;//GPSと1番目の測位システムとの時計差
	init[5] = 0.0;//GPSと2番目の測位システムとの時計差
	init[6] = 0.0;//GPSと3番目の測位システムとの時計差
	LS_update[0]=10;LS_update[1]=10;LS_update[2]=10;

	if(SATn[rcvn]>=MinSvNum){
		while(fabs(LS_update[0])>0.001 || fabs(LS_update[1])>0.001 || fabs(LS_update[2])>0.001)
//		while(iter<=6)
		{
			least_square(rcvn, init, iter, pos);
			iter = iter + 1;
		}
	}
	else{
		printf("%f,***4衛星未満***\n",GPSTIME);
	}

	//測位演算結果を取得
	POSx[rcvn] = init[0];
	POSy[rcvn] = init[1];
	POSz[rcvn] = init[2];
	Clock[rcvn] = init[3];

	//受信機時計の誤差分を補正する
	for(i=0;i<SATn[rcvn];i++){
		j=SVn[rcvn][i];
		Pr1[rcvn][j]=Pr1[rcvn][j]-Clock[rcvn];

		if(rcvn==1){//基準局と移動局の遅延考慮（特にRTK）で利用
			Pr1ref[j] = Pr1[rcvn][j];
			Cp1ref[j] = Cp1[rcvn][j];
		}
	}

//3次元地心座標から楕円体座標への変換//
	double a,b,e,f,n,oldt,p,pi;
	double s[3],t[3];
	
	s[0] = POSx[rcvn];
	s[1] = POSy[rcvn];
	s[2] = POSz[rcvn];

	a = 6378137;
	pi     = 3.1415926535898;
	f = 1/298.257223563;

	b = a*(1-f);
	p = sqrt(s[0]*s[0]+s[1]*s[1]);
	e = f*(2-f);
	t[0] = atan(s[2]/((1-e)*p));
	oldt = 1;
	t[2]=0;

	while(fabs(t[2]-oldt) > 0.000001)
	{
		oldt = t[2];
		n = a/sqrt(1-(e*sin(t[0])*sin(t[0])));
		t[2] = p/cos(t[0]) - n;
		t[0] = atan((s[2]/p)/((1-e*n/(n+t[2]))));
	}
	t[1] = 2*atan(s[1]/(s[0]+sqrt(s[0]*s[0]+s[1]*s[1])));
//////////////////////////////////////////////////////////////
	t[0] = t[0]*180/pi;
	t[1] = t[1]*180/pi;

	double x,y;//順番に緯度方向、経度方向の誤差
	
	//基準側の単独測位結果ファイル出力
	if(rcvn==1 && POS==1){

		x = (t[0]-POSrcvlat[1])*111319.49;
		y = (t[1]-POSrcvlon[1])*cos(t[0]*PI/180.0)*111319.49;
/*
		fprintf(fp[3],"%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f",
		GPSTIME,SATn[rcvn],y,x,t[2]-POSrcvhgt[1],t[0],t[1],t[2],HDOP,VDOP,MinSvNum,
		Clock_ext[rcvn],Clock5[rcvn],Clock6[rcvn],Clock7[rcvn]);

		for(i=0;i<SATn[rcvn];i++){
			fprintf(fp[3],",%d",SVn[rcvn][i]);
		}
		fprintf(fp[3],"\n");
		*/
	}

	//移動側の単独測位結果ファイル出力
	if(rcvn==0 && POS==1){
		x = (t[0]-POSrcvlat[1])*111319.49;
		y = (t[1]-POSrcvlon[1])*cos(t[0]*PI/180.0)*111319.49;
///*
		fprintf(fp[3],"%f,%d,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%f,%f",
		GPSTIME,SATn[rcvn],y,x,t[2]-POSrcvhgt[1],t[0],t[1],t[2],HDOP,VDOP,MinSvNum,
		Clock_ext[rcvn],Clock5[rcvn],Clock6[rcvn],Clock7[rcvn]);

		for(i=0;i<SATn[rcvn];i++){
			fprintf(fp[3],",%d",SVn[rcvn][i]);
		}
		fprintf(fp[3],"\n");
//		*/

		if(t[0]>-90 && t[0]<90){//緯度だけチェック
			User_lat = t[0];//電離層等の計算用
			User_lon = t[1];//電離層等の計算用
			User_hgt = t[2];
			User_pos[0] = POSx[rcvn];
			User_pos[1] = POSy[rcvn];
			User_pos[2] = POSz[rcvn];
		}

	}

}//calc_pos2の最後
