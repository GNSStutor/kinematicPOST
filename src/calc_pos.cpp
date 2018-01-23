//////////////////////////////////////////////////////////////////////////
//
// 単独測位計算
//
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <stdio.h>
#include <math.h>
#include "global_extern.h"

void trans_coordinates(int);
void least_square(int,double[],int,int);
#define PI 3.1415926535897932384626433832795

const double cs     = 299792458.0;//光速
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//円周率
const double f      = -4.442807633e-10;//
const double F1     = 1575420000.0;
const double B1     = 1561098000.0;

//測位演算部分のトップ
void calc_pos(int rcvn, int iter_main, int pos)
{
	int i,j;
	int iter = 0;
	static int iter_pos = 1;
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
			iter = iter+1;
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
	Clock5[rcvn] = init[4];
	Clock6[rcvn] = init[5];
	Clock7[rcvn] = init[6];
	Clock_ext[rcvn] = Clock[rcvn];

	//受信機時計の誤差分を補正する
	for(i=0;i<SATn[rcvn];i++){
		j=SVn[rcvn][i];

		Pr1[rcvn][j]=Pr1[rcvn][j]-Clock[rcvn];

		//ここではRTK用に衛星時計側のクロック変動の影響を補正する
		//基準局と移動局で補正データの遅延（AGE）があるので、それを考慮
		if(rcvn==0){
			Pr1[1][j] = Pr1ref[j] - cs*(DGPSTIME-GPSTIME)*Ephe.af1[j];

			if(j>=1 && j<=70)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.af1[j]*F1/cs);
			if(j>=71 && j<=100)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.af1[j]*B1/cs);
			if(j>=101 && j<=135)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.GammaN[j]*GF1[j]/cs);

		}

		if(j<=70 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/F1);
		if(j<=70 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/F1);

		if(j>=71 && j<=100 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/B1);
		if(j>=71 && j<=100 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/B1);

		if(j>=101 && j<=130 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/GF1[j]);
		if(j>=101 && j<=130 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/GF1[j]);
	}

	if(rcvn==1)
		GPSTIME=GPSTIME-Clock[rcvn]/cs;
	if(rcvn==0)
		DGPSTIME=DGPSTIME-Clock[rcvn]/cs;

}
