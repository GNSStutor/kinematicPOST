//////////////////////////////////////////////////
//
//	NAVFILE と OBSFILE の時刻を比較する
//　エフェメリスを使用するタイミングを決定
//
//////////////////////////////////////////////////
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "global_extern.h"

extern double sub_gpstime[PRN][256];


int r_time_count = 0;
int ED_count[PRN] = {1};

void rinex_time(int rcvn)
{
	int i,j;
	int count = 0;
	double stock_toe[PRN]={0};
	
	if(r_time_count == 1){
		for(i=1;i<=100;i++){//GPS/QZS/GALILEO/BeiDou
			j = ED_count[i];
			if(sub_gpstime[i][j] == 0.0){;}
			else if(sub_gpstime[i][j] <= GPSTIME){//GPS時刻が航法メッセージ内の対象衛星の時刻に達した時点
				Sub_E.af0[i][0] = Sub_E.af0[i][j];
				Sub_E.af1[i][0] = Sub_E.af1[i][j];
				Sub_E.af2[i][0] = Sub_E.af2[i][j];
				Sub_E.accuracy[i][0] = Sub_E.accuracy[i][j];
				Sub_E.cic[i][0] = Sub_E.cic[i][j];
				Sub_E.cis[i][0] = Sub_E.cis[i][j];
				Sub_E.code[i][0] = Sub_E.code[i][j];
				Sub_E.crc[i][0] = Sub_E.crc[i][j];
				Sub_E.crs[i][0] = Sub_E.crs[i][j];
				Sub_E.cuc[i][0] = Sub_E.cuc[i][j];
				Sub_E.cus[i][0] = Sub_E.cus[i][j];
				Sub_E.di0[i][0] = Sub_E.di0[i][j];
				Sub_E.dn[i][0] = Sub_E.dn[i][j];
				Sub_E.domega0[i][0] = Sub_E.domega0[i][j];
				Sub_E.e[i][0] = Sub_E.e[i][j];
				Sub_E.health[i][0] = Sub_E.health[i][j];
				Sub_E.i0[i][0] = Sub_E.i0[i][j];
				Sub_E.interval[i][0] = Sub_E.interval[i][j];
				Sub_E.iodc[i][0] = Sub_E.iodc[i][j];
				Sub_E.iode[i][0] = Sub_E.iode[i][j];
				Sub_E.m0[i][0] = Sub_E.m0[i][j];
				Sub_E.omega[i][0] = Sub_E.omega[i][j];
				Sub_E.omega0[i][0] = Sub_E.omega0[i][j];
				Sub_E.pflag[i][0] = Sub_E.pflag[i][j];
				Sub_E.roota[i][0] = Sub_E.roota[i][j];
				Sub_E.tgd[i][0] = Sub_E.tgd[i][j];
				Sub_E.toc[i][0] = Sub_E.toc[i][j];
				Sub_E.toe[i][0] = Sub_E.toe[i][j];
				ED_count[i]++;
			}				
		}

	}
	if(r_time_count == 1){	
		int k;
		double tk[PRN][200];
		for(i=101;i<=130;i++){//GLONASS
			for(k=0;k<200;k++){
			
				tk[i][k] =GPSTIME - sub_gpstime[i][k]-LeapSecond;
				if(tk[i][k]<=900||sub_gpstime[i][k]==0)
				break;
			}
						
			tk[i][k+1] = GPSTIME-sub_gpstime[i][k+1]-LeapSecond;

			if(k!=0 && fabs(tk[i][k])>fabs(tk[i][k+1]))
				k=k+1;
			
			Sub_E.TauN[i][0] = Sub_E.TauN[i][k];
			Sub_E.GammaN[i][0] = Sub_E.GammaN[i][k];
			Sub_E.tk[i][0] = Sub_E.tk[i][k];
			Sub_E.Xp[i][0] = Sub_E.Xp[i][k];
			Sub_E.Xv[i][0] = Sub_E.Xv[i][k];
			Sub_E.Xa[i][0] = Sub_E.Xa[i][k];
			Sub_E.Yp[i][0] = Sub_E.Yp[i][k];
			Sub_E.Yv[i][0] = Sub_E.Yv[i][k];
			Sub_E.Ya[i][0] = Sub_E.Ya[i][k];
			Sub_E.Zp[i][0] = Sub_E.Zp[i][k];
			Sub_E.Zv[i][0] = Sub_E.Zv[i][k];
			Sub_E.Za[i][0] = Sub_E.Za[i][k];
			Sub_E.toc[i][0] = Sub_E.toc[i][k];
			Ephe.roota[i]=1.0;
		}
	}

	//実際に上記で使用タイミングを決定した衛星をそのまま変数に代入して計算できるようにする
	for(i=1;i<=100;i++){//GPS/QZS/GALILEO/BeiDou
		Ephe.af0[i] = Sub_E.af0[i][0];
		Ephe.af1[i] = Sub_E.af1[i][0];
		Ephe.af2[i] = Sub_E.af2[i][0];
		Ephe.accuracy[i] = Sub_E.accuracy[i][0];
		Ephe.cic[i] = Sub_E.cic[i][0];
		Ephe.cis[i] = Sub_E.cis[i][0];
		Ephe.code[i] = Sub_E.code[i][0];
		Ephe.crc[i] = Sub_E.crc[i][0];
		Ephe.crs[i] = Sub_E.crs[i][0];
		Ephe.cuc[i] = Sub_E.cuc[i][0];
		Ephe.cus[i] = Sub_E.cus[i][0];
		Ephe.di0[i] = Sub_E.di0[i][0];
		Ephe.dn[i] = Sub_E.dn[i][0];
		Ephe.domega0[i] = Sub_E.domega0[i][0];
		Ephe.e[i] = Sub_E.e[i][0];
		Ephe.health[i] = Sub_E.health[i][0];
		Ephe.i0[i] = Sub_E.i0[i][0];
		Ephe.interval[i] = Sub_E.interval[i][0];
		Ephe.iodc[i] = Sub_E.iodc[i][0];
		Ephe.iode[i] = Sub_E.iode[i][0];
		Ephe.m0[i] = Sub_E.m0[i][0];
		Ephe.omega[i] = Sub_E.omega[i][0];
		Ephe.omega0[i] = Sub_E.omega0[i][0];
		Ephe.pflag[i] = Sub_E.pflag[i][0];
		Ephe.roota[i] = Sub_E.roota[i][0];
		Ephe.tgd[i] = Sub_E.tgd[i][0];
		Ephe.toc[i] = Sub_E.toc[i][0];
		Ephe.toe[i] = Sub_E.toe[i][0];
		stock_toe[i] = Ephe.toe[i];
		if(fabs(stock_toe[i])<0.1)
			stock_toe[i] = 604800;
		r_time_count=1;
		count++;

	}
	for(i=101;i<=130;i++){//GLONASS
			Ephe.TauN[i] = Sub_E.TauN[i][0];
			Ephe.GammaN[i] = Sub_E.GammaN[i][0];
			Ephe.tk[i] = Sub_E.tk[i][0];
			Ephe.Xp[i] = Sub_E.Xp[i][0];
			Ephe.Xv[i] = Sub_E.Xv[i][0];
			Ephe.Xa[i] = Sub_E.Xa[i][0];
			Ephe.Yp[i] = Sub_E.Yp[i][0];
			Ephe.Yv[i] = Sub_E.Yv[i][0];
			Ephe.Ya[i] = Sub_E.Ya[i][0];
			Ephe.Zp[i] = Sub_E.Zp[i][0];
			Ephe.Zv[i] = Sub_E.Zv[i][0];
			Ephe.Za[i] = Sub_E.Za[i][0];
			Ephe.toc[i] = Sub_E.toc[i][0];
			Ephe.roota[i]=1.0;
	}
	
	
	//Ephemeris Life Determination
	for(i=1;i<=39;i++){
		if(stock_toe[i]+7200 < GPSTIME)
			Ephe.roota[i]=0.0;//not used in positioing (calc_sat_pos.cpp)
	}
	for(i=101;i<=130;i++){
		if(fabs(Ephe.toc[i]-GPSTIME)>7200)
			Ephe.roota[i]=0.0;//not used in positioing (calc_sat_pos.cpp)
	}

}
	