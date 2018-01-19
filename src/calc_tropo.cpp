///////////////////////////////////////////////////////////////////////
//
// 対流圏遅延量の計算（Saastamoinen Model）
//
///////////////////////////////////////////////////////////////////////
 
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <process.h>

#include "global_extern.h"

#define rad(x) ((x)*3.14159265359/180)

void calc_tropo(int rcvn)
{
	//代表的な気象条件
	double ps1=1000.0;
	double temp=15.0 + 273.16;
	double ps2=38.03*0.77;
	
	double tropo_d[RCVN][PRN],tropo_w[RCVN][PRN];
	double map_d[PRN],map_w[PRN];
	double pi = 3.14159265359;
	int j,k;

	for(k=1;k<=PRN-1;k++)
		Tropo[rcvn][k]=0.0;

	for(k=0;k<SATn[rcvn];k++){
		j=SVn[rcvn][k];
		
		//Saastamoinen modelを利用
		if(rcvn==1)
			tropo_d[rcvn][j] = 0.002277*(1.0+0.0026*cos(2*POSrcvlat[rcvn]*pi/180.0)+0.00028*POSrcvhgt[rcvn])*ps1;
		else
			tropo_d[rcvn][j] = 0.002277*(1.0+0.0026*cos(2*User_lat*pi/180.0)+0.00028*User_hgt)*ps1;
		tropo_w[rcvn][j] = 0.002277*(1255/temp+0.05)*ps2;

		map_d[j] = 0.00143/(tan(Elevation[rcvn][j]*pi/180)+0.0445);
		map_d[j] = 1/(sin(Elevation[rcvn][j]*pi/180)+map_d[j]);
		map_w[j] = 0.00035/(tan(Elevation[rcvn][j]*pi/180)+0.017);
		map_w[j] = 1/(sin(Elevation[rcvn][j]*pi/180)+map_d[j]);

		Tropo[rcvn][j] = tropo_d[rcvn][j]*map_d[j] + tropo_w[rcvn][j]*map_w[j];
	}
}