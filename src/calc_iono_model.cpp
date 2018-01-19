///////////////////////////////////////////////////////////////////////
//
// クロバッチャモデルを利用して電離層遅延量を計算
// 他国のシステムも周波数を考慮することでそのまま計算
//
///////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include "global_extern.h"

#define rad(x) ((x)*3.14159265359/180)

const double F1     = 1575420000.0;
const double B1     = 1561098000.0;

void calc_iono_model(int rcvn)
{
	double earth_centered_angle[PRN];//((degree)
	double subionospheric_latitud[PRN];//((degree)
	double subionospheric_longitude[PRN];//(degree)
	double geomagnetic_latitud[PRN];//(degree)
	double local_time[PRN];
	double slant_factor[PRN];
	double first_computing[PRN];
	double ionospheric_time[PRN];
	double geodetic_latitud;////(degree)
	double geodetic_longitud;////(degree)
	double pai=3.14159265359;
	double siguma_b;	
	double siguma_a;
	double abc;
	int j,i;
	int k;
	double deg_c= pai /180;

	if(rcvn==1){//基準局
		geodetic_latitud = POSrcvlat[rcvn];
		geodetic_longitud = POSrcvlon[rcvn];
	}
	else{//移動側
		geodetic_latitud = User_lat;
		geodetic_longitud = User_lon;
	}

	for(i=0;i<SATn[rcvn];i++){
		
		j=SVn[rcvn][i];
		siguma_b=0.0;
		siguma_a=0.0;

		earth_centered_angle[j]=(0.0137/((Elevation[rcvn][j]/180)+0.11))-0.022;//(degree)

		subionospheric_latitud[j] = geodetic_latitud /180+ 
			earth_centered_angle[j] * cos(Azimuth[rcvn][j] * deg_c);

				if(subionospheric_latitud[j] > 0.416){
					subionospheric_latitud[j] = 0.416;
				}
				if(subionospheric_latitud[j] < -0.416){
					subionospheric_latitud[j] = -0.416;
				}

		subionospheric_longitude[j] = geodetic_longitud/180 + 
			((earth_centered_angle[j] * sin(Azimuth[rcvn][j] * deg_c)) 
			/ (cos(subionospheric_latitud[j] * 180*pai/180)));

		geomagnetic_latitud[j] = subionospheric_latitud[j] + 
			(0.064 * cos((subionospheric_longitude[j]*180*pai/180 - 1.617)));

		local_time[j] = 4.32 * 10000 * subionospheric_longitude[j] 
			+ GPSTIME;
			//+ GPStime;
		while(local_time[j] > 86400){ 		
			local_time[j] = local_time[j] -86400;
		}
		while(local_time[j] < 0){
			local_time[j] = local_time[j] +86400;
		}

//slant factor
		slant_factor[j] = 1.0 + 16.0 * pow((0.53- (Elevation[rcvn][j]/ 180) ),3.0);
	
		for(k=0;k<=3;k++){
			siguma_b = siguma_b + Bet[k] * pow(geomagnetic_latitud[j],k);
						

		}
		if(siguma_b < 72000){
			siguma_b = 72000;
		}

		first_computing[j] = (2*pai*(local_time[j]-50400)) / siguma_b;

		if(fabs(first_computing[j]) >= 1.57 ){
			ionospheric_time[j] = slant_factor[j] * (5 * pow(10.0,-9.0));
		}
		else{
			for(k=0;k<=3;k++){
				siguma_a = siguma_a + Arp[k] *
					pow(geomagnetic_latitud[j],k);
			}
			if(siguma_a < 0){
				siguma_a=0;
			}

			abc = (1.0 - (pow(first_computing[j],2.0) / 2.0)
				+ (pow(first_computing[j],4.0) / 24.0));

			ionospheric_time[j] = slant_factor[j] * ((5 * pow(10.0,-9.0))
				+ siguma_a * abc);
		}

		//GPS/QZS/GALILEOは同じ1575.42MHz
		Iono[rcvn][j] = ionospheric_time[j] * 299792458;

		if(j>=71 && j<=100)//for BeiDou
			Iono[rcvn][j] = Iono[rcvn][j]*(F1*F1)/(B1*B1);

		if(j>=101 && j<=130)//for GLONASS
			Iono[rcvn][j] = Iono[rcvn][j]*(F1*F1)/(GF1[j]*GF1[j]);

	}
}