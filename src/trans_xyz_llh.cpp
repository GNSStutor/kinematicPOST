/////////////////////////////////////////////////////
//
// 地球中心座標系のx,y,zを緯度、経度、楕円体高度に変換
//
/////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>

#define pi 3.1415926535898

void trans_xyz_llh(double a1,double a2,double a3, double xx[])
{
	double p[3],f,a,b;
	double b1;
	double stock_p = 1000000000;

	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);

	p[0] = atan2(a2,a1);
	b1 = sqrt(a1*a1+a2*a2);
	p[1] = atan2(a3,b1);
	p[1] = atan2(a3+a*f*(2-f)*sin(p[1])/
				sqrt(1-f*(2-f)*sin(p[1])*sin(p[1])),b1);

	while(fabs(stock_p-p[1]) >= 0.000000001){
		stock_p = p[1];
		p[1] = atan2(a3+a*f*(2-f)*sin(stock_p)/
			sqrt(1-f*(2-f)*sin(stock_p)*sin(stock_p)),b1);
			//cout << stock_p << " " << p[1] << endl;
	}

	xx[2] = b1/cos(p[1])-a/sqrt(1-f*(2-f)*sin(p[1])*sin(p[1]));
	xx[1] = p[0]*180/pi;
	xx[0] = p[1]*180/pi;
}