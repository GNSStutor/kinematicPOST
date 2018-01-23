////////////////////////////////////////////////////////////////
//
// 緯度、経度、楕円体高度を地球中心座標系のx,y,zに変換
//
////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>

#include "global_extern.h"

void trans_coordinates(int rcvn)
{
    double a,b,f,n,pi;
	double rcv_pos[3];
	pi = 3.14159265359;	
	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);
	n = a*a/sqrt(a*a*cos(POSrcvlat[rcvn]*pi/180)*
		cos(POSrcvlat[rcvn]*pi/180)+b*b*sin(POSrcvlat[rcvn]*pi/180)
		*sin(POSrcvlat[rcvn]*pi/180));

	rcv_pos[0] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*cos(POSrcvlon[rcvn]*pi/180);
	rcv_pos[1] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*sin(POSrcvlon[rcvn]*pi/180);
	rcv_pos[2] = (n*b*b/a/a+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		sin(POSrcvlat[rcvn]*pi/180);
}