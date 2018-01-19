//////////////////////////////////////////////////////////////////
//
// 最小二乗法による単独測位
//
/////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stdlib.h>

#include "global_extern.h"

using namespace std;

int minv(double[], double, int);
void trans_xyz_llh(double,double,double,double[]);

const double pi = 3.1415926535898;//円周率
const double cs = 299792458.0;//光速

							  //測位計算用の反復//
							  //--- iteration ---//
void least_square(int rcvn, double init[4], int iter, int pos)
{
	int i, prn;
	int stock_sat[PRN] = { 0 };
	static int j = 0, m = 0;
	double r2[PRN], r3[PRN], a[PRN] = { 0.0 }, b[PRN] = { 0.0 }, c[PRN] = { 0.0 }, e[PRN] = { 0.0 }, f[PRN] = { 0.0 }, g[PRN] = { 0.0 }, a2[PRN], a3[PRN], a4[PRN],
		aa[PRN], bb[PRN], cc[PRN], dd[PRN], ee[PRN], ff[PRN], gg[PRN], delta[8], b3[PRN] = { 0.0 }, v[5] = { 0.0 }, v2[5] = { 0.0 };
	double lat, lon;
	double m1[9] = { 0 }, m2[9] = { 0 };
	double xx[64] = { 0 }, xx2[64] = { 0 };



	//DOP計算用に自身の位置（緯度、経度）を入力　地心座標系から測地座標で利用
	lat = (POSrcvlat[rcvn])*2.0*pi / 360.0;
	lon = POSrcvlon[rcvn] * 2.0*pi / 360.0;



	for (i = 0; i<SATn[rcvn]; i++) {
		prn = SVn[rcvn][i];

		r2[i] = sqrt((SVx[rcvn][prn] - init[0])*(SVx[rcvn][prn] - init[0])
			+ (SVy[rcvn][prn] - init[1])*(SVy[rcvn][prn] - init[1])
			+ (SVz[rcvn][prn] - init[2])*(SVz[rcvn][prn] - init[2]));

		//単独測位
		if (POS == 1) {
			r3[i] = Pr1[rcvn][prn] + SV_corrtime[rcvn][prn] -
				Iono[rcvn][prn] - Tropo[rcvn][prn] - r2[i];
		}

		//DGNSS
//		if (POS == 2) {
//			r3[i] = Pr1[rcvn][prn] + SV_corrtime[rcvn][prn] +
//				Correction[prn] - r2[i];
//		}

		//行列Gの生成(視線方向ベクトル:-a[i],-b[i],-c[i])
		a[i] = -(SVx[rcvn][prn] - init[0]) / r2[i];
		b[i] = -(SVy[rcvn][prn] - init[1]) / r2[i];
		c[i] = -(SVz[rcvn][prn] - init[2]) / r2[i];

		//最小二乗法の重みを生成　基本仰角依存のシンプルなもの
		SV_omomi[rcvn][prn] = 1.0*1.0*(90.0 / Elevation[rcvn][prn]);

		aa[i] = a[i] / SV_omomi[rcvn][prn];
		bb[i] = b[i] / SV_omomi[rcvn][prn];
		cc[i] = c[i] / SV_omomi[rcvn][prn];
		dd[i] = 1 / SV_omomi[rcvn][prn];

		if (MinSvNum == 5) {//GPS/QZSS + 他の測位システム1つ
			if (prn >= 41) {
				ee[i] = 1 / SV_omomi[rcvn][prn];
				e[i] = 1;
			}
		}

		if (MinSvNum == 6) {//GPS/QZSS + 他の測位システム2つ
							//測位システムごと　Systemflag[2]:GLONASS　Systemflag[3]:BeiDou　Systemflag[4]:GALILEO
			if (Systemflag[4] == 1 && Systemflag[3] == 1) {
				if (prn >= 41 && prn <= 70) {
					ff[i] = 1 / SV_omomi[rcvn][prn];//system_time//GALILLEO
					f[i] = 1;
				}
				else if (prn >= 71 && prn <= 100) {
					ee[i] = 1 / SV_omomi[rcvn][prn];//system_time;//BeiDou
					e[i] = 1;
				}
			}
			else if (Systemflag[4] == 1 && Systemflag[2] == 1) {
				if (prn >= 41 && prn <= 70) {
					ff[i] = 1 / SV_omomi[rcvn][prn];//system_time//GALILLEO
					f[i] = 1;
				}
				else if (prn >= 101 && prn <= 130) {
					ee[i] = 1 / SV_omomi[rcvn][prn];//system_time;//GLONASS
					e[i] = 1;
				}
			}
			else if (Systemflag[3] == 1 && Systemflag[2] == 1) {//GPS
				if (prn >= 71 && prn <= 100) {
					ff[i] = 1 / SV_omomi[rcvn][prn];//system_time/BeiDo
					f[i] = 1;
				}
				else if (prn >= 101 && prn <= 130) {
					ee[i] = 1 / SV_omomi[rcvn][prn];//system_time;//GLONASS
					e[i] = 1;
				}
			}
		}

		if (MinSvNum == 7) {//GPS/QZSS + 他の測位システム3つ
							//測位システムごと　Systemflag[2]:GLONASS　Systemflag[3]:BeiDou　Systemflag[4]:GALILEO
			if (prn >= 41 && prn <= 70) {
				gg[i] = 1 / SV_omomi[rcvn][prn];//system_time//GALILLEO
				g[i] = 1;
			}
			else if (prn >= 71 && prn <= 100) {
				ff[i] = 1 / SV_omomi[rcvn][prn];//system_time;//BeiDo
				f[i] = 1;
			}
			else if (prn >= 101 && prn <= 130) {
				ee[i] = 1 / SV_omomi[rcvn][prn];//system_time;//GLONASS
				e[i] = 1;
			}
		}
	}

	for (i = 0; i<64; i++) {
		a2[i] = 0.0;
		a4[i] = 0.0;
	}




	//GPS/QZSS　＋　他の測位衛星システム毎　最小二乗法の未知数において、GPS時刻との差分を求める必要がある

	//GPS/QZSS
	if (MinSvNum == 4) {

		//行列Gの転置*Gの計算(重み考慮している場合はGの転置とGの間に重み行列)
		for (i = 0; i<SATn[rcvn]; i++) {
			a2[0] += aa[i] * a[i];
			a2[1] += aa[i] * b[i];
			a2[2] += aa[i] * c[i];
			a2[3] += aa[i];
			a2[4] += bb[i] * a[i];
			a2[5] += bb[i] * b[i];
			a2[6] += bb[i] * c[i];
			a2[7] += bb[i];
			a2[8] += cc[i] * a[i];
			a2[9] += cc[i] * b[i];
			a2[10] += cc[i] * c[i];
			a2[11] += cc[i];
			a2[12] += dd[i] * a[i];
			a2[13] += dd[i] * b[i];
			a2[14] += dd[i] * c[i];
			a2[15] += dd[i];

			a4[0] += a[i] * a[i];
			a4[1] += a[i] * b[i];
			a4[2] += a[i] * c[i];
			a4[3] += a[i];
			a4[4] += b[i] * a[i];
			a4[5] += b[i] * b[i];
			a4[6] += b[i] * c[i];
			a4[7] += b[i];
			a4[8] += c[i] * a[i];
			a4[9] += c[i] * b[i];
			a4[10] += c[i] * c[i];
			a4[11] += c[i];
			a4[12] += a[i];
			a4[13] += b[i];
			a4[14] += c[i];
			a4[15] += 1;
		}

		double eps = 1.0e-25;

		//DOP計算用　重みをすべて同じにする（HDOP,VDOP,PDOP,GDOP,TDOPの計算）
		for (i = 0; i <= 15; i++)
			xx2[i] = a4[i];
		//DOPのための逆行列計算//
		minv(xx2, eps, MinSvNum);

		//ここは地心座標系での値なので、ローカル座標系への変換　
		m2[0] = (-sin(lat)*cos(lon))*xx2[0] + (-sin(lat)*sin(lon))*xx2[4] + cos(lat)*xx2[8];
		m2[1] = (-sin(lat)*cos(lon))*xx2[1] + (-sin(lat)*sin(lon))*xx2[5] + cos(lat)*xx2[9];
		m2[2] = (-sin(lat)*cos(lon))*xx2[2] + (-sin(lat)*sin(lon))*xx2[6] + cos(lat)*xx2[10];
		m2[3] = (-sin(lon))*xx2[0] + cos(lon)*xx2[4];
		m2[4] = (-sin(lon))*xx2[1] + cos(lon)*xx2[5];
		m2[5] = (-sin(lon))*xx2[2] + cos(lon)*xx2[6];
		m2[6] = cos(lat)*cos(lon)*xx2[0] + cos(lat)*sin(lon)*xx2[4] + sin(lat)*xx2[8];
		m2[7] = cos(lat)*cos(lon)*xx2[1] + cos(lat)*sin(lon)*xx2[5] + sin(lat)*xx2[9];
		m2[8] = cos(lat)*cos(lon)*xx2[2] + cos(lat)*sin(lon)*xx2[6] + sin(lat)*xx2[10];

		m1[0] = m2[0] * (-sin(lat)*cos(lon)) + m2[1] * (-sin(lat)*sin(lon)) + m2[2] * cos(lat);
		m1[1] = m2[0] * (-sin(lon)) + m2[1] * cos(lon);
		m1[2] = m2[0] * cos(lat)*cos(lon) + m2[1] * cos(lat)*sin(lon) + m2[2] * sin(lat);
		m1[3] = m2[3] * (-sin(lat)*cos(lon)) + m2[4] * (-sin(lat)*sin(lon)) + m2[5] * cos(lat);
		m1[4] = m2[3] * (-sin(lon)) + m2[4] * cos(lon);
		m1[5] = m2[3] * cos(lat)*cos(lon) + m2[4] * cos(lat)*sin(lon) + m2[5] * sin(lat);
		m1[6] = m2[6] * (-sin(lat)*cos(lon)) + m2[7] * (-sin(lat)*sin(lon)) + m2[8] * cos(lat);
		m1[7] = m2[6] * (-sin(lon)) + m2[7] * cos(lon);
		m1[8] = m2[6] * cos(lat)*cos(lon) + m2[7] * cos(lat)*sin(lon) + m2[8] * sin(lat);

		//DOP値をアウトプット
		HDOP = sqrt(m1[0] + m1[4]);
		VDOP = sqrt(m1[8]);
		PDOP = sqrt(xx2[0] + xx2[5] + xx2[10]);
		TDOP = sqrt(xx2[15]);

		//単独測位計算用　重みを考慮する
		for (i = 0; i <= 15; i++)
			xx[i] = a2[i];
		//単独測位のための逆行列計算//
		minv(xx, eps, MinSvNum);

		//推定される値：X,Y,Z軸方向に移動させる量（dx,dy,dz）とクロック誤差
		delta[0] = 0.0;
		delta[1] = 0.0;
		delta[2] = 0.0;
		delta[3] = 0.0;

		//最後にGの転置（重み含む）と上で求めたr3[]をかける
		for (i = 0; i<SATn[rcvn]; i++) {
			a3[i] = xx[0] * aa[i] + xx[1] * bb[i] + xx[2] * cc[i] + xx[3] * dd[i];
			a3[i + SATn[rcvn]] = xx[4] * aa[i] + xx[5] * bb[i] + xx[6] * cc[i] + xx[7] * dd[i];
			a3[i + SATn[rcvn] * 2] = xx[8] * aa[i] + xx[9] * bb[i] + xx[10] * cc[i] + xx[11] * dd[i];
			a3[i + SATn[rcvn] * 3] = xx[12] * aa[i] + xx[13] * bb[i] + xx[14] * cc[i] + xx[15] * dd[i];
			delta[0] += a3[i] * r3[i];
			delta[1] += a3[i + SATn[rcvn]] * r3[i];
			delta[2] += a3[i + SATn[rcvn] * 2] * r3[i];
			delta[3] += a3[i + SATn[rcvn] * 3] * r3[i];
		}

		//初期値に上記で計算されたdx,dy,dzを更新
		for (i = 0; i<3; i++) {
			init[i] += delta[i];
		}
		init[3] = delta[3];//時計誤差は毎回新たに計算される
		LS_update[0]=delta[0];LS_update[1]=delta[1];LS_update[2]=delta[2];//最小二乗法の修正量が十分ちいさいか

	}//GPS/QZSS

	 //GPS/QZSS+1つの他の測位衛星
	if (MinSvNum == 5) {
		//行列Gの転置*Gの計算
		for (i = 0; i<SATn[rcvn]; i++) {
			prn = SVn[rcvn][i];

			a2[0] += aa[i] * a[i];
			a2[1] += aa[i] * b[i];
			a2[2] += aa[i] * c[i];
			a2[3] += aa[i] * 1;//clock
			a2[4] += aa[i] * e[i];

			a2[5] += bb[i] * a[i];
			a2[6] += bb[i] * b[i];
			a2[7] += bb[i] * c[i];
			a2[8] += bb[i];
			a2[9] += bb[i] * e[i];

			a2[10] += cc[i] * a[i];
			a2[11] += cc[i] * b[i];
			a2[12] += cc[i] * c[i];
			a2[13] += cc[i];
			a2[14] += cc[i] * e[i];

			a2[15] += dd[i] * a[i];
			a2[16] += dd[i] * b[i];
			a2[17] += dd[i] * c[i];
			a2[18] += dd[i];
			a2[19] += dd[i] * e[i];

			a2[20] += ee[i] * a[i];
			a2[21] += ee[i] * b[i];
			a2[22] += ee[i] * c[i];
			a2[23] += ee[i];
			a2[24] += ee[i] * e[i];

			a4[0] += a[i] * a[i];
			a4[1] += a[i] * b[i];
			a4[2] += a[i] * c[i];
			a4[3] += a[i];
			a4[4] += a[i] * e[i];
			a4[5] += b[i] * a[i];
			a4[6] += b[i] * b[i];
			a4[7] += b[i] * c[i];
			a4[8] += b[i];
			a4[9] += b[i] * e[i];
			a4[10] += c[i] * a[i];
			a4[11] += c[i] * b[i];
			a4[12] += c[i] * c[i];
			a4[13] += c[i];
			a4[14] += c[i] * e[i];
			a4[15] += a[i];
			a4[16] += b[i];
			a4[17] += c[i];
			a4[18] += 1;
			a4[19] += e[i];
			a4[20] += e[i] * a[i];
			a4[21] += e[i] * b[i];
			a4[22] += e[i] * c[i];
			a4[23] += e[i];
			a4[24] += e[i] * e[i];
		}

		double eps = 1.0e-25;

		//DOP計算用　重みをすべて同じにする（HDOP,VDOP,PDOP,GDOP,TDOPの計算）
		for (i = 0; i <= 24; i++)
			xx2[i] = a4[i];
		//DOPのための逆行列計算//
		minv(xx2, eps, MinSvNum);

		//ここは地心座標系での値なので、ローカル座標系への変換　
		m2[0] = (-sin(lat)*cos(lon))*xx2[0] + (-sin(lat)*sin(lon))*xx2[5] + cos(lat)*xx2[10];
		m2[1] = (-sin(lat)*cos(lon))*xx2[1] + (-sin(lat)*sin(lon))*xx2[6] + cos(lat)*xx2[11];
		m2[2] = (-sin(lat)*cos(lon))*xx2[2] + (-sin(lat)*sin(lon))*xx2[7] + cos(lat)*xx2[12];
		m2[3] = (-sin(lon))*xx2[0] + cos(lon)*xx2[5];
		m2[4] = (-sin(lon))*xx2[1] + cos(lon)*xx2[6];
		m2[5] = (-sin(lon))*xx2[2] + cos(lon)*xx2[7];
		m2[6] = cos(lat)*cos(lon)*xx2[0] + cos(lat)*sin(lon)*xx2[5] + sin(lat)*xx2[10];
		m2[7] = cos(lat)*cos(lon)*xx2[1] + cos(lat)*sin(lon)*xx2[6] + sin(lat)*xx2[11];
		m2[8] = cos(lat)*cos(lon)*xx2[2] + cos(lat)*sin(lon)*xx2[7] + sin(lat)*xx2[12];

		m1[0] = m2[0] * (-sin(lat)*cos(lon)) + m2[1] * (-sin(lat)*sin(lon)) + m2[2] * cos(lat);
		m1[1] = m2[0] * (-sin(lon)) + m2[1] * cos(lon);
		m1[2] = m2[0] * cos(lat)*cos(lon) + m2[1] * cos(lat)*sin(lon) + m2[2] * sin(lat);
		m1[3] = m2[3] * (-sin(lat)*cos(lon)) + m2[4] * (-sin(lat)*sin(lon)) + m2[5] * cos(lat);
		m1[4] = m2[3] * (-sin(lon)) + m2[4] * cos(lon);
		m1[5] = m2[3] * cos(lat)*cos(lon) + m2[4] * cos(lat)*sin(lon) + m2[5] * sin(lat);
		m1[6] = m2[6] * (-sin(lat)*cos(lon)) + m2[7] * (-sin(lat)*sin(lon)) + m2[8] * cos(lat);
		m1[7] = m2[6] * (-sin(lon)) + m2[7] * cos(lon);
		m1[8] = m2[6] * cos(lat)*cos(lon) + m2[7] * cos(lat)*sin(lon) + m2[8] * sin(lat);

		//DOP値をアウトプット
		HDOP = sqrt(m1[0] + m1[4]);
		VDOP = sqrt(m1[8]);
		PDOP = sqrt(xx2[0] + xx2[6] + xx2[12]);
		TDOP = sqrt(xx2[18]);



		for (i = 0; i <= 24; i++)
			xx[i] = a2[i];
		//逆行列計算//
		minv(xx, eps, MinSvNum);

		delta[0] = 0.0;
		delta[1] = 0.0;
		delta[2] = 0.0;
		delta[3] = 0.0;
		delta[4] = 0.0;

		for (i = 0; i<SATn[rcvn]; i++) {
			a3[i] = xx[0] * aa[i] + xx[1] * bb[i] + xx[2] * cc[i] + xx[3] * dd[i] + xx[4] * ee[i];
			a3[i + SATn[rcvn]] = xx[5] * aa[i] + xx[6] * bb[i] + xx[7] * cc[i] + xx[8] * dd[i] + xx[9] * ee[i];
			a3[i + SATn[rcvn] * 2] = xx[10] * aa[i] + xx[11] * bb[i] + xx[12] * cc[i] + xx[13] * dd[i] + xx[14] * ee[i];
			a3[i + SATn[rcvn] * 3] = xx[15] * aa[i] + xx[16] * bb[i] + xx[17] * cc[i] + xx[18] * dd[i] + xx[19] * ee[i];
			a3[i + SATn[rcvn] * 4] = xx[20] * aa[i] + xx[21] * bb[i] + xx[22] * cc[i] + xx[23] * dd[i] + xx[24] * ee[i];
			delta[0] += a3[i] * r3[i];
			delta[1] += a3[i + SATn[rcvn]] * r3[i];
			delta[2] += a3[i + SATn[rcvn] * 2] * r3[i];
			delta[3] += a3[i + SATn[rcvn] * 3] * r3[i];
			delta[4] += a3[i + SATn[rcvn] * 4] * r3[i];
		}

		for (i = 0; i<3; i++) {
			init[i] += delta[i];
		}
		init[3] = delta[3];
		init[4] = delta[4];
		LS_update[0]=delta[0];LS_update[1]=delta[1];LS_update[2]=delta[2];//最小二乗法の修正量が十分ちいさいか

	}//GPS/QZSS+1つの他の測位衛星

	 //GPS/QZSS+2つの他の測位衛星
	if (MinSvNum == 6) {
		//行列Gの転置*Gの計算
		for (i = 0; i<SATn[rcvn]; i++) {
			prn = SVn[rcvn][i];

			a2[0] += aa[i] * a[i];
			a2[1] += aa[i] * b[i];
			a2[2] += aa[i] * c[i];
			a2[3] += aa[i] * 1;//clock
			a2[4] += aa[i] * e[i];
			a2[5] += aa[i] * f[i];

			a2[6] += bb[i] * a[i];
			a2[7] += bb[i] * b[i];
			a2[8] += bb[i] * c[i];
			a2[9] += bb[i];
			a2[10] += bb[i] * e[i];
			a2[11] += bb[i] * f[i];

			a2[12] += cc[i] * a[i];
			a2[13] += cc[i] * b[i];
			a2[14] += cc[i] * c[i];
			a2[15] += cc[i];
			a2[16] += cc[i] * e[i];
			a2[17] += cc[i] * f[i];

			a2[18] += dd[i] * a[i];
			a2[19] += dd[i] * b[i];
			a2[20] += dd[i] * c[i];
			a2[21] += dd[i];
			a2[22] += dd[i] * e[i];
			a2[23] += dd[i] * f[i];

			a2[24] += ee[i] * a[i];
			a2[25] += ee[i] * b[i];
			a2[26] += ee[i] * c[i];
			a2[27] += ee[i];
			a2[28] += ee[i] * e[i];
			a2[29] += ee[i] * f[i];

			a2[30] += ff[i] * a[i];
			a2[31] += ff[i] * b[i];
			a2[32] += ff[i] * c[i];
			a2[33] += ff[i];
			a2[34] += ff[i] * e[i];
			a2[35] += ff[i] * f[i];

			a4[0] += a[i] * a[i];
			a4[1] += a[i] * b[i];
			a4[2] += a[i] * c[i];
			a4[3] += a[i];
			a4[4] += a[i] * e[i];
			a4[5] += a[i] * f[i];

			a4[6] += b[i] * a[i];
			a4[7] += b[i] * b[i];
			a4[8] += b[i] * c[i];
			a4[9] += b[i];
			a4[10] += b[i] * e[i];
			a4[11] += b[i] * f[i];

			a4[12] += c[i] * a[i];
			a4[13] += c[i] * b[i];
			a4[14] += c[i] * c[i];
			a4[15] += c[i];
			a4[16] += c[i] * e[i];
			a4[17] += c[i] * f[i];

			a4[18] += a[i];
			a4[19] += b[i];
			a4[20] += c[i];
			a4[21] += 1.0;
			a4[22] += e[i];
			a4[23] += f[i];

			a4[24] += e[i] * a[i];
			a4[25] += e[i] * b[i];
			a4[26] += e[i] * c[i];
			a4[27] += e[i];
			a4[28] += e[i] * e[i];
			a4[29] += e[i] * f[i];

			a4[30] += f[i] * a[i];
			a4[31] += f[i] * b[i];
			a4[32] += f[i] * c[i];
			a4[33] += f[i];
			a4[34] += f[i] * e[i];
			a4[35] += f[i] * f[i];
		}

		double eps = 1.0e-25;

		//DOP計算用　重みをすべて同じにする（HDOP,VDOP,PDOP,GDOP,TDOPの計算）
		for (i = 0; i <= 35; i++)
			xx2[i] = a4[i];
		//DOPのための逆行列計算//
		minv(xx2, eps, MinSvNum);

		//ここは地心座標系での値なので、ローカル座標系への変換　
		m2[0] = (-sin(lat)*cos(lon))*xx2[0] + (-sin(lat)*sin(lon))*xx2[6] + cos(lat)*xx2[12];
		m2[1] = (-sin(lat)*cos(lon))*xx2[1] + (-sin(lat)*sin(lon))*xx2[7] + cos(lat)*xx2[13];
		m2[2] = (-sin(lat)*cos(lon))*xx2[2] + (-sin(lat)*sin(lon))*xx2[8] + cos(lat)*xx2[14];
		m2[3] = (-sin(lon))*xx2[0] + cos(lon)*xx2[6];
		m2[4] = (-sin(lon))*xx2[1] + cos(lon)*xx2[7];
		m2[5] = (-sin(lon))*xx2[2] + cos(lon)*xx2[8];
		m2[6] = cos(lat)*cos(lon)*xx2[0] + cos(lat)*sin(lon)*xx2[6] + sin(lat)*xx2[12];
		m2[7] = cos(lat)*cos(lon)*xx2[1] + cos(lat)*sin(lon)*xx2[7] + sin(lat)*xx2[13];
		m2[8] = cos(lat)*cos(lon)*xx2[2] + cos(lat)*sin(lon)*xx2[8] + sin(lat)*xx2[14];

		m1[0] = m2[0] * (-sin(lat)*cos(lon)) + m2[1] * (-sin(lat)*sin(lon)) + m2[2] * cos(lat);
		m1[1] = m2[0] * (-sin(lon)) + m2[1] * cos(lon);
		m1[2] = m2[0] * cos(lat)*cos(lon) + m2[1] * cos(lat)*sin(lon) + m2[2] * sin(lat);
		m1[3] = m2[3] * (-sin(lat)*cos(lon)) + m2[4] * (-sin(lat)*sin(lon)) + m2[5] * cos(lat);
		m1[4] = m2[3] * (-sin(lon)) + m2[4] * cos(lon);
		m1[5] = m2[3] * cos(lat)*cos(lon) + m2[4] * cos(lat)*sin(lon) + m2[5] * sin(lat);
		m1[6] = m2[6] * (-sin(lat)*cos(lon)) + m2[7] * (-sin(lat)*sin(lon)) + m2[8] * cos(lat);
		m1[7] = m2[6] * (-sin(lon)) + m2[7] * cos(lon);
		m1[8] = m2[6] * cos(lat)*cos(lon) + m2[7] * cos(lat)*sin(lon) + m2[8] * sin(lat);

		//DOP値をアウトプット
		HDOP = sqrt(m1[0] + m1[4]);
		VDOP = sqrt(m1[8]);
		PDOP = sqrt(xx2[0] + xx2[7] + xx2[14]);
		TDOP = sqrt(xx2[21]);



		for (i = 0; i <= 35; i++)
			xx[i] = a2[i];
		//逆行列計算//
		minv(xx, eps, MinSvNum);

		delta[0] = 0.0;
		delta[1] = 0.0;
		delta[2] = 0.0;
		delta[3] = 0.0;
		delta[4] = 0.0;
		delta[5] = 0.0;

		for (i = 0; i<SATn[rcvn]; i++) {
			a3[i] = xx[0] * aa[i] + xx[1] * bb[i] + xx[2] * cc[i] + xx[3] * dd[i] + xx[4] * ee[i] + xx[5] * ff[i];
			a3[i + SATn[rcvn]] = xx[6] * aa[i] + xx[7] * bb[i] + xx[8] * cc[i] + xx[9] * dd[i] + xx[10] * ee[i] + xx[11] * ff[i];
			a3[i + SATn[rcvn] * 2] = xx[12] * aa[i] + xx[13] * bb[i] + xx[14] * cc[i] + xx[15] * dd[i] + xx[16] * ee[i] + xx[17] * ff[i];
			a3[i + SATn[rcvn] * 3] = xx[18] * aa[i] + xx[19] * bb[i] + xx[20] * cc[i] + xx[21] * dd[i] + xx[22] * ee[i] + xx[23] * ff[i];
			a3[i + SATn[rcvn] * 4] = xx[24] * aa[i] + xx[25] * bb[i] + xx[26] * cc[i] + xx[27] * dd[i] + xx[28] * ee[i] + xx[29] * ff[i];
			a3[i + SATn[rcvn] * 5] = xx[30] * aa[i] + xx[31] * bb[i] + xx[32] * cc[i] + xx[33] * dd[i] + xx[34] * ee[i] + xx[35] * ff[i];

			delta[0] += a3[i] * r3[i];
			delta[1] += a3[i + SATn[rcvn]] * r3[i];
			delta[2] += a3[i + SATn[rcvn] * 2] * r3[i];
			delta[3] += a3[i + SATn[rcvn] * 3] * r3[i];
			delta[4] += a3[i + SATn[rcvn] * 4] * r3[i];
			delta[5] += a3[i + SATn[rcvn] * 5] * r3[i];
		}

		for (i = 0; i<3; i++) {
			init[i] += delta[i];
		}
		init[3] = delta[3];
		init[4] = delta[4];
		init[5] = delta[5];
		LS_update[0]=delta[0];LS_update[1]=delta[1];LS_update[2]=delta[2];//最小二乗法の修正量が十分ちいさいか

	}//GPS/QZSS+2つの他の測位衛星

	 //GPS/QZSS+3つの他の測位衛星
	if (MinSvNum == 7) {
		//行列Gの転置*Gの計算
		for (i = 0; i<SATn[rcvn]; i++) {
			prn = SVn[rcvn][i];

			a2[0] += aa[i] * a[i];
			a2[1] += aa[i] * b[i];
			a2[2] += aa[i] * c[i];
			a2[3] += aa[i] * 1;//clock
			a2[4] += aa[i] * e[i];
			a2[5] += aa[i] * f[i];
			a2[6] += aa[i] * g[i];

			a2[7] += bb[i] * a[i];
			a2[8] += bb[i] * b[i];
			a2[9] += bb[i] * c[i];
			a2[10] += bb[i];
			a2[11] += bb[i] * e[i];
			a2[12] += bb[i] * f[i];
			a2[13] += bb[i] * g[i];

			a2[14] += cc[i] * a[i];
			a2[15] += cc[i] * b[i];
			a2[16] += cc[i] * c[i];
			a2[17] += cc[i];
			a2[18] += cc[i] * e[i];
			a2[19] += cc[i] * f[i];
			a2[20] += cc[i] * g[i];

			a2[21] += dd[i] * a[i];
			a2[22] += dd[i] * b[i];
			a2[23] += dd[i] * c[i];
			a2[24] += dd[i];
			a2[25] += dd[i] * e[i];
			a2[26] += dd[i] * f[i];
			a2[27] += dd[i] * g[i];

			a2[28] += ee[i] * a[i];
			a2[29] += ee[i] * b[i];
			a2[30] += ee[i] * c[i];
			a2[31] += ee[i];
			a2[32] += ee[i] * e[i];
			a2[33] += ee[i] * f[i];
			a2[34] += ee[i] * g[i];

			a2[35] += ff[i] * a[i];
			a2[36] += ff[i] * b[i];
			a2[37] += ff[i] * c[i];
			a2[38] += ff[i];
			a2[39] += ff[i] * e[i];
			a2[40] += ff[i] * f[i];
			a2[41] += ff[i] * g[i];

			a2[42] += gg[i] * a[i];
			a2[43] += gg[i] * b[i];
			a2[44] += gg[i] * c[i];
			a2[45] += gg[i];
			a2[46] += gg[i] * e[i];
			a2[47] += gg[i] * f[i];
			a2[48] += gg[i] * g[i];

			a4[0] += a[i] * a[i];
			a4[1] += a[i] * b[i];
			a4[2] += a[i] * c[i];
			a4[3] += a[i];
			a4[4] += a[i] * e[i];
			a4[5] += a[i] * f[i];
			a4[6] += a[i] * g[i];

			a4[7] += b[i] * a[i];
			a4[8] += b[i] * b[i];
			a4[9] += b[i] * c[i];
			a4[10] += b[i];
			a4[11] += b[i] * e[i];
			a4[12] += b[i] * f[i];
			a4[13] += b[i] * g[i];

			a4[14] += c[i] * a[i];
			a4[15] += c[i] * b[i];
			a4[16] += c[i] * c[i];
			a4[17] += c[i];
			a4[18] += c[i] * e[i];
			a4[19] += c[i] * f[i];
			a4[20] += c[i] * g[i];

			a4[21] += a[i];
			a4[22] += b[i];
			a4[23] += c[i];
			a4[24] += 1.0;
			a4[25] += e[i];
			a4[26] += f[i];
			a4[27] += g[i];

			a4[28] += e[i] * a[i];
			a4[29] += e[i] * b[i];
			a4[30] += e[i] * c[i];
			a4[31] += e[i];
			a4[32] += e[i] * e[i];
			a4[33] += e[i] * f[i];
			a4[34] += e[i] * g[i];

			a4[35] += f[i] * a[i];
			a4[36] += f[i] * b[i];
			a4[37] += f[i] * c[i];
			a4[38] += f[i];
			a4[39] += f[i] * e[i];
			a4[40] += f[i] * f[i];
			a4[41] += f[i] * g[i];

			a4[42] += g[i] * a[i];
			a4[43] += g[i] * b[i];
			a4[44] += g[i] * c[i];
			a4[45] += g[i];
			a4[46] += g[i] * e[i];
			a4[47] += g[i] * f[i];
			a4[48] += g[i] * g[i];
		}

		double eps = 1.0e-25;

		//DOP計算用　重みをすべて同じにする（HDOP,VDOP,PDOP,GDOP,TDOPの計算）
		for (i = 0; i <= 48; i++)
			xx2[i] = a4[i];
		//DOPのための逆行列計算//
		minv(xx2, eps, MinSvNum);

		//ここは地心座標系での値なので、ローカル座標系への変換　
		m2[0] = (-sin(lat)*cos(lon))*xx2[0] + (-sin(lat)*sin(lon))*xx2[7] + cos(lat)*xx2[14];
		m2[1] = (-sin(lat)*cos(lon))*xx2[1] + (-sin(lat)*sin(lon))*xx2[8] + cos(lat)*xx2[15];
		m2[2] = (-sin(lat)*cos(lon))*xx2[2] + (-sin(lat)*sin(lon))*xx2[9] + cos(lat)*xx2[16];
		m2[3] = (-sin(lon))*xx2[0] + cos(lon)*xx2[7];
		m2[4] = (-sin(lon))*xx2[1] + cos(lon)*xx2[8];
		m2[5] = (-sin(lon))*xx2[2] + cos(lon)*xx2[9];
		m2[6] = cos(lat)*cos(lon)*xx2[0] + cos(lat)*sin(lon)*xx2[7] + sin(lat)*xx2[14];
		m2[7] = cos(lat)*cos(lon)*xx2[1] + cos(lat)*sin(lon)*xx2[8] + sin(lat)*xx2[15];
		m2[8] = cos(lat)*cos(lon)*xx2[2] + cos(lat)*sin(lon)*xx2[9] + sin(lat)*xx2[16];

		m1[0] = m2[0] * (-sin(lat)*cos(lon)) + m2[1] * (-sin(lat)*sin(lon)) + m2[2] * cos(lat);
		m1[1] = m2[0] * (-sin(lon)) + m2[1] * cos(lon);
		m1[2] = m2[0] * cos(lat)*cos(lon) + m2[1] * cos(lat)*sin(lon) + m2[2] * sin(lat);
		m1[3] = m2[3] * (-sin(lat)*cos(lon)) + m2[4] * (-sin(lat)*sin(lon)) + m2[5] * cos(lat);
		m1[4] = m2[3] * (-sin(lon)) + m2[4] * cos(lon);
		m1[5] = m2[3] * cos(lat)*cos(lon) + m2[4] * cos(lat)*sin(lon) + m2[5] * sin(lat);
		m1[6] = m2[6] * (-sin(lat)*cos(lon)) + m2[7] * (-sin(lat)*sin(lon)) + m2[8] * cos(lat);
		m1[7] = m2[6] * (-sin(lon)) + m2[7] * cos(lon);
		m1[8] = m2[6] * cos(lat)*cos(lon) + m2[7] * cos(lat)*sin(lon) + m2[8] * sin(lat);

		//DOP値をアウトプット
		HDOP = sqrt(m1[0] + m1[4]);
		VDOP = sqrt(m1[8]);
		PDOP = sqrt(xx2[0] + xx2[8] + xx2[16]);
		TDOP = sqrt(xx2[24]);

		for (i = 0; i <= 48; i++)
			xx[i] = a2[i];
		//逆行列計算//
		minv(xx, eps, MinSvNum);

		delta[0] = 0.0;
		delta[1] = 0.0;
		delta[2] = 0.0;
		delta[3] = 0.0;
		delta[4] = 0.0;
		delta[5] = 0.0;
		delta[6] = 0.0;

		for (i = 0; i<SATn[rcvn]; i++) {
			a3[i] = xx[0] * aa[i] + xx[1] * bb[i] + xx[2] * cc[i] + xx[3] * dd[i] + xx[4] * ee[i] + xx[5] * ff[i] + xx[6] * gg[i];
			a3[i + SATn[rcvn]] = xx[7] * aa[i] + xx[8] * bb[i] + xx[9] * cc[i] + xx[10] * dd[i] + xx[11] * ee[i] + xx[12] * ff[i] + xx[13] * gg[i];
			a3[i + SATn[rcvn] * 2] = xx[14] * aa[i] + xx[15] * bb[i] + xx[16] * cc[i] + xx[17] * dd[i] + xx[18] * ee[i] + xx[19] * ff[i] + xx[20] * gg[i];
			a3[i + SATn[rcvn] * 3] = xx[21] * aa[i] + xx[22] * bb[i] + xx[23] * cc[i] + xx[24] * dd[i] + xx[25] * ee[i] + xx[26] * ff[i] + xx[27] * gg[i];
			a3[i + SATn[rcvn] * 4] = xx[28] * aa[i] + xx[29] * bb[i] + xx[30] * cc[i] + xx[31] * dd[i] + xx[32] * ee[i] + xx[33] * ff[i] + xx[34] * gg[i];
			a3[i + SATn[rcvn] * 5] = xx[35] * aa[i] + xx[36] * bb[i] + xx[37] * cc[i] + xx[38] * dd[i] + xx[39] * ee[i] + xx[40] * ff[i] + xx[41] * gg[i];
			a3[i + SATn[rcvn] * 6] = xx[42] * aa[i] + xx[43] * bb[i] + xx[44] * cc[i] + xx[45] * dd[i] + xx[46] * ee[i] + xx[47] * ff[i] + xx[48] * gg[i];

			delta[0] += a3[i] * r3[i];
			delta[1] += a3[i + SATn[rcvn]] * r3[i];
			delta[2] += a3[i + SATn[rcvn] * 2] * r3[i];
			delta[3] += a3[i + SATn[rcvn] * 3] * r3[i];
			delta[4] += a3[i + SATn[rcvn] * 4] * r3[i];
			delta[5] += a3[i + SATn[rcvn] * 5] * r3[i];
			delta[6] += a3[i + SATn[rcvn] * 6] * r3[i];
		}

		for (i = 0; i<3; i++) {
			init[i] += delta[i];
		}
		init[3] = delta[3];
		init[4] = delta[4];
		init[5] = delta[5];
		init[6] = delta[6];
		LS_update[0]=delta[0];LS_update[1]=delta[1];LS_update[2]=delta[2];//最小二乗法の修正量が十分ちいさいか

	}//GPS/QZSS+3つの他の測位衛星


	double xxx[3]={0},dis;
	trans_xyz_llh(init[0],init[1],init[2],xxx);
	dis = sqrt((init[0]-Ref_pos[0][1])*(init[0]-Ref_pos[0][1])+(init[1]-Ref_pos[1][1])*(init[1]-Ref_pos[1][1])+(init[2]-Ref_pos[2][1])*(init[2]-Ref_pos[2][1]));
}
