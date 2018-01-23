//////////////////////////////////////////////////////////////////////
//
// GPS/QZS/GALILEOでのRTK
//
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "global_extern.h"

using namespace std;

const double cs     = 299792458.0;//光速
const double F1     = 1575420000.0;
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//円周率
const double f      = -4.442807633e-10;//

int minv(double[],double,int);
void select_sat(int[],int);
void w_inv(double[],int,int,int,int);
void trans_xyz_llh(double,double,double, double []);
extern "C" {int lambda(int, int, const double *, const double *, double *, double *); }

//int lambda(int, int, const double[], const double[], double[], double[]);

void calc_rtk_GQE(int rcvn)
{
	int i,j,k,ndim;
	int prn;
	int max_prn_gps,max_prn=1;
	int kaisu=0;
	int novalid=0,miss=0;
	int l=0;
	int flag=1,gpsnum=0,sv[PRN]={0};
	
	double lrl[PRN],lrk[PRN][PRN],lrln[PRN],lrkn[PRN][PRN];
	double lrlnb[PRN];
	double dx,dy,dz;
	double y[PRN]={0};
	double eps;
	static double g[PRN*PRN]={0},sg[PRN*PRN]={0},gtw[PRN*PRN]={0},g2[PRN*PRN]={0},g3[PRN*PRN]={0};
	static double gt[PRN*PRN]={0},w[PRN*PRN]={0},gtwg[PRN*PRN]={0},inv_gtwg[PRN*PRN]={0},gtwg_gtw[PRN*PRN]={0};
	double q[PRN*PRN]={0};
	double ref_diff[PRN]={0},rov_diff[PRN]={0};
	double max_ele=0;
	double dd_code[PRN],dd_carrier[PRN];
	double n1[PRN]={0};
	double n[PRN]={0};
	double dx_float,dy_float,dz_float;
	double omomi[PRN]={0},distance;
	double rov_diff_l1[PRN]={0},ref_diff_l1[PRN]={0},n_init[PRN]={0};

	static double s_factor=0.5,pre_ratio=0,Pre_F[256]={0};
	double threshold_ratio = Ratio_limit;
	int No_dp_flag=1;

	static int Pre_SATn,Pre_SVn[PRN] = {0},Pre_nsvg,Pre_nsvb,SKIP_flag=0,Pre_max_gps=0,Pre_max_bei=0;
	int stock_sat[PRN]={0},gps1=0,bei1=0;
	
	if(No_rtk>=1)
		s_factor=0.5;

	int nsv = SATn[rcvn]-1;//TOTAL
	int nsvg;//GPS
	int nsvb = 0;

	dx = 0;
	dy = 0;
	dz = 0;

	for(i=0;i<PRN;i++){
		lrl[i]=0,lrln[i]=0,lrlnb[i]=0;
		for(j=0;j<PRN;j++){
			lrk[i][j]=0,lrkn[i][j]=0;
		}
	}

	//reference satellite setting
	max_ele=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		
		if(prn<=70)
			gpsnum++;

		if(max_ele < Elevation[rcvn][prn] && prn<=70){
			max_ele = Elevation[rcvn][prn];
			max_prn_gps = prn;
		}
	}
	max_prn = max_prn_gps;

	nsvg = gpsnum-1;

	j=0;

	for(i=0;i<PRN*PRN;i++)
		w[i]=0;

	j=0;
	//可視衛星の2重位相差の算出
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=70){
			//搬送波の符号に注意
			dd_carrier[prn] = cs/F1*((Cp1[0][prn] - Cp1[1][prn])
				-(Cp1[0][max_prn_gps] - Cp1[1][max_prn_gps]));

			dd_code[prn] = (Pr1[0][prn] - Pr1[1][prn])
				-(Pr1[0][max_prn_gps] - Pr1[1][max_prn_gps]);

			n1[prn] = dd_carrier[prn]*F1/cs - dd_code[prn]*F1/cs;

			OmomiP[j] = Code_noise+Code_noise/sin(Elevation[rcvn][prn]*pi/180.0);
			OmomiC[j] = Carrier_noise+Carrier_noise/sin(Elevation[rcvn][prn]*pi/180.0);

//			OmomiC[j] = sqrt(90.0/Elevation[rcvn][prn]);
//			OmomiP[j] = sqrt(1.0/);  
			OmomiP[j] = OmomiP[j]*OmomiP[j];
			OmomiC[j] = OmomiC[j]*OmomiC[j];
			j++;

		}//最大仰角衛星をはずす
	}//衛星数

	w_inv(w,SATn[rcvn]-1,flag,nsvg,nsvb);

	j=0;//gps l1 bei b1 + gps l1 bei b1
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=70)
			y[j++] = dd_carrier[prn];
	}
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=70)
			y[j++] = dd_code[prn];
	}

	distance=0;

////////////////////////////////////////////////////////////////////////////////////////////
	lrl[0] = SVx[0][max_prn_gps]-Ref_pos[0][1];
	lrl[1] = SVy[0][max_prn_gps]-Ref_pos[1][1];
	lrl[2] = SVz[0][max_prn_gps]-Ref_pos[2][1];
	lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_gps && prn<=70){
			lrk[0][j] = SVx[0][prn]-Ref_pos[0][1];
			lrk[1][j] = SVy[0][prn]-Ref_pos[1][1];
			lrk[2][j] = SVz[0][prn]-Ref_pos[2][1];
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////

///*

	int iter;

	for(iter=1;iter<=4;iter++){

		//G行列の生成　上からL1搬送波, L1コード
		for(i=0;i<nsv+3;i++){
			for(j=0;j<nsv*2;j++){
				g[i+j*(nsv+3)]=-1.0*(lrkn[i][j%nsv]-lrln[i]);
			}
		}
		for(i=0;i<(nsv+3)*nsv;i++){
			if((i-3)%(nsv+3+1)==0)
				g[i] = cs/F1;
		}

		//Gの転置行列の生成
		for(i=0;i<nsv*2;i++){
			for(j=0;j<nsv+3;j++){
				gt[i+j*nsv*2]=g[i*(nsv+3)+j];
			}
		}

		for(i=0;i<PRN*PRN;i++){
			gtw[i]=0;
			gtwg[i]=0;
			gtwg_gtw[i]=0;
		}

		//G(t)*W
		for(i=0;i<nsv+3;i++){
			for(j=0;j<nsv*2;j++){
				for(k=0;k<nsv*2;k++){
					gtw[i*nsv*2+j]=gtw[i*nsv*2+j]+gt[k+i*nsv*2]*w[k*nsv*2+j];
				}
			}
		}

		//G(t)*W*G
		for(i=0;i<nsv+3;i++){
			for(j=0;j<nsv+3;j++){
				for(k=0;k<nsv*2;k++){
					gtwg[i*(nsv+3)+j]=gtwg[i*(nsv+3)+j]+gtw[k+i*nsv*2]*g[k*(nsv+3)+j];
				}
			}
		}

		//G(t)*W*Gの逆行列を求める
		eps=1.0e-25;
		minv(gtwg,eps,nsv+3);

		//inv((G(t)*W*G))*G(t)W
		for(i=0;i<nsv+3;i++){
			for(j=0;j<nsv*2;j++){
				for(k=0;k<=nsv+3;k++){
					gtwg_gtw[i*nsv*2+j]=gtwg_gtw[i*nsv*2+j]+gtwg[k+i*(nsv+3)]*gtw[k*nsv*2+j];
				}
			}
		}

		if(iter==1){
			j=0;
			for(i=0;i<SATn[rcvn];i++){
				prn = SVn[rcvn][i];
				if(prn!=max_prn){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))
						-sqrt((SVx[1][max_prn]-Ref_pos[0][1])*(SVx[1][max_prn]-Ref_pos[0][1])
						+(SVy[1][max_prn]-Ref_pos[1][1])*(SVy[1][max_prn]-Ref_pos[1][1])
						+(SVz[1][max_prn]-Ref_pos[2][1])*(SVz[1][max_prn]-Ref_pos[2][1]));
					j++;
				}
			}
			for(i=0;i<j*2;i++){
				ref_diff[i] = ref_diff[i%j];
			}
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
				j++;
			}
		}
		for(i=0;i<j*2;i++){
			rov_diff[i] = rov_diff[i%j];
		}

		//dx,dy,dz,float_integerの計算
		for(i=0;i<nsv*2;i++){
			dx = dx+gtwg_gtw[i]*(y[i]+ref_diff[i]-rov_diff[i]);
			dy = dy+gtwg_gtw[i+nsv*2]*(y[i]+ref_diff[i]-rov_diff[i]);
			dz = dz+gtwg_gtw[i+nsv*2*2]*(y[i]+ref_diff[i]-rov_diff[i]);

			if(iter==4){
				for(j=3;j<nsv+3;j++)
					n[j-3] = n[j-3]+gtwg_gtw[i+nsv*2*j]*(y[i]+ref_diff[i]-rov_diff[i]);
			}
		}

		//衛星の方向の単位ベクトルを算出
		//基準局と基準衛星の単位ベクトル
		lrl[0] = SVx[0][max_prn]-(Ref_pos[0][1]+dx);
		lrl[1] = SVy[0][max_prn]-(Ref_pos[1][1]+dy);
		lrl[2] = SVz[0][max_prn]-(Ref_pos[2][1]+dz);
		lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);

		//移動局と選択衛星の単位ベクトル
		//衛星数分まわす
		//最大仰角の衛星はすでに選択されている

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}
	}//iter

	//float解の抽出
	dx_float = dx;
	dy_float = dy;
	dz_float = dz;

	double llh[3]={0};
	trans_xyz_llh(dx+Ref_pos[0][1],dy+Ref_pos[1][1],dz+Ref_pos[2][1],llh);

	if(llh[0]>-90 && llh[0]<90){//緯度だけチェック
		User_lat = llh[0];//電離層等の計算用
		User_lon = llh[1];//電離層等の計算用
		User_hgt = llh[2];
		User_pos[0] = dx+Ref_pos[0][1];
		User_pos[1] = dx+Ref_pos[1][1];
		User_pos[2] = dx+Ref_pos[2][1];
	}


	//LAMBDA法の適用
	double F[PRN]={0},s[PRN]={0};

	j=0;
	for(i=0;i<=(nsv+3)*(nsv+3);i++){
		if(i>=(nsv+3)*3){
			q[j] = gtwg[i];
			j++;
			if(i%(nsv+3)==0)
				j--;
			if((i-1)%(nsv+3)==0)
				j--;
			if((i-2)%(nsv+3)==0)
				j--;
		}
	}

	k=0;
	ndim=nsv;
	double det=1.0,buf;
	double aa[PRN][PRN]={0};
	for(i=0;i<ndim;i++){
		for(j=0;j<nsv;j++){
			aa[i][j] = q[k++];
		}
	}
	//三角行列を作成
	for(i=0;i<ndim;i++){
		for(j=0;j<ndim;j++){
			if(i<j){
				buf=aa[j][i]/aa[i][i];
				for(k=0;k<ndim;k++){
					aa[j][k]-=aa[i][k]*buf;
				}
			}
		}
	}
	//対角部分の積
	for(i=0;i<ndim;i++){
		det*=aa[i][i];
	}
	det = pow(det,1.0/(nsv*2));//ADOP

	//LAMBDA法の呼び出し
	lambda(nsv,2,n,q,F,s);

///////////////////////////////////////////////////////////////////////////////////////////////////////
	//LAMBDA法で求めたアンビギュイテイを利用して測位計算を行う
///////////////////////////////////////////////////////////////////////////////////////////////////////

	//各衛星の仰角による重みの生成
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn != max_prn_gps && prn<=70)
			omomi[j++] = 1.0*1.0*sqrt(90.0/Elevation[rcvn][prn]);
	}

	dx=0;dy=0;dz=0;

	for(i=nsv;i<nsv*2;i++)//残差最小の候補（アンビギュイティそのもの）のみ残す
		F[i]=0;

/////////////////////////////////////////////////////////////////////////////////////////////
	lrl[0] = SVx[0][max_prn_gps]-(Ref_pos[0][1]+dx);
	lrl[1] = SVy[0][max_prn_gps]-(Ref_pos[1][1]+dy);
	lrl[2] = SVz[0][max_prn_gps]-(Ref_pos[2][1]+dz);
	lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
	j=0;
	for(i=0;i<SATn[rcvn];i++){
		prn = SVn[rcvn][i];
		if(prn!=max_prn_gps && prn<=70){
			lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
			lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
			lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
			lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
			j++;
		}
	}
/////////////////////////////////////////////////////////////////////////////////////////////

	for(iter=1;iter<=3;iter++){
		//dd_carrierをｍから波数にすべき？
		if(iter==1){
			for(i=0;i<nsv;i++){
				y[i] = y[i]*F1/cs;
			}
		}

		//G行列の生成
		for(i=0;i<=2;i++){
			for(j=0;j<nsv;j++){
				g[i+j*3]=-1.0*(lrkn[i][j]-lrln[i])*F1/cs;
			}
		}

		//Gの転置行列の生成
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=g[i*3+j];
			}
		}

		//Gの転置行列とWの掛け算　これ以降はgtはgtwとなっている
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				gt[i+j*nsv]=gt[i+j*nsv]*omomi[i];
			}
		}

		for(i=0;i<=8;i++)
			g3[i]=0;

		//Gの転置行列とG行列を掛ける
		for(i=0;i<nsv;i++){
			g3[0]=g3[0]+gt[i]*g[i*3];
			g3[1]=g3[1]+gt[i]*g[i*3+1];
			g3[2]=g3[2]+gt[i]*g[i*3+2];
			g3[3]=g3[3]+gt[i+nsv]*g[i*3];
			g3[4]=g3[4]+gt[i+nsv]*g[i*3+1];
			g3[5]=g3[5]+gt[i+nsv]*g[i*3+2];
			g3[6]=g3[6]+gt[i+nsv*2]*g[i*3];
			g3[7]=g3[7]+gt[i+nsv*2]*g[i*3+1];
			g3[8]=g3[8]+gt[i+nsv*2]*g[i*3+2];
		}

		//G3行列の逆行列を求める
		//逆行列計算//
		eps = 1.0e-25;
		minv(g3,eps,3);

		for(i=0;i<nsv*3;i++)
			g2[i]=0;
		//逆行列にGの転置行列を掛ける
		for(i=0;i<nsv;i++){
			for(j=0;j<3;j++){
				g2[i]      =g2[i]+g3[j]*gt[i+nsv*j];
				g2[nsv+i]  =g2[nsv+i]+g3[j+3]*gt[i+nsv*j];
				g2[nsv*2+i]=g2[nsv*2+i]+g3[j+6]*gt[i+nsv*j];
			}
		}

		if(iter==1){
			j=0;
			for(i=0;i<SATn[rcvn];i++){
				prn = SVn[rcvn][i];
				if(prn!=max_prn){
					ref_diff[j] = sqrt((SVx[1][prn]-Ref_pos[0][1])*(SVx[1][prn]-Ref_pos[0][1])
						+(SVy[1][prn]-Ref_pos[1][1])*(SVy[1][prn]-Ref_pos[1][1])
						+(SVz[1][prn]-Ref_pos[2][1])*(SVz[1][prn]-Ref_pos[2][1]))
						-sqrt((SVx[1][max_prn]-Ref_pos[0][1])*(SVx[1][max_prn]-Ref_pos[0][1])
						+(SVy[1][max_prn]-Ref_pos[1][1])*(SVy[1][max_prn]-Ref_pos[1][1])
						+(SVz[1][max_prn]-Ref_pos[2][1])*(SVz[1][max_prn]-Ref_pos[2][1]));
					ref_diff[j] = ref_diff[j]*F1/cs;
					j++;
				}
			}
		}

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn){
				rov_diff[j] = sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j])
					- sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
				rov_diff[j] = rov_diff[j]*F1/cs;
				j++;
			}
		}

		//dx,dy,dzの計算
		for(i=0;i<nsv;i++){
			dx = dx+g2[i]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dy = dy+g2[i+nsv]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
			dz = dz+g2[i+nsv*2]*(y[i]+ref_diff[i]-rov_diff[i]-F[i]);
		}

		//衛星の方向の単位ベクトルを算出
		//基準局と基準衛星の単位ベクトル
		lrl[0] = SVx[0][max_prn]-(Ref_pos[0][1]+dx);
		lrl[1] = SVy[0][max_prn]-(Ref_pos[1][1]+dy);
		lrl[2] = SVz[0][max_prn]-(Ref_pos[2][1]+dz);
		lrln[0] = lrl[0]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[1] = lrl[1]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);
		lrln[2] = lrl[2]/sqrt(lrl[0]*lrl[0]+lrl[1]*lrl[1]+lrl[2]*lrl[2]);

		//移動局と選択衛星の単位ベクトル
		//衛星数分まわす
		//最大仰角の衛星はすでに選択されている

		j=0;
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(prn!=max_prn){
				lrk[0][j] = SVx[0][prn]-(Ref_pos[0][1]+dx);
				lrk[1][j] = SVy[0][prn]-(Ref_pos[1][1]+dy);
				lrk[2][j] = SVz[0][prn]-(Ref_pos[2][1]+dz);
				lrkn[0][j] = lrk[0][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[1][j] = lrk[1][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				lrkn[2][j] = lrk[2][j]/sqrt(lrk[0][j]*lrk[0][j]+lrk[1][j]*lrk[1][j]+lrk[2][j]*lrk[2][j]);
				j++;
			}
		}
	}//iter

	//DGNSSとRTKの結果を出力

	//x,y,z -> lat,lon,hgt
	double ecef[3]={0},t[3]={0};
	double a,f,b,p,e,oldt,nn;
	double lat_diff,lon_diff,hgt_diff;

/////////////////////////////////////////////////////////////////////////////
	//DGNSS
	ecef[0]=dx_float+Ref_pos[0][1];
	ecef[1]=dy_float+Ref_pos[1][1];
	ecef[2]=dz_float+Ref_pos[2][1];

	a = 6378137;
	f = 1/298.257223563;
	b = a*(1-f);
	p = sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1]);
	e = f*(2-f);
	t[0] = atan(ecef[2]/((1-e)*p));
	oldt = 1;
	while(fabs(t[2]-oldt) > 0.000001)
	{
		oldt = t[2];
		nn = a/sqrt(1-(e*sin(t[0])*sin(t[0])));
		t[2] = p/cos(t[0]) - nn;
		t[0] = atan((ecef[2]/p)/((1-e*nn/(nn+t[2]))));
	}
	t[1] = 2*atan(ecef[1]/(ecef[0]+sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1])));
	t[0] = t[0]*180/pi;
	t[1] = t[1]*180/pi;

	lat_diff = (t[0]-POSrcvlat[1])*110947.0;
	lon_diff = (t[1]-POSrcvlon[1])*cos(t[0]*pi/180.0)*111319.0;
	hgt_diff = t[2] - POSrcvhgt[1];

	fprintf(fp[4],"%f,%f,%f,%f,%15.10f,%15.10f,%f,%d\n",
		DGPSTIME,lon_diff,lat_diff,hgt_diff,t[0],t[1],t[2],SATn[rcvn]);
	
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
	//RTK
	ecef[0]=dx+Ref_pos[0][1];
	ecef[1]=dy+Ref_pos[1][1];
	ecef[2]=dz+Ref_pos[2][1];

	a = 6378137;
	f = 1/298.257223563;
	b = a*(1-f);
	p = sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1]);
	e = f*(2-f);
	t[0] = atan(ecef[2]/((1-e)*p));
	oldt = 1;
	while(fabs(t[2]-oldt) > 0.000001)
	{
		oldt = t[2];
		nn = a/sqrt(1-(e*sin(t[0])*sin(t[0])));
		t[2] = p/cos(t[0]) - nn;
		t[0] = atan((ecef[2]/p)/((1-e*nn/(nn+t[2]))));
	}
	t[1] = 2*atan(ecef[1]/(ecef[0]+sqrt(ecef[0]*ecef[0]+ecef[1]*ecef[1])));
	t[0] = t[0]*180/pi;
	t[1] = t[1]*180/pi;

	lat_diff = (t[0]-POSrcvlat[1])*110947.0;
	lon_diff = (t[1]-POSrcvlon[1])*cos(t[0]*pi/180.0)*111319.0;
	hgt_diff = t[2] - POSrcvhgt[1];

	s_factor = s[1]/s[0];
	//ADOPでのチェック
	if(det>=1.00)
		s_factor = 1.0;
	Ratio = s_factor;
	distance = sqrt((lon_diff)*(lon_diff)+(lat_diff)*(lat_diff));


	if(s_factor>=threshold_ratio){
		Sol_flag[2]++;//FIX回数カウント
		fprintf(fp[5],"%f,%f,%f,%f,%15.10f,%15.10f,%f,%f,%d\n",
			DGPSTIME,lon_diff,lat_diff,hgt_diff,t[0],t[1],t[2],Ratio,SATn[rcvn]);
	
//		for(i=0;i<SATn[rcvn];i++){
//			fprintf(fp[5],",%d",SVn[rcvn][i]);
//		}
//		fprintf(fp[5],"\n");
	}
/////////////////////////////////////////////////////////////////////////////

}