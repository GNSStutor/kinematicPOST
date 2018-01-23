//////////////////////////////////////////////////////////////
//
// 全ての衛星の位置を計算する（衛星時計の補正も実施）
// GPS/QZS/GALILEO/BeiDouはほぼ同じ　GLONASSは別
//
//////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "global_extern.h"

using namespace std;

double newton(double,double,int);

const double cs     = 299792458.0;//
const double myu    = 3.986004418e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//
const double f      = -4.442807633e-10;//

#define J2_GLO   1.0826257E-3
#define OMGE_GLO 7.292115E-5
#define omg2     5.317494117322500E-009
#define MU_GLO   3.9860044E+14
#define RE_GLO   6378136.0

static double calc_glosvpos( double y0[], double t1, double *yout);
static double RK4( double y0[], double t1, double *yout);
static double satfunc( double y0[], double t1, double k[],double kout[]);

void calc_satpos(int rcvn)
{
	int i,j;
	int prn;
	double range[PRN];
	double tt[PRN],tk[PRN],t_toc,a[PRN],n[PRN];
	double n0[PRN],mk[PRN],e[PRN],ek[PRN],sinvk[PRN],cosvk[PRN];
	double vk[PRN],dtr[PRN]={0},dk[PRN],duk[PRN],uk[PRN],drk[PRN],rk[PRN];
	double dik[PRN],ik[PRN],omegak[PRN],ypos2[PRN],xpos2[PRN],zpos2[PRN];
	double xpos[PRN],ypos[PRN],zpos[PRN],xposold[PRN]={0},yposold[PRN]={0};
	double xposold1[PRN],yposold1[PRN];
	double system_offset=0,xxx[9]={0};

	j = 0;
	for(i=1;i<=PRN-1;i++)
		SVpos_flag[rcvn][i]=0;

	for(i=0;i<SATn[rcvn];i++){ 
		prn = SVn[rcvn][i];//

		if(prn==0)
			continue;
		if(prn<=100){
			if(prn>=71)
				system_offset = -14.0;
			else if(prn>=41 && prn<=70)
				system_offset = 0.0000000;
			else
				system_offset = 0;

			if(rcvn==1)
				t_toc=GPSTIME + system_offset - Ephe.toe[prn];
			else if(rcvn==0)
				t_toc=DGPSTIME + system_offset - Ephe.toe[prn];

			if(t_toc>302400)
				t_toc=t_toc-604800;
			if(t_toc<-302400)
				t_toc=t_toc+604800;

			//satellite time correction
			range[prn] = Pr1[rcvn][prn] + cs*(Ephe.af0[prn]+Ephe.af1[prn]*t_toc
		 			+Ephe.af2[prn]*t_toc*t_toc);

			//satellite time correction
			range[prn] = range[prn] - cs*Ephe.tgd[prn];

			//transmit time estimation
			if(rcvn==1)
				tt[prn] = GPSTIME + system_offset - range[prn]/cs;
			else if(rcvn==0)
				tt[prn] = DGPSTIME + system_offset - range[prn]/cs;

			tk[prn] = tt[prn] - Ephe.toe[prn];

			if(tk[prn]>302400)
				tk[prn]=tk[prn]-604800;
			if(tk[prn]<-302400)
				tk[prn]=tk[prn]+604800;

			a[prn]     = Ephe.roota[prn]*Ephe.roota[prn];
			n0[prn]    = sqrt(myu/(a[prn]*a[prn]*a[prn]));
			n[prn]     = n0[prn] + Ephe.dn[prn];
			mk[prn]    = Ephe.m0[prn] + n[prn]*tk[prn];
			e[prn]     = Ephe.e[prn];

			//エフェメリスのヘルス情報より利用可否判断
			//ついでにエフェメリスが存在するかどうかを軌道半径でチェック
			if(Ephe.roota[prn]<=0.1 && Ephe.roota[prn]>=-0.1 || Ephe.health[prn]>0.1){
//				fprintf(fp[8],"No Ephemeris   PRN=%d   time=%7.1f\n",prn,GPSTIME);
				j = 0;
			}
			else
				SVpos_flag[rcvn][prn]=1;

			ek[prn]    = newton(mk[prn],e[prn],prn);

			sinvk[prn] = sqrt(1.0-e[prn]*e[prn])*sin(ek[prn])/(1.0-e[prn]*cos(ek[prn]));
			cosvk[prn] = (cos(ek[prn])-e[prn])/(1.0-e[prn]*cos(ek[prn]));		
			if(cosvk[prn] == 0 || e[prn]*cos(vk[prn]) == -1.0){
				printf("satpos error\n");
				exit(1);
			}		
			vk[prn]     = atan2(sinvk[prn],cosvk[prn]);
		
			if(ek[prn]>=0)
				ek[prn] = acos((e[prn]+cos(vk[prn]))/(1.0+e[prn]*cos(vk[prn])));
			else
				ek[prn] = -acos((e[prn]+cos(vk[prn]))/(1.0+e[prn]*cos(vk[prn])));
		
			dtr[prn]  = f*e[prn]*Ephe.roota[prn]*sin(ek[prn]);

			range[prn]  = range[prn] + cs*dtr[prn];
	       			
			dk[prn]     = vk[prn] + Ephe.omega[prn];
			duk[prn]    = Ephe.cus[prn]*sin(2*dk[prn])+Ephe.cuc[prn]*cos(2*dk[prn]);
			uk[prn]     = dk[prn] + duk[prn];
			drk[prn]    = Ephe.crc[prn]*cos(2*dk[prn])+Ephe.crs[prn]*sin(2*dk[prn]);
			rk[prn]     = a[prn]*(1.0-e[prn]*cos(ek[prn]))+drk[prn];
			dik[prn]    = Ephe.cic[prn]*cos(2*dk[prn])+Ephe.cis[prn]*sin(2*dk[prn]);
			ik[prn]     = Ephe.i0[prn] + dik[prn] + Ephe.di0[prn]*tk[prn];
			omegak[prn] = Ephe.omega0[prn] + (Ephe.domega0[prn] - omegae)*tk[prn]-omegae*Ephe.toe[prn];

	//		omegak[prn] = Ephe.omega0[prn] + (Ephe.domega0[prn] - omegae)*(tk[prn]+range[prn]/cs)-omegae*Ephe.toe[prn];

			xpos2[prn]  = rk[prn]*cos(uk[prn]);
			ypos2[prn]  = rk[prn]*sin(uk[prn]);
			xpos[prn]   = xpos2[prn]*cos(omegak[prn]) - ypos2[prn]*cos(ik[prn])*sin(omegak[prn]);
			ypos[prn]   = xpos2[prn]*sin(omegak[prn]) + ypos2[prn]*cos(ik[prn])*cos(omegak[prn]);
			zpos[prn]   = ypos2[prn]*sin(ik[prn]);		
			xposold1[prn] = xpos[prn];
			yposold1[prn] = ypos[prn];
	      
			for(j=1;j<=2;j++){
				xposold[prn] = +omegae*ypos[prn]*range[prn]/cs;
				yposold[prn] = -omegae*xpos[prn]*range[prn]/cs;
				xpos[prn] = xposold1[prn] + xposold[prn];
				ypos[prn] = yposold1[prn] + yposold[prn];
			}

			//BEIDOU GEO
			if(prn>=71 && prn<=75){
				xposold[prn]=0;
				yposold[prn]=0;

				omegak[prn] = Ephe.omega0[prn] + Ephe.domega0[prn]*(tk[prn])-omegae*Ephe.toe[prn];

				xpos[prn]   = xpos2[prn]*cos(omegak[prn]) - ypos2[prn]*cos(ik[prn])*sin(omegak[prn]);
				ypos[prn]   = xpos2[prn]*sin(omegak[prn]) + ypos2[prn]*cos(ik[prn])*cos(omegak[prn]);
				zpos[prn]   = ypos2[prn]*sin(ik[prn]);

				xpos2[prn] = xpos[prn];
				ypos2[prn] = ypos[prn]*cos(-5.0*pi/180.0) + zpos[prn]*sin(-5.0*pi/180.0);
				zpos2[prn] = -1.0*ypos[prn]*sin(-5.0*pi/180.0) + zpos[prn]*cos(-5.0*pi/180.0);

				xpos[prn] = xpos2[prn]*cos(omegae*tk[prn]) + ypos2[prn]*sin(omegae*tk[prn]);
				ypos[prn] = -1.0*xpos2[prn]*sin(omegae*tk[prn]) + ypos2[prn]*cos(omegae*tk[prn]);
				zpos[prn] = zpos2[prn];

				xposold1[prn] = xpos[prn];
				yposold1[prn] = ypos[prn];

				for(j=1;j<=2;j++){
					xposold[prn] = +omegae*ypos[prn]*range[prn]/cs;
					yposold[prn] = -omegae*xpos[prn]*range[prn]/cs;
					xpos[prn] = xposold1[prn] + xposold[prn];
					ypos[prn] = yposold1[prn] + yposold[prn];
				}
			}

			SVx[rcvn][prn] = xpos[prn];
			SVy[rcvn][prn] = ypos[prn];
			SVz[rcvn][prn] = zpos[prn];

//			if(Rinex==1){
				SV_corrtime[rcvn][prn]=range[prn] - Pr1[rcvn][prn];
//			}

		}//prn<=100
		else {
			int TSTEP = 180; //step幅
			int i;
			int count = 0;
			double stock_toe[PRN]={0};
			double G_dt[70]={0.0};
			double G_tkmin[70] = {0.0};
			double y0[9] = {0.0};
			double tkmin,*yout=0, yy0[9] = {0.0};
			yout = yy0;
			double xutc,temptime,t1,leapsecond=LeapSecond;
	 
		   if(rcvn==1)      temptime = GPSTIME;
		   else if(rcvn==0) temptime = DGPSTIME;

			xutc = temptime - leapsecond;
		    G_tkmin[prn] = xutc- Ephe.toc[prn];
						
			//衛星が電波を発射したGPSシステム時刻を求める
			//衛星時刻補正を行う
			G_dt[prn] = Ephe.TauN[prn] + Ephe.GammaN[prn] * G_tkmin[prn];// - TauC;
		//	G_dt[prn] = G_TauN[prn] + G_GammaN[prn] * G_tkmin[prn];
			
			range[prn] = Pr1[rcvn][prn] + cs*G_dt[prn];
      
			//衛星電波発射時刻の予測
			G_tkmin[prn] = xutc - range[prn]/cs;
			
			//衛星が軌道の元期からどれだけ進んでいる？
			 G_tkmin[prn] =  G_tkmin[prn] -  Ephe.toc[prn];
			
			tkmin = G_tkmin[prn] ; 

			y0[0] = Ephe.Xp[prn];
			y0[1] = Ephe.Yp[prn];
			y0[2] = Ephe.Zp[prn];
			
			y0[3] = Ephe.Xv[prn];
			y0[4] = Ephe.Yv[prn];
			y0[5] = Ephe.Zv[prn];
			
			y0[6] = Ephe.Xa[prn];
			y0[7] = Ephe.Ya[prn];
			y0[8] = Ephe.Za[prn];
			
			if(Ephe.roota[prn]<0.1)//エフェメリスのない衛星は計算しない
				tkmin=0;

			while((fabs(tkmin))>=TSTEP){
				if(tkmin>0){
				 	t1 = TSTEP;
				}else{
					t1 = -TSTEP;			
				}
				calc_glosvpos(y0, t1 ,yout);

				for(i=0;i<9;i++)y0[i] = yout[i];
				tkmin = tkmin - t1;  
			}
			t1=tkmin;
			calc_glosvpos( y0, t1 ,yout);

		for(i=0;i<9;i++)y0[i] = yout[i];

			SVx[rcvn][prn] = y0[0];
			SVy[rcvn][prn] = y0[1];
			SVz[rcvn][prn] = y0[2];
						
			omegak[prn] = omegae*range[prn]/cs;

			SVx[rcvn][prn]   = SVx[rcvn][prn]*cos(omegak[prn]) + SVy[rcvn][prn]*sin(omegak[prn]);
			SVy[rcvn][prn]   = -SVx[rcvn][prn]*sin(omegak[prn]) + SVy[rcvn][prn]*cos(omegak[prn]);

			SV_corrtime[rcvn][prn] = range[prn] - Pr1[rcvn][prn];

				
			//衛星位置計算ができたかどうか
			if( G_tkmin[prn]>7200 || SVx[rcvn][prn]>10000000000000000 || Ephe.health[prn]>0.1 || Ephe.roota[prn]<=0.1)
				SVpos_flag[rcvn][prn]=0;
			else
				SVpos_flag[rcvn][prn]=1;
		}
		//GLONASS
///////////////////////////////////////////////////////////
		//
	}//for
}

double calc_glosvpos( double y0[], double t1, double *yout){
    //アドレス渡しだと,同じアドレスのものは,関数外でも保持及び連動しているyoutとy0
       int i = 0;
       double b[9] = {0.0};
       for ( i = 0;i < 9; i++) b[i] = y0[i];//y0ではなく,bを使う理由は,bを変更してもy0の値には影響しないため
       RK4( b, t1, yout);
       return 0;
}
double RK4( double b[], double t1, double *yout){


       int i = 0;
       double dummy[9] = {0.0}, k1[9] = {0.0}, k2[9] = {0.0}, k3[9] = {0.0}, k4[9] = {0.0}, kout[9] = {0.0};
       double y[9] = {0.0};
		satfunc(b,    1, dummy, kout);for ( i = 0; i < 9; i++) k1[i] = kout[i];
		satfunc(b, t1/2,    k1, kout);for ( i = 0; i < 9; i++) k2[i] = kout[i];
		satfunc(b, t1/2,    k2, kout);for ( i = 0; i < 9; i++) k3[i] = kout[i];
		satfunc(b,   t1,    k3, kout);for ( i = 0; i < 9; i++) k4[i] = kout[i];
       
       for ( i = 0; i < 9; i++) y[i] = b[i] + t1/6*(k1[i] + 2*k2[i]+2*k3[i] + k4[i]);
       for ( i = 0; i < 9; i++) yout[i] = y[i];
       return 0;
}
double satfunc( double b[], double t1, double k[],double *kout ){
       //配列は戻り値にできない(returnで返せない) ので,引数で返す
       int i = 0;
       double r2 = 0.0, r3 = 0.0, a = 0.0, bb = 0.0, acc[3] = {0.0}, ydot[9] = {0.0};
       double yy0[9] = {0.0};
       for ( i = 0; i < 9; i++) yy0[i] = b[i];
       

       for ( i = 0; i < 9; i++) {
              yy0[i] = yy0[i] +  t1*k[i];
       }
       for ( i = 0; i < 3; i++)
       r2 = r2 + yy0[i]*yy0[i];
       r3 = r2*sqrt(r2);
       a = 1.5*J2_GLO*MU_GLO*RE_GLO*RE_GLO/r2/r3;
	   bb = 5.0*yy0[2]*yy0[2]/r2;                 

	   for (i = 0;i < 3;i++) acc[i] = yy0[i+6];
       for (i = 0;i < 3;i++) ydot[i] = yy0[i+3];

       ydot[3] = -MU_GLO*yy0[0]/r3-a*yy0[0]*(1.0-bb) + omg2*yy0[0] + 2.0*OMGE_GLO*yy0[4] + acc[0];
       ydot[4] = -MU_GLO*yy0[1]/r3-a*yy0[1]*(1.0-bb) + omg2*yy0[1] - 2.0*OMGE_GLO*yy0[3] + acc[1];
       ydot[5] = -MU_GLO*yy0[2]/r3-a*yy0[2]*(3.0-bb) + acc[2];

       //kout = ydot;//kout[]のアドレスを,ydotのアドレスにするこれだと,関数を出たら,ydotは消えるので,koutは更新されない
       for ( i = 0; i < 9; i++) kout[i] = ydot[i];
       return 0;
}

double newton(double mk, double e, int i)
{
	double ek;
	double a1,a2,a3;
	double tol = 1.0e-15;
	int iter, lmax = 4;
	int j=0;

	ek = 0.0;
	if(fabs(mk) > 0.0)
	{
		a1 = mk + e*sin(mk);
		iter = 1;
		while(iter <= lmax)
		{
			a2 = a1 - e*sin(a1) - mk;
			a3 = a1 - a2/2.0;
			a3 = 1.0 - e*cos(a3);
			ek = a1 - a2/a3;
			a2 = a1 - ek;
			if(fabs(a2) < tol){break ;}
			++iter;
			if(iter <= lmax)
			{
				a1 = ek;
			}
			else{
//				fprintf(fp[8],"Newton Error   PRN=%d\n",i);
				j = 0;
			}
		}
	}
	return(ek);
}
