///////////////////////////////////////////////////////////////////////////////
//
//　外部変数の宣言
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>

#define RCVN 21
#define PRN  140
extern FILE *fp[RCVN];

extern int Rinex,prn;
extern int SATn[RCVN],SVn[RCVN][PRN];//使用衛星数、使用衛星番号
extern int SVn_sat[RCVN][PRN];
extern int Iteration;
extern int Sunday,End_flag,POS;
extern int No_rtk,Max_prn;
extern int LLI[RCVN][PRN];
extern int Sol_flag[3],GPSWEEK;
extern int Gal_Num,Glo_Num,Gps_Num,Bei_Num;
extern int MinSvNum;
extern int Systemflag[5];
extern int Gflag,Rflag,Jflag,Cflag,Eflag;
extern int RTK;

extern double Elevation_mask1;
extern double Threshold_cn,HDOP,VDOP,PDOP,GDOP,TDOP,Ratio_limit;
extern double Pr1[RCVN][PRN];//擬似距離
extern double Cp1[RCVN][PRN];//搬送波位相
extern double Cn1[RCVN][PRN],Dp1[RCVN][PRN];
extern double SVx[RCVN][PRN],SVy[RCVN][PRN],SVz[RCVN][PRN];//衛星位置
extern double SV_corrtime[RCVN][PRN],Iono[RCVN][PRN],Tropo[RCVN][PRN];
extern double SV_omomi[RCVN][PRN];
extern double GPSTIME,POSx[RCVN],POSy[RCVN],POSz[RCVN];
extern double POSrcvlat[RCVN],POSrcvlon[RCVN],POSrcvhgt[RCVN],POSrcvund[RCVN];
extern double DGPSTIME;
extern double Elevation[RCVN][PRN],Azimuth[RCVN][PRN];
extern double Ref_pos[3][4];//アンテナの精密位置
extern double Arp[4],Bet[4],Clock[RCVN],Clock5[RCVN],Clock6[RCVN],Clock7[RCVN];
extern double Code_noise,Carrier_noise;
extern double Interval,Ref_lon,Ref_lat,Ref_hgt;
extern double Ratio;
extern double Clock_ext[RCVN];
extern double LS_update[3];
extern double OmomiP[PRN],OmomiC[PRN],OmomiPB[PRN],OmomiCB[PRN];
extern double GF1[PRN];
extern double Pr1ref[PRN],Cp1ref[PRN];
extern double User_lat,User_lon,User_hgt,User_pos[3];
extern double LeapSecond;

typedef struct peph{
	int gpsweek[PRN],iode[PRN],accuracy[PRN],prn[PRN],
		health[PRN],code[PRN],pflag[PRN];
	double interval[PRN],gpstime[PRN],toc[PRN],af0[PRN],
		af1[PRN],af2[PRN],crs[PRN],dn[PRN],
		m0[PRN],cuc[PRN],e[PRN],cus[PRN],
		roota[PRN],toe[PRN],cic[PRN],omega0[PRN],
		cis[PRN],i0[PRN],crc[PRN],omega[PRN],
		domega0[PRN],di0[PRN],tgd[PRN],iodc[PRN];		
		double	TauN[PRN],GammaN[PRN],tk[PRN],
			Xp[PRN],Xv[PRN],Xa[PRN],
			Yp[PRN],Yv[PRN],Ya[PRN],
			Zp[PRN],Zv[PRN],Za[PRN];
}peph_t;
extern peph_t Ephe;

typedef struct nav{
	int		gpsweek[PRN][128],iode[PRN][128],accuracy[PRN][128],prn[PRN][128],
			health[PRN][128],code[PRN][128],pflag[PRN][128];
	double	interval[PRN][128],gpstime[PRN][128],toc[PRN][128],af0[PRN][128],
			af1[PRN][128],af2[PRN][128],crs[PRN][128],dn[PRN][128],
			m0[PRN][128],cuc[PRN][128],e[PRN][128],cus[PRN][128],
			roota[PRN][128],toe[PRN][128],cic[PRN][128],omega0[PRN][128],
			cis[PRN][128],i0[PRN][128],crc[PRN][128],omega[PRN][128],
			domega0[PRN][128],di0[PRN][128],tgd[PRN][128],iodc[PRN][128];	
	//for GLONASS
	double	TauN[PRN][128],GammaN[PRN][128],tk[PRN][128],
			Xp[PRN][128],Xv[PRN][128],Xa[PRN][128],
			Yp[PRN][128],Yv[PRN][128],Ya[PRN][128],
			Zp[PRN][128],Zv[PRN][128],Za[PRN][128];

}nav_t;
extern nav_t Sub_E;