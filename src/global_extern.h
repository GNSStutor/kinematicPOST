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
extern int SVpos_flag[RCVN][PRN];
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
	int		gpsweek[PRN][256],iode[PRN][256],accuracy[PRN][256],prn[PRN][256],
			health[PRN][256],code[PRN][256],pflag[PRN][256];
	double	interval[PRN][256],gpstime[PRN][256],toc[PRN][256],af0[PRN][256],
			af1[PRN][256],af2[PRN][256],crs[PRN][256],dn[PRN][256],
			m0[PRN][256],cuc[PRN][256],e[PRN][256],cus[PRN][256],
			roota[PRN][256],toe[PRN][256],cic[PRN][256],omega0[PRN][256],
			cis[PRN][256],i0[PRN][256],crc[PRN][256],omega[PRN][256],
			domega0[PRN][256],di0[PRN][256],tgd[PRN][256],iodc[PRN][256];	
	//for GLONASS
	double	TauN[PRN][256],GammaN[PRN][256],tk[PRN][256],
			Xp[PRN][256],Xv[PRN][256],Xa[PRN][256],
			Yp[PRN][256],Yv[PRN][256],Ya[PRN][256],
			Zp[PRN][256],Zv[PRN][256],Za[PRN][256];

}nav_t;
extern nav_t Sub_E;