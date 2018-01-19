///////////////////////////////////////////////////////////////////////////////
//
//　外部変数の宣言
//
///////////////////////////////////////////////////////////////////////////////

void read_data(int);
void calc_satpos(int);
void calc_pos(int,int,int);
void calc_tropo(int);
void calc_direction(int,int);
void choose_sat(int,int);
void set_initial_value();
void file_close(int);
void calc_iono_model(int);
void calc_pos2(int,int,int);
void calc_rtk_GQER(int);
void calc_rtk_GQEB(int);
void calc_rtk_GQE(int);
void trans_xyz_llh(double,double,double, double []);
int lambda(int, int, const double[], const double[], double[], double[]);

#define RCVN 21
#define PRN  140
FILE *fp[RCVN];

int Rinex,prn;
int SATn[RCVN],SVn[RCVN][PRN];//使用衛星数、使用衛星番号
int SVn_sat[RCVN][PRN];
int Iteration;
int Sunday,End_flag,POS;
int No_rtk,Max_prn;
int LLI[RCVN][PRN];
int Sol_flag[3],GPSWEEK;
int Gal_Num,Glo_Num,Gps_Num,Bei_Num;
int MinSvNum;
int Systemflag[5];
int Gflag,Rflag,Jflag,Cflag,Eflag;
int RTK;

double Elevation_mask1;
double Threshold_cn,HDOP,VDOP,PDOP,GDOP,TDOP,Ratio_limit;
double Pr1[RCVN][PRN];//擬似距離
double Cp1[RCVN][PRN];//搬送波位相
double Cn1[RCVN][PRN],Dp1[RCVN][PRN];//信号強度、ドップラ周波数
double SVx[RCVN][PRN],SVy[RCVN][PRN],SVz[RCVN][PRN];//衛星位置
double SV_corrtime[RCVN][PRN],Iono[RCVN][PRN],Tropo[RCVN][PRN];
double SV_omomi[RCVN][PRN];
double GPSTIME,POSx[RCVN],POSy[RCVN],POSz[RCVN];
double POSrcvlat[RCVN],POSrcvlon[RCVN],POSrcvhgt[RCVN],POSrcvund[RCVN];
double DGPSTIME;
double Elevation[RCVN][PRN],Azimuth[RCVN][PRN];
double Ref_pos[3][4];//アンテナの精密位置
double Arp[4],Bet[4],Clock[RCVN],Clock5[RCVN],Clock6[RCVN],Clock7[RCVN];
double Code_noise,Carrier_noise;
double Interval,Ref_lon,Ref_lat,Ref_hgt;
double Ratio;
double Clock_ext[RCVN];
double LS_update[3];
double OmomiP[PRN],OmomiC[PRN],OmomiPB[PRN],OmomiCB[PRN];
double GF1[PRN];
double Pr1ref[PRN],Cp1ref[PRN];
double User_lat,User_lon,User_hgt,User_pos[3];
double LeapSecond;

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
peph_t Ephe;

typedef struct nav{
	int		gpsweek[PRN][128],iode[PRN][128],accuracy[PRN][128],prn[PRN][128],
			health[PRN][128],code[PRN][128],pflag[PRN][128];
	double	interval[PRN][128],gpstime[PRN][128],toc[PRN][128],af0[PRN][128],
			af1[PRN][128],af2[PRN][128],crs[PRN][128],dn[PRN][128],
			m0[PRN][128],cuc[PRN][128],e[PRN][128],cus[PRN][128],
			roota[PRN][128],toe[PRN][128],cic[PRN][128],omega0[PRN][128],
			cis[PRN][128],i0[PRN][128],crc[PRN][128],omega[PRN][128],
			domega0[PRN][128],di0[PRN][128],tgd[PRN][128],iodc[PRN][128];
	//for Glonass
	double	TauN[PRN][128],GammaN[PRN][128],tk[PRN][128],
			Xp[PRN][128],Xv[PRN][128],Xa[PRN][128],
			Yp[PRN][128],Yv[PRN][128],Ya[PRN][128],
			Zp[PRN][128],Zv[PRN][128],Za[PRN][128];

}nav_t;
nav_t Sub_E;
