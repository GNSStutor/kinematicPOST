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
int SVpos_flag[RCVN][PRN];
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
	int		gpsweek[PRN][256],iode[PRN][256],accuracy[PRN][256],prn[PRN][256],
			health[PRN][256],code[PRN][256],pflag[PRN][256];
	double	interval[PRN][256],gpstime[PRN][256],toc[PRN][256],af0[PRN][256],
			af1[PRN][256],af2[PRN][256],crs[PRN][256],dn[PRN][256],
			m0[PRN][256],cuc[PRN][256],e[PRN][256],cus[PRN][256],
			roota[PRN][256],toe[PRN][256],cic[PRN][256],omega0[PRN][256],
			cis[PRN][256],i0[PRN][256],crc[PRN][256],omega[PRN][256],
			domega0[PRN][256],di0[PRN][256],tgd[PRN][256],iodc[PRN][256];
	//for Glonass
	double	TauN[PRN][256],GammaN[PRN][256],tk[PRN][256],
			Xp[PRN][256],Xv[PRN][256],Xa[PRN][256],
			Yp[PRN][256],Yv[PRN][256],Ya[PRN][256],
			Zp[PRN][256],Zv[PRN][256],Za[PRN][256];

}nav_t;
nav_t Sub_E;
