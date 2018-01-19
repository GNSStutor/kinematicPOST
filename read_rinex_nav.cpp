//////////////////////////////////////////////////////////////////////////////
//
//	read navigation file
//
//////////////////////////////////////////////////////////////////////////////
#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include "global_extern.h"

int header_flag = 0;
double sub_gpstime[PRN][256];

//—j“ú‚ğ‹‚ß‚é
int dayofweek(int y, int m, int d) /* 0 = Sunday , 1 <= m <= 12,  y > 1752 or so */
{
	return y -= m < 3, (y + y/4 - y/100 + y/400 + "032503514624"[m-1]-'0' + d) % 7;
}

const double f1  = 1602e+6;
const double df1 = 0.5625e+6;
const double f2  = 1246e+6;
const double df2 = 0.4375e+6;

void read_rinex_nav(int rcvn)
{
	char a,aa[512],buff[512],buff2[512],b;
	int i=0;
	int return_flag=0;
	int type_sat[5]={8,8,8,8,4},num=0;
	int sv,offset,type,dotw;
	int year,month,day,hour,minute;
	double second,kawari,check[PRN];
	static int t_number=0;
	static int prn_count[PRN]={0};

	Sunday = 0;//‰Šú’lİ’è

	while((a = fgetc(fp[2])) != EOF){

		if(header_flag<3){
			if(a == ' ' || a == '\n'){
				aa[i] ='\0';
				i=0;
			}//if
			else {
				aa[i] = a;
				i++;
			}//else

			//read header part
			if((strcmp(aa,"END")) == 0)
				header_flag = 1;
			if(header_flag == 1 && (strcmp(aa,"OF")) == 0)
				header_flag = 2;
			if(header_flag == 2 && (strcmp(aa,"HEADER")) == 0){
				header_flag = 3;
			}
		}
		else{
			i=0;
			buff[i++]=a;
			return_flag=0;
			while(return_flag != type_sat[num]){
					a = fgetc(fp[2]);
				buff[i++]=a;
				while(a != '\n'){
					a = fgetc(fp[2]);
					buff[i++]=a;
				}
				i=0;
				return_flag++;

				if(return_flag==1){//read line1
					b = buff[0];//type of sat
					if(b=='G'){//GPS     0
						num=0;offset=0;type=0;}
					if(b=='R'){//GLONASS 1
						num=4;offset=100;type=1;}
					if(b=='E'){//GALILEO 2
						num=1;offset=40;type=2;}
					if(b=='J'){//QZS     3
						num=2;offset=32;type=3;}
					if(b=='C'){//BEIDOU  4
						num=3;offset=70;type=4;}
				}

				if(type==0 || type==3 || type==4){//GPS, QZS, BEIDOU
					if(return_flag==1){
						for(i=1;i<=2;i++)//sv number
							buff2[i-1]=buff[i];
						buff2[2]='\0';sv=atoi(buff2)+offset;

						t_number=prn_count[sv];
						prn_count[sv]++;
						t_number = prn_count[sv]-1;

						Sub_E.prn[sv][t_number]=sv;

						for(i=4;i<=7;i++)//year
							buff2[i-4] = buff[i];
						buff2[4] = '\0';year = atoi(buff2);

						for(i=9;i<=10;i++)//month
							buff2[i-9] = buff[i];
						buff2[2] = '\0';month = atoi(buff2);

						for(i=12;i<=13;i++)//day
							buff2[i-12] = buff[i];
						buff2[2] = '\0';day = atoi(buff2);

						for(i=15;i<=16;i++)//hour
							buff2[i-15] = buff[i];
						buff2[2] = '\0';hour = atoi(buff2);

						for(i=18;i<=19;i++)//minute
							buff2[i-18] = buff[i];
						buff2[2] = '\0';minute = atoi(buff2);

						for(i=21;i<=22;i++)//second
							buff2[i-21] = buff[i];
						buff2[2] = '\0';second = atof(buff2);

						for(i=24;i<=41;i++)//af0
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.af0[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//af1
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.af1[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//af2
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.af2[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==2){//read line2
						for(i=5;i<=22;i++)//iode
							buff2[i-5] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);Sub_E.iode[sv][t_number]=(int)kawari;

						for(i=24;i<=41;i++)//crs
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.crs[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//dn
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.dn[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//m0
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.m0[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==3){//read line3
						for(i=5;i<=22;i++)//cuc
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.cuc[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)//e
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.e[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//cus
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.cus[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//roota
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.roota[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==4){//read line4
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.toe[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.cic[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.omega0[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.cis[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==5){//read line5
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.i0[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.crc[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.omega[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.domega0[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==6){//read line6
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.di0[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);i=0;
					}

					if(return_flag==7){//read line7
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.health[sv][t_number] = atoi(buff2);i=0;

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.tgd[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.iodc[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==8){//read line8
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';check[sv] = atof(buff2);
						dotw = dayofweek(2000+year,month,day);
						Sunday = dotw;
						
						Sub_E.toc[sv][t_number] =  
							(Sunday*86400)+
							((double)hour * 3600.0)+
							((double)minute * 60.0)+second;

						sub_gpstime[sv][t_number] = Sub_E.toc[sv][t_number];

						if(check[sv]<0)
							check[sv] = 604800+check[sv];
						sub_gpstime[sv][t_number] = check[sv];i=0;
										
						}
				}//GPS

				if(type==2){//GALILEO
					if(return_flag==1){
						for(i=1;i<=2;i++)//sv number
							buff2[i-1]=buff[i];
						buff2[2]='\0';sv=atoi(buff2)+offset;

						t_number=prn_count[sv];
						prn_count[sv]++;
						t_number = prn_count[sv]-1;

						Sub_E.prn[sv][t_number]=sv;

						for(i=4;i<=7;i++)//year
							buff2[i-4] = buff[i];
						buff2[4] = '\0';year = atoi(buff2);

						for(i=9;i<=10;i++)//month
							buff2[i-9] = buff[i];
						buff2[2] = '\0';month = atoi(buff2);

						for(i=12;i<=13;i++)//day
							buff2[i-12] = buff[i];
						buff2[2] = '\0';day = atoi(buff2);

						for(i=15;i<=16;i++)//hour
							buff2[i-15] = buff[i];
						buff2[2] = '\0';hour = atoi(buff2);

						for(i=18;i<=19;i++)//minute
							buff2[i-18] = buff[i];
						buff2[2] = '\0';minute = atoi(buff2);

						for(i=21;i<=22;i++)//second
							buff2[i-21] = buff[i];
						buff2[2] = '\0';second = atof(buff2);

						for(i=24;i<=41;i++)//af0
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.af0[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//af1
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.af1[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//af2
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.af2[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==2){//read line2
						for(i=5;i<=22;i++)//iode
							buff2[i-5] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);Sub_E.iode[sv][t_number]=(int)kawari;

						for(i=24;i<=41;i++)//crs
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.crs[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//dn
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.dn[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//m0
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.m0[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==3){//read line3
						for(i=5;i<=22;i++)//cuc
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.cuc[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)//e
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.e[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//cus
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.cus[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//roota
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.roota[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==4){//read line4
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.toe[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.cic[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.omega0[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.cis[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==5){//read line5
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.i0[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.crc[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.omega[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.domega0[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==6){//read line6
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.di0[sv][t_number] = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);i=0;
					}

					if(return_flag==7){//read line7
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						for(i=24;i<=41;i++)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.health[sv][t_number] = atoi(buff2);i=0;

						for(i=43;i<=60;i++)
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.tgd[sv][t_number] = atof(buff2);i=0;
					}

					if(return_flag==8){//read line8
						for(i=5;i<=22;i++)
							buff2[i-5] = buff[i];
						buff2[18] = '\0';check[sv] = atof(buff2);
						dotw = dayofweek(2000+year,month,day);
						Sunday = dotw;
						
						Sub_E.toc[sv][t_number] =  
							(Sunday*86400)+
							((double)hour * 3600.0)+
							((double)minute * 60.0)+second;

						sub_gpstime[sv][t_number] = Sub_E.toc[sv][t_number];

						if(check[sv]<0)
							check[sv] = 604800+check[sv];
						sub_gpstime[sv][t_number] = check[sv];
						for(i=23;i<=79;i++)
							buff[i]='\0';
							i=0;
					}//GALILEO
				}
				if(type==1){//Glonass
					if(return_flag==1){
						for(i=1;i<=2;i++)//sv number
							buff2[i-1]=buff[i];
						buff2[2]='\0';sv=atoi(buff2)+offset;

						t_number=prn_count[sv];
						prn_count[sv]++;
						t_number = prn_count[sv]-1;

						Sub_E.prn[sv][t_number]=sv;

						for(i=4;i<=7;i++)//year
							buff2[i-4] = buff[i];
						buff2[4] = '\0';year = atoi(buff2);

						for(i=9;i<=10;i++)//month
							buff2[i-9] = buff[i];
						buff2[2] = '\0';month = atoi(buff2);

						for(i=12;i<=13;i++)//day
							buff2[i-12] = buff[i];
						buff2[2] = '\0';day = atoi(buff2);

						for(i=15;i<=16;i++)//hour
							buff2[i-15] = buff[i];
						buff2[2] = '\0';hour = atoi(buff2);

						for(i=18;i<=19;i++)//minute	
							buff2[i-18] = buff[i];
						buff2[2] = '\0';minute = atoi(buff2);

						for(i=21;i<=22;i++)//second
							buff2[i-21] = buff[i];
						buff2[2] = '\0';second = atof(buff2);

						for(i=24;i<=41;i++)// SV clock bias (sec)  (-TauN)
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.TauN[sv][t_number] = atof(buff2);

						for(i=43;i<=60;i++)//SV relative frequency bias
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.GammaN[sv][t_number] = atof(buff2);

						for(i=62;i<=79;i++)//message frame time(sec of day UTC) 
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.tk[sv][t_number] = atof(buff2);i=0;
					}
					if(return_flag==2){//read line2
						for(i=5;i<=22;i++)//X
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.Xp[sv][t_number]=atof(buff2)*1000;

						for(i=24;i<=41;i++)//x
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.Xv[sv][t_number] = atof(buff2)*1000;

						for(i=43;i<=60;i++)//x
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.Xa[sv][t_number] = atof(buff2)*1000;

						for(i=62;i<=79;i++)//health
							buff2[i-62] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);i=0;
					}if(return_flag==3){//read line3
						for(i=5;i<=22;i++)//y
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.Yp[sv][t_number]=atof(buff2)*1000;

						for(i=24;i<=41;i++)//crs
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.Yv[sv][t_number] = atof(buff2)*1000;

						for(i=43;i<=60;i++)//dn
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.Ya[sv][t_number] = atof(buff2)*1000;

						for(i=62;i<=79;i++)//m0
							buff2[i-62] = buff[i];
						buff2[18] = '\0';kawari = atof(buff2);

						GF1[sv]=f1+kawari*df1;
//						GF2[sv]=f2+kawari*df2;i=0;

					}if(return_flag==4){//read line4
						for(i=5;i<=22;i++)//y
							buff2[i-5] = buff[i];
						buff2[18] = '\0';Sub_E.Zp[sv][t_number]=atof(buff2)*1000;

						for(i=24;i<=41;i++)//crs
							buff2[i-24] = buff[i];
						buff2[18] = '\0';Sub_E.Zv[sv][t_number] = atof(buff2)*1000;

						for(i=43;i<=60;i++)//dn
							buff2[i-43] = buff[i];
						buff2[18] = '\0';Sub_E.Za[sv][t_number] = atof(buff2)*1000;


						for(i=62;i<=79;i++)//health
							buff2[i-62] = buff[i];
						buff2[18] = '\0';Sub_E.health[sv][t_number] = atoi(buff2);i=0;

//						for(i=62;i<=79;i++)//Age of oper. information  (days) 
//							buff2[i-62] = buff[i];
//						buff2[18] = '\0';kawari = atof(buff2);
					
						dotw = dayofweek(2000+year,month,day);
						Sunday = dotw;
						
						Sub_E.toc[sv][t_number] =  
							(Sunday*86400)+
							((double)hour * 3600.0)+
							((double)minute * 60.0)+second;

						sub_gpstime[sv][t_number] = Sub_E.toc[sv][t_number];

		//				if(check[sv]<0)
			//				check[sv] = 604800+check[sv];
					//	sub_gpstime[sv][t_number] = check[sv];
						i=0;
					}//line4
					}//glonass
			}//while
			return_flag = 0;
		}
	}
}