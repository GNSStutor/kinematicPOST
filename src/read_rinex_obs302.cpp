////////////////////////////////////////////////////////////////////////////////
//
//@observation
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>

#include "global_extern.h"
using namespace std;
const double cs     = 299792458.0;//
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;

int dayofweek(int y, int m, int d);


void read_rinex_obs302(int rcvn)
{
	int i=0,j,k=0,amount=0,num=0,num2=0,flag,total=0;
	int read_finish_flag=0,return_flag=0,satn_count=0;
	int sv[PRN]={0};
	int no_blank=0;
	int year,month,day,hour,minute,satn=0;
	int first_return=0;
	int oem4_flag=0,ii=0;
	int pre_amount=0,amount1=0;
	int type_sat=0;
	int set_inc=32;
	int count_return=0;
	int sv_number=0,sv_offset=0,gps_l2c_flag=0;
	double quality=0;
	double second;
	static int forward[2]={1,1};
	char a,buff[512],buff2[512];

	for(j=0;j<=PRN-1;j++){
		Cp1[rcvn][j]=0;
		Pr1[rcvn][j]=0;
		Dp1[rcvn][j]=0;
		Cn1[rcvn][j]=0;
	}

	//1:UBLOX 8T RINEX300
		if(rcvn==1)
			flag=26;
		else
			flag=26;


	while(read_finish_flag != 1){
		if(forward[rcvn]==1){//‘O’u‚«‚Ì•”•ªiend of header‚Ü‚Åj
			a = fgetc(fp[rcvn]);
			if( a == ' ' || a == '\n'){
				buff[i]='\0';
				i=0;
			}
			else{
				buff[i++]=a;
			}
			if((strcmp(buff,"HEADER")) == 0){
				forward[rcvn]=0;
				while(a != '\n')
					a = fgetc(fp[rcvn]);
				i=0;
			}
		}
		else{//end of headerˆÈ~
			a = fgetc(fp[rcvn]);

			if(rcvn==1 && a==EOF){
				cout << "Read File End" << endl;
				read_finish_flag = 1;
				End_flag = 1;
				exit(1);
			}
			if(rcvn==0 && a==EOF){
				cout << "rover file finish" << endl;
				read_finish_flag = 1;
				End_flag = 1;
				exit(1);
			}

			buff[i++]=a;
			amount++;
			amount1++;

			switch(amount){
				case 6: buff[0]=' ';buff[i]='\0';i=0;year=atoi(buff);break;
				case 9: buff[i]='\0';i=0;month=atoi(buff);break;
				case 12: buff[i]='\0';i=0;day=atoi(buff);break;
				case 15: buff[i]='\0';i=0;hour=atoi(buff);break;
				case 18: buff[i]='\0';i=0;minute=atoi(buff);break;
				case 29: buff[i]='\0';i=0;second=atof(buff);break;
				case 32: buff[i]='\0';i=0;break;
				case 35: buff[i]='\0';i=0;satn=atoi(buff);
						while(a != '\n'){
							a = fgetc(fp[rcvn]);
						}break;
				default: i=i;break;
			}

			int dotw;
			//GPSTIME‚Ìİ’è
			if(amount==29 && rcvn==1){
				dotw = dayofweek(year,month,day);
				GPSTIME=(dotw*86400)+
						((double)hour * 3600.0)+
						((double)minute * 60.0)+second;
			}
			//DGPSTIME‚Ìİ’è
			if(amount==29 && rcvn==0){
				dotw = dayofweek(year,month,day);
				DGPSTIME=(dotw*86400)+
							((double)hour * 3600.0)+
							((double)minute * 60.0)+second;
			}


			if(flag==26 && amount>=35){//UBLOX M8 rinex3.02

				//ÀÛ‚Ì‹[—‹——£‚â”À‘—”g‚Ì“Ç‚İ‚İ
				while(return_flag != satn){
					i=0;

					a = fgetc(fp[rcvn]);//G R E S J C
					if(a == 'G' || a == 'R' || a=='E' || a=='J' || a=='W' || a=='S' || a=='C'){
						if(a=='G')
							type_sat=0;
						else if(a=='R')
							type_sat=1;
						else if(a=='E')
							type_sat=2;
						else if(a=='J')
							type_sat=3;
						else if(a=='W')
							type_sat=4;
						else if(a=='S')
							type_sat=5;
						else if(a=='C')
							type_sat=6;
						else
							cout << "error in read_rinex_obs" << endl;
					}

					i=0;
					a = fgetc(fp[rcvn]);buff[i++]=a;
					a = fgetc(fp[rcvn]);buff[i++]=a;buff[2]='\0';
					sv_number = atoi(buff);

					i=0;
					while(a != '\n'){
						a = fgetc(fp[rcvn]);
						buff[i++]=a;
					}
					return_flag++;

					//GPS QZS GLO BEI
					if(type_sat==0 || type_sat==3 || type_sat==1 || type_sat==6 || type_sat==2){
						sv_offset = 0;
						if(type_sat==1)
							sv_offset = 100;
						else if(type_sat==3)
							sv_offset = 32;//+1‚È‚Ì‚Å‚±‚±‚Í32‚»‚µ‚Ä‰º‚Å33‚É‚È‚é
						else if(type_sat==6)
							sv_offset = 70;
						else if(type_sat==2)
							sv_offset = 40;

						SVn[rcvn][return_flag-1]=sv_number+sv_offset;

						for(i=0;i<=15;i++)
							buff2[i]=buff[i];
						buff2[15]='\0';
						Pr1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16;i<=16+15;i++)
							buff2[i-16]=buff[i];
						buff2[15]='\0';
						Cp1[rcvn][sv_number+sv_offset]=atof(buff2);


						LLI[rcvn][sv_number+sv_offset] = atoi(&buff2[14]);

 //|             | Loss of lock indicator (LLI). Range: 0-7        |            |
 //|             |  0 or blank: OK or not known                    |            |
 //|             |  Bit 0 set : lost lock between previous and     |            |
 //|             |     current observation: cycle slip possible    |            |
 //|             |  Bit 1 set : Inverse wavelength factor to       |            |
 //|             |     default (does NOT change default)           |            |
 //|             |  Bit 2 set : observation under Antispoofing     |            |
 //|             |              (may suffer from increased noise) 


						for(i=16*2;i<=16*2+15;i++)
							buff2[i-16*2]=buff[i];
						buff2[15]='\0';
						Dp1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*3;i<=16*3+15;i++)
							buff2[i-16*3]=buff[i];
						buff2[15]='\0';
						Cn1[rcvn][sv_number+sv_offset]=atof(buff2);
					}

				}//while

				read_finish_flag=1;
				SATn[rcvn]=satn;
			}//flag==26

		}//else
	}//while
}
