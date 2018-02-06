////////////////////////////////////////////////////////////////////////////////
//
//　観測データの読み込み
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
	double kawari;
	int netr9_satn=0;

	for(j=0;j<=PRN-1;j++){
		Cp1[rcvn][j]=0;
		Pr1[rcvn][j]=0;
		Dp1[rcvn][j]=0;
		Cn1[rcvn][j]=0;
	}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
	//読み込みRINEX情報(まずRTKLIBのRTKCONV2.4.2でUBLOXのUBXファイルを変換)

	//RINEX versionは3.02
	//必要なSatellite Systemsをチェック
	//Observation Typesは全てチェック(C,L,D,S)
	//FrequenciesはL1のみ！

	//以上で、以下のflag=1で読み込む
		if(rcvn==1)
			flag=1;//上記で変換したタイプのRINEX
		else
			flag=1;//上記で変換したタイプのRINEX
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

	while(read_finish_flag != 1){
		if(forward[rcvn]==1){//前置きの部分（end of headerまで）
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
		else{//end of header以降
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
				getchar();				//////test
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
			//GPSTIMEの設定
			if(amount==29 && rcvn==1){
				dotw = dayofweek(year,month,day);
				GPSTIME=(dotw*86400)+
						((double)hour * 3600.0)+
						((double)minute * 60.0)+second;
			}
			//DGPSTIMEの設定
			if(amount==29 && rcvn==0){
				dotw = dayofweek(year,month,day);
				DGPSTIME=(dotw*86400)+
							((double)hour * 3600.0)+
							((double)minute * 60.0)+second;
			}


			if(flag==1 && amount>=35){//UBLOX M8 rinex3.02

				//実際の擬似距離や搬送波の読み込み
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
							sv_offset = 32;//+1なのでここは32そして下で33になる
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
			}//flag==1

			if(flag==22 && amount>=35){//

//R   16 C1C C1P C2C C2P D1C D1P D2C D2P L1C L1P L2C L2P S1C  SYS / # / OBS TYPES 
//       S1P S2C S2P                                          SYS / # / OBS TYPES 
//G   16 C1C C2W C2X C5X D1C D2W D2X D5X L1C L2W L2X L5X S1C  SYS / # / OBS TYPES 
//       S2W S2X S5X                                          SYS / # / OBS TYPES 
//J   20 C1  C1C C2X C5X C6X D1  D1C D2X D5X D6X L1  L1C L2X  SYS / # / OBS TYPES 
//       L5X L6X S1  S1C S2X S5X S6X                          SYS / # / OBS TYPES 
//C   12 C2I C6I C7I D2I D6I D7I L2I L6I L7I S2I S6I S7I      SYS / # / OBS TYPES 

				//実際の擬似距離や搬送波の読み込み
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

					//GPS
					if(type_sat==0){
						sv_offset = 0;
						SVn[rcvn][netr9_satn]=sv_number+sv_offset;
						netr9_satn++;

						for(i=0;i<=15;i++)//CODE
							buff2[i]=buff[i];
						buff2[15]='\0';
						Pr1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16;i<=16+15;i++)
							buff2[i-16]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*2;i<=16*2+15;i++)
							buff2[i-16*2]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*3;i<=16*3+15;i++)
							buff2[i-16*3]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*4;i<=16*4+15;i++)//DOPPLER
							buff2[i-16*4]=buff[i];
						buff2[15]='\0';
						Dp1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*5;i<=16*5+15;i++)
							buff2[i-16*5]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*6;i<=16*6+15;i++)
							buff2[i-16*6]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*7;i<=16*7+15;i++)
							buff2[i-16*7]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*8;i<=16*8+15;i++)//CARRIER
							buff2[i-16*8]=buff[i];
						buff2[15]='\0';
						Cp1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*9;i<=16*9+15;i++)
							buff2[i-16*9]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*10;i<=16*10+15;i++)
							buff2[i-16*10]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*11;i<=16*11+15;i++)
							buff2[i-16*11]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*12;i<=16*12+15;i++)//SN
							buff2[i-16*12]=buff[i];
						buff2[15]='\0';
						Cn1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*131;i<=16*13+15;i++)
							buff2[i-16*13]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*14;i<=16*14+15;i++)
							buff2[i-16*14]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*15;i<=16*15+15;i++)
							buff2[i-16*15]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
					}


//					/*
					//QZS    24
					if(type_sat==3){

						sv_offset = 32;

						SVn[rcvn][netr9_satn]=sv_number+sv_offset;
						netr9_satn++;


						for(i=0;i<=15;i++)//CODE
							buff2[i]=buff[i];
						buff2[15]='\0';
						Pr1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16;i<=16+15;i++)
							buff2[i-16]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*2;i<=16*2+15;i++)
							buff2[i-16*2]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*3;i<=16*3+15;i++)
							buff2[i-16*3]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*4;i<=16*4+15;i++)
							buff2[i-16*4]=buff[i];
						buff2[15]='\0';
	//					[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*6;i<=16*6+15;i++)//DOPPLER
							buff2[i-16*6]=buff[i];
						buff2[15]='\0';
						Dp1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*7;i<=16*7+15;i++)
							buff2[i-16*7]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*8;i<=16*8+15;i++)
							buff2[i-16*8]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*9;i<=16*9+15;i++)
							buff2[i-16*9]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*10;i<=16*10+15;i++)
							buff2[i-16*10]=buff[i];
						buff2[15]='\0';
//						Cp5[rcvn][sv_number]=atof(buff2);

						for(i=16*12;i<=16*12+15;i++)//CARRIER
							buff2[i-16*12]=buff[i];
						buff2[15]='\0';
						Cp1[rcvn][sv_number+sv_offset]=atof(buff2);
						for(i=16*13;i<=16*13+15;i++)//CARRIER
							buff2[i-16*13]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*14;i<=16*14+15;i++)//CARRIER
							buff2[i-16*14]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*15;i<=16*15+15;i++)//CARRIER
							buff2[i-16*15]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*16;i<=16*16+15;i++)//CARRIER
							buff2[i-16*16]=buff[i];
						buff2[15]='\0';
//						Cp1[rcvn][sv_number]=atof(buff2);

						for(i=16*18;i<=16*18+15;i++)//SN
							buff2[i-16*18]=buff[i];
						buff2[15]='\0';
						Cn1[rcvn][sv_number+sv_offset]=atof(buff2);
						for(i=16*19;i<=16*19+15;i++)//SN
							buff2[i-16*19]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*20;i<=16*20+15;i++)//SN
							buff2[i-16*20]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*21;i<=16*21+15;i++)//SN
							buff2[i-16*21]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*22;i<=16*22+15;i++)//SN
							buff2[i-16*22]=buff[i];
						buff2[15]='\0';
//						Cn1[rcvn][sv_number]=atof(buff2);

					}
//					*/
//					/*
					//BEIDOU no B3 data
					if(type_sat==6){

						sv_offset = 70;
						SVn[rcvn][netr9_satn]=sv_number+sv_offset;
						netr9_satn++;

						for(i=0;i<=15;i++)//CODE
							buff2[i]=buff[i];
						buff2[15]='\0';
						Pr1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16;i<=16+15;i++)
							buff2[i-16]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*2;i<=16*2+15;i++)
							buff2[i-16*2]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*3;i<=16*3+15;i++)//DOPPLER
							buff2[i-16*3]=buff[i];
						buff2[15]='\0';
						Dp1[rcvn][sv_number+sv_offset]=atof(buff2);

						for(i=16*4;i<=16*4+15;i++)
							buff2[i-16*4]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*5;i<=16*5+15;i++)
							buff2[i-16*5]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*6;i<=16*6+15;i++)//CARRIER
							buff2[i-16*6]=buff[i];
						buff2[15]='\0';
						Cp1[rcvn][sv_number+sv_offset]=atof(buff2);
						for(i=16*7;i<=16*7+15;i++)
							buff2[i-16*7]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*8;i<=16*8+15;i++)
							buff2[i-16*8]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);

						for(i=16*9;i<=16*9+15;i++)//SN
							buff2[i-16*9]=buff[i];
						buff2[15]='\0';
						Cn1[rcvn][sv_number+sv_offset]=atof(buff2);
						for(i=16*10;i<=16*10+15;i++)
							buff2[i-16*10]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
						for(i=16*11;i<=16*11+15;i++)
							buff2[i-16*11]=buff[i];
						buff2[15]='\0';
						kawari=atof(buff2);
					}
//					*/

				}//while

				read_finish_flag=1;
				SATn[rcvn]=netr9_satn;
			}//flag==22

		}//else
	}//while
}
