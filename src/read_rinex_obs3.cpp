////////////////////////////////////////////////////////////////////////////////
//
//Å@observation
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <string.h>

#include "global_extern.h"
#include "rtklib.h"
using namespace std;
const double cs     = 299792458.0;//
const double F1     = 1575420000.0;
const double F2     = 1227600000.0;

int dayofweek(int y, int m, int d);

/* set string without tail space ---------------------------------------------*/
static void setstr(char *dst, const char *src, int n)
{
	char *p = dst;
	const char *q = src;
	while (*q&&q<src + n) *p++ = *q++;
	*p-- = '\0';
	while (p >= dst&&*p == ' ') *p-- = '\0';
}


void read_rinex_obs3(int rcvn)
{
	int i=0,j,k=0,num=0,flag;
	int read_finish_flag=0;
	int sv[PRN]={0};
	int year,month,day,hour,minute,satn=0;
	int type_sat=0;
	int sv_number=0;
	double second;
	static int forward[2]={1,1};
	char a, buff[512],buff2[82];
	char obstype[4];
	int types=0;

	for(j=0;j<=PRN-1;j++){
		Cp1[rcvn][j]=0;   // L1: Carrier phase
		Pr1[rcvn][j]=0;   // L1: Pseudo-range
		Dp1[rcvn][j]=0;  // L1: Doppler
		Cn1[rcvn][j]=0;   // L1: SNR
		LLI1[rcvn][j]=0;   // L1: LLI
		Cp2[rcvn][j]=0;   // L2: Carrier phase
		Pr2[rcvn][j]=0;   // L2: Pseudo-range
		Dp2[rcvn][j]=0;  // L2: Doppler
		Cn2[rcvn][j]=0;   // L2: SNR
		LLI2[rcvn][j]=0;   // L2: LLI
		Cp5[rcvn][j]=0;   // L5: Carrier phase
		Pr5[rcvn][j]=0;   // L5: Pseudo-range
		Dp5[rcvn][j]=0;  // L5: Doppler
		Cn5[rcvn][j]=0;   // L5: SNR
		LLI5[rcvn][j]=0;   // L5: LLI
		SVn[rcvn][j]=0;   // SV number
	}

	while (read_finish_flag !=1) {
		if(forward[rcvn]==1){   // read RINEX header
			fgets(buff, MAXCHAR, fp[rcvn]);
			if (strstr(buff,"RINEX VERSION / TYPE")){
				ObsHead[rcvn].version=(int)str2num(buff, 0, 6);
				for (i=0;i<15;i++) {
					for (j=0;j<5;j++)
						ObsHead[rcvn].type[j][i]=-1;  // Initilization
				}
			}
			else if (strstr(buff,"END OF HEADER"))
				forward[rcvn]=0;    // exit read RINEX header
			else if (strstr(buff,"# / TYPES OF OBSERV")){  // RINEX 2.XX
				types=(int)str2num(buff, 0, 6);
				if (types>30) {
					cout<<buff<<endl;
					cout<<"***ERROR: Types of observation exceeds 30."<<endl;
					exit(1);
				}
				ObsHead[rcvn].types=types;
				for (j=0,k=10;j<types;j++,k+=6) {
					if (k>58) {
						fgets(buff, MAXCHAR, fp[rcvn]);
						k=10;
					}
					setstr(obstype,buff+k,2);
					if (strstr(obstype,"C1"))  // L1 C/A
						ObsHead[rcvn].type[0][0]=j;
					else if (strstr(obstype,"L1"))   // L1
						ObsHead[rcvn].type[0][1]=j;
					else if (strstr(obstype,"D1"))  // D1
						ObsHead[rcvn].type[0][2]=j;
					else if (strstr(obstype,"S1"))  // S1
						ObsHead[rcvn].type[0][3]=j;
					else if (strstr(obstype,"P2"))  // P2(Y)
						ObsHead[rcvn].type[0][4]=j;
					else if (strstr(obstype,"L2"))   // L2
						ObsHead[rcvn].type[0][5]=j;
					else if (strstr(obstype,"D2"))  // D2
						ObsHead[rcvn].type[0][6]=j;
					else if (strstr(obstype,"S2"))  // S2
						ObsHead[rcvn].type[0][7]=j;
					else if (strstr(obstype,"C5"))  // P5
						ObsHead[rcvn].type[0][8]=j;
					else if (strstr(obstype,"L5"))   // L5
						ObsHead[rcvn].type[0][9]=j;
					else if (strstr(obstype,"D5"))  // D5
						ObsHead[rcvn].type[0][10]=j;
					else if (strstr(obstype,"S5"))  // S5
						ObsHead[rcvn].type[0][11]=j;
				}
			}
			else if (strstr(buff,"SYS / # / OBS TYPES")){  // RINEX 3.XX
				types=(int)str2num(buff, 1, 5);   // GJECR
				for (j=0,k=7;j<types;j++,k+=4) {
					if (k>58) {
						fgets(buff, MAXCHAR, fp[rcvn]);
						k=7;
					}
					setstr(obstype,buff+k,3);
					if (buff[0]=='G') {
						if (strstr(obstype,"C1C"))  // L1 C/A
							ObsHead[rcvn].type[0][0]=j;
						else if (strstr(obstype,"L1C"))   // L1
							ObsHead[rcvn].type[0][1]=j;
						else if (strstr(obstype,"D1C"))  // D1
							ObsHead[rcvn].type[0][2]=j;
						else if (strstr(obstype,"S1C"))  // S1
							ObsHead[rcvn].type[0][3]=j;
						else if (strstr(obstype,"C2W"))  // P2(Y)
							ObsHead[rcvn].type[0][4]=j;
						else if (strstr(obstype,"L2W"))   // L2
							ObsHead[rcvn].type[0][5]=j;
						else if (strstr(obstype,"D2W"))  // D2
							ObsHead[rcvn].type[0][6]=j;
						else if (strstr(obstype,"S2W"))  // S2
							ObsHead[rcvn].type[0][7]=j;
						else if (strstr(obstype,"C5Q"))  // P5
							ObsHead[rcvn].type[0][8]=j;
						else if (strstr(obstype,"L5Q"))   // L5
							ObsHead[rcvn].type[0][9]=j;
						else if (strstr(obstype,"D5W"))  // D5
							ObsHead[rcvn].type[0][10]=j;
						else if (strstr(obstype,"S5Q"))  // S5
							ObsHead[rcvn].type[0][11]=j;
					}
					else if (buff[0]=='J') {
						if (strstr(obstype,"C1C"))  // P1
							ObsHead[rcvn].type[1][0]=j;
						else if (strstr(obstype,"L1C"))   // L1
							ObsHead[rcvn].type[1][1]=j;
						else if (strstr(obstype,"D1C"))  // D1
							ObsHead[rcvn].type[1][2]=j;
						else if (strstr(obstype,"S1C"))  // S1
							ObsHead[rcvn].type[1][3]=j;
						else if (strstr(obstype,"C2X"))  // P2
							ObsHead[rcvn].type[1][4]=j;
						else if (strstr(obstype,"L2X"))   // L2
							ObsHead[rcvn].type[1][5]=j;
						else if (strstr(obstype,"D2X"))  // D2
							ObsHead[rcvn].type[1][6]=j;
						else if (strstr(obstype,"S2X"))  // S2
							ObsHead[rcvn].type[1][7]=j;
						else if (strstr(obstype,"C5X"))  // P5
							ObsHead[rcvn].type[1][8]=j;
						else if (strstr(obstype,"L5X"))   // L5
							ObsHead[rcvn].type[1][9]=j;
						else if (strstr(obstype,"D5X"))  // D5
							ObsHead[rcvn].type[1][10]=j;
						else if (strstr(obstype,"S5X"))  // S5
							ObsHead[rcvn].type[1][11]=j;
					}
					else if (buff[0]=='E') {
						if (strstr(obstype,"C1C"))  // P1
							ObsHead[rcvn].type[2][0]=j;
						else if (strstr(obstype,"L1C"))   // L1
							ObsHead[rcvn].type[2][1]=j;
						else if (strstr(obstype,"D1C"))  // D1
							ObsHead[rcvn].type[2][2]=j;
						else if (strstr(obstype,"S1C"))  // S1
							ObsHead[rcvn].type[2][3]=j;
						else if (strstr(obstype,"C5X"))  // P5a
							ObsHead[rcvn].type[2][4]=j;
						else if (strstr(obstype,"L5X"))   // L5a
							ObsHead[rcvn].type[2][5]=j;
						else if (strstr(obstype,"D5X"))  // D5a
							ObsHead[rcvn].type[2][6]=j;
						else if (strstr(obstype,"S5X"))  // S5a
							ObsHead[rcvn].type[2][7]=j;
						else if (strstr(obstype,"C7X"))  // P5b
							ObsHead[rcvn].type[2][8]=j;
						else if (strstr(obstype,"L7X"))   // L5b
							ObsHead[rcvn].type[2][9]=j;
						else if (strstr(obstype,"D7X"))  // D5b
							ObsHead[rcvn].type[2][10]=j;
						else if (strstr(obstype,"S7X"))  // S5b
							ObsHead[rcvn].type[2][11]=j;
					}
					else if (buff[0]=='C') {
						if ((strstr(obstype,"C1I")) || (strstr(obstype,"C1X")) || (strstr(obstype,"C2I")) || (strstr(obstype,"C2X"))) //  B1
							ObsHead[rcvn].type[3][0]=j;  // C2 is for RINEX 3.01/3.03 and C1 is for RINEX 3.02
						else if ((strstr(obstype,"L1I")) || (strstr(obstype,"L1X")) || (strstr(obstype,"L2I")) || (strstr(obstype,"L2X"))) 
							ObsHead[rcvn].type[3][1]=j;
						else if ((strstr(obstype,"D1I")) || (strstr(obstype,"D1X")) || (strstr(obstype,"D2I")) || (strstr(obstype,"D2X"))) 
							ObsHead[rcvn].type[3][2]=j;
						else if ((strstr(obstype,"S1I")) || (strstr(obstype,"S1X")) || (strstr(obstype,"S2I")) || (strstr(obstype,"S2X"))) 
							ObsHead[rcvn].type[3][3]=j;
						else if ((strstr(obstype,"C7I")) || (strstr(obstype,"C7X"))) // B3
							ObsHead[rcvn].type[3][4]=j;
						else if ((strstr(obstype,"L7I")) || (strstr(obstype,"L7X")))
							ObsHead[rcvn].type[3][5]=j;
						else if ((strstr(obstype,"D7I")) || (strstr(obstype,"D7X")))
							ObsHead[rcvn].type[3][6]=j;
						else if ((strstr(obstype,"S7I")) || (strstr(obstype,"S7X")))
							ObsHead[rcvn].type[3][7]=j;
						else if ((strstr(obstype,"C6I")) || (strstr(obstype,"C6X"))) // B3
							ObsHead[rcvn].type[3][8]=j;
						else if ((strstr(obstype,"L6I")) || (strstr(obstype,"L6X")))
							ObsHead[rcvn].type[3][9]=j;
						else if ((strstr(obstype,"D6I")) || (strstr(obstype,"D6X")))
							ObsHead[rcvn].type[3][10]=j;
						else if ((strstr(obstype,"S6I")) || (strstr(obstype,"S6X")))
							ObsHead[rcvn].type[3][11]=j;
					}
					else if (buff[0]=='R') {
						if (strstr(obstype,"C1C"))  // P1
							ObsHead[rcvn].type[4][0]=j;
						else if (strstr(obstype,"L1C"))   // L1
							ObsHead[rcvn].type[4][1]=j;
						else if (strstr(obstype,"D1C"))  // D1
							ObsHead[rcvn].type[4][2]=j;
						else if (strstr(obstype,"S1C"))  // S1
							ObsHead[rcvn].type[4][3]=j;
						else if (strstr(obstype,"C2C"))  // P2  C2P???
							ObsHead[rcvn].type[4][4]=j;
						else if (strstr(obstype,"L2C"))   // L2
							ObsHead[rcvn].type[4][5]=j;
						else if (strstr(obstype,"D2C"))  // D2
							ObsHead[rcvn].type[4][6]=j;
						else if (strstr(obstype,"S2C"))  // S2
							ObsHead[rcvn].type[4][7]=j;
						else if (strstr(obstype,"C3I"))  // P3
							ObsHead[rcvn].type[4][8]=j;
						else if (strstr(obstype,"L3I"))   // L3
							ObsHead[rcvn].type[4][9]=j;
						else if (strstr(obstype,"D3I"))  // D3
							ObsHead[rcvn].type[4][10]=j;
						else if (strstr(obstype,"S3I"))  // S3
							ObsHead[rcvn].type[4][11]=j;
					}
				}  // for (j=0,k=7;j<types;j++,k+=4) {
			}  // else if (strstr(buff,"SYS / # / OBS TYPES")){ 
		}  // if(forward[rcvn]==1)

		  // *****************start to read RINEX body***********************
		else{
			a = fgetc(fp[rcvn]);
			fgets(buff, MAXCHAR, fp[rcvn]);
			if(rcvn==1 && a==EOF){
				cout << "Read File End" << endl;
				read_finish_flag = 1;
				exit(1);
			}
			if(rcvn==0 && a==EOF){
				cout << "rover file finish" << endl;
				read_finish_flag = 1;
				exit(1);
			}
			if (ObsHead[rcvn].version==2) {
				year=(int)str2num(buff, 0, 2)+2000;
				month=(int)str2num(buff, 2, 3);
				day=(int)str2num(buff, 5, 3);
				hour=(int)str2num(buff, 8, 3);
				minute=(int)str2num(buff, 11, 3);
				second=str2num(buff, 14, 11);
				flag=(int)str2num(buff, 25, 3);
				satn=(int)str2num(buff, 28, 3);
				if (flag>1) {  // Special event 
					for (i=0;i<flag;i++)
						fgets(buff, MAXCHAR, fp[rcvn]);
				}
				
				int dotw;
				//GPSTIME/DGPSTIME
				dotw = dayofweek(year,month,day);
				if(rcvn==1)
					GPSTIME=(dotw*86400)+((double)hour * 3600.0)+((double)minute * 60.0)+second;
				else if(rcvn==0)
					DGPSTIME=(dotw*86400)+((double)hour * 3600.0)+	((double)minute * 60.0)+second;
								
				for (j=0,k=31;j<satn;j++,k+=3) {    // read satellite SVN
					if (k>68) {
						fgets(buff, MAXCHAR, fp[rcvn]);
						k=32;
					}
					sv_number=(int)str2num(buff, k+1, 2);
					if (buff[k]=='G')
						SVn[rcvn][j]=sv_number;
					else if (buff[k]=='J')
						SVn[rcvn][j]=sv_number+32;
					else if (buff[k]=='E')
						SVn[rcvn][j]=sv_number+40;
					else if (buff[k]=='C')
						SVn[rcvn][j]=sv_number+70;
					else if (buff[k]=='R')
						SVn[rcvn][j]=sv_number+100;
				}
				for (i=0,k=0;i<satn;i++,k+=16) {   // read observations of each satellite
					fgets(buff, MAXCHAR, fp[rcvn]);
					for (j=strlen(buff);j<=80;j++)
						buff[j]=' ';
					buff[80]='\0';
					if (ObsHead[rcvn].types>5) {
						fgets(buff2, MAXCHAR, fp[rcvn]);
						for (j=strlen(buff2);j<=79;j++)
							buff2[j]=' ';
						buff2[80]='\0';
						strcat(buff, buff2);
						if (ObsHead[rcvn].types>10) {
							fgets(buff2, MAXCHAR, fp[rcvn]);
							for (j=strlen(buff2);j<=80;j++)
								buff2[j]=' ';
							buff2[80]='\0';
							strcat(buff,buff2);
							if (ObsHead[rcvn].types>15) {
								fgets(buff2, MAXCHAR, fp[rcvn]);
								for (j=strlen(buff2);j<=80;j++)
									buff2[j]=' ';
								buff2[80]='\0';
								strcat(buff,buff2);
								if (ObsHead[rcvn].types>20) {
									fgets(buff2, MAXCHAR, fp[rcvn]);
									for (j=strlen(buff2);j<=80;j++)
										buff2[j]=' ';
									buff2[80]='\0';
									strcat(buff,buff2);
									if (ObsHead[rcvn].types>25) {
										fgets(buff2, MAXCHAR, fp[rcvn]);
										for (j=strlen(buff2);j<=80;j++)
											buff2[j]=' ';
										buff2[80]='\0';
										strcat(buff,buff2); 
									}
								}
							}
						}
					}
					if (ObsHead[rcvn].type[0][0]!=-1)
						Pr1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][0]*16, 14);   // L1: Pseudo-range
					if (ObsHead[rcvn].type[0][1]!=-1)
						Cp1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][1]*16, 14);;   // L1: Carrier phase
					if (ObsHead[rcvn].type[0][1]!=-1)
						LLI1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][1]*16+14, 1);;   // L1: LLI
					if (ObsHead[rcvn].type[0][2]!=-1)
						Dp1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][2]*16, 14);;  // L1: Doppler
					if (ObsHead[rcvn].type[0][3]!=-1)
						Cn1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][3]*16, 14);;   // L1: SNR

					if (ObsHead[rcvn].type[0][4]!=-1)
						Pr2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][4]*16, 14);   // L2: Pseudo-range
					if (ObsHead[rcvn].type[0][5]!=-1)
						Cp2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][5]*16, 14);;   // L2: Carrier phase
					if (ObsHead[rcvn].type[0][5]!=-1)
						LLI2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][5]*16+14, 1);;   // L2: LLI
					if (ObsHead[rcvn].type[0][6]!=-1)
						Dp2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][6]*16, 14);;  // L2: Doppler
					if (ObsHead[rcvn].type[0][7]!=-1)
						Cn2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][7]*16, 14);;   // L2: SNR

					if (ObsHead[rcvn].type[0][8]!=-1)
						Pr5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][8]*16, 14);   // L5: Pseudo-range
					if (ObsHead[rcvn].type[0][9]!=-1)
						Cp5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][9]*16, 14);;   // L5: Carrier phase
					if (ObsHead[rcvn].type[0][9]!=-1)
						LLI5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][9]*16+14, 1);;   // L5: LLI
					if (ObsHead[rcvn].type[0][10]!=-1)
						Dp5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][10]*16, 14);;  // L5: Doppler
					if (ObsHead[rcvn].type[0][11]!=-1)
						Cn5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[0][11]*16, 14);;   // L5: SNR	

				}
			} //if (ObsHead[rcvn].version==2)
			else if (ObsHead[rcvn].version==3) {
				year=(int)str2num(buff, 1, 4);
				month=(int)str2num(buff, 5, 3);
				day=(int)str2num(buff, 8, 3);
				hour=(int)str2num(buff, 11, 3);
				minute=(int)str2num(buff, 14, 3);
				second=str2num(buff, 17, 11);
				flag=(int)str2num(buff, 28, 3);
				satn=(int)str2num(buff, 31, 3);
				if (flag>1) {  // Special event 
					for (i=0;i<flag;i++)
						fgets(buff, MAXCHAR, fp[rcvn]);
				}
				
				int dotw;
				//GPSTIME/DGPSTIME
				dotw = dayofweek(year,month,day);
				if(rcvn==1)
					GPSTIME=(dotw*86400)+((double)hour * 3600.0)+((double)minute * 60.0)+second;
				else if(rcvn==0)
					DGPSTIME=(dotw*86400)+((double)hour * 3600.0)+	((double)minute * 60.0)+second;

				for (i=0;i<satn;i++) {   // read observations of each satellite
					fgets(buff, MAXCHAR, fp[rcvn]);
					sv_number=(int)str2num(buff, 1, 2);
					if (buff[0]=='G') {
						type_sat=0;
						SVn[rcvn][i]=sv_number;
					}
					else if (buff[0]=='J'){
						type_sat=1;
						SVn[rcvn][i]=sv_number+32;
					}
					else if (buff[0]=='E'){
						type_sat=2;
						SVn[rcvn][i]=sv_number+40;
					}
					else if (buff[0]=='C'){
						type_sat=3;
						SVn[rcvn][i]=sv_number+70;
					}
					else if (buff[0]=='R'){
						type_sat=4;
						SVn[rcvn][i]=sv_number+100;
					}
					
					if (ObsHead[rcvn].type[type_sat][0]!=-1)
						Pr1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][0]*16+3, 14);   // L1: Pseudo-range
					if (ObsHead[rcvn].type[type_sat][1]!=-1)
						Cp1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][1]*16+3, 14);;   // L1: Carrier phase
					if (ObsHead[rcvn].type[type_sat][1]!=-1)
						LLI1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][1]*16+17, 1);;   // L1: LLI
					if (ObsHead[rcvn].type[type_sat][2]!=-1)
						Dp1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][2]*16+3, 14);;  // L1: Doppler
					if (ObsHead[rcvn].type[type_sat][3]!=-1)
						Cn1[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][3]*16+3, 14);;   // L1: SNR
					
					if (ObsHead[rcvn].type[type_sat][4]!=-1)
						Pr2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][4]*16+3, 14);   // L2: Pseudo-range
					if (ObsHead[rcvn].type[type_sat][5]!=-1)
						Cp2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][5]*16+3, 14);;   // L2: Carrier phase
					if (ObsHead[rcvn].type[type_sat][5]!=-1)
						LLI2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][5]*16+17, 1);;   // L2: LLI
					if (ObsHead[rcvn].type[type_sat][6]!=-1)
						Dp2[rcvn][i]=str2num(buff, ObsHead[rcvn].type[type_sat][6]*16+3, 14);;  // L2: Doppler
					if (ObsHead[rcvn].type[type_sat][7]!=-1)
						Cn2[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][7]*16+3, 14);;   // L2: SNR
					
					if (ObsHead[rcvn].type[type_sat][8]!=-1)
						Pr5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][8]*16+3, 14);   // L5: Pseudo-range
					if (ObsHead[rcvn].type[type_sat][9]!=-1)
						Cp5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][9]*16+3, 14);;   // L5: Carrier phase
					if (ObsHead[rcvn].type[type_sat][9]!=-1)
						LLI5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][9]*16+17, 1);;   // L5: LLI
					if (ObsHead[rcvn].type[type_sat][10]!=-1)
						Dp5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][10]*16+3, 14);;  // L5: Doppler
					if (ObsHead[rcvn].type[type_sat][11]!=-1)
						Cn5[rcvn][SVn[rcvn][i]]=str2num(buff, ObsHead[rcvn].type[type_sat][11]*16+3, 14);;   // L5: SNR					
				} //for (i=0;i<satn;i++) {  
			} //else if (ObsHead[rcvn].version==3) {			
			read_finish_flag=1;
			SATn[rcvn]=satn;
		}  //end of read rinex body
	}

}