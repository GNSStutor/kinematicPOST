///////////////////////////////////////////////////////////////////
//  
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include "global_extern.h"

void choose_sat(int rcvn, int iter)
{
	int i,j,k,stock_sat[PRN],sat=0,prn,stock[PRN]={0};
////////  selection of satellite  ////////////
	//reject code from receiver
	//minimum elevation angle
	//minimum carrier to noise ratio
	//SVn_sat[rcvn][prn]
//////////////////////////////////////////////

	if(rcvn == 1){
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if(Elevation[rcvn][prn] <= Elevation_mask1 ||
				Cn1[rcvn][prn] <= Threshold_cn ||
				fabs(Pr1[rcvn][prn]) <= 1.0 ||//
				fabs(Cp1[rcvn][prn]) <= 1.0 ||//
				LLI[rcvn][prn] >= 1 ||
				SVn_sat[rcvn][prn] == 0 ||//
				Ephe.health[prn] > 0.9)//
				i=i;
			else{
				stock_sat[sat] = SVn[rcvn][i];
				sat++;

				if(prn>=101 && prn<=130)
					Glo_Num++;
				if(prn>=71 && prn<=100)
					Bei_Num++;
				if(prn>=1 && prn<=39)
					Gps_Num++;
				if(prn>=41 && prn<=70)
					Gal_Num++;

			}
		}
		for(i=0;i<sat;i++){
			SVn[rcvn][i] = stock_sat[i];
		}
		SATn[rcvn] = sat;
	}

	sat = 0;
	//�ړ��ǁiDGPS���j�̃f�[�^�ŏd�ݕW���΍����ُ�l�̏ꍇ�͉q���r��
	//SATA��reject_code������
	if(rcvn == 0){
		for(i=0;i<SATn[rcvn];i++){
			prn = SVn[rcvn][i];
			if (
				Elevation[rcvn][prn] <= Elevation_mask1 ||
				Cn1[rcvn][prn] <= Threshold_cn ||
				SVn_sat[rcvn][prn] == 0 ||

				fabs(Pr1[rcvn][prn]) <= 1.0 ||//�[�������̃`�F�b�N]
				fabs(Dp1[rcvn][prn]) <= 1.0 ||//�����g�̃`�F�b�N
				fabs(Cp1[rcvn][prn]) <= 1.0 ||//�����g�̃`�F�b�N
				LLI[rcvn][prn] >= 1 ||
				prn > 1000)//�ǂ��ł��������
				i = i;
			else{
				stock_sat[sat] = SVn[rcvn][i];
				sat++;
				fprintf(fp[7], ",%d", prn);

				
				//if(prn>=101 && prn<=130)
				//	Glo_Num++;
				//if(prn>=71 && prn<=100)
				//	Bei_Num++;
				//if(prn>=1 && prn<=33)
				//	Gps_Num++;
				//if(prn>=41 && prn<=70)
				//	Gal_Num++;

			}
		}
		for(i=0;i<sat;i++){
			SVn[rcvn][i] = stock_sat[i];
		}
		SATn[rcvn] = sat;
	}


	//DGPS�̏ꍇ�q������v�����Ȃ���΂Ȃ�Ȃ�
	if(rcvn == 0){
		k = 0;
		for(i=0;i<SATn[1];i++){
			for(j=0;j<SATn[0];j++){
				if(SVn[1][i] == SVn[0][j]){
					stock_sat[k++] = SVn[0][j];
					prn = SVn[0][j];

					if(prn>=101 && prn<=130)
						Glo_Num++;
					if(prn>=71 && prn<=100)
						Bei_Num++;
					if(prn>=1 && prn<=39)
						Gps_Num++;
					if(prn>=41 && prn<=70)
						Gal_Num++;
				}
			}
		}
		for(i=0;i<k;i++){
			SVn[0][i] = stock_sat[i];
//			prn=SVn[0][i];
//			stock[prn]=1;
//			if(SLIP[rcvn][prn]==1)
//				stock[prn]=2;
		}
		SATn[0] = k;
	}


	if(rcvn==1){//����̗��p�q���V�X�e��
		MinSvNum=3;
		for(i=0;i<5;i++)
			Systemflag[i]=0;
		
		for(i=0;i<SATn[rcvn];i++){	
				prn=SVn[rcvn][i];
			if(prn<=39)
				Systemflag[0]=1;//GPS+QZSS
			else if(prn>=41&&prn<=70)
				Systemflag[4]=1;//GALILEO
			else if(prn>=71&&prn<=100)
				Systemflag[3]=1;//BeiDou
			else if(prn>=101&&prn<=130)
				Systemflag[2]=1;//GLONASS
		}
		for(i=0;i<5;i++)
			MinSvNum+=Systemflag[i];
	}

	if(rcvn==0){//�ړ����̗��p�q���V�X�e���i�{����Ƃ̋��ʉq���Ƃ���j
		MinSvNum=3;
		for(i=0;i<5;i++)
			Systemflag[i]=0;
		
		for(i=0;i<SATn[rcvn];i++){	
				prn=SVn[rcvn][i];
			if(prn<=39)
				Systemflag[0]=1;//GPS+QZSS
			else if(prn>=41&&prn<=70)
				Systemflag[4]=1;//GALILEO
			else if(prn>=71&&prn<=100)
				Systemflag[3]=1;//BeiDou
			else if(prn>=101&&prn<=130)
				Systemflag[2]=1;//GLONASS
		}
		for(i=0;i<5;i++)
			MinSvNum+=Systemflag[i];
	}

}
