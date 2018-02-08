//////////////////////////////////////////////////////////////////////////
//
// �P�Ƒ��ʌv�Z
//
//////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <stdio.h>
#include <math.h>
#include "global_extern.h"

void trans_coordinates(int);
void least_square(int,double[],int,int);
#define PI 3.1415926535897932384626433832795

const double cs     = 299792458.0;//����
const double myu    = 3.986005e+14;//
const double omegae = 7.2921151467e-5;//
const double pi     = 3.1415926535898;//�~����
const double f      = -4.442807633e-10;//
const double F1     = 1575420000.0;
const double B1     = 1561098000.0;

//���ʉ��Z�����̃g�b�v
void calc_pos(int rcvn, int iter_main, int pos)
{
	int i,j;
	int iter = 0;
	static int iter_pos = 1;
	double init[7];
	double ref[3]={0};
	
	init[0] = 0;//X�iECEF�j�̏����l
	init[1] = 0;//Y�iECEF�j�̏����l
	init[2] = 0;//Z�iECEF�j�̏����l
//	init[0] = -3960000.0;//X�iECEF�j�̏����l
//	init[1] = 3350000.0;//Y�iECEF�j�̏����l
//	init[2] = 3698000.0;//Z�iECEF�j�̏����l
	init[3] = 0.0;//��M�@�N���b�N�덷
	init[4] = 0.0;//GPS��1�Ԗڂ̑��ʃV�X�e���Ƃ̎��v��
	init[5] = 0.0;//GPS��2�Ԗڂ̑��ʃV�X�e���Ƃ̎��v��
	init[6] = 0.0;//GPS��3�Ԗڂ̑��ʃV�X�e���Ƃ̎��v��
	LS_update[0]=10;LS_update[1]=10;LS_update[2]=10;

	if(SATn[rcvn]>=MinSvNum){
		while(fabs(LS_update[0])>0.001 || fabs(LS_update[1])>0.001 || fabs(LS_update[2])>0.001)
//		while(iter<=6)
		{
			least_square(rcvn, init, iter, pos);
			iter = iter+1;
		}
	}
	else{
		printf("%f,***4�q������***\n",GPSTIME);
	}

	//���ʉ��Z���ʂ��擾
	POSx[rcvn] = init[0];
	POSy[rcvn] = init[1];
	POSz[rcvn] = init[2];
	Clock[rcvn] = init[3];
	Clock5[rcvn] = init[4];
	Clock6[rcvn] = init[5];
	Clock7[rcvn] = init[6];
	Clock_ext[rcvn] = Clock[rcvn];

	//��M�@���v�̌덷����␳����
	for(i=0;i<SATn[rcvn];i++){
		j=SVn[rcvn][i];

		Pr1[rcvn][j]=Pr1[rcvn][j]-Clock[rcvn];

		//�����ł�RTK�p�ɉq�����v���̃N���b�N�ϓ��̉e����␳����
		//��ǂƈړ��ǂŕ␳�f�[�^�̒x���iAGE�j������̂ŁA������l��
		if(rcvn==0){
			Pr1[1][j] = Pr1ref[j] - cs*(DGPSTIME-GPSTIME)*Ephe.af1[j];

			if(j>=1 && j<=70)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.af1[j]*F1/cs);
			if(j>=71 && j<=100)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.af1[j]*B1/cs);
			if(j>=101 && j<=135)
				Cp1[1][j] = Cp1ref[j] - cs*(DGPSTIME-GPSTIME)*(Ephe.GammaN[j]*GF1[j]/cs);

		}

		if(j<=70 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/F1);
		if(j<=70 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/F1);

		if(j>=71 && j<=100 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/B1);
		if(j>=71 && j<=100 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/B1);

		if(j>=101 && j<=130 && rcvn==1 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/GF1[j]);
		if(j>=101 && j<=130 && rcvn==0 && POS==1)
			Cp1[rcvn][j]=Cp1[rcvn][j]-Clock[rcvn]/(cs/GF1[j]);
	}

	if(rcvn==1)
		GPSTIME=GPSTIME-Clock[rcvn]/cs;
	if(rcvn==0)
		DGPSTIME=DGPSTIME-Clock[rcvn]/cs;

}
