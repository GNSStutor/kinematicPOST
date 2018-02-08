/////////////////////////////////////////////////////////////////////////////////
//
// UBLOX(M8P,M8T)�̎�M�@�ϑ��f�[�^�𗘗p����RTK�̉��Z
// GPS/QZSS L1, GALILEO E1, BeiDou B1, GLONASS G1
//
//		ver.1.3
//
/////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "global.h"

using namespace std;

int main(){
	int iter;//�ǂݍ��݉�
	int rcvn;//rcvn=1�i����j�@rcvn=0�i�ړ����j
	int pos=1;//�H
	int rover_kaisu=0;
	End_flag = 0;

	cout.precision(7);
	Sol_flag[0]=0;Sol_flag[1]=0;Sol_flag[2]=0;//�ړ����̓ǂݍ��݉񐔁A�P�Ƒ��ʉ񐔁ARTK��FIX��

	set_initial_value();//�����ݒ�{�t�@�C���ݒ�
//	cout << "GPSTIME" << " " << "Ref" << " " << "Rov" << " " << "�S��" << " " << "���ʉ�" << " " << "FIX��" << endl;
	printf("  GPSTIME    RefSAT RovSAT      ALLcount     POSTcount     FIXcount\n");

	for(iter=1;iter<=Iteration;iter++){//�ݒ�񐔕��ǂݍ���

		rcvn = 1;//1:��� 0:�ړ���
		read_data(rcvn);//�ϑ��f�[�^�A�q�@���b�Z�[�W�̓ǂݍ���

		if(End_flag==1){//�t�@�C���ǂݍ��ݍŌ�ŋ����I��
			iter = Iteration;continue;
		}

		calc_satpos(rcvn);//�G�t�F�����X����̉q���ʒu�v�Z
		calc_direction(rcvn,iter);//�p�E���ʊp�v�Z�i����Ȃ̂Ń��[�U�ʒu�͌Œ�̊�ʒu�Ōv�Z�j
		calc_iono_model(rcvn);//�d���w���f���ł̓d���w�x���ʌv�Z
		calc_tropo(rcvn);//�Η����x���ʌv�Z

		Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;//�e�q���V�X�e���̐���������
		choose_sat(rcvn,iter);//�q���I��
		
		//�P�Ƒ��ʂ��s���iGPS�q����4�@�ȏ゠�邱�Ƃ��O��j
		POS=1;//�P�Ƒ���
		if(SATn[rcvn]>=MinSvNum && Gps_Num>=4){
			calc_pos(rcvn,iter,pos);//�ŏ����@�ł̒P�Ƒ���
			calc_satpos(rcvn);//��M�@�N���b�N�덷�ɂ��[���������ς�邽�߁A�ēx���ˎ������Čv�Z���A�q���ʒu���Čv�Z����
			calc_pos2(rcvn,iter,pos);//�����C�������q���ʒu�ł̒P�Ƒ��ʍČv�Z
		}

		//RTK�ɓ���ꍇ�ړ����̊ϑ��f�[�^��ǂݍ��݁A�������r�i����ƈړ����j����RTK�̉��Z�֐i��
		if(RTK == 1){
			rcvn=0;//1:��� 0:�ړ���
			if(GPSTIME>=604000)
				GPSTIME=GPSTIME-604800;
			
			//��Ǌϑ��f�[�^�̎����ƈړ��Ǌϑ��f�[�^�̎������r���ǂݍ���
			while((GPSTIME-DGPSTIME)>0.05 || rover_kaisu==0){
				read_data(0);//�ړ����ϑ��f�[�^�̓ǂݍ���
				rover_kaisu++;
			}

			//�������f��������
//			if(fabs(GPSTIME - DGPSTIME) < 0.1 && SATn[1]>=0){//����GPS����
			while(DGPSTIME-GPSTIME >= -0.005 && DGPSTIME-GPSTIME <= 0.995){//�����1Hz�A�ړ�����5Hz�Ȃǂ̂Ƃ�
				read_data(0);//while�̂Ƃ�����
				int i;

				rover_kaisu++;
				Sol_flag[0]++;//�ړ����̑S�񐔃J�E���g

				calc_satpos(0);///�G�t�F�����X����̉q���ʒu�v�Z
				calc_direction(0,iter);//�p�E���ʊp�v�Z�i�ړ����Ȃ̂Ń��[�U�ʒu����j
				calc_iono_model(rcvn);//�d���w���f���ł̓d���w�x���ʌv�Z
				calc_tropo(0);//�Η����x���ʂ̌v�Z

				Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
				choose_sat(0,iter);//�q���I��

				//�ړ����̒P�Ƒ��ʂ��s���iGPS�q����4�@�ȏ゠�邱�Ƃ��O��j
				POS=1;//�P�Ƒ���
				if(SATn[rcvn]>=MinSvNum && Gps_Num>=4){
					calc_pos(rcvn,iter,pos);//�ŏ����@�ł̒P�Ƒ���
					calc_satpos(rcvn);//��M�@�N���b�N�덷�ɂ��[���������ς�邽�߁A�ēx���ˎ������Čv�Z���A�q���ʒu���Čv�Z����
					calc_pos2(rcvn,iter,pos);//�����C�������q���ʒu�ł̒P�Ƒ��ʍČv�Z
					Sol_flag[1]++;//�P�Ƒ��ʂ̉񐔃J�E���g
				}

				//RTK���s������(GPS/QZS/GALILEO�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=5 && Bei_Num==0 && Glo_Num==0)
					calc_rtk_GQE(rcvn);

				//RTK���s������(GPS/QZS/GALILEO+BeiDou�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=6 && Gps_Num+Gal_Num>=2 && Bei_Num>=2 && Gps_Num>=2)
					calc_rtk_GQEB(rcvn);

				//RTK���s������(GPS/QZS/GALILEO+GLONASS�̏ꍇ)�@���킹��DGNSS���s��
				if(SATn[rcvn]>=6 && Gps_Num+Gal_Num>=2 && Glo_Num>=2 && Gps_Num>=2)
					calc_rtk_GQER(rcvn);

			}//��ǎ����ɂ��킹���ړ����̃^�C�~���O�ł�if�܂���while���̏I���

		}//DGNSS�{RTK

		if (((int)iter % 1) == 0 && GPSTIME >= -0.1 && DGPSTIME >= -0.1)//�r���o�߂̏����o��
	//		cout << GPSTIME << " " << SATn[1] << " " << SATn[0] << " " << Sol_flag[0] << " " << Sol_flag[1] << " " << Sol_flag[2] << endl;
			printf("%10.4f     %3d    %3d    %10d    %10d   %10d\n", GPSTIME, SATn[1], SATn[0], Sol_flag[0], Sol_flag[1], Sol_flag[2]);

		
	}//�����܂ł��J��Ԃ��v�Z��


	file_close(RTK);//�t�@�C�������
	return(0);
}