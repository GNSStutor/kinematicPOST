/////////////////////////////////////////////////////////////////////////////
//
// �q�@���b�Z�[�W�A�ϑ��f�[�^��ǂݍ���
//
/////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include "global_extern.h"

void read_rinex_nav(int);
void read_rinex_obs302(int);
void rinex_time(int);

void read_data(int rcvn)
{
	static int nav_i = 0;

	if(nav_i == 0){
		read_rinex_nav(rcvn);//�q�@���b�Z�[�W�ǂݍ���(�ŏ��ɑS��1���)
		nav_i = 1;
	}
	if(nav_i == 1){
		read_rinex_obs302(rcvn);//�ϑ��f�[�^�ǂݍ��݁@���G�|�b�N
		if(rcvn==1){
			rinex_time(rcvn);//�q�@���b�Z�[�W�̎g�p�^�C�~���O��GPS������茈�߂�
		}
	}
}
