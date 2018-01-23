/////////////////////////////////////////////////////////////////////////////
//
// 航法メッセージ、観測データを読み込み
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
		read_rinex_nav(rcvn);//航法メッセージ読み込み(最初に全て1回で)
		nav_i = 1;
	}
	if(nav_i == 1){
		read_rinex_obs302(rcvn);//観測データ読み込み　毎エポック
		if(rcvn==1){
			rinex_time(rcvn);//航法メッセージの使用タイミングをGPS時刻より決める
		}
	}
}
