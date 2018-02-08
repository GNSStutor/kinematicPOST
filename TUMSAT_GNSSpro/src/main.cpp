/////////////////////////////////////////////////////////////////////////////////
//
// UBLOX(M8P,M8T)の受信機観測データを利用したRTKの演算
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
	int iter;//読み込み回数
	int rcvn;//rcvn=1（基準側）　rcvn=0（移動側）
	int pos=1;//？
	int rover_kaisu=0;
	End_flag = 0;

	cout.precision(7);
	Sol_flag[0]=0;Sol_flag[1]=0;Sol_flag[2]=0;//移動側の読み込み回数、単独測位回数、RTKのFIX回数

	set_initial_value();//初期設定＋ファイル設定
//	cout << "GPSTIME" << " " << "Ref" << " " << "Rov" << " " << "全回数" << " " << "測位回数" << " " << "FIX回数" << endl;
	printf("  GPSTIME    RefSAT RovSAT      ALLcount     POSTcount     FIXcount\n");

	for(iter=1;iter<=Iteration;iter++){//設定回数分読み込む

		rcvn = 1;//1:基準側 0:移動側
		read_data(rcvn);//観測データ、航法メッセージの読み込み

		if(End_flag==1){//ファイル読み込み最後で強制終了
			iter = Iteration;continue;
		}

		calc_satpos(rcvn);//エフェメリスからの衛星位置計算
		calc_direction(rcvn,iter);//仰角・方位角計算（基準側なのでユーザ位置は固定の基準位置で計算）
		calc_iono_model(rcvn);//電離層モデルでの電離層遅延量計算
		calc_tropo(rcvn);//対流圏遅延量計算

		Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;//各衛星システムの数を初期化
		choose_sat(rcvn,iter);//衛星選択
		
		//単独測位を行う（GPS衛星が4機以上あることが前提）
		POS=1;//単独測位
		if(SATn[rcvn]>=MinSvNum && Gps_Num>=4){
			calc_pos(rcvn,iter,pos);//最小二乗法での単独測位
			calc_satpos(rcvn);//受信機クロック誤差により擬似距離が変わるため、再度発射時刻を再計算し、衛星位置も再計算する
			calc_pos2(rcvn,iter,pos);//少し修正した衛星位置での単独測位再計算
		}

		//RTKに入る場合移動側の観測データを読み込み、時刻を比較（基準側と移動側）してRTKの演算へ進む
		if(RTK == 1){
			rcvn=0;//1:基準側 0:移動側
			if(GPSTIME>=604000)
				GPSTIME=GPSTIME-604800;
			
			//基準局観測データの時刻と移動局観測データの時刻を比較し読み込む
			while((GPSTIME-DGPSTIME)>0.05 || rover_kaisu==0){
				read_data(0);//移動側観測データの読み込み
				rover_kaisu++;
			}

			//時刻判断処理部分
//			if(fabs(GPSTIME - DGPSTIME) < 0.1 && SATn[1]>=0){//同じGPS時刻
			while(DGPSTIME-GPSTIME >= -0.005 && DGPSTIME-GPSTIME <= 0.995){//基準側が1Hz、移動側が5Hzなどのとき
				read_data(0);//whileのときだけ
				int i;

				rover_kaisu++;
				Sol_flag[0]++;//移動側の全回数カウント

				calc_satpos(0);///エフェメリスからの衛星位置計算
				calc_direction(0,iter);//仰角・方位角計算（移動側なのでユーザ位置を基準）
				calc_iono_model(rcvn);//電離層モデルでの電離層遅延量計算
				calc_tropo(0);//対流圏遅延量の計算

				Glo_Num=0;Gps_Num=0;Gal_Num=0;Bei_Num=0;
				choose_sat(0,iter);//衛星選択

				//移動側の単独測位を行う（GPS衛星が4機以上あることが前提）
				POS=1;//単独測位
				if(SATn[rcvn]>=MinSvNum && Gps_Num>=4){
					calc_pos(rcvn,iter,pos);//最小二乗法での単独測位
					calc_satpos(rcvn);//受信機クロック誤差により擬似距離が変わるため、再度発射時刻を再計算し、衛星位置も再計算する
					calc_pos2(rcvn,iter,pos);//少し修正した衛星位置での単独測位再計算
					Sol_flag[1]++;//単独測位の回数カウント
				}

				//RTKを行う部分(GPS/QZS/GALILEOの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=5 && Bei_Num==0 && Glo_Num==0)
					calc_rtk_GQE(rcvn);

				//RTKを行う部分(GPS/QZS/GALILEO+BeiDouの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=6 && Gps_Num+Gal_Num>=2 && Bei_Num>=2 && Gps_Num>=2)
					calc_rtk_GQEB(rcvn);

				//RTKを行う部分(GPS/QZS/GALILEO+GLONASSの場合)　あわせてDGNSSも行う
				if(SATn[rcvn]>=6 && Gps_Num+Gal_Num>=2 && Glo_Num>=2 && Gps_Num>=2)
					calc_rtk_GQER(rcvn);

			}//基準局時刻にあわせた移動側のタイミングでのifまたはwhile文の終わり

		}//DGNSS＋RTK

		if (((int)iter % 1) == 0 && GPSTIME >= -0.1 && DGPSTIME >= -0.1)//途中経過の書き出し
	//		cout << GPSTIME << " " << SATn[1] << " " << SATn[0] << " " << Sol_flag[0] << " " << Sol_flag[1] << " " << Sol_flag[2] << endl;
			printf("%10.4f     %3d    %3d    %10d    %10d   %10d\n", GPSTIME, SATn[1], SATn[0], Sol_flag[0], Sol_flag[1], Sol_flag[2]);

		
	}//ここまでが繰り返し計算部


	file_close(RTK);//ファイルを閉じる
	return(0);
}