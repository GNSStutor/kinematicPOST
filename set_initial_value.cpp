////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////
#include <math.h>
#include <iostream>

#include "global_extern.h"

using namespace std;

void set_initial_value()
{
	int i,j,k;
	char c,buff[128];
	char file1[128],file2[128],file3[128];

	if((fp[8] = fopen("initial_setting.txt","r")) == NULL){
		cout << "ファイルを開けない initial_setting.txt" << endl;exit(1);
	}

	for(k=1; k<=18; k++){
		for(j=0; j<100; j++){
			i = 0;
			while((c = fgetc(fp[8])) != ',' && c != '\n' && c != EOF && c != ' ')
			{
				buff[i++] = c;
			}
			buff[i] = '\0';

			if(c=='\n'){
				switch(k){
					case 1: Elevation_mask1=atof(buff); break;
					case 2: strcpy(file1,buff); break;
					case 3: strcpy(file2,buff); break;
					case 4: strcpy(file3,buff); break;
					case 5: POSrcvlat[1]=atof(buff); break;
					case 6: POSrcvlon[1]=atof(buff); break;
					case 7: POSrcvhgt[1]=atof(buff); break;
					case 8: Code_noise=atof(buff); break;
					case 9: Carrier_noise=atof(buff); break;
					case 10: Iteration=atoi(buff); break;
					case 11: Threshold_cn=atof(buff); break;
					case 12: RTK=atoi(buff); break;
					case 13: Ratio_limit=atof(buff); break;
					case 14: Gflag = atoi(buff);  break;
					case 15: Jflag = atoi(buff);  break;
					case 16: Eflag = atoi(buff);  break;
					case 17: Cflag = atoi(buff);  break;
					case 18: Rflag = atoi(buff);  break;
					default: cout << "error in initial setting" << endl;break;
				}
				j=100;
			}
		}
	}

/////////////////////////////////////////////////////////////////////////////////
	//////////////// 緯度、経度、高度を地心座標系に変換しておく /////////////////
	int rcvn=1;//reference
    double a,b,f,n,pi;
	pi = 3.14159265359;	
	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);
	n = a*a/sqrt(a*a*cos(POSrcvlat[rcvn]*pi/180)*
		cos(POSrcvlat[rcvn]*pi/180)+b*b*sin(POSrcvlat[rcvn]*pi/180)
		*sin(POSrcvlat[rcvn]*pi/180));

	Ref_pos[0][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*cos(POSrcvlon[rcvn]*pi/180);
	Ref_pos[1][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*sin(POSrcvlon[rcvn]*pi/180);
	Ref_pos[2][rcvn] = (n*b*b/a/a+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		sin(POSrcvlat[rcvn]*pi/180);

	rcvn=0;//rover
	pi = 3.14159265359;	
	a = 6378137.0;
	f = 1/298.257223563;
	b = a*(1-f);
	n = a*a/sqrt(a*a*cos(POSrcvlat[rcvn]*pi/180)*
		cos(POSrcvlat[rcvn]*pi/180)+b*b*sin(POSrcvlat[rcvn]*pi/180)
		*sin(POSrcvlat[rcvn]*pi/180));

	Ref_pos[0][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*cos(POSrcvlon[rcvn]*pi/180);
	Ref_pos[1][rcvn] = (n+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		cos(POSrcvlat[rcvn]*pi/180)*sin(POSrcvlon[rcvn]*pi/180);
	Ref_pos[2][rcvn] = (n*b*b/a/a+POSrcvhgt[rcvn]+POSrcvund[rcvn])*
		sin(POSrcvlat[rcvn]*pi/180);
/////////////////////////////////////////////////////////////////////////////////

	//電離層モデルの係数を設定する（手入力、全ての衛星はこれで計算）
	Arp[0]=1.118E-08;
	Arp[1]=1.49E-08;
	Arp[2]=-5.96E-08;
	Arp[3]=-5.96E-08;
	Bet[0]=8.806E+04;
	Bet[1]=1.638E+04;
	Bet[2]=-1.966E+05;
	Bet[3]=-1.311E+05;

	User_lat = POSrcvlat[1];//移動側の初期緯度
	User_lon = POSrcvlon[1];//移動側の初期経度
	User_pos[0] = Ref_pos[0][1];//移動側の初期X
	User_pos[1] = Ref_pos[1][1];//移動側の初期Y
	User_pos[2] = Ref_pos[2][1];//移動側の初期Z

	LeapSecond = 18.0;//うるう秒の設定

////////////////////// 観測データの読み込み　ファイル指定 /////////////////////////
///*
	if((fp[0] = fopen(file1,"r")) == NULL){
		cout << "ファイルを開けない0" << endl;exit(1);
	}
	//base
	if((fp[1] = fopen(file2,"r")) == NULL){
		cout << "ファイルを開けない1" << endl;exit(1);
	}

	if((fp[2] = fopen(file3,"r")) == NULL){
		cout << "ファイルを開けない2" << endl;exit(1);
	}
//*/

/*
	if((fp[0] = fopen("0623\\rovR.obs","r")) == NULL){
		cout << "ファイルを開けない0" << endl;exit(1);
	}
	//base
	if((fp[1] = fopen("0623\\refR.obs","r")) == NULL){//11140924test.14o
		cout << "ファイルを開けない1" << endl;exit(1);
	}

	if((fp[2] = fopen("0623\\nav1.rnx","r")) == NULL){
		cout << "ファイルを開けない19" << endl;exit(1);
	}
*/
////////////////////// 観測データの読み込み　ファイル指定 /////////////////////////

////////////////////// 出力する結果のためのファイル設定   /////////////////////////

	//単独測位結果
	if((fp[3] = fopen("pos.csv","w")) == NULL){
		cout << "出力ファイルを開けない3" << endl;exit(1);
	}
	//DGNSS結果
	if((fp[4] = fopen("dgnss.csv","w")) == NULL){
		cout << "出力ファイルを開けない4" << endl;exit(1);
	}
	//RTK結果
	if((fp[5] = fopen("rtk.csv","w")) == NULL){
		cout << "出力ファイルを開けない5"<< endl;exit(1);
	}
	//利用衛星情報
	if((fp[6] = fopen("sv.csv","w")) == NULL){
		cout << "出力ファイルを開けない6" << endl;exit(1);
	}
	//テスト用
	if((fp[7] = fopen("test.csv","w")) == NULL){
		cout << "出力ファイルを開けない7" << endl;exit(1);
	}

////////////////////// 出力する結果のためのファイル設定   /////////////////////////

}
