/////////////////////////////////////////
// FILEを閉じる
/////////////////////////////////////////
#include <stdio.h>

#include "global_extern.h"

void file_close(int RTK)
{
	if(RTK == 1){
		fclose(fp[0]);//rover obs
	}
	fclose(fp[1]);//ref obs
	fclose(fp[3]);//単独測位ファイル
	fclose(fp[4]);//DGNSSファイル
	fclose(fp[5]);//RTKファイル
	fclose(fp[6]);//衛星番号ファイル
	fclose(fp[7]);//テスト用ファイル
	fclose(fp[8]);//コメント用ファイル
	fclose(fp[2]);//navigation file

}