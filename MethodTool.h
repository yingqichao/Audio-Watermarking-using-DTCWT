#pragma once
#include <stdlib.h> 
#include "DTCWT_BAND.h"
#include "WaterEncode.h"

class MethodTool
{
public:
	MethodTool();
	~MethodTool();
	void get16arr(int mlen, int flen);    //产生16个序列   并根据传入数据 确定嵌入哪一个矩阵
	void getsrand(vector<double>& in, int num, int m_wwidth, int m_wheight);

	//16个随机序列,嵌入字符用 
	vector<vector<double>> m_seq;
	vector<vector<double>>  f_seq;
	WaterEncode* encode;
	vector<vector<double>> m_seq_mid;
	vector<vector<double>> f_seq_mid;
	double energyCalc(vector<double>& v, int index, int len);
	double eng_ave;

};

