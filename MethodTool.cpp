
#include "MethodTool.h"
#include "pch.h"
#include <iostream>
#include <algorithm>

//using std::vector;
//using std::begin;
//using std::end;
//using std::cout;
using namespace std;

MethodTool::MethodTool()
{
}


MethodTool::~MethodTool()
{

}



/*
获得16个序列,并确定嵌入的字符水印
*/
void MethodTool::get16arr(int mlen, int flen)
{
	m_seq.resize(16); f_seq.resize(2);
	m_seq_mid.resize(16); f_seq_mid.resize(2);
	encode = new WaterEncode();
	cout << "Generating Mid-Freq for Pseudo-Random Sequences...\r\n## Note that this initialization is only required once! ##" << endl;

	int num = 216;
	for (int i = 0; i < 16; i++) {
		m_seq[i].resize(mlen); m_seq_mid[i].resize(mlen);
		getsrand(m_seq[i], i, mlen, 1);
		/*auto maxPosition = max_element(m_seq[i].begin(), m_seq[i].end());
		cout << *maxPosition << " at the postion of " << maxPosition - m_seq[i].begin() << endl;*/

		DTCWT_BANDS* dtcwt = encode->dtwavexfm_1d(m_seq[i], 0, mlen, 2);
		vector<complex<double>> ilist(dtcwt->BAND[1].size());
		dtcwt->BAND[1] = ilist;
		vector<double> ilist2(dtcwt->Lo.size()); 
		dtcwt->Lo = ilist2;
		encode->dtwaveifm_1d(dtcwt, mlen, 2, m_seq_mid[i], 0);
		if (i == 0)
			eng_ave = energyCalc(m_seq_mid[i], 0, m_seq_mid[i].size());
	}
	for (int i = 0; i < 2; i++) {
		f_seq[i].resize(flen); f_seq_mid[i].resize(flen);
		getsrand(f_seq[i], 16+i, flen, 1);
		DTCWT_BANDS* dtcwt = encode->dtwavexfm_1d(f_seq[i], 0, flen, 2);
		vector<complex<double>> ilist(dtcwt->BAND[1].size());
		dtcwt->BAND[1] = ilist;
		vector<double> ilist2(dtcwt->Lo.size());
		dtcwt->Lo = ilist2;
		encode->dtwaveifm_1d(dtcwt, 1024, 2, f_seq_mid[i], 0);
	}

}

double MethodTool::energyCalc(vector<double>& v, int index, int len) {
	double sum = 0;
	for (int i = index; i < index + len; i++)
		sum += v[i] * v[i];

	return sum;
}


/*
产生随机矩阵
*/
void MethodTool::getsrand(vector<double>& in, int num, int m_wwidth, int m_wheight)
{
	srand(num);
	for (int i = 0; i < m_wwidth*m_wheight; i++)
	{
		//int temp = rand() % 10;
		//double temp1 = (double)temp / 10;
		/*if (temp1 <= 0.5)
			temp1 = 0;
		else
			temp1 = 1;*/
		in[i] = (rand() % (99999 + 1) / (float)(99999 + 1))-0.1637;//
	}
}
