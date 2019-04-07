// audioTest.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

//Created by Qichao Ying on 2019-3-30
//Final Version Completed on 2019-4-5
//Git repository:yingqichao/audio-watermark-based-on-DTCWT latest commit 2019-4-6

#include "pch.h"
#include <iostream>
#include "AudioFile.h"
#include <string>
#include "MethodTool.h"
#include <complex>
#include <io.h>
#include <fstream>
#include "WaterEncode.h"
#include "DTCWT_BAND.h"
#include <map>
#include<ctime>
#include"maxheap.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "matlibXcorr/rt_nonfinite.h"
#include "matlibXcorr/rtwtypes.h"
#include "omp.h"
#include "matlibXcorr/matlibXcorr_types.h"
#include "entity.h"

//Public Methods:
void EmbedController(string a,string txt);
string ExtractController(string a);
//Private Methods:
void computecorrelation(vector<double> xnum, vector<double> ynum, double meanxnum, double meanynum, double stdeviationx, double stdeviationy, int const size, double& precorrelation, double& correlation, int ybias);
double Corrcoef(vector<double> xnum, vector<double> ynum, int ybias);
void computemean(double totalxnum, int const size, vector<double> xnum, double& meanxnum, vector<double> ynum, double totalynum, double& meanynum, int ybias);
void computevariance(double varianceaddy, double varianceaddx, vector<double> xnum, double meanxnum, double meanynum, double& variancey, double& variancex, int const size, vector<double> ynum, int ybias);
void computestandarddeviation(double& stdeviationx, double& stdeviationy, double variancey, double variancex,int ybias);

string dataToString(vector<int> data);
bool detectBark(vector<int> retrieved, int bark[]);
void genInfo(string a, vector<int>& info);
string read(string name);
void Xcorr(vector<double>& r, vector<double>& x, vector<double>& y, int startIndex, int xlen, int ylen);
int getInfo(int index, vector<int> info);
void saveFile(string retrieved, string foldpath);


using namespace std;
int main()
{
	clock_t startTime, endTime;
	startTime = clock();//计时开始

	//Complex numbers: complex<double> num(3.5, 4.5);
	EmbedController("ZOUHECAN.wav", "toBeEmbedded.txt");
	string extracted = ExtractController("audioFile.wav");
	
	cout << "Extracted Message: "<<extracted << endl << endl;
	cout << endl<<endl<< "Embedding Process Finished Successfully..." << endl<<endl;
	endTime = clock();//计时结束
	cout << "The run time is: " << (double)(endTime - startTime)/1000 << "s" << endl;
	saveFile(extracted, "Extracted.txt");

}

void EmbedController(string a,string txtname) {
	//Settings
	//string txtname = "toBeEmbedded.txt";
	int bark[] = { 1,1,1,0,0,1,0 };	
	int outputChannels = 2;
	double alpha = 0.2; int mlen = 1024; int flen = 1024;//mlen是信息帧的长度，flen是同步帧的长度
	//Initialization
	cout << "########### Initialization ################" << endl << endl;
	MethodTool *methodTool = new MethodTool();
	methodTool->get16arr(mlen, flen);
	WaterEncode *encode = new WaterEncode();
	string txt = read(txtname);
	cout << "Length of Texts: " << txt.length() << endl;
	vector<int> info;
	genInfo(txt,info);
	AudioFile<double> audioFile;
	audioFile.load(a);
	int sampleRate = audioFile.getSampleRate();
	int bitDepth = audioFile.getBitDepth();
	int numSamples = audioFile.getNumSamplesPerChannel();
	double lengthInSeconds = audioFile.getLengthInSeconds();
	int numChannels = audioFile.getNumChannels();
	bool isMono = audioFile.isMono();
	bool isStereo = audioFile.isStereo();
	//Print a summary to the console
	audioFile.printSummary();
	int numFrame = numSamples / mlen;
	cout << "Number of Frames: " << numFrame << endl;
	int reqFrame = 2*txt.length()+7;
	cout << "Frames Required: " << 3 * txt.length() + 7 << endl;
	int numSegment = numFrame / (3 * txt.length() + 7);
	cout << "Number of Segments: " << numSegment << endl;
	AudioFile<double>::AudioBuffer origin = audioFile.samples;
	//Create a new buffer for the output audio
	AudioFile<double>::AudioBuffer watered; 
	watered.resize(outputChannels);
	for (int i = 0; i < outputChannels; i++) {
		watered[i].resize(numSamples);
	}

	//###############Embedding###############
	int i = 0; int i_frame = 0; vector<double> wAdd = methodTool->f_seq_mid[0];
	for (; i < numFrame;i++) {
		if ((i + 1) % 3 == 0) {
			wAdd = methodTool->f_seq_mid[0];
		
		}
		else {
			int code = 0;
			//Get Code
			if (i_frame >= reqFrame - 7) { 
				code = bark[i_frame+7 - reqFrame]; i_frame++; }
			else { code = getInfo(i_frame * 4, info); i_frame ++; }
			cout << "Embedded codeword " << code << " at Frame: " << i << endl;
			//Select Watermark
			wAdd = methodTool->m_seq_mid[code];
			
		}

		vector<double> fraction(origin[0].begin()+i*mlen, origin[0].begin()+ (i+1)*mlen);
		DTCWT_BANDS* dtcwt = encode->dtwavexfm_1d(fraction, 0, mlen, 2);
		vector<complex<double>> ilist(dtcwt->BAND[1].size());
		dtcwt->BAND[1] = ilist;
		vector<double> ilist2(dtcwt->Lo.size());
		dtcwt->Lo = ilist2;
		encode->dtwaveifm_1d(dtcwt, mlen, 2, fraction, 0);

		double eng = methodTool->energyCalc(fraction, 0, mlen);
		/*if (numChannels == 2) {
			double eng2 = methodTool->energyCalc(origin[1], i* mlen, mlen);
			eng = (eng + eng2) / 2;
		}*/
		for (int j = 0; j < mlen; j++) {
			for (int c = 0; c < outputChannels; c++) {
				int ch = (c == 1 && numChannels == 2) ? 1 : 0;
				
				watered[c][i*mlen + j] = origin[ch][i*mlen+j]+max(0.2,eng)/methodTool->eng_ave*wAdd[j];
				
			}
		}

		if (i_frame >= reqFrame)	
			i_frame = 0;
		
	}

	//Fill the rest of the audio with the original one
	for (int j = i*mlen; j < numSamples; j++)
	{
		for (int c = 0; c < outputChannels; c++) {
			int ch = (c == 1 && numChannels == 2) ? 1 : 0;
			watered[c][j] = origin[ch][j];
		}
	}

	//Save and Quit

	bool ok = audioFile.setAudioBuffer(watered);
	cout <<endl<<endl<< "Saving File...."<< endl;
	// Wave file (implicit)
	audioFile.save("audioFile.wav");


}

string ExtractController(string a) {
	int outputChannels = 2;
	double alpha = 0.2; int mlen = 1024; int flen = 1024;//mlen是信息帧的长度，flen是同步帧的长度
	//Initialization
	int bark[] = { 1,1,1,0,0,1,0 };
	cout << "########### Initialization ################" << endl << endl;
	MethodTool *methodTool = new MethodTool();
	methodTool->get16arr(mlen, flen);
	WaterEncode *encode = new WaterEncode();
	
	AudioFile<double> audioFile;
	audioFile.load(a);
	int sampleRate = audioFile.getSampleRate();
	int bitDepth = audioFile.getBitDepth();
	int numSamples = audioFile.getNumSamplesPerChannel();
	double lengthInSeconds = audioFile.getLengthInSeconds();
	int numChannels = audioFile.getNumChannels();
	bool isMono = audioFile.isMono();
	bool isStereo = audioFile.isStereo();
	//Print a summary to the console
	audioFile.printSummary();
	int numFrame = numSamples / mlen;
	cout << "Number of Frames: " << numFrame << endl;
	
	AudioFile<double>::AudioBuffer origin = audioFile.samples;
	//Create a new buffer for the output audio
	AudioFile<double>::AudioBuffer watered;

	//同步
	vector<int> tongbu; 
	double xlen = mlen; double ylen = 2 * mlen+2*flen;
	vector<double> r; r.resize(ylen + xlen - 1);
	Maxheap<entity> H(ylen + xlen - 1);
	for (int i = 0; i < 5; i++) {
		H.clear();

		Xcorr(r, methodTool->f_seq_mid[0], origin[0], i*(2*mlen+flen), xlen, ylen);
		double Max = 0; int MaxIndex = 0;
		for (int i = 0; i < xlen + ylen - 1; i++) {
			entity e(i, r[i]);
			H.push(e);
		}
		vector<int> Ind; Ind.resize(5);
		for (int i = 0; i < 5; i++) {
			Ind[i] = H.top().index;
			H.pop();
		}
		cout << "同步检测 Round： " << i << endl;
		if (i==0) tongbu = Ind;
		else {
			auto it = tongbu.begin();
			while ( it != tongbu.end()) { 
				bool remove = true;
				for (int k = 0; k < 5; k++) {
					if (abs(Ind[k] - *it)<10) {
						remove = false; break;
					}
				}
				if (remove) {
					it = tongbu.erase(it);
					if (it == tongbu.end())
						break;
				}
				else
				{
					it++;
				}
			}
		
		}
		if (i>3 && tongbu.size() <=1 ) {
			if (tongbu.size() == 1) {
				cout << "同步检测 已完成： 位置 " << tongbu[0] << endl;
				break;
			}
			else {
				cout << "同步信号检测失败，此音频被认为不含有水印信息 " << endl;
				return "";
			
			}
		}
	}
	int begin = tongbu[0] + 1; int frameNum = 1; vector<vector<int>> data;
	//提取信息
	vector<int> retrieved; bool startOfData = false;
	for (int i = begin; i < origin[0].size(); i += mlen) {
		if (frameNum%3!=0) {
			vector<double> debug;
			int codeword = 0; double max = 0;
			for (int k = 0; k < 16; k++) {
				double corr = Corrcoef(methodTool->m_seq_mid[k], origin[0], i);
				//cout << k << " : " << corr << endl;
				debug.push_back(corr);
				if (corr > max) {
					codeword = k; max = corr;
				}
			}
			cout << codeword << endl;
			retrieved.push_back(codeword);
			if (detectBark(retrieved, bark)) {
				//检测到bark码
				if (!startOfData)	startOfData = true;
				else {
					vector<int> tmp(retrieved);
					data.push_back(tmp);
				}
				retrieved.clear();
			}
		}

		frameNum++;
	}
	cout << "Finished Retriving All Data..." << endl << endl;
	vector<int> finalData;
	for (int i = 0; i < data[0].size(); i++) {
		int a[16] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0 };//记录频数
		for (int j = 0; j < data.size(); j++) {
			a[data[j][i]]++;
		}
		int max = a[0], flag = 0;
		for (int k = 0; k < 16; k++)
		{
			if (max < a[k])//找出频数最大的一位
			{
				max = a[k];
				flag = k;
			}
		}
		finalData.push_back(flag);
	}

	return dataToString(finalData);

}

string dataToString(vector<int> data) {
	string a = "";
	for (int i = 0; i < data.size()-7; i+=2) {
		//int tmp1 = (data[i]/8)+(data[i])
		int tmp = data[i]+(data[i+1]<<4);
		a = a + (char)tmp;
	}
	return a;
}

void genInfo(string a, vector<int>& info) {
	//注意：低位在前！
	for (int i = 0; i < a.length(); i++) {
		int num = (int)a[i];
		for (int j = 0; j < 8; j++) {
			info.push_back(num % 2);
			num >>= 1;
		}
	}
}

void saveFile(string retrieved, string foldpath) {

	ofstream   ofresult(foldpath);

	if (!ofresult.is_open()) {
		cout << "File open fail!" << endl;
		return;
	}

	ofresult << retrieved << endl;

	ofresult.close();
}

bool detectBark(vector<int> retrieved, int bark[]) {
	int len = 7;
	if (retrieved.size() < len)	return false;
	int same = 0;
	for (int i = 0; i < len;i++) {
		if (retrieved[retrieved.size() - len + i] == bark[i])	same++;
	}
	return same >= 5;
}


int getInfo(int index,vector<int> info) {
	return	info[index] + (info[index + 1] << 1) + (info[index + 2] << 2) + (info[index + 3] << 3);
}

string read(string name) {
	fstream  reader(name);
	string line = "";
	
		getline(reader, line);
		reader.close();

	return line;
}

void Xcorr(vector<double>& r, vector<double>& x, vector<double>& y,int startIndex,int xlen,int ylen)
{
	//x:原始信号，y:随机信号
	double sxy = 0;
	int delay, i, j;

	for (delay = - xlen + 1; delay < ylen; delay++)
	{
		sxy = 0;
		for (i = 0; i < xlen; i++)
		{
			j = i + delay;
			if ((j < 0) || (j >= ylen))//ignore
				continue;
			else
				sxy += (x[i] * y[j+startIndex]);
		}
		r[delay + xlen - 1] = sxy;
	}
}

double Corrcoef(vector<double> xnum, vector<double> ynum,int ybias) {
	int size = xnum.size();
	//double xnum[size] = { 56, 65, 47, 57, 62, 48, 68, 75, 79, 49 };
	//double ynum[size] = { 45, 49, 35, 44, 45, 40, 52, 57, 62, 39 };

	double totalxnum = 0;
	double meanxnum = 0;
	double totalynum = 0;
	double meanynum = 0;
	double varianceaddy = 0;
	double varianceaddx = 0;
	double variancex = 0;
	double variancey = 0;
	double stdeviationx = 0;
	double stdeviationy = 0;
	double precorrelation = 0;
	double correlation = 0;

	computemean(totalxnum, size, xnum, meanxnum, ynum, totalynum, meanynum,ybias);
	computevariance(varianceaddy, varianceaddx, xnum, meanxnum, meanynum, variancey, variancex, size, ynum,ybias);
	computestandarddeviation(stdeviationx, stdeviationy, variancey, variancex,ybias);
	computecorrelation(xnum, ynum, meanxnum, meanynum, stdeviationx, stdeviationy, size, precorrelation, correlation,ybias);

	return correlation;

}

void computemean(double totalxnum, int const size, vector<double> xnum, double& meanxnum, vector<double> ynum, double totalynum, double& meanynum,int ybias) {
	for (int y = 0; y < size; y++) {
		totalxnum += xnum[y];
		totalynum += ynum[y+ybias];
	}
	meanxnum = totalxnum / size;
	meanynum = totalynum / size;
}

void computevariance(double varianceaddy, double varianceaddx, vector<double> xnum, double meanxnum, double meanynum, double& variancey, double& variancex, int const size, vector<double> ynum,int ybias) {
	for (int i = 0; i < size; i++) {
		varianceaddx = varianceaddx + pow(xnum[i] - meanxnum, 2);
		varianceaddy = varianceaddy + pow(ynum[i+ybias] - meanynum, 2);
	}
	variancex = varianceaddx / 9;
	variancey = varianceaddy / 9;
}

void computestandarddeviation(double& stdeviationx, double& stdeviationy, double variancey, double variancex, int ybias) {
	stdeviationx = sqrt(variancex);
	stdeviationy = sqrt(variancey);
}



void computecorrelation(vector<double> xnum, vector<double> ynum, double meanxnum, double meanynum, double stdeviationx, double stdeviationy, int const size, double& precorrelation, double& correlation,int ybias) {
	for (int i = 0; i < size; i++) {
		precorrelation = precorrelation + (xnum[i] - meanxnum)*(ynum[i+ybias] - meanynum);
	}
	correlation = (precorrelation / 9) / (stdeviationx*stdeviationy);
}


//README.md
//
//Created by Qichao Ying 2019-3-30
//Final Version Completed on 2019-4-5
//###########Default Settings#############
//1.字符存储ascll码，7位
//
//################About DTCWT#############
//1.Yh{2},也就是BAND[1]是中频系数，Yh{1}是高频，Lo是低频
//
//################Utils###################
//1.copy(a.begin(), a.end(), b.begin() + 1); 
//把a中的从a.begin()（包括它）到a.end()（不包括它）的元素复制到b中，从b.begin()+1的位置（包括它）开始复制，覆盖掉原有元素
//2.vector<double> ivec(m_seq[i], m_seq[i]+mlen);
//这里用了begin和end函数;//使用初始化函数直接拷贝
//################Tests###################
//MethodTool *methodTool = new MethodTool();
//methodTool->get16arr(mlen, flen);
//WaterEncode *encode = new WaterEncode();
//DTCWT_BANDS* dtcwt = encode->dtwavexfm_1d(audioFile.samples, 0, 1024, 2);
//encode->dtwaveifm_1d(dtcwt, 1024, 2, out, 0);
//for (int i = 0; i < 1024; i++)
//	cout << out[0][i] << endl;
//cout << "Embedding Process Finished Successfully..." << endl;
//
//############About AudioFile##############
//Create an AudioFile object :
//#include "AudioFile.h"
//AudioFile<double> audioFile;
//Load an audio file :
//audioFile.load("/path/to/my/audiofile.wav");
//Get some information about the loaded audio :
//int sampleRate = audioFile.getSampleRate();
//int bitDepth = audioFile.getBitDepth();
//
//int numSamples = audioFile.getNumSamplesPerChannel();
//double lengthInSeconds = audioFile.getLengthInSeconds();
//
//int numChannels = audioFile.getNumChannels();
//bool isMono = audioFile.isMono();
//bool isStereo = audioFile.isStereo();
//
//// or, just use this quick shortcut to print a summary to the console
//audioFile.printSummary();
//Access the samples directly :
//int channel = 0;
//int numSamples = audioFile.getNumSamplesPerChannel();
//
//for (int i = 0; i < numSamples; i++)
//{
//	double currentSample = audioFile.samples[channel][i];
//}
//Replace the AudioFile audio buffer with another
//// 1. Create an AudioBuffer 
//// (BTW, AudioBuffer is just a vector of vectors)
//
//AudioFile<double>::AudioBuffer buffer;
//
//// 2. Set to (e.g.) two channels
//buffer.resize(2);
//
//// 3. Set number of samples per channel
//buffer[0].resize(100000);
//buffer[1].resize(100000);
//
//// 4. do something here to fill the buffer with samples, e.g.
//
//#include <math.h> // somewhere earler (for M_PI and sinf())
//
//// then...
//
//int numChannels = 2;
//int numSamplesPerChannel = 100000;
//float sampleRate = 44100.f;
//float frequency = 440.f;
//
//for (int i = 0; i < numSamplesPerChannel; i++)
//{
//	float sample = sinf(2. * M_PI * ((float)i / sampleRate) * frequency);
//
//	for (int channel = 0; channel < numChannels; channel++)
//		buffer[channel][i] = sample * 0.5;
//}
//
//// 5. Put into the AudioFile object
//bool ok = audioFile.setAudioBuffer(buffer);
//Resize the audio buffer
//// Set both the number of channels and number of samples per channel
//audioFile.setAudioBufferSize(numChannels, numSamples);
//
//// Set the number of samples per channel
//audioFile.setNumSamplesPerChannel(numSamples);
//
//// Set the number of channels
//audioFile.setNumChannels(int numChannels);
//Set bit depth and sample rate
//audioFile.setBitDepth(24);
//audioFile.setSampleRate(44100);
//Save the audio file to disk
//// Wave file (implicit)
//audioFile.save("path/to/desired/audioFile.wav");
//
//// Wave file (explicit)
//audioFile.save("path/to/desired/audioFile.wav", AudioFileFormat::Wave);
//
//// Aiff file
//audioFile.save("path/to/desired/audioFile.aif", AudioFileFormat::Aiff);
//A Note On Types
//AudioFile is a template class and so it can be instantiated using floating point precision :
//AudioFile<float> audioFile;
//... or double precision :
//AudioFile<double> audioFile;
//This simply reflects the data type you would like to use to store the underlying audio samples.
//You can still read or write 8, 16 or 24 - bit audio files, regardless of the type that you use
//(unless your system uses a precision for floats less than your desired bit depth).
//
