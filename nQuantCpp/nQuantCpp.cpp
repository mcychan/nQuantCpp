// nQuantCpp.cpp
//

#include "stdafx.h"
#include <algorithm>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;

#include "nQuantCpp.h"

#include "PnnQuantizer.h"
#include "NeuQuantizer.h"
#include "WuQuantizer.h"
#include "PnnLABQuantizer.h"
#include "EdgeAwareSQuantizer.h"
#include "SpatialQuantizer.h"
#include "DivQuantizer.h"
#include "MoDEQuantizer.h"
#include "MedianCut.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#ifdef UNICODE
wostream& tcout = wcout;
#else
ostream& tcout = cout;
#endif

GdiplusStartupInput  m_gdiplusStartupInput;
ULONG_PTR m_gdiplusToken;

string algs[] = { "PNN", "PNNLAB", "NEU", "WU", "EAS", "SPA", "DIV", "MODE", "MMC" };

void PrintUsage()
{
    cout << endl;
    cout << "usage: nQuantCpp <input image path> [options]" << endl;
    cout << endl;
    cout << "Valid options:" << endl;
	cout << "  /a : Algorithm used - Choose one of them, otherwise give you the defaults from [";
	int i = 0;	
	for(; i<sizeof(algs)/sizeof(string) - 1; ++i)
		cout << algs[i] << ", ";
	cout << algs[i] << "] ." << endl;
    cout << "  /m : Max Colors (pixel-depth) - Maximum number of colors for the output format to support. The default is 256 (8-bit)." << endl;
	cout << "  /d : Dithering or not? y or n." << endl;
    cout << "  /o : Output image file dir. The default is <source image path directory>" << endl;
}

bool isdigit(const char* string) {
	const int string_len = strlen(string);
	for (int i = 0; i < string_len; ++i) {
		if (!isdigit(string[i]))
			return false;
	}
	return true;
}

bool isAlgo(const string& alg) {
	for (const auto& algo : algs) {
		if (alg == algo)
			return true;
	}		
	return false;
}

bool ProcessArgs(int argc, string& algo, UINT& nMaxColors, bool& dither, wstring& targetPath, char** argv)
{
	for (int index = 1; index < argc; ++index) {
		string currentArg(argv[index]);
		transform(currentArg.begin(), currentArg.end(), currentArg.begin(), ::toupper);

		auto currentCmd = currentArg[0];
		if (currentArg.length() > 1 && 
			(currentCmd == '-' || currentCmd == '–' || currentCmd == '/')) {
			if (currentArg[1] == 'A') {
				string strAlgo = argv[index + 1];
				transform(strAlgo.begin(), strAlgo.end(), strAlgo.begin(), ::toupper);
				if (index >= argc - 1 || !isAlgo(strAlgo)) {
					PrintUsage();
					return false;
				}
				algo = strAlgo;
			}
			else if (currentArg[1] == 'M') {
				if (index >= argc - 1 || !isdigit(argv[index + 1])) {
					PrintUsage();
					return false;
				}
				nMaxColors = atoi(argv[index + 1]);
				if (nMaxColors < 2)
					nMaxColors = 2;
				else if (nMaxColors > 65536)
					nMaxColors = 65536;
			}
			else if (currentArg[1] == 'D') {
				string strDither = argv[index + 1];
				transform(strDither.begin(), strDither.end(), strDither.begin(), ::toupper);
				if (index >= argc - 1 || !(strDither == "Y" || strDither == "N")) {
					PrintUsage();
					return false;
				}
				dither = strDither == "Y";
			}
			else if (currentArg[1] == 'O') {
				if (index >= argc - 1) {
					PrintUsage();
					return false;
				}
				wstring tmpPath(argv[index + 1], argv[index + 1] + strlen(argv[index + 1]));
				targetPath = tmpPath;
			}
			else {
				PrintUsage();
				return false;
			}
		}
	}
	return true;
}

inline bool fileExists(const wstring& path)
{
	return fs::exists(fs::path(path));
}

bool QuantizeImage(const string& algorithm, const wstring& sourceFile, wstring& targetDir, Bitmap* pSource, UINT nMaxColors, bool dither)
{
	// Create 8 bpp indexed bitmap of the same size
	auto pDest = make_unique<Bitmap>(pSource->GetWidth(), pSource->GetHeight(), (nMaxColors > 256) ? PixelFormat16bppARGB1555 : (nMaxColors > 16) ? PixelFormat8bppIndexed : (nMaxColors > 2) ? PixelFormat4bppIndexed : PixelFormat1bppIndexed);

	bool bSucceeded = false;
	if (algorithm == "PNN") {
		PnnQuant::PnnQuantizer pnnQuantizer;
		bSucceeded = pnnQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "PNNLAB") {
		PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
		bSucceeded = pnnLABQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither ? 1 : -1);
	}
	else if (algorithm == "NEU") {
		NeuralNet::NeuQuantizer neuQuantizer;
		bSucceeded = neuQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "WU") {
		nQuant::WuQuantizer wuQuantizer;
		bSucceeded = wuQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "EAS") {
		EdgeAwareSQuant::EdgeAwareSQuantizer easQuantizer;
		bSucceeded = easQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "SPA") {
		SpatialQuant::SpatialQuantizer spaQuantizer;
		bSucceeded = spaQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "DIV") {
		DivQuant::DivQuantizer divQuantizer;
		bSucceeded = divQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "MODE") {
		MoDEQuant::MoDEQuantizer moDEQuantizer;
		bSucceeded = moDEQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == "MMC") {
		MedianCutQuant::MedianCut mmcQuantizer;
		bSucceeded = mmcQuantizer.QuantizeImage(pSource, pDest.get(), nMaxColors, dither);
	}

	if (!bSucceeded)
		return bSucceeded;

	auto sourcePath = fs::canonical(fs::path(sourceFile));
	auto fileName = sourcePath.filename().wstring();
	fileName = fileName.substr(0, fileName.find_last_of(L'.'));

	targetDir = fileExists(targetDir) ? fs::canonical(fs::path(targetDir)) : fs::current_path();	
	auto destPath = targetDir + L"/" + fileName + L"-";
	wstring algo(algorithm.begin(), algorithm.end());
	destPath += algo + L"quant";
	destPath += std::to_wstring(nMaxColors) + L".png";
	
	// image/bmp  : {557cf400-1a04-11d3-9a73-0000f81ef32e}
	const CLSID bmpEncoderClsId = { 0x557cf400, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };

	// image/gif  : {557cf402-1a04-11d3-9a73-0000f81ef32e}
	const CLSID gifEncoderClsId = { 0x557cf402, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };

	// image/png  : {557cf406-1a04-11d3-9a73-0000f81ef32e}
	const CLSID pngEncoderClsId = { 0x557cf406, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
	Status status = pDest->Save(destPath.c_str(), &pngEncoderClsId);
	if (status == Status::Ok)
		tcout << "Converted image: " << destPath << endl;
	else
		tcout << "Failed to save image in '" << destPath << "' file" << endl;

	return status == Status::Ok;
}

int main(int argc, char** argv)
{
	if (argc <= 1) {
#ifndef _DEBUG
		PrintUsage();
		return 0;
#endif
	}
	
	auto szDir = fs::current_path().wstring();
	
	bool dither = true;
	UINT nMaxColors = 256;
	string algo = "";
	wstring targetDir = L"";
#ifdef _DEBUG
	wstring sourcePath = szDir + L"/../ImgV64.gif";
	nMaxColors = 1024;
#else
	if (!ProcessArgs(argc, algo, nMaxColors, dither, targetDir, argv))
		return 0;

	wstring sourcePath(argv[1], argv[1] + strlen(argv[1]));
	if (!fileExists(sourcePath) && sourcePath.find_first_of(L"\\/") != wstring::npos)
		sourcePath = szDir + L"/" + sourcePath;
#endif	
	
	if(!fileExists(sourcePath)) {
		cout << "The source file you specified does not exist." << endl;
		return 0;
	}		

	if(GdiplusStartup(&m_gdiplusToken, &m_gdiplusStartupInput, NULL) == Ok) {
		auto pSource = unique_ptr<Bitmap>(Bitmap::FromFile(sourcePath.c_str()));
		Status status = pSource->GetLastStatus();
		if (status == Ok) {
			if (!fileExists(targetDir))
				targetDir = fs::path(sourcePath).parent_path().wstring();

			auto sourceFile = (sourcePath[sourcePath.length() - 1] != L'/' && sourcePath[sourcePath.length() - 1] != L'\\') ? sourcePath : sourcePath.substr(0, sourcePath.find_last_of(L"\\/"));
			if (algo == "") {
				//QuantizeImage(_T("MMC"), sourceFile, targetDir, pSource.get(), nMaxColors, dither);
				QuantizeImage("DIV", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
				if (nMaxColors > 32) {
					QuantizeImage("PNN", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
					QuantizeImage("WU", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
					//QuantizeImage("MODE", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
					QuantizeImage("NEU", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
				}
				else {
					QuantizeImage("PNNLAB", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
					QuantizeImage("EAS", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
					QuantizeImage("SPA", sourceFile, targetDir, pSource.get(), nMaxColors, dither);
				}
			}
			else
				QuantizeImage(algo, sourceFile, targetDir, pSource.get(), nMaxColors, dither);
		}
		else
			tcout << "Failed to read image in '" << sourcePath.c_str() << "' file";
	}
	GdiplusShutdown(m_gdiplusToken);
    return 0;
}
