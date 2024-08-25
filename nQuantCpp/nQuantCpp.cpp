// nQuantCpp.cpp
//

#include "stdafx.h"
#include <tchar.h>
#include <algorithm>
#include <chrono>
#include <io.h>
#define tcslen wcslen
#include <iostream>
#include <fcntl.h>
#include <filesystem>
#include <limits>
#include <sstream>
namespace fs = std::filesystem;

#include "nQuantCpp.h"

#include "PnnQuantizer.h"
#include "NeuQuantizer.h"
#include "WuQuantizer.h"
#include "PnnLABQuantizer.h"
#include "PnnLABGAQuantizer.h"
#include "EdgeAwareSQuantizer.h"
#include "SpatialQuantizer.h"
#include "DivQuantizer.h"
#include "Dl3Quantizer.h"
#include "MedianCut.h"
#include "Otsu.h"
#include "GifWriter.h"
#include <unordered_map>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#ifdef UNICODE
	auto& tcout = std::wcout;
#else
	auto& tcout = std::cout;
#endif

GdiplusStartupInput  m_gdiplusStartupInput;
ULONG_PTR m_gdiplusToken;

wstring algs[] = { L"PNN", L"PNNLAB", L"PNNLAB+", L"NEU", L"WU", L"EAS", L"SPA", L"DIV", L"DL3", L"MMC", L"OTSU" };
unordered_map<LPCTSTR, CLSID> extensionMap;

void PrintUsage()
{
	tcout << endl;
	tcout << "usage: nQuantCpp <input image path> [options]" << endl;
	tcout << endl;
	tcout << "Valid options:" << endl;
	tcout << "  /a : Algorithm used - Choose one of them, otherwise give you the defaults from [";
	int i = 0;	
	for(; i<sizeof(algs)/sizeof(string) - 1; ++i)
		tcout << algs[i] << ", ";
	tcout << algs[i] << "] ." << endl;
	tcout << "  /m : Max Colors (pixel-depth) - Maximum number of colors for the output format to support. The default is 256 (8-bit)." << endl;
	tcout << "  /d : Dithering or not? y or n." << endl;
	tcout << "  /f : Frame delay in milliseconds for PNNLAB+ only." << endl;
	tcout << "  /o : Output image file dir. The default is <source image path directory>" << endl;
}

bool isdigit(const wchar_t* chars, const bool positiveOnly = true) {
	const int string_len = wcslen(chars);
	for (int i = 0; i < string_len; ++i) {
		if (!isdigit(chars[i]))
			return positiveOnly ? false : (chars[0] == L'-');
	}
	return true;
}

bool isAlgo(const wstring& alg) {
	for (const auto& algo : algs) {
		if (alg == algo)
			return true;
	}		
	return false;
}

bool ProcessArgs(int argc, wstring& algo, UINT& nMaxColors, bool& dither, wstring& targetPath, wstring* argv, long& delay)
{
	for (int index = 1; index < argc; ++index) {
		auto currentArg = argv[index];
		transform(currentArg.begin(), currentArg.end(), currentArg.begin(), ::toupper);

		auto currentCmd = currentArg[0];
		if (currentArg.length() > 1 && 
			(currentCmd == L'–' || currentCmd == L'/')) {
			if (index >= argc - 1) {
				PrintUsage();
				return false;
			}

			if (currentArg[1] == L'A') {
				auto strAlgo = argv[index + 1];
				transform(strAlgo.begin(), strAlgo.end(), strAlgo.begin(), ::toupper);
				if (!isAlgo(strAlgo)) {
					PrintUsage();
					return false;
				}
				algo = strAlgo;
			}
			else if (currentArg[1] == L'M') {
				if (!isdigit(argv[index + 1].c_str())) {
					PrintUsage();
					return false;
				}
				nMaxColors = stoi(argv[index + 1].c_str());
				if (nMaxColors < 2)
					nMaxColors = 2;
				else if (nMaxColors > 65536)
					nMaxColors = 65536;
			}
			else if (currentArg[1] == L'D') {
				auto strDither = argv[index + 1];
				transform(strDither.begin(), strDither.end(), strDither.begin(), ::toupper);
				if (!(strDither == L"Y" || strDither == L"N")) {
					PrintUsage();
					return false;
				}
				dither = strDither == L"Y";
			}
			else if (currentArg[1] == L'F') {
				if (!isdigit(argv[index + 1].c_str(), false)) {
					PrintUsage();
					return false;
				}
				delay = stol(argv[index + 1].c_str());
			}
			else if (currentArg[1] == L'O') {
				auto szPath = argv[index + 1].c_str();
				wstring tmpPath(szPath, szPath + wcslen(szPath));
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

bool OutputImage(const fs::path& sourcePath, const wstring& algorithm, const UINT& nMaxColors, wstring& targetDir, Bitmap* pDest, LPCTSTR defaultExtension = L".png")
{
	auto fileName = sourcePath.filename().wstring();
	fileName = fileName.substr(0, fileName.find_last_of(L'.'));

	targetDir = fileExists(targetDir) ? fs::canonical(fs::path(targetDir)) : fs::current_path();
	auto destPath = targetDir + L"/" + fileName + L"-";
	wstring algo(algorithm.begin(), algorithm.end());
	destPath += algo + L"quant";
	
	// image/bmp  : {557cf400-1a04-11d3-9a73-0000f81ef32e}
	const CLSID bmpEncoderClsId = { 0x557cf400, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
	extensionMap.emplace(L".bmp", bmpEncoderClsId);

	// image/gif  : {557cf402-1a04-11d3-9a73-0000f81ef32e}
	const CLSID gifEncoderClsId = { 0x557cf402, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
	extensionMap.emplace(L".gif", gifEncoderClsId);

	// image/png  : {557cf406-1a04-11d3-9a73-0000f81ef32e}
	const CLSID pngEncoderClsId = { 0x557cf406, 0x1a04, 0x11d3,{ 0x9a,0x73,0x00,0x00,0xf8,0x1e,0xf3,0x2e } };
	extensionMap.emplace(L".png", pngEncoderClsId);

	auto targetExtension = (pDest->GetPixelFormat() < PixelFormat16bppARGB1555 && nMaxColors > 256) ? L".bmp" : defaultExtension;
	destPath += std::to_wstring(nMaxColors) + targetExtension;
	auto status = pDest->Save(destPath.c_str(), &extensionMap[targetExtension]);
	if (status == Status::Ok)
		tcout << L"Converted image: " << destPath << endl;
	else
		tcout << L"Failed to save image in '" << destPath << L"' file" << endl;

	return status == Status::Ok;
}

bool QuantizeImage(const wstring& algorithm, const wstring& sourceFile, wstring& targetDir, shared_ptr<Bitmap> pSource, UINT nMaxColors, bool dither)
{
	// Create 8 bpp indexed bitmap of the same size
	auto pDest = make_shared<Bitmap>(pSource->GetWidth(), pSource->GetHeight(), (nMaxColors > 256) ? PixelFormat16bppARGB1555 : (nMaxColors > 16) ? PixelFormat8bppIndexed : (nMaxColors > 2) ? PixelFormat4bppIndexed : PixelFormat1bppIndexed);

	bool bSucceeded = false;
	if (algorithm == L"PNN") {
		PnnQuant::PnnQuantizer pnnQuantizer;
		bSucceeded = pnnQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"PNNLAB") {
		PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
		bSucceeded = pnnLABQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"PNNLAB+") {
		PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
		vector<shared_ptr<Bitmap> > sources(1, pSource);
		PnnLABQuant::PnnLABGAQuantizer pnnLABGAQuantizer(pnnLABQuantizer, sources, nMaxColors);
		nQuantGA::APNsgaIII alg(pnnLABGAQuantizer);
		alg.run(9999, -numeric_limits<double>::epsilon());
		auto pGAq = alg.getResult();
		tcout << L"\n" << pGAq->getResult().c_str() << endl;
		vector<shared_ptr<Bitmap> > dests;
		dests.emplace_back(pDest);
		bSucceeded = pGAq->QuantizeImage(dests, dither);
	}
	else if (algorithm == L"NEU") {
		NeuralNet::NeuQuantizer neuQuantizer;
		bSucceeded = neuQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"WU") {
		nQuant::WuQuantizer wuQuantizer;
		bSucceeded = wuQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"EAS") {
		EdgeAwareSQuant::EdgeAwareSQuantizer easQuantizer;
		bSucceeded = easQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"SPA") {
		SpatialQuant::SpatialQuantizer spaQuantizer;
		bSucceeded = spaQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"DIV") {
		DivQuant::DivQuantizer divQuantizer;
		bSucceeded = divQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"DL3") {
		Dl3Quant::Dl3Quantizer dl3Quantizer;
		bSucceeded = dl3Quantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"MMC") {
		MedianCutQuant::MedianCut mmcQuantizer;
		bSucceeded = mmcQuantizer.QuantizeImage(pSource.get(), pDest.get(), nMaxColors, dither);
	}
	else if (algorithm == L"OTSU") {
		nMaxColors = 2;
		OtsuThreshold::Otsu otsu;
		bSucceeded = otsu.ConvertGrayScaleToBinary(pSource.get(), pDest.get());
	}

	if (!bSucceeded)
		return bSucceeded;

	auto sourcePath = fs::canonical(fs::path(sourceFile));
	return OutputImage(sourcePath, algorithm, nMaxColors, targetDir, pDest.get());
}

static void OutputImages(const fs::path& sourceDir, wstring& targetDir, const UINT& nMaxColors, const bool dither, const wstring& algo, const long& delay)
{
	auto start = chrono::steady_clock::now();

	vector<fs::path> sourcePaths;
	vector<shared_ptr<Bitmap> > pSources, pDests;
	for (const auto& entry : fs::recursive_directory_iterator(sourceDir)) {
		if (entry.is_regular_file() && !entry.is_symlink()) {
			auto pSource = shared_ptr<Bitmap>(Bitmap::FromFile(entry.path().wstring().c_str()));
			auto status = pSource->GetLastStatus();
			if (status != Ok)
				continue;
			sourcePaths.emplace_back(entry.path());
			pSources.emplace_back(pSource);
			pDests.emplace_back(make_shared<Bitmap>(pSource->GetWidth(), pSource->GetHeight(), (nMaxColors > 256) ? PixelFormat16bppARGB1555
			: (nMaxColors > 16) ? PixelFormat8bppIndexed : (nMaxColors > 2) ? PixelFormat4bppIndexed : PixelFormat1bppIndexed));
		}
	}

	if (algo == _T("PNNLAB+")) {
		PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
		PnnLABQuant::PnnLABGAQuantizer pnnLABGAQuantizer(pnnLABQuantizer, pSources, nMaxColors);
		nQuantGA::APNsgaIII alg(pnnLABGAQuantizer);
		alg.run(9999, -numeric_limits<double>::epsilon());
		auto pGAq = alg.getResult();
		tcout << L"\n" << pGAq->getResult().c_str() << endl;
		if (pGAq->QuantizeImage(pDests, dither)) {
			if (nMaxColors > 256 || delay < 0) {
				int i = 0;
				for (auto& sourcePath : sourcePaths)
					OutputImage(sourcePath, L"PNNLAB+", nMaxColors, targetDir, pDests[i++].get(), nMaxColors > 256 || delay > -2 ? L".png" : L".gif");
			}
			else {
				auto fileName = sourcePaths[0].filename().wstring();
				fileName = fileName.substr(0, fileName.find_last_of(L'.'));

				targetDir = fileExists(targetDir) ? fs::canonical(fs::path(targetDir)) : fs::current_path();
				auto destPath = targetDir + L"/" + fileName + L"-";
				destPath += L"PNNLAB+quant.gif";
				GifEncode::GifWriter gifWriter(destPath, pGAq->hasAlpha(), abs(delay));
				auto status = gifWriter.AddImages(pDests);
				if (status == Status::Ok)
					tcout << L"Converted image: " << destPath << endl;
				else
					tcout << L"Failed to save image in '" << destPath << L"' file" << endl;
			}
		}
	}
	else {
		if (nMaxColors > 256 || delay < 0) {
			int i = 0;
			for (auto& sourcePath : sourcePaths)
				OutputImage(sourcePath, L"PNN", nMaxColors, targetDir, pDests[i++].get(), nMaxColors > 256 || delay > -2 ? L".png" : L".gif");
		}
		else {
			ostringstream ss;
			auto fileName = sourcePaths[0].filename().wstring();
			fileName = fileName.substr(0, fileName.find_last_of(L'.'));

			targetDir = fileExists(targetDir) ? fs::canonical(fs::path(targetDir)) : fs::current_path();
			auto destPath = targetDir + L"/" + fileName + L"-";
			if (algo == _T("PNNLAB")) {
				destPath += L"PNNLABquant.gif";

				UINT maxColors = nMaxColors;
				for (int i = 0; i < pSources.size(); ++i) {
					ss << "\r" << i << " of " << pSources.size() << " completed." << showpoint;
					tcout << ss.str().c_str();

					PnnLABQuant::PnnLABQuantizer pnnLABQuantizer;
					pnnLABQuantizer.QuantizeImage(pSources[i].get(), pDests[i].get(), maxColors, dither);					
				}				
			}
			else {
				destPath += L"PNNquant.gif";

				UINT maxColors = nMaxColors;
				for (int i = 0; i < pSources.size(); ++i) {
					ss << "\r" << i << " of " << pSources.size() << " completed." << showpoint;
					tcout << ss.str().c_str();

					PnnQuant::PnnQuantizer pnnQuantizer;
					pnnQuantizer.QuantizeImage(pSources[i].get(), pDests[i].get(), maxColors, dither);
				}
			}
			tcout << L"\rWell done!!!                             " << endl;

			GifEncode::GifWriter gifWriter(destPath, false, abs(delay));
			auto status = gifWriter.AddImages(pDests);
			if (status == Status::Ok)
				tcout << L"Converted image: " << destPath << endl;
			else
				tcout << L"Failed to save image in '" << destPath << L"' file" << endl;
		}
	}

	auto dur = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() / 1000000.0;
	tcout << "Completed in " << dur << " secs." << endl;
}

int _tmain(int argc, _TCHAR** argv)
{
	_setmode(_fileno(stdout), _O_U16TEXT);

#ifndef _DEBUG
	if (argc <= 1) {
		PrintUsage();
		return 0;
	}
#endif
	
	auto szDir = fs::current_path().wstring();
	
	bool dither = true;
	UINT nMaxColors = 256;
	long delay = -1;
	wstring algo = L"";
	wstring targetDir = L"";

	vector<wstring> argList(argc);
	for (int i = 1; i < argc; ++i)
		argList[i] = wstring(argv[i], argv[i] + tcslen(argv[i]));

#ifdef _DEBUG
	wstring sourceFile = szDir + L"/../ImgV64.gif";
	nMaxColors = 1024;
#else
	if (!ProcessArgs(argc, algo, nMaxColors, dither, targetDir, argList.data(), delay))
		return 0;

	wstring sourceFile(argv[1], argv[1] + tcslen(argv[1]));
	if (!fileExists(sourceFile) && sourceFile.find_first_of(L"\\/") != wstring::npos)
		sourceFile = szDir + L"/" + sourceFile;
#endif
	
	if(!fileExists(sourceFile)) {
		tcout << "The source file you specified does not exist." << endl;
		return 0;
	}		

	if(GdiplusStartup(&m_gdiplusToken, &m_gdiplusStartupInput, NULL) == Ok) {
		if(fs::is_directory(fs::status(sourceFile.c_str())) ) {
			if (!targetDir.empty() && !fileExists(targetDir))
				fs::create_directories(targetDir);
			OutputImages(sourceFile, targetDir, nMaxColors, dither, algo, delay);
			GdiplusShutdown(m_gdiplusToken);
			return 0;
		}
		
		auto pSource = shared_ptr<Bitmap>(Bitmap::FromFile(sourceFile.c_str()));
		auto status = pSource->GetLastStatus();
		if (status == Ok) {
			auto start = chrono::steady_clock::now();
			if (!fileExists(targetDir))
				targetDir = fs::path(sourceFile).parent_path().wstring();

	#ifdef _DEBUG
			auto sourcePath = fs::canonical(fs::path(sourceFile));
			sourceFile = sourcePath;
	#endif

			sourceFile = (sourceFile[sourceFile.length() - 1] != L'/' && sourceFile[sourceFile.length() - 1] != L'\\') ? sourceFile : sourceFile.substr(0, sourceFile.find_last_of(L"\\/"));
			if (algo == L"") {
				//QuantizeImage(L"MMC", sourceFile, targetDir, pSource, nMaxColors, dither);
				QuantizeImage(L"DIV", sourceFile, targetDir, pSource, nMaxColors, dither);
				if (nMaxColors > 32) {
					QuantizeImage(L"PNN", sourceFile, targetDir, pSource, nMaxColors, dither);
					QuantizeImage(L"WU", sourceFile, targetDir, pSource, nMaxColors, dither);
					QuantizeImage(L"NEU", sourceFile, targetDir, pSource, nMaxColors, dither);
				}
				else {
					QuantizeImage(L"PNNLAB", sourceFile, targetDir, pSource, nMaxColors, dither);
					QuantizeImage(L"EAS", sourceFile, targetDir, pSource, nMaxColors, dither);
					QuantizeImage(L"SPA", sourceFile, targetDir, pSource, nMaxColors, dither);
				}
			}
			else
				QuantizeImage(algo, sourceFile, targetDir, pSource, nMaxColors, dither);

			auto dur = chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count() / 1000000.0;
			tcout << "Completed in " << dur << " secs." << endl;
		}
		else
			tcout << "Failed to read image in '" << sourceFile.c_str() << "' file";
	}
	GdiplusShutdown(m_gdiplusToken);
	return 0;
}
