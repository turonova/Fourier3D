#pragma once
#include <vector>
#include "common.h"
#include <string>
using namespace std;

class ParameterSetup
{

public:
	ParameterSetup(std::vector<string> argList);
	ParameterSetup();

	~ParameterSetup();

	string InputFilename();
	string OutputFilename();
	
	fourier::Vec3i BinFactor();
	fourier::Vec3i NewDimensions();

	unsigned int ContinueFFT();
	unsigned int MemoryLimit();

	bool DeleteFiles();
private:
	
	void parseCommandLine(std::vector<string> argList);
	void storeValues(string paramName, string paramValue, char separator);
	void parseIntegers(std::string& aParamValue,std::vector<int>& dim, const char separator);

	string inputFilename;
	string outputFilename;
	
	fourier::Vec3i binFactor;
	fourier::Vec3i newDimensions;

	bool binFactorInitialized;
	bool newDimensionsInitialized;

	unsigned int continueFFT;
	unsigned int memoryLimit;

	bool deleteFiles;

};
