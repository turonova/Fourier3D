#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "parameterSetup.h"

ParameterSetup::ParameterSetup()
{

}

void ParameterSetup::parseIntegers(std::string& aParamValue,std::vector<int>& dim, const char separator)
{
	unsigned int start=0;

	for(unsigned int i = 0; i < aParamValue.length(); i++)
	    {
	        if(aParamValue[i]==separator)
	        {	            
	            string tmpValue = aParamValue.substr(start, i-start);
				
				//do not want the space neither in name or value
				while(isspace(aParamValue[i+1]))
				{
					i++;
				}
				start=i+1;				
				dim.push_back(atoi(tmpValue.c_str()));
	        }
	    }
	 string tmpValue = aParamValue.substr(start, aParamValue.length());
	 dim.push_back(atoi(tmpValue.c_str()));
}


void ParameterSetup::parseCommandLine(std::vector<string> argList)
{
	for(std::vector<string>::iterator it=argList.begin(); it!=argList.end();it++)
	{
		string paramName=*it;
		string paramValue;
		std::vector<string>::iterator nextArgument = it+1;
		if(nextArgument!=argList.end())
		{
			paramValue=*nextArgument;
			it++;
		}
		else
		{
			cout << "Missing value of the following parameter: " << *it << endl;
			return;
		}


		paramName.erase(0,1); //erase "-" at the beginning of command line arguments
		storeValues(paramName,paramValue,',');
	}
}

void ParameterSetup::storeValues(string paramName, string paramValue, char separator)
{
		if(paramName == "OutputFile")
		{
			outputFilename = paramValue;
		}
		else if(paramName == "InputFile")
		{
			inputFilename = paramValue;
		}
		else if(paramName == "ContinueFFT")
		{
			continueFFT = atoi(paramValue.c_str());;
		}
		else if(paramName == "DeleteFiles")
		{
			deleteFiles = atoi(paramValue.c_str());;
		}
		else if(paramName == "MemoryLimit")
		{
			memoryLimit = atoi(paramValue.c_str());;
		}
		else if(paramName == "BinFactor")
		{	
			std::vector<int> tempValues;
			
			if(paramValue.find(separator)!=string::npos)
			{
				parseIntegers(paramValue, tempValues, separator);
				binFactor.x = tempValues[0];
				binFactor.y = tempValues[1];
				binFactor.z = tempValues[2];
			}
			else
			{
				binFactor.x = atoi(paramValue.c_str());
				binFactor.y = binFactor.x;
				binFactor.z = binFactor.x;
			}

			binFactorInitialized=true;

			if(newDimensionsInitialized)
			{
				cout << "BinFactor was specified after NewDimensions and will be used for cropping." << endl;
				newDimensionsInitialized = false;

				newDimensions.x = 0;
				newDimensions.y = 0;
				newDimensions.z = 0;
			}
		}
		else if(paramName == "NewDimensions")
		{
			std::vector<int> tempValues;

			if(paramValue.find(separator)!=string::npos)
			{
				parseIntegers(paramValue, tempValues, separator);
				newDimensions.x = tempValues[0];
				newDimensions.y = tempValues[1];
				newDimensions.z = tempValues[2];

				newDimensionsInitialized = true;

				if(binFactorInitialized)
				{
					cout << "NewDimensions were specified after BinFactor and will be used for cropping." << endl;
					binFactorInitialized = false;

					binFactor.x = 0;
					binFactor.y = 0;
					binFactor.z = 0;
				}
			}
			else
			{
				cout << "You have to specify all three dimensions and separate them by comma!!!" << endl;
			}
		}
		else
		{
			cout << "Ignoring following parameter: " << paramName << endl;
		}
}



ParameterSetup::ParameterSetup(std::vector<string> argList)
{
	binFactorInitialized = false;
	newDimensionsInitialized = false;

	continueFFT = 0;
	memoryLimit = 0;

	binFactor.x = 0;
	binFactor.y = 0;
	binFactor.z = 0;

	newDimensions.x = 0;
	newDimensions.y = 0;
	newDimensions.z = 0;

	deleteFiles = true;

	parseCommandLine(argList);

	if(!binFactorInitialized && !newDimensionsInitialized)
	{
		cout << "You have to specify either BinFactor or NewDimensions!!!" << endl;
		exit(EXIT_FAILURE);
	}

	if(binFactorInitialized)
	{
		if(binFactor.x==0 || binFactor.y==0 || binFactor.z==0)
		{
			cout << "BinFactor cannot be zero!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}
	else // if(newDimensionsInitialized)
	{
		if(newDimensions.x==0 || newDimensions.y==0 || newDimensions.z==0)
		{
			cout << "Neither of the dimensions can be 0!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

ParameterSetup::~ParameterSetup()
{}


string ParameterSetup::InputFilename()
{
	return inputFilename;
}

string ParameterSetup::OutputFilename()
{
	return outputFilename;
}

fourier::Vec3i ParameterSetup::BinFactor()
{
	return binFactor;
}

fourier::Vec3i ParameterSetup::NewDimensions()
{
	return newDimensions;
}

unsigned int ParameterSetup::ContinueFFT()
{
	return continueFFT;
}

unsigned int ParameterSetup::MemoryLimit()
{
	return memoryLimit;
}

bool ParameterSetup::DeleteFiles()
{
	return deleteFiles;
}

