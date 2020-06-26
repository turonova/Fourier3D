//============================================================================
// Name        : WBP.cpp
// Author      : 
// Version     :
//============================================================================

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include "parameterSetup.h"
#include "exception.h"
#include "volumeIO.h"

#include "fourier3D.h"

using namespace std;

string dateAndTimeStamp(time_t startTime)
{
	time_t endTime;
	//struct tm * timeinfo;
	//char buffer[80];

	time(&endTime);

	double seconds_total = difftime(endTime,startTime);
	unsigned int hours = ((unsigned int)seconds_total)/3600;
	unsigned int minutes = (((int)seconds_total)%3600)/60;
	unsigned int seconds = (((int)seconds_total)%3600)%60;
	//timeinfo = localtime(&rawtime);
	stringstream timeDifference;


	if(hours>0)
		timeDifference << hours << "h ";

	if(minutes > 0 || hours > 0)
		timeDifference << minutes << "m ";

	timeDifference << seconds << "s";

	//strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
	//std::string stamp(buffer);

	return timeDifference.str();
}

int main (int argc, const char* argv[]) {

	time_t startTime;
	time(&startTime);

	vector<std::string> argList;

	cout << "Reading parameters... " <<  endl;

	for(int i=1;i<argc;i++)
	{
		argList.push_back(argv[i]);
	}

    ParameterSetup params(argList);

    cout << "...done" <<  endl << endl;
	cout << "Starting 3D Fourier crop." << endl << endl;
	Fourier3D* newRec = new Fourier3D(params);
	newRec->run();
	cout << "3D Fourier crop finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
	delete newRec;


    return 0;
}
