#pragma once

#include <string>
#include "common.h"
#include "mrcHeader.h"
#include <fstream>
#include <vector>
#include <typeinfo>

using namespace std;

class MRCStack
{
public:

	MRCStack(string aFileName, bool getRange, bool swapYZ, bool keepOpen);

	~MRCStack();

	void readProjections(float* aData, unsigned int numberOfProjections);
	void readProjections(float* aData, unsigned int numberOfProjections, unsigned int aIndex);
	void readAllProjections(float* aData);

	void writeOutHeader();
	float getStackMin();
	float getStackMax();
	float getStackMean();
	float getInputMeanValue();
	unsigned int getNumberOfProjections();
	size_t getProjectionSize();
	fourier::Vec3ui getResolution();

	MRCHeader getStackHeader();

	ifstream mStackFile;

private:

	MRCHeader mHeaderMRC;
	size_t mProjectionSize;
	string mStackName;

	template <typename voxelType>
	void readSlices(float* aData, size_t elementsToRead);


	void prepareFilePosition(size_t offset = 0);

	template <typename voxelType>
	void getDataRange(ifstream& infile);

	float mDataMin;				//min density value in the stack
	float mDataMax;				//max denisty value in the stack
	float mDataMean;
	float mInputMeanValue;
	unsigned int mNumberOfProjections;
	size_t mSizeOfVoxelType;

	bool mKeepOpen;

};
