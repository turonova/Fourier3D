#pragma once
#include "parameterSetup.h"
#include "mrcStack.h"


class Fourier3D
{
public:
	Fourier3D(ParameterSetup& aParams);
	~Fourier3D();

	void run();

private:

	void generate2DFFTVolume();
	void generate3DFFT();
	void computeCroppedIndices2D(std::vector<size_t>& indices, size_t oldDimX, size_t oldDimY, size_t newDimX, size_t newDimY);
	void computeCroppedIndices1D(std::vector<size_t>& indices, size_t oldDim, size_t newDim);
	void initParameters();
	void inverse1DFFT();
	void transposeInPlace(std::vector<float>& data, size_t sizeX, size_t sizeY);
	void transposeInPlaceDiffSizes(std::vector<float>& data, size_t originalSizeX, size_t sizeX, size_t sizeY);
	void prepareFiles();
	void deleteFiles();
	void openFile(ifstream& file, string filename);
	void computeBatchSize();

	string generateTempName(std::string base, std::string extension, int index = -1);

	MRCStack* inputVolume;
	ParameterSetup params;

	fourier::Vec3t dimensions;
	fourier::Vec3t newDimensions;
	//Vec3t newDimHalf;
	
	std::vector<size_t> indicesToUse2D;
	std::vector<size_t> indicesToUse1Dz;
	std::vector<size_t> indicesToUse1Dx;
	std::vector<size_t> indicesToUse1Dy;

	unsigned int batchXY;
	unsigned int batchZY;
	unsigned int batchInvXY;

	unsigned int reportXY;
	unsigned int reportZY;
	unsigned int reportInvXY;

	unsigned int zyNumberOfFiles;
	unsigned int xyNumberOfFiles;
};
