#include "fourier3D.h"
#include "common.h"
#include "volumeIO.h"
#include <math.h>
#include <algorithm>
#include <float.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include "exception.h"
#include <fstream>
#include "fftRoutines.h"

using namespace fourier;
using namespace std;


Fourier3D::Fourier3D(ParameterSetup& aParams)
{
	params=aParams;
	inputVolume = new MRCStack(params.InputFilename(),false,false,1);
}

Fourier3D::~Fourier3D()
{
	delete inputVolume;
}

void Fourier3D::run()
{
	initParameters();
	prepareFiles();
	computeBatchSize();

	if(params.ContinueFFT()==0)
	{
		generate2DFFTVolume();
	}
	else
	{
		cout << "Skipping the FFT in X and Y directions! \n" << std:: endl;
	}

	if(params.ContinueFFT()<=1)
	{
		generate3DFFT();
	}

	if(params.ContinueFFT()==2)
	{
		cout << "Skipping the FFT in Z direction and inverse FFT in Z and Y directions! \n" << std:: endl;
	}

	inverse1DFFT();

	deleteFiles();

}

void Fourier3D::initParameters()
{
	dimensions = inputVolume->getResolution();

	if(params.NewDimensions().x!=0)
	{
		newDimensions.x=params.NewDimensions().x;
		newDimensions.y=params.NewDimensions().y;
		newDimensions.z=params.NewDimensions().z;
	}
	else if(params.BinFactor().x!=0)
	{
		newDimensions.x=dimensions.x/params.BinFactor().x;
		newDimensions.y=dimensions.y/params.BinFactor().y;
		newDimensions.z=dimensions.z/params.BinFactor().z;
	}

	if(newDimensions.x>dimensions.x || newDimensions.y>dimensions.y || newDimensions.z>dimensions.z)
	{
		cout << "At least one of the new dimensions is larger than the original one!" << std::endl;
		exit(EXIT_FAILURE);
	}

	computeCroppedIndices1D(indicesToUse1Dx,dimensions.x,newDimensions.x);
	computeCroppedIndices1D(indicesToUse1Dy,dimensions.y,newDimensions.y);
	computeCroppedIndices1D(indicesToUse1Dz,dimensions.z,newDimensions.z);
}


void Fourier3D::prepareFiles()
{
	VolumeIO::writeHeader(params.InputFilename(),params.OutputFilename(),newDimensions,2,fourier::VolumeRotation::ALONG_XY);

	/*if(params.ContinueFFT()==2)
	{
		//all files are there already, no need to do anything
		return;
	}
	else if(params.ContinueFFT()==0)
	{
		VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D_real.rec"),Vec3t(newDimensions.x,newDimensions.y,dimensions.z),2,fourier::VolumeRotation::ALONG_XY,true);
		VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D_complex.rec"),Vec3t(newDimensions.x,newDimensions.y,dimensions.z),2,fourier::VolumeRotation::ALONG_XY,true);
	}

	VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D_real_inv.rec"),Vec3t(newDimensions.y,newDimensions.z,newDimensions.x),2,fourier::VolumeRotation::ALONG_XY,true);
	VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D_complex_inv.rec"),Vec3t(newDimensions.y,newDimensions.z,newDimensions.x),2,fourier::VolumeRotation::ALONG_XY,true);
	*/
}

void Fourier3D::deleteFiles()
{
	if(!params.DeleteFiles())
		return;

	for(unsigned int i=0; i<xyNumberOfFiles; i++)
	{
		remove(generateTempName(params.OutputFilename(),"fft2D_inv",i*batchZY).c_str());
	}

	for(size_t i=0; i<zyNumberOfFiles; i++)
	{
		remove(generateTempName(params.OutputFilename(),"fft2D",i*batchXY).c_str());
	}

}

void Fourier3D::transposeInPlaceDiffSizes(std::vector<float>& data, size_t originalSizeX, size_t sizeX, size_t sizeY)
{
	std::vector<float> temp;
	temp.resize(originalSizeX*sizeY);

	std::copy(data.begin(),data.begin()+ originalSizeX*sizeY, temp.begin());

	for(size_t x=0; x<sizeX; x++)
	{
		for(size_t y=0; y<sizeY; y++)
		{
			data[y+x*sizeY]=temp[x+y*originalSizeX];
		}
	}
}

void Fourier3D::transposeInPlace(std::vector<float>& data, size_t sizeX, size_t sizeY)
{
	std::vector<float> temp;
	temp.resize(sizeX*sizeY);

	std::copy(data.begin(),data.begin()+ sizeX*sizeY, temp.begin());

	for(size_t x=0; x<sizeX; x++)
	{
		for(size_t y=0; y<sizeY; y++)
		{
			data[y+x*sizeY]=temp[x+y*sizeX];
		}
	}
}

void Fourier3D::generate2DFFTVolume()
{
	cout << "Starter xy processing! \n" << std::endl;

	std::vector<float> xySlice;
	xySlice.resize(dimensions.x*dimensions.y*batchXY);

	std::vector<float> fft_real;
	std::vector<float> fft_complex;

	fft_real.resize(newDimensions.x*dimensions.y);
	fft_complex.resize(newDimensions.x*dimensions.y);

	size_t slizeSize = dimensions.x*dimensions.y;

	std::vector<float> partialVolume;
	partialVolume.resize(newDimensions.x*newDimensions.y*batchXY*2);

	for(size_t z=0; z<dimensions.z; z=z+batchXY)
	{
		if(z%reportXY==0)
		{
			cout << "Processing of xy slice number " << z << " out of " << dimensions.z << " started."<< endl;
		}

		inputVolume->readProjections(&xySlice[0],batchXY,z);

		size_t zDoubleIndex=0;
		for(size_t bz=0; bz < batchXY; bz++)
		{
			for(size_t y=0; y<dimensions.y; y++)
				FFTRoutines::realToComplex1DTransform(dimensions.x,&xySlice[y*dimensions.x+slizeSize*bz],&fft_real[y*newDimensions.x], &fft_complex[y*newDimensions.x],indicesToUse1Dx);

			transposeInPlace(fft_real, newDimensions.x,dimensions.y);
			transposeInPlace(fft_complex, newDimensions.x,dimensions.y);

			for(size_t x=0; x<newDimensions.x; x++)
				FFTRoutines::complex1DTransformForward(dimensions.y,&fft_real[x*dimensions.y], &fft_complex[x*dimensions.y],indicesToUse1Dy);

			transposeInPlaceDiffSizes(fft_real, dimensions.y, newDimensions.y, newDimensions.x);
			transposeInPlaceDiffSizes(fft_complex, dimensions.y, newDimensions.y, newDimensions.x);

			for(size_t px=0; px < newDimensions.x; px++)
				for(size_t py=0; py < newDimensions.y; py++)
				{
					partialVolume[zDoubleIndex+py*batchXY*2+px*batchXY*2*newDimensions.y] = fft_real[px+py*newDimensions.x];
					partialVolume[zDoubleIndex+1+py*batchXY*2+px*batchXY*2*newDimensions.y] = fft_complex[px+py*newDimensions.x];
				}

			zDoubleIndex=zDoubleIndex+2;

		}

		VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D",z),Vec3t((size_t)batchXY*2,newDimensions.y,newDimensions.x),2,fourier::VolumeRotation::ALONG_XY,true);
		VolumeIO::writeVolumeSliceInFloat(generateTempName(params.OutputFilename(),"fft2D",z),partialVolume,newDimensions.x*newDimensions.y*batchXY*2,2);
	}

	cout << "Finished FFT of all XY slices!!! \n" << endl;
}

void Fourier3D::generate3DFFT()
{
	std::vector<ifstream*> zyFiles;

	for(size_t i=0; i<zyNumberOfFiles; i++)
	{
		ifstream* temp = new ifstream(generateTempName(params.OutputFilename(),"fft2D",i*batchXY).c_str(),std::ios::binary);
		(*temp).seekg(sizeof(MRCHeader));
		zyFiles.push_back(temp);
	}

	std::vector<float> realData;
	std::vector<float> complexData;

	realData.resize(newDimensions.y*newDimensions.z);
	complexData.resize(newDimensions.y*newDimensions.z);

	std::vector<float> partialYZSlice;

	partialYZSlice.resize(newDimensions.y*batchXY*2);

	std::vector<std::vector<float>> batchXYSlices_real;
	std::vector<std::vector<float>> batchXYSlices_complex;

	batchXYSlices_real.resize(batchZY);
	batchXYSlices_complex.resize(batchZY);
	for(unsigned int i=0; i<batchZY; i++)
	{
		batchXYSlices_real[i].resize(dimensions.z*newDimensions.y);
		batchXYSlices_complex[i].resize(dimensions.z*newDimensions.y);
	}


	std::vector<float> partialVolume;
	partialVolume.resize(newDimensions.y*newDimensions.z*batchZY*2);


	for(size_t x=0; x <newDimensions.x; x=x+batchZY)
	{
		if(x%reportZY==0)
		{
			cout << "Processing of yz slice number " << x << " out of " << newDimensions.x << " started."<< endl;
		}

		for(unsigned int b=0; b<batchZY; b++)
		{
			for(size_t i=0; i<zyNumberOfFiles; i++)
			{
				(*zyFiles[i]).read((char*)&partialYZSlice[0],newDimensions.y*batchXY*2*sizeof(float));
				for(size_t y =0 ; y < newDimensions.y; y++)
				{
					size_t incr=0;
					for(size_t z=0; z< batchXY*2; z=z+2)
					{
						batchXYSlices_real[b][incr+i*batchXY+y*dimensions.z]=partialYZSlice[z+y*batchXY*2];
						batchXYSlices_complex[b][incr+i*batchXY+y*dimensions.z]=partialYZSlice[z+1+y*batchXY*2];
						incr++;
					}
				}
			}
		}

		size_t xDoubleIndex=0;
		for(unsigned int i=0; i<batchZY; i++)
		{
			for(size_t y=0; y<newDimensions.y; y++)
			{
				FFTRoutines::complex1DTransformForward(dimensions.z,&batchXYSlices_real[i][y*dimensions.z],&batchXYSlices_complex[i][y*dimensions.z],indicesToUse1Dz);
				FFTRoutines::complex1DTransformBackward(newDimensions.z,&batchXYSlices_real[i][y*dimensions.z],&batchXYSlices_complex[i][y*dimensions.z]);
				for(size_t j=0; j<newDimensions.z; j++)
				{
					realData[y+j*newDimensions.y] = batchXYSlices_real[i][j+y*dimensions.z]/newDimensions.z;
					complexData[y+j*newDimensions.y] = batchXYSlices_complex[i][j+y*dimensions.z]/newDimensions.z;
				}
			}

			for(size_t z=0; z<newDimensions.z; z++)
				FFTRoutines::complex1DTransformBackward(newDimensions.y,&realData[z*newDimensions.y],&complexData[z*newDimensions.y]);


			for(size_t z=0; z<newDimensions.z; z++)
				for(size_t y=0; y<newDimensions.y; y++)
				{
					partialVolume[xDoubleIndex+y*batchZY*2+z*batchZY*2*newDimensions.y]=realData[y+z*newDimensions.y];
					partialVolume[xDoubleIndex+1+y*batchZY*2+z*batchZY*2*newDimensions.y]=complexData[y+z*newDimensions.y];
				}

			xDoubleIndex=xDoubleIndex+2;
		}

		VolumeIO::writeHeader(params.InputFilename(),generateTempName(params.OutputFilename(),"fft2D_inv",x),Vec3t((size_t)batchZY*2,newDimensions.y,newDimensions.z),2,fourier::VolumeRotation::ALONG_XY,true);
		VolumeIO::writeVolumeSliceInFloat(generateTempName(params.OutputFilename(),"fft2D_inv",x),partialVolume,batchZY*2*newDimensions.z*newDimensions.y,2);

	}

	for(size_t i=0; i<zyNumberOfFiles; i++)
	{
		(*zyFiles[i]).close();
	}

	cout << "Finished FFT of all Z rows and inverse FFT of all ZY slices!!! \n" << endl;
}

void Fourier3D::inverse1DFFT()
{
	std::vector<ifstream*> xyFiles;

	for(size_t i=0; i<xyNumberOfFiles; i++)
	{
		ifstream* temp = new ifstream(generateTempName(params.OutputFilename(),"fft2D_inv",i*batchZY).c_str(),std::ios::binary);
		(*temp).seekg(sizeof(MRCHeader));
		xyFiles.push_back(temp);
	}

	std::vector<float> partialXYSlice;
	partialXYSlice.resize(newDimensions.y*batchZY*2);

	std::vector<std::vector<float>> batchXYSlices_real;
	std::vector<std::vector<float>> batchXYSlices_complex;

	batchXYSlices_real.resize(batchInvXY);
	batchXYSlices_complex.resize(batchInvXY);
	for(unsigned int i=0; i<batchInvXY; i++)
	{
		batchXYSlices_real[i].resize(newDimensions.x*newDimensions.y);
		batchXYSlices_complex[i].resize(newDimensions.x*newDimensions.y);
	}

	float min= FLT_MAX;
	float max= -FLT_MAX;
	float mean = 0.0f;

	std::vector<float> partialVolume;
	partialVolume.resize(newDimensions.x*newDimensions.y*batchInvXY);

	for(size_t z=0; z <newDimensions.z; z=z+batchInvXY)
	{
		if(z%reportInvXY==0)
		{
			cout << "Processing of xy slice number " << z << " out of " << newDimensions.z << " started."<< endl;
		}

		for(unsigned int b=0; b<batchInvXY; b++)
		{
			for(size_t i=0; i<xyNumberOfFiles; i++)
			{
				(*xyFiles[i]).read((char*)&partialXYSlice[0],newDimensions.y*batchZY*2*sizeof(float));
				for(size_t y =0 ; y < newDimensions.y; y++)
				{
					size_t incr=0;
					for(size_t x=0; x< batchZY*2; x=x+2)
					{
						batchXYSlices_real[b][incr+i*batchZY+y*newDimensions.x]=partialXYSlice[x+y*batchZY*2];
						batchXYSlices_complex[b][incr+i*batchZY+y*newDimensions.x]=partialXYSlice[x+1+y*batchZY*2];
						incr++;
					}
				}
			}

			for(size_t y=0; y<newDimensions.y; y++)
			{
				FFTRoutines::complex1DTransformBackward(newDimensions.x,&batchXYSlices_real[b][y*newDimensions.x],&batchXYSlices_complex[b][y*newDimensions.x]);

				for(size_t i=0; i<newDimensions.x; i++)
				{
					float value = batchXYSlices_real[b][i+y*newDimensions.x]/newDimensions.x;
					partialVolume[i+y*newDimensions.x+b*newDimensions.x*newDimensions.y] = value;
					max=std::max(value,max);
					min=std::min(value,min);
					mean+=value;
				}
			}

		}

		VolumeIO::writeVolumeSliceInFloat(params.OutputFilename(),partialVolume,newDimensions.x*newDimensions.y*batchInvXY,2);

	}

	cout << "Finished inverse FFT of all X rows!!! \n" << endl;

	VolumeIO::convertVolumeValues(params.OutputFilename(), min, max,mean/(float)(newDimensions.z*newDimensions.y*newDimensions.x),2);

	for(size_t i=0; i<xyNumberOfFiles; i++)
	{
		(*xyFiles[i]).close();
	}
}


void Fourier3D::computeBatchSize()
{
	if(params.MemoryLimit()==0)
	{
		batchXY = 1;
		batchZY = 1;
		reportXY = 200;
		reportZY = 200;
		zyNumberOfFiles = dimensions.z;
		xyNumberOfFiles = newDimensions.x;
		cout << "No memory limit was set. The batch size was set to 1!" << std::endl;
		return;
	}

	unsigned int memoryLimit = params.MemoryLimit()*1000000; //from MB to B

	// batchXY computation
	float xyAdditionalMemory = (2*newDimensions.x*dimensions.y+dimensions.x*dimensions.y)*sizeof(float);
	float xyBatchDependentMemory = (dimensions.x*dimensions.y+2*newDimensions.x*newDimensions.y)*sizeof(float);

	batchXY=dimensions.z;

	while((batchXY*xyBatchDependentMemory+xyAdditionalMemory)>=memoryLimit)
	{
		batchXY--;
		if(batchXY==0)
		{
			batchXY=1;
			break;
		}
	}

	while(dimensions.z%batchXY!=0)
	{
		batchXY--;
		if(batchXY==0)
		{
			batchXY=1;
			break;
		}
	}

	zyNumberOfFiles = dimensions.z/batchXY;

	// batchZY computation
	float zyAdditionalMemory = 2*newDimensions.y*(newDimensions.z + batchXY)*sizeof(float);
	float zyBatchDependentMemory = 2*newDimensions.y*(dimensions.z+newDimensions.z)*sizeof(float);

	batchZY=newDimensions.x;

	while((batchZY*zyBatchDependentMemory+zyAdditionalMemory)>=memoryLimit)
	{
		batchZY--;
		if(batchZY==0)
		{
			batchZY=1;
			break;
		}
	}

	while(newDimensions.x%batchZY!=0)
	{
		batchZY--;
		if(batchZY==0)
		{
			batchZY=1;
			break;
		}
	}

	xyNumberOfFiles = newDimensions.x/batchZY;

	// batchInvXY computation
	float xyInvAdditionalMemory = 2*newDimensions.y*batchZY*sizeof(float);
	float xvInvBatchDependentMemory = (3*newDimensions.x*newDimensions.y)*sizeof(float);

	batchInvXY = newDimensions.z;
	while((batchInvXY*xvInvBatchDependentMemory+xyInvAdditionalMemory)>=memoryLimit)
	{
		batchInvXY--;
		if(batchInvXY==0)
		{
			batchInvXY=1;
			break;
		}
	}

	while(newDimensions.z%batchInvXY!=0)
	{
		batchInvXY--;
		if(batchInvXY==0)
		{
			batchInvXY=1;
			break;
		}
	}


	cout << "The batch size for XY was set to " << batchXY << " resulting in " << zyNumberOfFiles << " files." << std::endl;
	cout << "The batch size for YZ was set to " << batchZY << " resulting in " << xyNumberOfFiles << " files." <<  std::endl;
	cout << "The batch size for inverse XY was set to " << batchInvXY << std::endl << std::endl;

	reportXY=batchXY;
	reportZY=batchZY;
	reportInvXY = batchInvXY;

	while(reportXY<200)
	{
		reportXY = reportXY*2;
	}

	while(reportZY<200)
	{
		reportZY = reportZY*2;
	}

	while(reportInvXY<200)
	{
		reportInvXY = reportInvXY*2;
	}


}




void Fourier3D::openFile(ifstream& file, string filename)
{
	file.open(filename.c_str(), std::ios::binary);

	file.seekg(sizeof(MRCHeader));

	if (!file.good())
	{
		throw ExceptionFileOpen(filename);
	}
}

string Fourier3D::generateTempName(std::string base, std::string extension, int index)
{
	stringstream tempName;
	tempName << base << "_" << extension;

	if(index!=-1)
	{
		tempName << "_part" << (unsigned int) index;
	}

	tempName << ".rec" ;

	return tempName.str();
}


void Fourier3D::computeCroppedIndices2D(std::vector<size_t>& indices, size_t oldDimX, size_t oldDimY, size_t newDimX, size_t newDimY)
{
	indices.reserve(newDimX*newDimY);

	size_t newDimHalfX=newDimX/2;
	size_t newDimHalfY=newDimY/2;

	for(size_t y=0; y<oldDimY; y++)
	{
		if(y<newDimHalfY || y>=(oldDimY-newDimHalfY))
		{
			for(size_t x=0; x < newDimHalfX; x++)
			{
				indices.push_back(x+y*oldDimX);
			}

			for(size_t x=(oldDimX-newDimHalfX); x<oldDimX; x++)
			{
				indices.push_back(x+y*oldDimX);
			}
		}
	}

	if(indices.size()!=newDimX*newDimY)
	{
		cout << "Indices to keep do not correspond to the new dimensions" << endl;
		exit(EXIT_FAILURE);
	}
}

void Fourier3D::computeCroppedIndices1D(std::vector<size_t>& indices, size_t oldDim, size_t newDim)
{
	indices.reserve(newDim);
	size_t newDimHalf = newDim / 2;
	size_t add = 0;

	if((newDim%2)!=0)
	{
		add=1;
	}


	for(size_t i=0; i < newDimHalf+add; i++)
	{
		indices.push_back(i);
	}

	for(size_t i=(oldDim-newDimHalf); i<oldDim; i++)
	{
		indices.push_back(i);
	}

	if(indices.size()!=newDim)
	{
		cout << "Indices to keep do not correspond to the new dimensions" << endl;
		exit(EXIT_FAILURE);
	}
}


