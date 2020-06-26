#include "volumeIO.h"
#include "parameterSetup.h"
#include "mrcHeader.h"
#include <fstream>
#include <iostream>
#include "exception.h"
#include <stdio.h>
#include <cstring>
#include <cstdio>

using namespace std;


void VolumeIO::writeMRCStack(MRCHeader header, std::vector<std::vector<float>>& data, string outputName)
{
	std::vector<float> linearData;
	linearData.resize(data.size()*data[0].size());
	size_t k=0;

	for(size_t i=0; i<data.size(); i++)
	{
		for (size_t j=0; j<data[0].size(); j++)
		{
			linearData[k]=data[i][j];
			k++;
		}
	}

	writeMRCStack(header, linearData, outputName);
}

void VolumeIO::writeMRCStack(MRCHeader header, std::vector<float>& data, string outputName)
{
	std::ofstream outfile;

	try
	{
		outfile.open(outputName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

		if(!outfile.good())
		{
			throw new ExceptionFileOpen(outputName);
		}


		outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));

		writeProjections(outfile,data,header.nx*header.ny, header.nz);

		outfile.close();

	}
	catch(ExceptionFileOpen& e)
	{
			cout << e.what() << endl;
	}
	catch(ExceptionFileFormat& e)
	{
				cout << e.what() << endl;
	}
	catch(...)
	{
		cout << "Error while writing volume to disk!" << endl;
	}
}

void VolumeIO::writeProjections(std::ofstream& outfile, std::vector<float>& data, size_t projectionSize, unsigned int numberOfProjections)
{
	float* currentProjection = new float[projectionSize];

	for(unsigned int z = 0; z < numberOfProjections; z++)
	{
		for(size_t i= 0; i < projectionSize; i++)
		{
			currentProjection[i] = data[i+z*projectionSize];
		}

		outfile.write((const char*)currentProjection, projectionSize*sizeof(float));
	}

	delete[] currentProjection;

}


void VolumeIO::write(std::vector<float>& data, fourier::Vec3ui resolution, std::string outputVolumeFileName, int mode)
{
	fourier::VolumeRotation rotation = fourier::VolumeRotation::ALONG_XZ;
    write( data, resolution, outputVolumeFileName, mode, rotation );
}


void VolumeIO::writeHeader(std::string inputStackFileName, std::string outputVolumeFileName, fourier::Vec3t volumeDimensions, int mode, fourier::VolumeRotation rotation, bool changeMode)
{
	fourier::Vec3ui dims;
	dims.x=(unsigned int)volumeDimensions.x;
	dims.y=(unsigned int)volumeDimensions.y;
	dims.z=(unsigned int)volumeDimensions.z;

	writeHeader(inputStackFileName, outputVolumeFileName, dims, mode, rotation,changeMode);
}

void VolumeIO::writeHeader(std::string inputStackFileName, std::string outputVolumeFileName, fourier::Vec3ui volumeDimensions, int mode, fourier::VolumeRotation rotation, bool changeMode)
{
	ifstream infile;
	infile.open(inputStackFileName.c_str(), std::ios::binary);

	if (!infile.good())
	{
		 throw ExceptionFileOpen(inputStackFileName);
	}

	MRCHeader header;

	infile.read((char*)&header, sizeof(MRCHeader));

	if(rotation==fourier::VolumeRotation::ALONG_XZ)
	{
		header.nx = volumeDimensions.x;
		header.ny = volumeDimensions.z;
		header.nz = volumeDimensions.y;

		header.mx = volumeDimensions.x;
		header.my = volumeDimensions.z;
		header.mz = volumeDimensions.y;

		float cellValue = header.my;
		header.cellDimY = header.cellDimZ;
		header.cellDimZ = cellValue;
	}
	else //if (rotation==ALONG_XY)
	{
		header.nx = volumeDimensions.x;
		header.ny = volumeDimensions.y;
		header.nz = volumeDimensions.z;
	}

	if(changeMode)
		header.mode = mode;

	std::ofstream outfile;
	outfile.open(generateTempName(outputVolumeFileName, mode).c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

	if(!outfile.good())
	{
		throw new ExceptionFileOpen(generateTempName(outputVolumeFileName, mode));
	}

	outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));
	outfile.close();
	infile.close();
}

std::string VolumeIO::generateTempName(std::string filename, int mode)
{
	std::stringstream newFilename;

	newFilename << filename;

	if(mode!=2)
		newFilename << "_temp";

	return newFilename.str();

}

void VolumeIO::convertVolumeValues(std::string outputVolumeFileName, float min, float max, float mean, int mode)
{
	float newMin;
	float newMax;

	switch(mode)
	{
		case 0: newMin = 0.f;
				newMax = 255.f;
				break;
		case 1: newMin = -32768.f;
				newMax = 32767.f;
				break;
		case 2: newMin = min;
				newMax = max;
				break;
		case 6: newMin = 0.f;
				newMax = 65535.f;
				break;
		default:	throw new ExceptionFileFormat();
				break;
	}

	float newMean = fourier::convertRange(mean,min, max,newMin, newMax);

	MRCHeader header;
	correctHeader(outputVolumeFileName,header,newMin,newMax,newMean,mode);

	size_t sliceSize = header.nx * header.ny;

	switch(mode)
	{
		case 0: convertSliceData<unsigned char>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
				break;
		case 1: convertSliceData<signed short>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
				break;
		case 2: break;
		case 6: convertSliceData<unsigned short>(outputVolumeFileName, sliceSize, header.nz, min, max, newMin, newMax);
				break;
		default:	throw new ExceptionFileFormat();
				break;
	}
}

template <typename voxelType>
void VolumeIO::convertSliceData(std::string outputVolumeFileName, size_t sliceSize, unsigned int sliceNumber, float min, float max, float newMin, float newMax)
{
	ifstream originalFile;
	ofstream newFile;

	originalFile.open(generateTempName(outputVolumeFileName, 1).c_str(), std::ios_base::in | std::ios_base::binary);
	newFile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::app);

	if (!newFile.good())
	{
		throw ExceptionFileOpen(outputVolumeFileName);
	}

	if (!originalFile.good())
	{
		throw ExceptionFileOpen(generateTempName(outputVolumeFileName, 1).c_str());
	}

	originalFile.seekg(sizeof(MRCHeader));
	newFile.seekp(sizeof(MRCHeader));

	voxelType* convertedData = new voxelType[sliceSize];
	float* floatData = new float[sliceSize];

	for(unsigned int slice = 0; slice <sliceNumber; slice++)
	{
		originalFile.read((char*)(floatData), sliceSize * sizeof(float));

		for(size_t i =0; i<sliceSize;i++)
		{
			convertedData[i]=static_cast<voxelType>(fourier::convertRange(floatData[i],min,max,newMin, newMax));
		}

		newFile.write((const char*)convertedData, sliceSize*sizeof(voxelType));
	}

	delete[] floatData;
	delete[] convertedData;

	newFile.close();
	originalFile.close();

	std::remove(generateTempName(outputVolumeFileName, 1).c_str());
}

void VolumeIO::correctHeader(std::string outputVolumeFileName, MRCHeader& header, float min, float max, float mean, int mode)
{
	fstream file;

	file.open(generateTempName(outputVolumeFileName, mode).c_str(), std::ios_base::in | std::ios_base::out | std::ios::binary);

	if (!file.good())
	{
		 throw ExceptionFileOpen(generateTempName(outputVolumeFileName, mode));
	}

	file.seekg(0);
	file.read((char*)&header, sizeof(MRCHeader));

	header.dMin = min;
	header.dMax = max;
	header.dMean = mean;

	header.mode = mode;
	header.extra = 0;

	if(mode==2)
	{
		file.seekp(0);
		file.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));
	}
	else
	{
		ofstream newFile;
		newFile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

		if (!newFile.good())
		{
			throw ExceptionFileOpen(outputVolumeFileName);
		}

		newFile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));
		newFile.close();
	}

	file.close();
}

void VolumeIO::writeVolumeSliceInFloat(std::string outputVolumeFileName, std::vector<float> data,  size_t sliceSize, int mode)
{
	std::ofstream outfile;
	outfile.open(generateTempName(outputVolumeFileName.c_str(),mode), fstream::out | fstream::binary | fstream::app);

	if(!outfile.good())
	{
		throw new ExceptionFileOpen(outputVolumeFileName);
	}

	outfile.write((const char*)&data[0], sliceSize*sizeof(float));
	outfile.close();
}

void VolumeIO::writeVolumeSliceInFloat(std::string outputVolumeFileName, float& data,  size_t sliceSize, int mode)
{
	std::ofstream outfile;
	outfile.open(generateTempName(outputVolumeFileName.c_str(),mode), fstream::out | fstream::binary | fstream::app);

	if(!outfile.good())
	{
		throw new ExceptionFileOpen(outputVolumeFileName);
	}

	outfile.write((const char*)(&data), sliceSize*sizeof(float));
	outfile.close();
}


void VolumeIO::computeMinMaxMean(std::vector<float>& data, size_t dataSize, float& min, float& max, float& mean)
{
	min=FLT_MAX;
	max=-FLT_MAX;
	mean=0.0f;

	for(size_t i=0; i<dataSize; i++)
	{
		min=std::min(min,data[i]);
		max=std::max(max,data[i]);
		mean+=data[i];
	}

	mean/=dataSize;
}

void VolumeIO::write(std::vector<float>& data, fourier::Vec3ui resolution, std::string outputVolumeFileName, int mode, fourier::VolumeRotation rotation)
{
	std::ofstream outfile;

    try
    {
        outfile.open(outputVolumeFileName.c_str(), std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

        if(!outfile.good())
        {
            throw new ExceptionFileOpen(outputVolumeFileName);
        }

        MRCHeader header;
        memset(&header, '\0', sizeof(MRCHeader));

        header.mode = mode;

        if(rotation==fourier::VolumeRotation::ALONG_XZ)
        {
            header.nx = resolution.x;
            header.ny = resolution.z;
            header.nz = resolution.y;

            header.mx = resolution.x;
            header.my = resolution.z;
            header.mz = resolution.y;

            header.cellDimX = (float) resolution.x;
            header.cellDimY = (float) resolution.z;
            header.cellDimZ = (float) resolution.y;
        }
        else //if (rotation==ALONG_XY)
        {
            header.nx = resolution.x;
            header.ny = resolution.y;
            header.nz = resolution.z;

            header.mx = resolution.x;
            header.my = resolution.y;
            header.mz = resolution.z;

            header.cellDimX = (float) resolution.x;
            header.cellDimY = (float) resolution.y;
            header.cellDimZ = (float) resolution.z;
        }

        header.map[0] = 'M';
        header.map[1] = 'A';
        header.map[2] = 'P';
        header.map[3] = ' ';

        header.mapC = 1;
        header.mapR = 2;
        header.mapS = 3;
        header.extra = 0;
        header.nint = 0;
        header.nreal = 0;

        float min, max, mean;

        computeMinMaxMean(data,(size_t)(resolution.x*resolution.y*resolution.z),min, max, mean);

        switch(header.mode)
		{
			case 0: header.dMin = 0.f;
					header.dMax = 255.f;
					header.dMean = fourier::convertRange(mean,min,max,header.dMin, header.dMax);
					break;
			case 1: header.dMin = -32768.f;
					header.dMax = 32767.f;
					header.dMean = fourier::convertRange(mean,min,max,header.dMin, header.dMax);
					break;
			case 2: header.dMin = min;
					header.dMax = max;
					header.dMean = mean;
					break;
			case 6: header.dMin = 0.f;
					header.dMax = 65535.f;
					header.dMean = fourier::convertRange(mean,min,max,header.dMin, header.dMax);
					break;
			default:	throw new ExceptionFileFormat();
					break;
		}

        outfile.write(reinterpret_cast<char*>(&header), sizeof(MRCHeader));

        switch( mode )
        {
        case 0: writeData<unsigned char>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        case 1: writeData<signed short>(data, resolution, outfile, rotation, min, max,  header.dMin, header.dMax);
            break;
        case 2: writeData<float>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        case 6: writeData<unsigned short>(data, resolution, outfile, rotation, min, max, header.dMin, header.dMax);
            break;
        default:
        	throw new ExceptionFileFormat();
            break;
        }

        outfile.close();

    }
    catch(ExceptionFileOpen& e)
    {
            cout << e.what() << endl;
    }
    catch(ExceptionFileFormat& e)
    {
                cout << e.what() << endl;
    }
    catch(...)
    {
        cout << "Error while writing volume to disk!" << endl;
    }
}


template <typename voxelType>
void VolumeIO::writeData(std::vector<float>& data, fourier::Vec3ui resolution, std::ofstream& outfile, fourier::VolumeRotation rotation, float oldMin, float oldMax,float newMin, float newMax )
{
	cout << "Writing the data..."<< endl;
	if(rotation==fourier::VolumeRotation::ALONG_XZ)
       writeDataAlongXZPlane<voxelType>(data, resolution, outfile, oldMin, oldMax, newMin, newMax);
    else
       writeDataAlongXYPlane<voxelType>(data, resolution, outfile, oldMin, oldMax, newMin, newMax);
}

// slices along xz plane - classic eTomo reconstruction output
template <typename voxelType>
void VolumeIO::writeDataAlongXZPlane(std::vector<float>& data, fourier::Vec3ui resolution, std::ofstream&  outfile ,float oldMin, float oldMax,float newMin, float newMax)
{

	voxelType* dataToWrite = new voxelType[resolution.x * resolution.z];

	for(unsigned int y = 0; y < resolution.y; y++)
	{
		size_t bufferIndex = 0;

		for(size_t z = 0; z < resolution.z; z++)
		{
			for (size_t x = 0; x < resolution.x; x++)
			{
				dataToWrite[bufferIndex] = static_cast<voxelType>(fourier::convertRange(data[x+y*resolution.x+z*resolution.x*resolution.y],oldMin,oldMax,newMin, newMax));
				++bufferIndex;
			}
		}

		outfile.write((const char*)dataToWrite, (size_t)resolution.x * (size_t)resolution.z*sizeof(voxelType));
	}

	delete[] dataToWrite;
}

template <typename voxelType>
void VolumeIO::writeDataAlongXYPlane(std::vector<float>& data, fourier::Vec3ui resolution, std::ofstream&  outfile, float oldMin, float oldMax, float newMin, float newMax)
{

	voxelType* dataToWrite = new voxelType[resolution.x * resolution.y];

	for(size_t z = 0; z < resolution.z; z++)
	{
		size_t bufferIndex = 0;

		for(size_t y= 0; y < resolution.y; y++)
		{
			for (size_t x = 0; x < resolution.x; x++)
			{
				dataToWrite[bufferIndex] = static_cast<voxelType>(fourier::convertRange(data[x+y*resolution.x+z*resolution.x*resolution.y],oldMin,oldMax,newMin, newMax));
				++bufferIndex;
			}
		}

		outfile.write((const char*)dataToWrite, (size_t)resolution.x * (size_t)resolution.y*sizeof(voxelType));
	}

	delete[] dataToWrite;
}

std::ostream& operator<<(std::ostream& oss, int mode)
{
    switch (mode)
    {
    case 2:
        oss << "gray scale (32 bit float)";
        break;
    case 0:
        oss << "gray scale (8 bit unsigned)";
        break;
    case 6:
        oss << "gray scale (16 bit unsigned)";
        break;
    case 1:
        oss << "gray scale (16 bit signed)";
        break;
    }
    return oss;
}
