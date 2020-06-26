#include "fftRoutines.h"
#include <fftw3.h>
#include "common.h"
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "volumeIO.h"

using namespace std;

void FFTRoutines::real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft)
{

		double* fftIn;
	  fftw_complex* fftOut;
	  fftw_plan plan_forward;

	  fftIn = (double*) fftw_malloc ( sizeof ( double ) * size );

	  for(unsigned int i=0; i<size;i++)
	  {
		  fftIn[i] = (double)data[i];
	  }

	  size_t fftSize = ( size / 2 ) + 1;

	  fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * fftSize );

	  plan_forward = fftw_plan_dft_r2c_1d( size, fftIn, fftOut, FFTW_ESTIMATE );

	  fftw_execute( plan_forward );

	  if (size % 2 == 0) /* N is even */
	 	  		  fft.push_back(fftOut[size/2][0]);

	  for (int k = (size+1)/2-1; k >=1 ; k--)  // (k < N/2 rounded up)
		  fft.push_back(fftOut[k][0]);

	  for (int k = 1; k < (int)(size+1)/2; k++)  // (k < N/2 rounded up)
		  fft.push_back(fftOut[k][0]);
	  if (size % 2 == 0) // N is even
		  fft.push_back(fftOut[size/2][0]);  // Nyquist freq.

	  fftw_destroy_plan ( plan_forward );

	  fftw_free ( fftOut );
	  fftw_free ( fftIn );

	  return;

}

double FFTRoutines::logarithmizeValue(double value, bool logarithmizeData)
{
	if(logarithmizeData && value>=1.0)
		return log(value);
	else
		return value;
}

void inversefftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)
{
	for(size_t i=0; i<ydim/2;i++)
	{
		for(size_t j=0; j<xdim/2;j++)
		{
			size_t outIndex1 = (j+xdim/2)+(i+ydim/2)*xdim;
			size_t inIndex1 = j+i*xdim;
			size_t outIndex2 = j+i*xdim;
			size_t inIndex2 = (j+xdim/2)+(i+ydim/2)*xdim;
			out[inIndex1]=in[outIndex1];	//1->4
			out[inIndex2]=in[outIndex2];	//4->1
		}
	}

	for(size_t i=0; i<ydim/2;i++)
	{
		for(size_t j=xdim/2; j<xdim;j++)
		{
			size_t outIndex1 = (j-xdim/2)+(i+ydim/2)*xdim;
			size_t inIndex1 = j+i*xdim;
			size_t outIndex2 =j+i*xdim;
			size_t inIndex2 = (j-xdim/2)+(i+ydim/2)*xdim;
			out[inIndex1]=in[outIndex1];	//2->3
			out[inIndex2]=in[outIndex2];	//3->2
		}
	}
}


void fftshift(std::vector<float>& out, std::vector<float>& in, size_t xdim, size_t ydim)
{
	for(size_t i=0; i<ydim/2;i++)
	{
		for(size_t j=0; j<xdim/2;j++)
		{
			size_t outIndex1 = (j+xdim/2)+(i+ydim/2)*xdim;
			size_t inIndex1 = j+i*xdim;
			size_t outIndex2 = j+i*xdim;
			size_t inIndex2 = (j+xdim/2)+(i+ydim/2)*xdim;
			out[outIndex1]=in[inIndex1];	//1->4
			out[outIndex2]=in[inIndex2];	//4->1
		}
	}

	for(size_t i=0; i<ydim/2;i++)
	{
		for(size_t j=xdim/2; j<xdim;j++)
		{
			size_t outIndex1 = (j-xdim/2)+(i+ydim/2)*xdim;
			size_t inIndex1 = j+i*xdim;
			size_t outIndex2 =j+i*xdim;
			size_t inIndex2 = (j-xdim/2)+(i+ydim/2)*xdim;
			out[outIndex1]=in[inIndex1];	//2->3
			out[outIndex2]=in[inIndex2];	//3->2
		}
	}
}

void FFTRoutines::computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData)
{
	size_t k=0;
	for(size_t j=0; j<nyh; j++)
	{
		for(size_t i = 0; i<sizeX; i++)
		{
			powerSpectrum[k] = logarithmizeValue(fftOut[j+i*nyh][0]*fftOut[j+i*nyh][0]+fftOut[j+i*nyh][1]*fftOut[j+i*nyh][1],logarithmizeData);
			k++;
		}
	}

	for(int j=nyh-2; j>0; j--)
	{
		powerSpectrum[k]= fftOut[j][0]*fftOut[j][0]+fftOut[j][1]*fftOut[j][1];
		k++;
		for(int i = sizeX-1; i>0; i--)
		{
			powerSpectrum[k] = logarithmizeValue(fftOut[j+i*nyh][0]*fftOut[j+i*nyh][0]+fftOut[j+i*nyh][1]*fftOut[j+i*nyh][1],logarithmizeData);
			k++;
		}
	}
}

void FFTRoutines::maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh)
{
	size_t k=0;
	for(size_t j=0; j<nyh; j++)
	{
		for(size_t i = 0; i<sizeX; i++)
		{
			if(mask[k]!=1.0f)
			{
				fftOut[j+i*nyh][0] = fftOut[j+i*nyh][0]*mask[k];
				fftOut[j+i*nyh][1] = fftOut[j+i*nyh][1]*mask[k];
			}
			k++;
		}
	}
}

void FFTRoutines::filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh)
{
	size_t k=0;
	for(size_t j=0; j<nyh; j++)
	{
		for(size_t i = 0; i<sizeX; i++)
		{
			if(filter[k]!=1.0f)
			{
				fftOut[j+i*nyh][0] = fftOut[j+i*nyh][0]*filter[k];
				fftOut[j+i*nyh][1] = fftOut[j+i*nyh][1]*filter[k];
			}
			k++;
		}
	}
}

void FFTRoutines::normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, fourier::DataStats& dataStats)
{
	for ( size_t i = 0; i < dataSize; i++ )
	{
		normalizedData[i]=originalData[i]/(dataSize);

		dataStats.mean+=normalizedData[i];
		dataStats.max = max(dataStats.max, normalizedData[i]);
		dataStats.min = min(dataStats.min, normalizedData[i]);
	}
}

void FFTRoutines::complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter)
{
	fftw_complex* fftOut;
	fftw_complex* fftIn;
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );
	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );

	for ( size_t i = 0; i < sizeX*sizeY; i++ )
	{
		fftIn[i][0] = (double)data[i];
		fftIn[i][1] = 0.0;
	}

	fftw_plan plan_forward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, -1, FFTW_ESTIMATE);
	fftw_execute ( plan_forward );

	std::vector<float> isFilter;
	isFilter.resize(filter.size());

	inversefftshift(isFilter,filter, sizeX, sizeY);
	filterFFT(fftOut, isFilter, sizeX, sizeY);

	fftw_plan plan_backward = fftw_plan_dft_2d ( sizeX, sizeY, fftOut, fftIn, 1, FFTW_ESTIMATE );
	fftw_execute ( plan_backward );

	fftw_destroy_plan ( plan_backward );
	fftw_destroy_plan ( plan_forward );
	fftw_free ( fftOut );

	for(size_t i= 0; i < data.size(); i++)
	{
		data[i]=fftIn[i][0]/data.size();
	}

	fftw_free(fftIn);

}

void FFTRoutines::complex2DTransformForward(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& fft_real, std::vector<float>& fft_complex, std::vector<size_t> indicesToCopy)
{
	fftw_complex* fftOut;
	fftw_complex* fftIn;
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );
	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );

	for ( size_t i = 0; i < sizeX*sizeY; i++ )
	{
		fftIn[i][0] = (double)data[i];
		fftIn[i][1] = 0.0;
	}

	fftw_plan plan_forward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, -1, FFTW_ESTIMATE);
	fftw_execute ( plan_forward );

	fftw_destroy_plan ( plan_forward );
	fftw_free(fftIn);

	for(size_t i=0; i<indicesToCopy.size();i++)
	{
		fft_real[i]=fftOut[indicesToCopy[i]][0];
		fft_complex[i]=fftOut[indicesToCopy[i]][1];
	}

	fftw_free ( fftOut );
}

void FFTRoutines::complex2DTransformBackward(size_t sizeX, size_t sizeY, std::vector<float>& fft_real, std::vector<float>& fft_complex)
{
	fftw_complex* fftOut;
	fftw_complex* fftIn;
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );
	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX * sizeY );

	for ( size_t i = 0; i < sizeX*sizeY; i++ )
	{
		fftIn[i][0] = fft_real[i];
		fftIn[i][1] = fft_complex[i];
	}

	fftw_plan plan_backward = fftw_plan_dft_2d(sizeX, sizeY, fftIn, fftOut, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute ( plan_backward );

	fftw_destroy_plan ( plan_backward );
	fftw_free(fftIn);

	for(size_t i=0; i<sizeX*sizeY;i++)
	{
		fft_real[i]=fftOut[i][0]/(sizeX*sizeY);
		fft_complex[i]=fftOut[i][1]/(sizeX*sizeY);
	}

	fftw_free ( fftOut );
}

void FFTRoutines::complex1DTransformBackward(size_t sizeX, std::vector<float>& fft_real, std::vector<float>& fft_complex)
{
	fftw_complex* fftIn;
	fftw_complex* fftOut;

	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );

	for(unsigned int i=0; i<sizeX;i++)
	{
		fftIn[i][0] = fft_real[i];
		fftIn[i][1] = fft_complex[i];
	}

	fftw_plan plan_backward = fftw_plan_dft_1d( sizeX, fftIn, fftOut, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute( plan_backward );
	fftw_free ( fftIn );


	for(size_t i=0; i<sizeX;i++)
	{
		fft_real[i]=fftOut[i][0];
		fft_complex[i]=fftOut[i][1];
	}

	fftw_destroy_plan ( plan_backward );
	fftw_free ( fftOut );

	return;
}

void FFTRoutines::complex1DTransformBackward(size_t sizeX, float* fft_real, float* fft_complex)
{
	fftw_complex* fftIn;
	fftw_complex* fftOut;

	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );

	for(unsigned int i=0; i<sizeX;i++)
	{
		fftIn[i][0] = fft_real[i];
		fftIn[i][1] = fft_complex[i];
	}

	fftw_plan plan_backward = fftw_plan_dft_1d( sizeX, fftIn, fftOut, FFTW_BACKWARD, FFTW_ESTIMATE );
	fftw_execute( plan_backward );
	fftw_free ( fftIn );


	for(size_t i=0; i<sizeX;i++)
	{
		fft_real[i]=fftOut[i][0];
		fft_complex[i]=fftOut[i][1];
	}

	fftw_destroy_plan ( plan_backward );
	fftw_free ( fftOut );

	return;
}

void FFTRoutines::complex1DTransformForward(size_t sizeX, std::vector<float>& fft_real, std::vector<float>& fft_complex, std::vector<size_t>& indicesToCopy)
{
	fftw_complex* fftIn;
	fftw_complex* fftOut;


	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );

	for(unsigned int i=0; i<sizeX;i++)
	{
		fftIn[i][0] = fft_real[i];
		fftIn[i][1] = fft_complex[i];
	}

	fftw_plan plan_forward = fftw_plan_dft_1d( sizeX, fftIn, fftOut, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute( plan_forward );
	fftw_free ( fftIn );


	for(size_t i=0; i<indicesToCopy.size();i++)
	{
		fft_real[i]=fftOut[indicesToCopy[i]][0];
		fft_complex[i]=fftOut[indicesToCopy[i]][1];
	}


	fftw_destroy_plan ( plan_forward );
	fftw_free ( fftOut );

	return;
}

void FFTRoutines::complex1DTransformForward(size_t sizeX, float* fft_real,float* fft_complex, std::vector<size_t>& indicesToCopy)
{
	fftw_complex* fftIn;
	fftw_complex* fftOut;


	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );

	for(unsigned int i=0; i<sizeX;i++)
	{
		fftIn[i][0] = fft_real[i];
		fftIn[i][1] = fft_complex[i];
	}

	fftw_plan plan_forward = fftw_plan_dft_1d( sizeX, fftIn, fftOut, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute( plan_forward );
	fftw_free ( fftIn );


	for(size_t i=0; i<indicesToCopy.size();i++)
	{
		fft_real[i]=fftOut[indicesToCopy[i]][0];
		fft_complex[i]=fftOut[indicesToCopy[i]][1];
	}


	fftw_destroy_plan ( plan_forward );
	fftw_free ( fftOut );

	return;
}

void FFTRoutines::realToComplex1DTransform(size_t sizeX, float* data, float* fft_real, float* fft_complex, std::vector<size_t>& indicesToCopy)
{
	fftw_complex* fftIn;
	fftw_complex* fftOut;

	fftIn = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );
	fftOut = (fftw_complex*) fftw_malloc ( sizeof ( fftw_complex ) * sizeX );

	for(unsigned int i=0; i<sizeX;i++)
	{
		fftIn[i][0] = data[i];
		fftIn[i][1] = 0.0;
	}

	fftw_plan plan_forward = fftw_plan_dft_1d( sizeX, fftIn, fftOut, FFTW_FORWARD, FFTW_ESTIMATE );
	fftw_execute( plan_forward );
	fftw_free ( fftIn );

	for(size_t i=0; i<indicesToCopy.size();i++)
	{
		fft_real[i]=fftOut[indicesToCopy[i]][0];
		fft_complex[i]=fftOut[indicesToCopy[i]][1];
	}


	fftw_destroy_plan ( plan_forward );
	fftw_free ( fftOut );

	return;
}


void FFTRoutines::odfft2(float* array, int nxp, int nyp, int idirp)
{
	int nx = nxp;
	int ny = nyp;
	int idir = idirp;
	int nxpad = 2 * (nx / 2 + 1);
	float invScale = (float)(1.0f / sqrt((double)nx));

	fftwf_plan plan;

	if (idir == 0)
	{
	  plan = fftwf_plan_many_dft_r2c(1, &nxp, ny, array, NULL, 1, nxpad, (fftwf_complex *)array, NULL, 1, nxpad / 2, FFTW_ESTIMATE);
	}
	else if (idir == 1)
	{
	  plan = fftwf_plan_many_dft_c2r(1, &nxp, ny, (fftwf_complex *)array, NULL, 1, nxpad / 2, array, NULL, 1, nxpad,FFTW_ESTIMATE);
	}

	fftwf_execute(plan);

	normalize(array, invScale, (size_t)nxpad * ny);

	fftwf_destroy_plan (plan);
}


/* Normalize the given number of real elements by the scaling factor */
void FFTRoutines::normalize(float *array, float scale, size_t dataSize)
{
  size_t i;
  for (i = 0; i < dataSize; i++)
    array[i] *= scale;
}

