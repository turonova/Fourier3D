#pragma once
#include <fftw3.h>
#include <vector>
#include "common.h"

class FFTRoutines
{
public:
	static void real1DTransform(size_t size, std::vector<float>& data, std::vector<float>& fft);
	static void realToComplex1DTransform(size_t sizeX, float* data, float* fft_real, float* fft_complex, std::vector<size_t>& indicesToCopy);
	static void complex1DTransformForward(size_t sizeX, std::vector<float>& fft_real, std::vector<float>& fft_complex, std::vector<size_t>& indicesToCopy);
	static void complex1DTransformForward(size_t sizeX, float*  fft_real, float*  fft_complex, std::vector<size_t>& indicesToCopy);
	static void complex1DTransformBackward(size_t sizeX, std::vector<float>& fft_real, std::vector<float>& fft_complex);
	static void complex1DTransformBackward(size_t sizeX, float* fft_real, float* fft_complex);

	static void real2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<double>& mask, fourier::DataStats& dataStats);
	static void complex2DTransformForward(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& fft_real, std::vector<float>& fft_complex, std::vector<size_t> indicesToCopy);
	static void complex2DTransformBackward(size_t sizeX, size_t sizeY, std::vector<float>& fft_real, std::vector<float>& fft_complex);
	static void computePowerSpectrum(std::vector<float>& powerSpectrum, fftw_complex* fftOut, size_t sizeX, size_t nyh, bool logarithmizeData);
	static void complex2DTransform(size_t sizeX, size_t sizeY, std::vector<float>& data, std::vector<float>& filter);
	static void odfft2(float* array, int nxp, int nyp, int idirp);
private:
	static double logarithmizeValue(double value, bool logarithmizeData);
	static void maskFFT(fftw_complex* fftOut, std::vector<double>& mask, size_t sizeX, size_t nyh);
	static void filterFFT(fftw_complex* fftOut, std::vector<float>& filter, size_t sizeX, size_t nyh);
	static void normalizeValues(std::vector<float>& normalizedData, std::vector<double>& originalData, size_t dataSize, fourier::DataStats& dataStats);
	static void normalize(float *array, float scale, size_t dataSize);

};
