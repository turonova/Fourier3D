# Fourier3D

A program for tomogram binning in Fourier space. It supports MRC format and can handle large tomograms. The amount of RAM used during the binning is specified by the "MemoryLimit" parameter (in MB) - the data will be split accordingly into chunks that will be temporarily stored in the folder from which the program is run. There is no parallelization in the code. 

The program was tested using foss 2017b. 

Install:

`make includepath="path_to_fftw_include_files" libpath="path_to_fftw_libraries"`

Example:

`make includepath="/path/to/fftw/fftw-3.3.4/include" libpath="/path/to/fftw/fftw-3.3.4/lib /path/to/fftw/fftw-3.1.2/lib64"`


Usage:

`./Fourier3D -BinFactor 2 -InputFile input_tomogram_name -OutputFile output_tomogram_name -MemoryLimit 8000`


The source code for Fourier3D is distributed under an GPL v.3 license. The code can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later version.

The code is distributed WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
