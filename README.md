# Installation
 Sys.setenv("PKG_CXXFLAGS"="-std=c++11")  
 library(devtools)  
 install_github("yiyangphd/HSVSM")  

Please install [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries) if the installation failed with warning messages like "ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/...". 
