# ---- Folder ----
MAFOT_DIR = /home/wingen/src/MAFOT
BIN_DIR = /home/wingen/bin
LIB_DIR = /home/wingen/lib/64

# ---- Support M3DC1 ----
# True or False
M3DC1 = True

# ---- Compiler ----
CXX = g++
CFLAGS = -O3 -fPIC
OMPFLAGS = -fopenmp -w

F90 = ifort
F90FLAGS = -O3 -fPIC

# ---- Linker ----
ARCH = ar crus

LDD = g++ --shared
LDFLAGS =

# ---- external Libraries ---- 
BLITZLIBS = -L/home/wingen/lib/64/blitz/lib -lblitz
FLIBS = -L/opt/intel/composerxe-2011.0.084/lib/intel64 -lsvml -lifcore -lirc -ldl -limf
OMPLIBS = -L/act/openmpi-1.6/gcc/lib -lmpi -lmpi_cxx
M3DC1LIBS = -L/home/ferraro/lib/_drop -lfusionio -lm3dc1
HDF5LIBS = -L/home/wingen/lib/64/hdf5/lib -lhdf5 -lz -Wl,-rpath,/home/wingen/lib/64/hdf5/lib -lstdc++

# ---- external Includes ----
BLITZINCLUDE = -I/home/wingen/lib/64/blitz/include
OMPINCLUDE = -I/act/openmpi-1.6/gcc/include -I/act/openmpi-1.6/gcc/include/openmpi
M3DC1INCLUDE = -I/home/ferraro/include
