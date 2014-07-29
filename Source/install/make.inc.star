# ---- Folder ----
MAFOT_DIR = /home/wingen/src/MAFOT
BIN_DIR = /home/wingen/bin
LIB_DIR = /home/wingen/lib/32

# ---- Support M3DC1 ----
# True or False
M3DC1 = False

# ---- Compiler ----
CXX = g++
CFLAGS = -m32 -O3 -fPIC -Dgenatomix
OMPFLAGS = -fopenmp -w

F90 = gfortran
F90FLAGS = -m32 -O3 -fPIC

# ---- Linker ----
ARCH = ar crus

LDD = g++ --shared
LDFLAGS = -m32

# ---- external Libraries ---- 
BLITZLIBS = -L/home/wingen/lib/32/blitz/lib -lblitz
FLIBS = -L/usr/lib -lgfortran
OMPLIBS = -L/home/wingen/lib/32/openmpi/lib -lmpi -lmpi_cxx
M3DC1LIBS = 
HDF5LIBS = 

# ---- external Includes ----
BLITZINCLUDE = -I/home/wingen/lib/32/blitz/include
OMPINCLUDE = -I/home/wingen/lib/32/openmpi/include -I/home/wingen/lib/32/openmpi/include/openmpi
M3DC1INCLUDE = 

# ---- build GUI with python (optional) ----
PYINSTALLER =
