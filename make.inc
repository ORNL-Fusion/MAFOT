# ---- Folder ----
MAFOT_DIR = /home/tom/source/MAFOT
#use this for local MAFOT builds:
BIN_DIR = /usr/bin
LIB_DIR = /usr/lib
#Use this for building appImage:
#BIN_DIR = /home/tom/source/HEAT/MAFOT/bin
#LIB_DIR = /home/tom/source/HEAT/MAFOT/lib



# ---- Support codes ----
# True or False
M3DC1 = False
VMEC = False
SIESTA = False

# ---- Compiler ----
CXX = g++
CFLAGS = -O3 -fPIC
OMPFLAGS = -fopenmp -w

F90 = gfortran
F90FLAGS = -O3 -fPIC

# ---- Linker ----
ARCH = ar crus

LDD = g++ --shared
LDFLAGS =

# ---- external Libraries ----
BLITZLIBS = -L/home/tom/source/MAFOT/blitz++/blitz-0.9/lib -lblitz
FLIBS = -L/usr/lib/x86_64-linux-gnu -lgfortran
OMPLIBS = -L/usr/lib/x86_64-linux-gnu/openmpi/lib -lmpi -lmpi_cxx
M3DC1LIBS =
HDF5LIBS =
NETCDFLIBS = -L/usr/lib/x86_64-linux-gnu -lnetcdf

# ---- external Includes ----
BLITZINCLUDE = -I/home/tom/source/MAFOT/blitz++/blitz-0.9/blitz
OMPINCLUDE = -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
M3DC1INCLUDE =
NETCDFINCLUDE = -I/usr/include

# ---- build GUI with python (optional) ----
#PYINSTALLER = /Library/Frameworks/EPD64.framework/Versions/7.2/bin/pyinstaller
