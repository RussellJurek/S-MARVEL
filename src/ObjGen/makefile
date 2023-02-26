# retrieve system architecture
ARCH := $(shell uname -s)
$(info Creating busy function C/C++ library and demo programs for ARCH = $(ARCH))

# specify compilation options
CXX = c++
CC = cc
CFLAGS = -O3 
LDFLAGS = -lm 
ifndef INSTALL_DIR
  INSTALL_DIR = /usr/local/lib
  $(info INSTALL_DIR not defined. Setting to default value. INSTALL_DIR = $(INSTALL_DIR).)
endif
ifndef INSTALL_INC
  INSTALL_INC = /usr/local/include
  $(info INSTALL_INC not defind. Setting to default value. INSTALL_INC = $(INSTALL_INC).)
endif

# set the library creation method
ifeq ($(ARCH),Darwin)
  LIB_METHOD = Libtool -static -o
else 
  LIB_METHOD = ar -crus
endif
$(info Using LIB_METHOD = $(LIB_METHOD))

# determine if cfitsio is installed 
THIS_FITS_DIR = $(THIS_FITS_LIB)$(THIS_FITS_INC)
ifeq ($(THIS_FITS_DIR), )
  $(info Looking for CFITSIO installation . . .)
  THIS_FITS_LIB = $(shell find / -name libcfitsio.* 2>/dev/null | xargs dirname | head -1)
  THIS_FITS_INC = $(shell find / -name *fitsio.h 2>/dev/null | xargs dirname | head -1)
  THIS_FITS_DIR = $(THIS_FITS_LIB)$(THIS_FITS_INC)
  #$(info THIS_FITS_DIR = '$(THIS_FITS_DIR)')
else
  THIS_FITS_LIB = $(FITS_LIB)
  THIS_FITS_INC = $(FITS_INC)
  ifeq ($(shell ls $(THIS_FITS_LIB)/libcfitsio.*), )
    THIS_FITS_LIB = 
    $(info ERROR --- $(THIS_FITS_LIB) doesn't contain a valid cfitsio library file.)
  endif
  ifeq ($(shell ls $(THIS_FITS_INC)/*fitsio.h), )
    THIS_FITS_INC = 
    $(info ERROR --- $(THIS_FITS_INC) doesn't contain a valid fitsio.h or cfitsio.h file.)
  endif
  THIS_FITS_DIR = $(THIS_FITS_LIBC)$(THIS_FITS_INC)
endif
ifeq ($(THIS_FITS_DIR), )
  $(info ERROR --- Complete cfitsio installation not found. Not using cfitsio. Change environment variable FITS_LIB and FITS_INC to \
  a valid installation location to use a cfitsio installation. Specify using, \
  FITS_INC=path-to-fitsio.h from FITS_DIR, FITS_LIB=path-to-libcfitsio.a/.so/.dylib library.)
else
  $(info cfitsio installation found at $(THIS_FITS_LIB) and $(THIS_FITS_INC). Using cfitsio.)
endif

# determine if pgplot is installed
THIS_PGPLOT_DIR = $(PGPLOT_DIR)
ifneq ($(THIS_PGPLOT_DIR), )
  $(info PGPLOT_DIR = $(THIS_PGPLOT_DIR), checking libraries can be found here . . . )
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/libpgplot.a), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- libpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/libcpgplot.a), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- libcpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_PGPLOT_DIR)/cpgplot.h), )
    THIS_PGPLOT_DIR = 
  $(info ERROR --- cpgplot.h not found.)
  endif
else
  THIS_PGPLOT_DIR =
endif
ifeq ($(THIS_PGPLOT_DIR), )
  $(info ERROR --- Complete pgplot installation not found. Not using pgplot. Change environment variable PGPLOT_DIR to \
  a valid installation directory to use a pgplot installation. Specify using, make PGPLOT_DIR=path-to-installation .)
else
  $(info Pgplot installation found at $(THIS_PGPLOT_DIR). Using pgplot.)
endif

# determine if wcslib is installed --- try /usr/local/lib & /usr/local/include/wcslib if nothing specified
# use libwcs.a as the default library
THIS_WCS_DIR = $(THIS_WCS_LIB)$(THIS_WCS_INC)
ifeq ($(WCS_LIB_NAME), )
  WCS_LIB_NAME = -lwcs
endif
ifeq ($(THIS_WCS_DIR), )
  $(info Looking for WCS installation . . .)
  THIS_WCS_LIB = $(shell find / -name libwcs.* 2>/dev/null | xargs dirname | sort -u | head -1)
  THIS_WCS_INC = $(shell find / -name wcs.h 2>/dev/null | xargs dirname | sort -u | head -1)
  THIS_WCS_DIR = $(THIS_WCS_LIB)$(THIS_WCS_INC)
else
  THIS_WCS_LIB = $(WCS_LIB)
  THIS_WCS_INC = $(WCS_INC)
  ifeq ($(shell ls $(THIS_WCS_LIB)/libwcs.*), )
    THIS_WCS_LIB = 
    $(info ERROR --- $(THIS_WCS_LIB) doesn't contain a valid wcs library file.)
  endif
  ifeq ($(shell ls $(THIS_WCS_INC)/*wcs.h), )
    THIS_WCS_INC = 
    $(info ERROR --- $(THIS_WCS_INC) doesn't contain a valid wcs.h file.)
  endif
  $(shell echo -n $(THIS_WCS_LIB) $(THIS_WCS_INC))
endif
ifeq ($(THIS_WCS_DIR), )
  $(info ERROR --- Complete wcs installation not found. Not using wcs. Change environment variable WCS_DIR to \
  a valid installation directory to use a wcs installation. Specify using, make WCS_DIR=path-to-installation , \
  WCS_INC=path-to-wcs.h from WCS_DIR, WCS_LIB=path-to-libwcs library from WCS_DIR .)
else
  $(info wcs installation found at $(THIS_WCS_LIB) & $(THIS_WCS_INC). Using wcs.)
endif

# define rules used to create BusyFunction library and test programs
all: ./librjj_objgen.a ./librjj_objgen.so ./librjj_objgen_plots.a ./librjj_objgen_plots.so ./librjj_objgen_wcs.a ./librjj_objgen_wcs.so ./create_catalog_NP ./create_catalog ./create_HUGE_catalog ./ProcessSourceOnlyCube ./Process_HUGE_SourceOnlyCube

ifneq ($(THIS_FITS_DIR), )

INCLUDES = RJJ_ObjGen.h RJJ_ObjGen_Plots.h RJJ_ObjGen_WCS.h
SOURCES = RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_MakeMask.cpp RJJ_ObjGen_Dmetric.cpp
OBJECTS = $(SOURCES:.cpp=.o)
./librjj_objgen.a: $(INCLUDES) $(SOURCES)
	@echo Creating C++ static library --- INCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -c $(SOURCES) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -lcfitsio
	$(LIB_METHOD) $@ $(OBJECTS) 

./librjj_objgen.so: $(INCLUDES) $(SOURCES)
	@echo Creating C++ shared library --- INCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -fPIC -c $(SOURCES) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -lcfitsio
	$(CXX) -shared -o $@ $(OBJECTS) -L$(THIS_FITS_LIB) 

./create_catalog_NP: ./create_catalog_NP.cpp ./librjj_objgen.a
	@echo Creating terminal application create_catalog_NP . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I$(THIS_FITS_INC) -L. -lrjj_objgen -L$(THIS_FITS_LIB) -lcfitsio 

else

INCLUDES = RJJ_ObjGen.h 
SOURCES = RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_Dmetric.cpp
OBJECTS = $(SOURCES:.cpp=.o)
./librjj_objgen.a: $(INCLUDES) $(SOURCES)
	@echo Creating C++ static library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -c $(SOURCES) -I. 
	$(LIB_METHOD) $@ $(OBJECTS)

./librjj_objgen.so: $(INCLUDES) $(SOURCES)
	@echo Creating C++ shared library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -fPIC -c $(SOURCES) -I. 
	$(CXX) -shared -o $@ $(OBJECTS)

./create_catalog_NP: ./create_catalog_NP.cpp ./librjj_objgen.a
	@echo NOT creating terminal application create_catalog_NP . . . Couldn\'t find cfitsio installation.

endif

ifneq ($(THIS_PGPLOT_DIR), )

./librjj_objgen_plots.a: RJJ_ObjGen_Plots.cpp RJJ_ObjGen_Plots.h 
	@echo Creating C++ static plotting library . . .
	@echo SOURCE = $(SOURCES) and OBJECTS = $(OBJECTS)
	$(CXX) -c RJJ_ObjGen_Plots.cpp $(CFLAGS) -I. -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L. -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -lX11
	$(LIB_METHOD) $@ RJJ_ObjGen_Plots.o 

./librjj_objgen_plots.so: RJJ_ObjGen_Plots.cpp RJJ_ObjGen_Plots.h 
	@echo Creating C++ shared plotting library . . .
	@echo SOURCE = $(SOURCES) and OBJECTS = $(OBJECTS)
	$(CXX) -fPIC -c RJJ_ObjGen_Plots.cpp $(CFLAGS) -I. -L. -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -lX11
	$(CXX) -shared -o $@ RJJ_ObjGen_Plots.o 

ifneq ($(THIS_FITS_DIR), )

./create_catalog: ./create_catalog.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application create_catalog . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) $(WCS_LIB_NAME) -L. -lrjj_objgen -lrjj_objgen_plots  -lcpgplot -lpgplot -lcfitsio -lX11 -lgfortran

./create_HUGE_catalog: ./create_HUGE_catalog.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application create_HUGE_catalog . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) $(WCS_LIB_NAME) -L. -lrjj_objgen -lrjj_objgen_plots -lcpgplot -lpgplot -lcfitsio -lX11 -lgfortran

./ProcessSourceOnlyCube: ./ProcessSourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application ProcessSourceOnlyCube . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -L. -lrjj_objgen -lrjj_objgen_plots -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) $(WCS_LIB_NAME) -lcpgplot -lpgplot -L$(THIS_FITS_LIB) -lcfitsio -lX11 -lgfortran

./Process_HUGE_SourceOnlyCube: ./Process_HUGE_SourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo Creating terminal application Process_HUGE_SourceOnlyCube . . . 
	$(CXX) $< -o $@ $(CFLAGS) -I. -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -L. -lrjj_objgen -lrjj_objgen_plots $(WCS_LIB_NAME) -lcpgplot -lpgplot -L$(THIS_FITS_LIB) -lcfitsio -lX11 -lgfortran

else

./create_catalog: ./create_catalog.cpp ./librjj_objgen.a
	@echo NOT creating terminal application create_catalog . . . Couldn\'t find cfitsio and pgplot installation.

./create_HUGE_catalog: ./create_HUGE_catalog.cpp ./librjj_objgen.a
	@echo NOT creating terminal application create_HUGE_catalog . . . Couldn\'t find cfitsio and pgplot installation.

./ProcessSourceOnlyCube: ./ProcessSourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo NOT creating terminal application ProcessSourceOnlyCube . . . Couldn\'t find cfitsio and pgplot installation.

./Process_HUGE_SourceOnlyCube: ./Process_HUGE_SourceOnlyCube.cpp ./librjj_objgen.a ./librjj_objgen_plots.a
	@echo NOT creating terminal application ProcessSource_HUGE_OnlyCube . . . Couldn\'t find cfitsio and pgplot installation.

endif

else

./librjj_objgen_plots.a: RJJ_ObjGen_Plots.cpp
	@echo NOT creating C++ plotting library . . . Couldn\'t find pgplot installation.

endif

ifneq ($(THIS_WCS_DIR), )

./librjj_objgen_wcs.a: RJJ_ObjGen_CatPrint_WCS.cpp RJJ_ObjGen_WCS.h
	@echo Creating C++ static WCS library . . .
	$(CXX) -c RJJ_ObjGen_CatPrint_WCS.cpp $(CFLAGS) -I. -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB)
	$(LIB_METHOD) $@ RJJ_ObjGen_CatPrint_WCS.o -L. -L$(THIS_WCS_LIB) $(WCS_LIB_NAME) -L$(THIS_FITS_LIB) -lcfitsio

./librjj_objgen_wcs.so: RJJ_ObjGen_CatPrint_WCS.cpp RJJ_ObjGen_WCS.h 
	@echo Creating C++ shared WCS library . . .
	$(CXX) -fPIC -c RJJ_ObjGen_CatPrint_WCS.cpp $(CFLAGS) -I. -I$(THIS_WCS_INC) -L$(THIS_WCS_LIB) -I$(THIS_FITS_INC) 
	$(CXX) -shared -o $@ RJJ_ObjGen_CatPrint_WCS.o -L. -lrjj_objgen -L$(THIS_WCS_LIB) $(WCS_LIB_NAME) -I$(THIS_FITS_INC) -L$(THIS_FITS_LIB) -lcfitsio

else

./librjj_objgen_wcs.a: RJJ_ObjGen_CatPrint_WCS.cpp RJJ_ObjGen_WCS.h
	@echo NOT creating C++ WCS library . . . Couldn\'t find WCS installation.
./librjj_objgen_wcs.so: RJJ_ObjGen_CatPrint_WCS.cpp RJJ_ObjGen_WCS.h 
	@echo NOT creating C++ shared WCS library . . . Couldn\'t find WCS installation.

endif

clean:
	@echo Removing all object files . . . 
	@rm -vf *.o

distclean:
	@echo Removing all object files, libraries and compiled programs . . . 
	@rm -vf *.o
	@rm -vf librjj_objgen.a
	@rm -vf librjj_objgen_plots.a
	@rm -vf librjj_objgen_wcs.a
	@rm -vf librjj_objgen.so
	@rm -vf librjj_objgen_plots.so
	@rm -vf librjj_objgen_wcs.so
	@rm -vf create_catalog
	@rm -vf create_HUGE_catalog
	@rm -vf create_catalog_NP
	@rm -vf ProcessSourceOnly
	@rm -vf Process_HUGE_SourceOnly

install:
	@echo Copying libraries and include files to $(INSTALL_DIR) and $(INSTALL_INC) . . .
	@cp librjj_objgen.a $(INSTALL_DIR)
	@cp librjj_objgen_plots.a $(INSTALL_DIR)
	@cp librjj_objgen_wcs.a $(INSTALL_DIR)
	@cp librjj_objgen.so $(INSTALL_DIR)
	@cp librjj_objgen_plots.so $(INSTALL_DIR)
	@cp librjj_objgen_wcs.so $(INSTALL_DIR)
	@cp RJJ_ObjGen.h $(INSTALL_INC)
	@cp RJJ_ObjGen_Plots.h $(INSTALL_INC)
	@cp RJJ_ObjGen_WCS.h $(INSTALL_INC)

