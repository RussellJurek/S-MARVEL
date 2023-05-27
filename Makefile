# retrieve system architecture
ARCH := $(shell uname -s)
CURR_DIR := $(shell pwd)
$(info Creating S-MARVEL for ARCH = $(ARCH) at base directory $(CURR_DIR))

# specify compilation options
CXX = c++
CC = cc
CFLAGS = -O3 
CXXFLAGS = -O3
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
  DYLIB_TYPE = dylib
else 
  LIB_METHOD = ar -crus
  DYLIB_TYPE = so
endif
$(info Using static library build mehod, LIB_METHOD = $(LIB_METHOD))

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

# define rules used to build S-MARVEL
# 1. Build ObjGen library
# 2. Build MathMorph library
# 3. Build S-MARVEL - MM removal of small-scale structure program
# 4. Build S-MARVEL - Find dots program

# define sources and includes for each library
OBJGEN_DIR = $(CURR_DIR)/src/ObjGen
OBJGEN_INCLUDES = $(addprefix $(OBJGEN_DIR)/,RJJ_Objgen.h) 
OBJGEN_SOURCES = $(addprefix $(OBJGEN_DIR)/,RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_Dmetric.cpp)
OBJGEN_OBJECTS = $(OBJGEN_SOURCES:.cpp=.o)

MM_DIR = $(CURR_DIR)/src/MathMorph
MM_INCLUDES = $(addprefix $(MM_DIR)/,MathMorph_RJJ.h)
MM_SOURCES = $(addprefix $(MM_DIR)/,MathMorph_RJJ.cpp)
MM_OBJECTS = $(MM_SOURCES:.cpp=.o)

# define rules used to create BusyFunction library and test programs
all: libs programs

$(OBJGEN_DIR)/%.o: $(OBJGEN_DIR)/%.cpp 
	$(CXX) -fPIC $(CXXFLAGS) -I$(OBJGEN_DIR) -c $< -o $@ 

$(OBJGEN_DIR)/librjj_objgen.a: $(OBJGEN_OBJECTS) 
	@echo Creating ObjGen C++ static library --- EXCLUDING cfitsio extensions . . . 
	$(LIB_METHOD) $@ $(OBJGEN_OBJECTS)

$(OBJGEN_DIR)/librjj_objgen.so: $(OBJGEN_SOURCES)
	@echo Creating ObjGen C++ shared library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) -shared -o $@ $(OBJGEN_OBJECTS)

$(MM_DIR)/%.o: $(MM_DIR)/%.cpp
	$(CXX) -fPIC $(CXXFLAGS) -I$(MM_DIR) -c $< -o $@ 

$(MM_DIR)/librjj_MM.a: $(MM_OBJECTS)
	@echo Creating MathMorph C++ static library . . . 
	$(LIB_METHOD) $@ $(MM_OBJECTS)

$(MM_DIR)/librjj_MM.so: $(MM_OBJECTS)
	@echo Creating MathMorph C++ dynamic library . . . 
	$(CXX) -shared -o $@ $(MM_OBJECTS)

bin/MMremove_ALL_Tiff_images: src/MMremove_ALL_Tiff_images.cpp $(MM_DIR)/librjj_MM.a $(MM_DIR)/librjj_MM.so 
	@echo Creating MM removal program for TIFF images . . . 
	$(CXX) $(CXXFLAGS) -o $@ $< -I$(MM_DIR) -L$(MM_DIR) -lrjj_MM -ltiff

bin/FindDots: src/FindDots.cpp $(OBJGEN_DIR)/librjj_objgen.a $(OBJGEN_DIR)/librjj_objgen.so
	@echo Creating FindDots program . . . 
	$(CXX) $(CXXFLAGS) -o $@ $< -ltiff -I$(OBJGEN_DIR) -L$(OBJGEN_DIR) -lrjj_objgen 

ifneq ($(THIS_PGPLOT_DIR), )

$(OBJGEN_DIR)/RJJ_ObjGen_Plots.o: $(OBJGEN_DIR)/RJJ_ObjGen_Plots.cpp
	$(CXX) -c $(OBJGEN_DIR)/RJJ_ObjGen_Plots.cpp -o $(OBJGEN_DIR)/RJJ_ObjGen_Plots.o -fPIC $(CXXFLAGS) -I$(OBJGEN_DIR) -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(OBJGEN_DIR) -lrjj_objgen -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -lX11

$(OBJGEN_DIR)/librjj_objgen_plots.a: src/ObjGen/RJJ_ObjGen_Plots.o
	@echo Creating C++ static plotting library . . .
	$(LIB_METHOD) $@ $(OBJGEN_DIR)/RJJ_ObjGen_Plots.o 

$(OBJGEN_DIR)/librjj_objgen_plots.so: src/ObjGen/RJJ_ObjGen_Plots.o
	@echo Creating C++ shared plotting library . . .
	$(CXX) -shared -o $@ $(OBJGEN_DIR)/RJJ_ObjGen_Plots.o 

bin/FindDots_wPlots: src/FindDots_wPlots.cpp $(OBJGEN_DIR)/librjj_objgen.a $(OBJGEN_DIR)/librjj_objgen.so $(OBJGEN_DIR)/librjj_objgen_plots.a $(OBJGEN_DIR)/librjj_objgen_plots.so
	@echo Creating FindDots_wPlots program . . . 
	$(CXX) $(CXXFLAGS) -o $@ $< -ltiff -I$(OBJGEN_DIR) -L$(OBJGEN_DIR) -lrjj_objgen -lrjj_objgen_plots -lcpgplot -lpgplot -lX11 -lgfortran

else

$(OBJGEN_DIR)/librjj_objgen_plots.a: $(OBJGEN_DIR)/RJJ_ObjGen_Plots.cpp 
	@echo NOT creating static C++ plotting library . . . Couldn\'t find pgplot installation.

$(OBJGEN_DIR)/librjj_objgen_plots.so: $(OBJGEN_DIR)/RJJ_ObjGen_Plots.cpp 
	@echo NOT creating dynamic C++ plotting library . . . Couldn\'t find pgplot installation.

bin/FindDots_wPlots: src/FindDots_wPlots.cpp $(OBJGEN_DIR)/librjj_objgen_plots.a $(OBJGEN_DIR)/librjj_objgen_plots.so
	@echo NOT creating FindDots_wPlots . . . Couldn\'t find pgplot installation.

endif

libs: $(OBJGEN_DIR)/librjj_objgen.a $(OBJGEN_DIR)/librjj_objgen.so $(OBJGEN_DIR)/librjj_objgen_plots.a $(OBJGEN_DIR)/librjj_objgen_plots.so $(MM_DIR)/librjj_MM.a $(MM_DIR)/librjj_MM.so 

programs: bin/MMremove_ALL_Tiff_images bin/FindDots bin/FindDots_wPlots

clean:
	@echo Removing all object files . . . 
	@rm -vf src/ObjGen/*.o src/MathMorph/*.o

distclean:
	@echo Removing all object files, libraries and compiled programs . . . 
	@rm -vf ObjGen/*.o ObjGen/*.a ObjGen/*.so
	@rm -vf MathMorph/*.o MathMorph/*.a MathMorph/*.so
	@rm -vf MMremove_ALL_Tiff_images FindDots FindDots_wPlots

binclean:
	@echo Removing binaries . . .
	@rm -vf bin/FindDots bin/FindDots_wPlots bin/MMremove_ALL_Tiff_images

install:
	@echo Copying libraries and include files to $(INSTALL_DIR) and $(INSTALL_INC) . . .
	@cp librjj_objgen.a $(INSTALL_DIR)
	@cp librjj_objgen_plots.a $(INSTALL_DIR)
	@cp librjj_objgen.so $(INSTALL_DIR)
	@cp librjj_objgen_plots.so $(INSTALL_DIR)
	@cp RJJ_ObjGen.h $(INSTALL_INC)
	@cp RJJ_ObjGen_Plots.h $(INSTALL_INC)
	@cp librjj_MM.a $(INSTALL_DIR)
	@cp librjj_MM.so $(INSTALL_DIR)
	@cp MathMorph_RJJ.h $(INSTALL_INC)


