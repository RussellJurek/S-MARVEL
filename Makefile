# retrieve system architecture
ARCH := $(shell uname -s)
$(info Creating S-MARVEL for ARCH = $(ARCH))

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

# determine if libtiff is installed
THIS_TIFF_DIR = $(PGPLOT_DIR)
ifneq ($(THIS_TIFF_DIR), )
  $(info TIFF_DIR = $(THIS_TIFF_DIR), checking libraries can be found here . . . )
  ifeq ($(shell ls $(THIS_TIFF_DIR)/libpgplot.a), )
    THIS_TIFF_DIR = 
  $(info ERROR --- libpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_TIFF_DIR)/libcpgplot.a), )
    THIS_TIFF_DIR = 
  $(info ERROR --- libcpgplot.a not found.)
  endif
  ifeq ($(shell ls $(THIS_TIFF_DIR)/cpgplot.h), )
    THIS_TIFF_DIR = 
  $(info ERROR --- cpgplot.h not found.)
  endif
else
  THIS_TIFF_DIR =
endif
ifeq ($(THIS_TIFF_DIR), )
  $(info ERROR --- Complete libtiff installation not found. Not using pgplot. Change environment variable TIFF_DIR to \
  a valid installation directory to use a pgplot installation. Specify using, make TIFF_DIR=path-to-installation .)
else
  $(info Pgplot installation found at $(THIS_TIFF_DIR). Using pgplot.)
endif

# define rules used to build S-MARVEL
# 1. Build ObjGen library
# 2. Build MathMorph library
# 3. Build S-MARVEL - MM removal of small-scale structure program
# 4. Build S-MARVEL - Find dots program

# define sources and includes for each library
OBJGEN_DIR = ObjGen
OBJGEN_INCLUDES = RJJ_ObjGen.h 
OBJGEN_INCLUDES = $(addprefix $(OBJGEN_DIR)/,$(OBJGEN_INCLUDES)) 
OBJGEN_SOURCES = RJJ_ObjGen_DetectDefn.cpp RJJ_ObjGen_CatPrint.cpp RJJ_ObjGen_PlotGlobal.cpp RJJ_ObjGen_CreateObjs.cpp RJJ_ObjGen_ThreshObjs.cpp RJJ_ObjGen_AddObjs.cpp RJJ_ObjGen_MemManage.cpp RJJ_ObjGen_Dmetric.cpp
OBJGEN_SOURCES = $(addprefix $(OBJGEN_DIR)/,$(OBJGEN_SOURCES)) 
OBJGEN_OBJGEN_OBJECTS = $(OBJGEN_SOURCES:.cpp=.o)

MM_DIR = MathMorph
MM_INCLUDES = MathMorph_RJJ.h 
MM_INCLUDES = $(addprefix $(MM_DIR)/,$(MM_INCLUDES)) 
MM_SOURCES = MathMorph_RJJ.cpp
MM_SOURCES = $(addprefix $(MM_DIR)/,$(MM_SOURCES)) 
MM_MM_OBJECTS = $(MM_SOURCES:.cpp=.o)

# Testing
#

# define rules used to create BusyFunction library and test programs
all: $(OBJGEN_DIR)/librjj_objgen.a $(OBJGEN_DIR)/librjj_objgen.so $(OBJGEN_DIR)/librjj_objgen_plots.a $(OBJGEN_DIR)/librjj_objgen_plots.so $(MM_DIR)/librjj_MM ./MMremove_ALL_TIF_images ./FindDots ./FindDots_wPlots

$(OBJGEN_DIR)/librjj_objgen.a: $(OBJGEN_INCLUDES) $(OBJGEN_SOURCES)
	@echo Creating ObjGen C++ static library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -c $(OBJGEN_SOURCES) -I$(OBJGEN_DIR)
	$(LIB_METHOD) $(OBJGEN_DIR)/$@ $(OBJGEN_OBJECTS)

$(OBJGEN_DIR)/librjj_objgen.so: $(OBJGEN_INCLUDES) $(OBJGEN_SOURCES)
	@echo Creating ObjGen C++ shared library --- EXCLUDING cfitsio extensions . . . 
	$(CXX) $(CFLAGS) -fPIC -c $(OBJGEN_SOURCES) -I$(OBJGEN_DIR)
	$(CXX) -shared -o $(OBJGEN_DIR)/$@ $(OBJGEN_OBJECTS)

MathMorph/librjj_MM.a: $(MM_INCLUDES) $(MM_SOURCES)
	@echo Creating MathMorph C++ static library . . . 
	$(CXX) $(CFLAGS) -c $(MM_SOURCES) -I$(MM_DIR)
	$(LIB_METHOD) $(MM_DIR)/$@ $(MM_OBJECTS)

MathMorph/librjj_MM.so: $(MM_INCLUDES) $(MM_SOURCES)
	@echo Creating MathMorph C++ dynamic library . . . 
	$(CXX) $(CFLAGS) -fPIC -c $(MM_SOURCES) -I$(MM_DIR)
	$(CXX) -shared -o $(MM_DIR)/$@ $(MM_OBJECTS)

MMremove_ALL_Tiff_images: MMremove_ALL_Tiff_images.cpp $(MM_DIR)/librjj_MM.a $(MM_DIR)/librjj_MM.so 
	@echo Creating MM removal program for TIFF images . . . 
	$(CXX) $(CFLAGS) -o $@ $< -I$(MM_DIR) -lrjj_MM

FindDots: FindDots.cpp ObjGen/librjj_objgen.a ObjGen/librjj_objgen.so
	@echo Creating FindDots program . . . 
	$(CXX) $(CFLAGS) -o $@ $< -ltiff -IObjGen -LObjGen -lrjj_objgen 

ifneq ($(THIS_PGPLOT_DIR), )

$(OBJGEN_DIR)/librjj_objgen_plots.a: ObjGen/RJJ_ObjGen_Plots.cpp ObjGen/RJJ_ObjGen_Plots.h 
	@echo Creating C++ static plotting library . . .
	@echo SOURCE = $(OBJGEN_SOURCES) and OBJGEN_OBJECTS = $(OBJGEN_OBJECTS)
	$(CXX) -c ObjGen/RJJ_ObjGen_Plots.cpp $(CFLAGS) -IObjGen -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -LObjGen -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -lX11
	$(LIB_METHOD) $@ ObjGen/RJJ_ObjGen_Plots.o 

$(OBJGEN_DIR)/librjj_objgen_plots.so: RJJ_ObjGen_Plots.cpp RJJ_ObjGen_Plots.h 
	@echo Creating C++ shared plotting library . . .
	@echo SOURCE = $(OBJGEN_SOURCES) and OBJGEN_OBJECTS = $(OBJGEN_OBJECTS)
	$(CXX) -fPIC -c RJJ_ObjGen_Plots.cpp $(CFLAGS) -IObjGen -LObjGen -I$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_INC) -L$(THIS_PGPLOT_DIR)/$(THIS_PGPLOT_LIB) -lX11
	$(CXX) -shared -o $@ ObjGen/RJJ_ObjGen_Plots.o 

FindDots_wPlots: FindDots_wPlots.cpp ObjGen/librjj_objgen.a ObjGen/librjj_objgen.so ObjGen/librjj_objgen_plots.a ObjGen/librjj_objgen_plots.so
	@echo Creating FindDots_wPlots program . . . 
	$(CXX) $(CFLAGS) -o $@ $< -ltiff -IObjGen -LObjGen -lrjj_objgen -lrjj_objgen_plots -lcpgplot -lpgplot -lX11 -lgfortran
	# Command used to build FindDots in my terminal 
	#c++ FindDots.cpp -o FindDots -O3 -I. -IObjGen -I/home/jurek83/src/Invadapodia_pipeline/GaussFit/latest/trunk -ltiff -L/home/jurek83/src/Invadapodia_pipeline/ObjGen -lrjj_objgen -lrjj_objgen_plots -lcpgplot -lpgplot -L/home/jurek83/src/Invadapodia_pipeline/GaussFit/latest/trunk -lGaussFit -lX11 -lgfortran -fopenmp -lgomp -I/usr/lib64 -lcfitsio

else

ObjGen/librjj_objgen_plots.a: RJJ_ObjGen_Plots.cpp RJJ_ObjGen_Plots.h
	@echo NOT creating static C++ plotting library . . . Couldn\'t find pgplot installation.

ObjGen/librjj_objgen_plots.so: RJJ_ObjGen_Plots.cpp RJJ_ObjGen_Plots.h
	@echo NOT creating dynamic C++ plotting library . . . Couldn\'t find pgplot installation.

FindDots_wPlots: FindDots_wPlots.cpp ObjGen/librjj_objgen_plots.a ObjGen/librjj_objgen_plots.so
	@echo NOT creating FindDots_wPlots . . . Couldn\'t find pgplot installation.

endif

clean:
	@echo Removing all object files . . . 
	@rm -vf *.o ObjGen/*.o MathMorph/*.o

distclean:
	@echo Removing all object files, libraries and compiled programs . . . 
	@rm -vf *.o ObjGen/*.o ObjGen/*.a ObjGen/*.so
	@rm -vf MathMorph/*.o MathMorph/*.a MathMorph/*.so
	@rm -vf MMremove_ALL_Tiff_images FindDots FindDots_wPlots

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


