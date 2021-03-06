CXX = $(shell root-config --cxx)
LD = $(shell root-config --ld)

INC = $(shell pwd)
Repo = $(shell git rev-parse --show-toplevel)
Aux = $(Repo)/Aux
HggRazor = $(Repo)/HggRazor
CommonTools = $(HggRazor)/CommonTools

CPPFLAGS := $(shell root-config --cflags) -I$(INC)/include -I$(CommonTools)/include -I$(Aux)/include
LDFLAGS := $(shell root-config --glibs) $(STDLIBDIR) -lRooFit -lRooFitCore -L$(INC)/include -lCustomPdfs

CPPFLAGS += -g -std=c++11

TARGET = MakeFit
TARGET1 = MakeFitRun2
TARGET2 = PlotAllBiasTest
TARGET4 = CreateDataCardCommands
TARGET5 = testCustomPdf

SRC = app/testapp.cc src/makeInvertedANFit_v2.cc $(Aux)/src/CommandLineInput.cc
SRC1 = app/MakeFitRun2.cc src/RunTwoFitMgg.cc src/DefinePdfs.cc src/PrintAIC.cc $(Aux)/src/CommandLineInput.cc
SRC2 = app/PlotAllBiasTest.cc src/PlotBias.cc $(Aux)/src/CommandLineInput.cc $(CommonTools)/src/HggAux.cc
SRC4 = app/CreateDataCardCommands.cc $(Aux)/src/CommandLineInput.cc
SRC5 = app/testCustomPdf.cc

OBJ = $(SRC:.cc=.o)
OBJ1 = $(SRC1:.cc=.o)
OBJ2 = $(SRC2:.cc=.o)
OBJ4 = $(SRC4:.cc=.o)
OBJ5 = $(SRC5:.cc=.o)

all : $(TARGET) $(TARGET1) $(TARGET2) $(TARGET4) $(TARGET5)

$(TARGET) : $(OBJ)
	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^


$(TARGET1) : $(OBJ1)
	$(LD) $(CPPFLAGS) -o $(TARGET1) $(OBJ1) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET2) : $(OBJ2)
	$(LD) $(CPPFLAGS) -o $(TARGET2) $(OBJ2) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^


$(TARGET4) : $(OBJ4)
	$(LD) $(CPPFLAGS) -o $(TARGET4) $(OBJ4) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

$(TARGET5) : $(OBJ5)
	$(LD) $(CPPFLAGS) -o $(TARGET5) $(OBJ5) $(LDFLAGS)
	@echo $@
	@echo $<
	@echo $^

%.o : %.cc
	$(CXX) $(CPPFLAGS) -o $@ -c $<
	@echo $@
	@echo $<


clean :
	rm -f *.o app/*.o src/*.o $(CommonTools)/src/*.o $(Aux)/src/*.o $(TARGET) $(TARGET1) $(TARGET2) $(TARGET4) $(TARGET5) *~
