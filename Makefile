EXEC = h5readout
CC = g++
OBJ = $(EXEC).o
# spdaq22
# DAQPATH=/usr/opt/daq/experimental/11.3-018
# DDASPATH=/usr/opt/ddas/3.4-002
# HDF5PATH=./hdf5-104
# CXXOPTPATH=./cxxopt_dd45a08/include
#
DAQPATH=/opt/nscldaq
DDASPATH=/opt/ddas
HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
CXXOPTPATH=./cxxopts

INC = -I$(DAQPATH)/include
INC += -I$(DDASPATH)/include
INC += -I$(HDF5PATH)/include
INC += -I$(CXXOPTPATH)/include
LIBS = -lhdf5 -lhdf5_cpp -L$(HDF5PATH)/lib
LIBS += -ldataformat -ldaqio -ldaqshm -lPortManager -lurl -lFragmentIndex -L$(DAQPATH)/lib
LIBS += -lddasformat -L$(DDASPATH)/lib
LIBS += -lException -L$(DAQPATH)/lib
LDFLAGS = -Wl,-rpath $(DAQPATH)/lib
LDFLAGS += -Wl,-rpath $(DDASPATH)/lib
LDFLAGS += -Wl,-rpath $(HDF5PATH)/lib
CXXFLAGS = -O2 -std=c++11

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $(OBJ) $(LIBS) $(LDFLAGS)
%.o: %.cpp
	$(CC) -c -o $@ $< $(CXXFLAGS) $(INC)

install:
	cp -a $(EXEC) /usr/local/bin

uninstall:
	rm -f /usr/local/bin/h5readout

testhdf: cleanout
	./h5readout -i run-0226-00.evt

test2:
	./h5readout -i run-8262-00.evt

test2-gdb:
	gdb --args ./h5readout -i run-8262-00.evt

repack:
	h5repack -v -f SHUF -f GZIP=9 run-0235-00.evt.h5 run-0235-00.evt1.h5

test1: test.cpp
	$(CC) $< -o $@ $(CXXFLAGS) -I/usr/include/hdf5/serial -lhdf5 -lhdf5_cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial

.PHONY: clean cleanout

cleanout:
	/bin/rm -rf run-0226-00.evt.h5

clean:
	rm -f $(OBJ)

run235:
	./h5readout -i ../../nscl_data/run-0235-00.evt
run226:
	./h5readout -i ../../nscl_data/run-0226-00.evt
run227:
	./h5readout -i ../../nscl_data/run-0227-00.evt
