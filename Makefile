EXEC = h5readout
CC = g++
OBJ = h5readout.o main.o
#
#spdaq22: Debian 8
# DAQPATH=/usr/opt/daq/experimental/11.3-018
# DDASPATH=/usr/opt/ddas/3.4-002
# HDF5PATH=./hdf5-104
#
#flagtail: Debian 10
# DAQPATH=/usr/opt/opt-buster/daq/11.3-018
# DDASPATH=/usr/opt/opt-buster/ddas/3.4-002
# HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial
#
#my vmphy0: Debian 10
DAQPATH=/usr/lib/nscldaq
DDASPATH=/usr/lib/opt/ddas
HDF5PATH=/usr/lib/x86_64-linux-gnu/hdf5/serial

#CXXOPTPATH=./cxxopts_dd45a08/include
CXXOPTPATH=./argh_d9964d4

INC = -I$(DAQPATH)/include
INC += -I$(DDASPATH)/include
INC += -I$(HDF5PATH)/include
INC += -I$(CXXOPTPATH)
LIBS = -lhdf5 -lhdf5_cpp -L$(HDF5PATH)/lib
LIBS += -lDataFlow -ldaqthreads
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

pack:
	tar cjvf h5readout-app_1.1.tar.bz2 h5readout hdf5-104

test3:
	./h5readout -i test/CCF-data/run-8262-00.evt
