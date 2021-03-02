EXEC = h5readout
CC = g++
OBJ = $(EXEC).o
INC = -I/opt/nscldaq/include
INC += -I/opt/ddas/include
INC += -I/usr/include/hdf5/serial
INC += -I../cs/cpp/cxxopts/include
LIBS = -lhdf5 -lhdf5_cpp -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LIBS += -ldataformat -ldaqio -ldaqshm -lPortManager -lurl -lFragmentIndex -L/opt/nscldaq/lib
LIBS += -lddasformat -L/opt/ddas/lib
LIBS += -lException -L/opt/libtclplus/lib
LDFLAGS = -Wl,-rpath /opt/nscldaq/lib
LDFLAGS += -Wl,-rpath /opt/ddas/lib
LDFLAGS += -Wl,-rpath /opt/libtclplus/lib
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
