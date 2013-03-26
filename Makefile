PROGRAM = compileme 

FILES = master.o fm.o env.o cat.o mix.o stripchannel.o
LOADARRAY = file2array.o loadarray.o

build: $(FILES)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(FILES) $(LIBS)
loadarray: $(LOADARRAY)
	$(CC) $(CFLAGS) -o loadarray $(LOADARRAY) $(LIBS)
clean:
	rm -f *.o core

rebuild: clean build
