CC=gcc 
OBJS=main.o sheet.o printfiles.o mainscripts.o grid.o int_weights.o sbtmodel.o nrutil.o systemsolver.o elliptic.o ludcmp.o lubksb.o
HDRS=sheet.h hsize.h printfiles.h mainscripts.h grid.h int_weights.h sbtmodel.h nr.h nrutil.h systemsolver.h elliptic.h ludcmp.h lubksb.h
EXE=prog.exe
CFLAGS=-Wall
LDFLAGS=-lm

all: $(EXE)

$(EXE) : $(OBJS)
	$(CC) $^ $(LDFLAGS)  -o $@

$(OBJS) : %.o : %.c $(HDRS) 
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: all clean spotless

clean:
	\rm -f $(OBJS)

spotless: clean
	\rm -f $(EXE)
