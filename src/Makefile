.PHONY: clean install check-syntax
CC	= gcc
#CFLAGS	= -Wall -O3 -fPIC #-fopenmp
CFLAGS	= -g
LDFLAGS	= -lm
INSTALL	= install
PREFIX	= ${HOME}

all: md2 mkconf granular1d

md2: md2.c
	$(CC) $(CFLAGS) -o md2 md2.c $(LDFLAGS)

mkconf: mkconf.c
	$(CC) $(CFLAGS) -o mkconf mkconf.c $(LDFLAGS)

granular1d: granular1d.c
	$(CC) $(CFLAGS) -o granular1d granular1d.c $(LDFLAGS)

install:
	$(INSTALL) -m 755 md2 ${PREFIX}/bin
	$(INSTALL) -m 755 mkconf ${PREFIX}/bin
	$(INSTALL) -m 755 granular1d ${PREFIX}/bin

clean:
	rm -f md2 mkconf granular1d *.o nul *~

check-syntax:
	$(CC) -o nul -S ${CHK_SOURCES}
