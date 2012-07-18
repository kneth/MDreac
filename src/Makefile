.PHONY: clean install check-syntax
CC	= gcc
CFLAGS	= -Wall -O3 -fPIC -fopenmp
LDFLAGS	= -lm
INSTALL	= install
PREFIX	= ${HOME}

all: md2 mkconf

md2: md2.c
	$(CC) $(CFLAGS) -o md2 md2.c $(LDFLAGS)

mkconf: mkconf.c
	$(CC) $(CFLAGS) -o mkconf mkconf.c $(LDFLAGS)

install:
	$(INSTALL) -m 755 md2 ${PREFIX}/bin
	$(INSTALL) -m 755 mkconf ${PREFIX}/bin

clean:
	rm -f md2 mkconf *.o nul *~

check-syntax:
	$(CC) -o nul -S ${CHK_SOURCES}
