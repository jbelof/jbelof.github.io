
CC=gcc
CFLAGS=-O3
DEFINES=-DDEBUG

.c.o:
	$(CC) -c $(CFLAGS) $(DEFINES) -I. $*.c

all:	astumian_game.o
	$(CC) $(CFLAGS) $(DEFINES) *.o -o ag

clean:
	rm -f *.o ag


