CC=mpicc
CFLAGS=-Wall -Wextra -Wconversion

main: main.c main.h
	$(CC) -o $@ $@.c $(CFLAGS)

	