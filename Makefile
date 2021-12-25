CC=mpicc
CFLAGS=-Wall -Wextra -Wconversion

%:
	$(CC) -o $@ $@.c $(CFLAGS)
	