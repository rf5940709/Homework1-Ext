cc ?= gcc
CFLAGS ?= -O0 -Wall -std=c99 -mavx
LIB ?= -lm
EXEC = compute_pi

all: $(EXEC)

astyle:
	astyle --style=kr --indent=spaces=4 --indent-switches --suffix=none *.[ch]

compute_pi: compute_pi.c
	$(CC) $(CFLAGS) -o $@ $@.c $(LIB)
clean:
	$(RM) $(EXEC)
