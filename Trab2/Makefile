CC     = gcc -g -std=c11
CFLAGS = -Wall -D_NO_LIKWID
LFLAGS = -lm

PROG = matrixInv 
OBJS = matrixLib.o utils.o

.PHONY: limpa faxina clean purge all

all: $(PROG)

%.o: %.c matrixLib.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(CFLAGS) $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

