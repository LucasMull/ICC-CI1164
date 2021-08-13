PROG  = matrixInv

CC   = gcc -std=c11 -g
OBJS = matrixLib.o

ifeq ($(build),debug)
	CFLAGS := -D_NO_LIKWID
	LFLAGS :=
else
	CFLAGS := -DLIKWID_PERFMON -I${LIKWID_INCLUDE} -O3 -mavx -march=native
	LFLAGS := -L${LIKWID_LIB} -llikwid -lm
endif

.PHONY: all debug clean limpa purge faxina

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

all: $(PRINT) $(PROG)

debug: CFLAGS += -DDEBUG
debug: $(PROG)

$(PROG):  $(PROG).o
$(PRINT): 
	$(CC) $(CFLAGS) -o $@ $@.c

$(PROG): $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

clean limpa:
	@echo "Limpando ...."
	@rm -f *~ *.bak *.tmp

purge faxina:   clean
	@echo "Faxina ...."
	@rm -f  $(PROG) $(PRINT) *.o core a.out
	@rm -f *.png marker.out *.log
