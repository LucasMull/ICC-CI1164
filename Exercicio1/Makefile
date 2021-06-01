# Lucas Müller GRR20197160

CC ?= gcc

CFLAGS := -Wall -g
LDFLAGS := -lm

all: run

run: analise_intervalar.c
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	@rm run

.PHONY: run clean
