VERSION=\"v0.0.1\"
CC = gcc
CFLAGS =-O3 -s -mtune=native -Wall -DVERSION=$(VERSION) -Wextra
CLIBS = libkmer.a

all: kmer_total_count

kmer_total_count: kmer_utils.o kmer_total_count.o
	$(CC) $(CFLAGS) -o $@ $^

clean:
	rm -vf kmer_total_count *.o

