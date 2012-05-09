
CFLAGS=-Wall -O4
LFLAGS=-lz
CC=gcc

objects=util.o err.o memutil.o

default: $(objects) check_fastq split_fastq

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

check_fastq: $(objects) check_fastq.c
	$(CC) $(LFLAGS) -o check_fastq check_fastq.c $(objects)

split_fastq: $(objects) split_fastq.c
	$(CC) $(LFLAGS) -o split_fastq split_fastq.c $(objects)

clean:
	rm -f $(objects) check_fastq
