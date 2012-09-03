
CFLAGS=-Wall -O4
LFLAGS=-lz
CC=gcc

objects=util.o err.o memutil.o

default: $(objects) check_fastq split_fastq filter_dup_reads

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

check_fastq: $(objects) check_fastq.c
	$(CC) $(LFLAGS) -o check_fastq check_fastq.c $(objects)

split_fastq: $(objects) split_fastq.c
	$(CC) $(LFLAGS) -o split_fastq split_fastq.c $(objects)


GLIB_FLAGS=`pkg-config --cflags --libs glib-2.0`

filter_dup_reads: $(objects) filter_dup_reads.c
	$(CC) $(LFLAGS) $(GLIB_FLAGS) -o filter_dup_reads filter_dup_reads.c $(objects)

clean:
	rm -f $(objects) check_fastq split_fastq filter_dup_reads
