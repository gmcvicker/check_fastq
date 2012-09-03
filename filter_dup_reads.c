#include <zlib.h>
#include <stdio.h>
#include <glib.h>

#include "util.h"
#include "memutil.h"

#define MAX_LINE 1024



int main(int argc, char **argv) {
  gzFile gzf;
  char *input_filename, *seq_str, *count_ptr;
  char line1[MAX_LINE], line2[MAX_LINE], line3[MAX_LINE], line4[MAX_LINE];
  GHashTable *hashtab;
  long *val, count;

  if(argc != 2) {
	fprintf(stderr, "usage: %s my_fastq_file.txt.gz | "
			"gzip > filtered_fastq.txt.gz\n", argv[0]);
	exit(2);
  }

  input_filename = argv[1];

  gzf = util_must_gzopen(input_filename, "rb");

  hashtab = g_hash_table_new_full(g_str_hash, g_str_equal,
								  free, free);

  count = 0;
  while(TRUE) {
	/* read all 4 fastq lines at once time */
	if((gzgets(gzf, line1, MAX_LINE) == NULL) ||
	   (gzgets(gzf, line2, MAX_LINE) == NULL) ||
	   (gzgets(gzf, line3, MAX_LINE) == NULL) ||
	   (gzgets(gzf, line4, MAX_LINE) == NULL)) {
	  break;
	}
	
	count += 1;
	if(count == 100000) {
	  count = 0;
	  fprintf(stderr, ".");
	}


	/* the second line is the sequence, we want to check if we've seen
	 * the same one before
	 */
	seq_str = &line2[1];
	
	val = g_hash_table_lookup(hashtab, seq_str);
	if(val == NULL) {
	  /* we have not seen a read with this sequence before */

	  /* record the read in the hash table */
	  val = my_new(long, 1);
	  *val = 1;
	  g_hash_table_insert(hashtab, util_str_dup(seq_str), val);
	  
	  /* print the read to stdout in fastq format */
	  fprintf(stdout, "%s%s%s%s", line1, line2, line3, line4);
	} else {
	  /* we have seen this sequence before, increment count and discard */
	  *val += 1;
	  /* fprintf(stderr, "discarding duplicate read (seen %d times) %s \n",
	   * seq_str, *val);
	   */
	}
  }  

  fprintf(stderr, "\n");

  gzclose(gzf);
  g_hash_table_destroy(hashtab);

  return 0;
}
