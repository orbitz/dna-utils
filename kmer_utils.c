// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "kmer_total_count.h"

#define ERROR 5
#define SPACE 6

const unsigned char kmer_alpha[256] =

//                             \n
{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

 //                                     A     C           G
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5,

 //                T                                      A     C           G
 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 5, 1, 5, 5, 5, 2,

 //                                  T
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};

static const char reverse_alpha[4] = { 'A', 'C', 'G', 'T' };

unsigned long long kmer_pow_four(unsigned long long x) {
  return (unsigned long long)1 << (x * 2);
}

// convert an index back into a kmer string
char *kmer_index_to_kmer(unsigned long long index, long kmer)  {

  int i = 0;
  int j = 0;
  char *num_array = calloc(kmer,  sizeof(char));
  char *ret = calloc(kmer + 1, sizeof(char));

  if(ret == NULL)
    exit(EXIT_FAILURE);

  // this is the core of the conversion. modulus 4 for base 4 conversion
  while (index != 0) {
    num_array[i] = index % 4;
    index /= 4;
    i++;
  }

  // for our first few nmers, like AAAAA, the output would only be "A" instead
  // of AAAAA so we prepend it
  for(j = 0; j < (kmer - i); j++)
    ret[j] = 'A';

  // our offset for how many chars we prepended
  int offset = j;
  // save i so we can print it
  int start = i ;

  // decrement i by 1 to reverse the last i++
  i--;
  j = 0;

  // reverse the array, as j increases, decrease i
  for(j = 0; j < start; j++, i--)
    ret[j + offset] = reverse_alpha[(int)num_array[i]];

  // set our last character to the null termination byte
  ret[kmer + 1] = '\0';

  free(num_array);
  return ret;
}

// convert a string into values from a lookup array
static void translate_nucleotides_to_numbers(char *str, size_t len, const unsigned char *lookup) {
  size_t i;

  for(i = 0; i < len; ++i)
    str[i] = lookup[(int)str[i]];
}

static size_t calculate_mer(const char *str, size_t *pos, size_t kmer_len, size_t error_mer) {
  size_t mer = 0;
  size_t multiply;
  size_t i;
  size_t end;

  // set our multiplier to be pow(4, mer) >> 2;
  multiply = (unsigned long long)1 << (kmer_len* 2);
  multiply = multiply >> 2;

  // for each char in the k-mer check if it is an error char
  end = *pos + kmer_len;
  for(i = *pos; i < end; ++i) {
    if(str[i] == SPACE) {
      end++;
      continue;
    }
    if(str[i] == ERROR) {
      *pos = i;
      return error_mer;
    }

    // multiply this char in the mer by the multiply
    // and bitshift the multiply for the next round
    mer += str[i] * multiply;
    multiply = multiply >> 2;
  }

  return mer;
}

unsigned long long *kmer_counts_from_file(FILE *fh, const unsigned int kmer) {
  char *line = NULL;
  char *seq = NULL;
  size_t len = 0;
  ssize_t read;

  size_t position = 0;

  // width is 4^kmer
  // there's a sneaky bitshift to avoid pow dependency
  const unsigned long width = kmer_pow_four(kmer);

  // malloc our return array
  unsigned long long *counts = calloc(width + 1, sizeof *counts);
  if(counts == NULL)  {
    fprintf(stderr, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while ((read = getdelim(&line, &len, '>', fh)) != -1) {

    // find our first \n, this should be the end of the header
    char *start = strchr(line, '\n');
    if(start == NULL)
      continue;

    // point to one past that.
    seq = start + 1;

    size_t seq_length = read - (seq - line);

    // relace A, C, G and T with 0, 1, 2, 3 respectively
    // unknowns are 5, newlines are 6
    translate_nucleotides_to_numbers(seq, seq_length, kmer_alpha);

    // loop through our string to process each k-mer
    for(position = 0; position < (seq_length - kmer + 1); position++) {
      size_t mer = calculate_mer(seq, &position, kmer, width);
      counts[mer]++;
    }
  }

  free(line);

  return counts;
}
