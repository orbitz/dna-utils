// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "kmer_total_count.h"

#define ERROR_CODE 5
#define SPACE_CODE 6

#define BUFFER_SIZE 256
const unsigned char kmer_alpha[256] =

//                             \n
{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

 //                            >        A     C           G
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
  return 1ULL << (x * 2);
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
size_t translate_nucleotides_to_numbers(char *str, size_t len, const unsigned char *lookup) {
  size_t i = 0;

  for(i = 0; i < len; ++i) {
    if(str[i] == '>') {
      return i;
    }
    else {
      str[i] = lookup[(int)str[i]];
    }
  }

  return i;
}

static int is_whitespace(char c) {
  return c == SPACE_CODE;
}

static int is_error_char(char c) {
  return c == ERROR_CODE;
}

static size_t calculate_mer(const char *str, size_t *pos, size_t kmer_len, size_t error_mer) {
  size_t mer = 0;
  size_t multiply = 1;
  size_t i;

  /* printf("-----\n"); */
  // for each char in the k-mer check if it is an error char
  for(i = *pos; i < *pos + kmer_len; ++i) {
    if(is_whitespace(str[i])) {
      continue;
    }

    if(is_error_char(str[i])) {
      mer = error_mer;
      *pos = i;
      return mer;
    }

    // multiply this char in the mer by the multiply
    // and bitshift the multiply for the next round
    /* printf("str[i] = %d\n", str[i]); */
    /* printf("i = %zu\n", i); */
    mer += str[i] * multiply;
    multiply = multiply << 2;
    /* printf("mer - %zu\n", mer); */
  }

  return mer;
}

size_t consume_header(const char *str, size_t len) {
  size_t i;

  for(i = 0; i < len; ++i) {
    if(str[i] == '\n') {
      ++i;
      return i;
    }
  }

  return len;
}

unsigned long long *kmer_counts_from_file(FILE *fh, const unsigned int kmer) {
  char buffer[BUFFER_SIZE] = { 0 };
  char *start = buffer;
  size_t read_len = BUFFER_SIZE;
  bool in_header = false;

  size_t len = 0;

  const unsigned long width = kmer_pow_four(kmer);

  // calloc our return array
  unsigned long long *counts = calloc(width + 1, sizeof *counts);

  if(counts == NULL)  {
    fprintf(stderr, strerror(errno));
    exit(EXIT_FAILURE);
  }

  while((len = fread(start, 1, read_len, fh)) != 0) {
    size_t total_len = (start + len) - buffer;
    size_t processed = 0;

    printf("total_len = %zu\n", total_len);
    while(processed != total_len) {
      printf("processed = %zu\n", processed);
      if(in_header) {
        processed += consume_header(buffer + processed, total_len - processed);
        in_header = (processed == total_len);
      }

      in_header = in_header || buffer[processed] == '>';

      if(!in_header) {
        size_t i;
        size_t to_process = processed;

        to_process += translate_nucleotides_to_numbers(buffer + processed, total_len, kmer_alpha);

        for(i = processed; i < (to_process - kmer + 1); ++i) {
          size_t mer = calculate_mer(buffer, &i, kmer, width);
          counts[mer]++;
        }

        processed = to_process;
        in_header = (buffer[processed] == '>');
      }
    }

    memcpy(buffer, buffer + processed, processed);
    start = buffer + (BUFFER_SIZE - processed);
    read_len = BUFFER_SIZE - (start - buffer);
  }

  return counts;
}
