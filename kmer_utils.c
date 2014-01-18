// Copyright 2013 Calvin Morrison
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#include "kmer_total_count.h"

#define ERROR 5
#define SPACE 6

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
bool translate_nucleotides_to_numbers(char *str, size_t len, const unsigned char *lookup, bool start_header) {

  bool header = start_header;
  size_t i;

  for(i = 0; i < len; ++i) {
    if(str[i] == '>') {
      header = true;
    } 
    if(header) {
      if(str[i] == '\n') {
        header = false;
      }
      str[i] = ERROR;
    } 
    else 
      str[i] = lookup[(int)str[i]];
  }

  return header;
}

static size_t calculate_mer(const char *str, size_t str_len, size_t *pos, size_t kmer_len, size_t error_mer) {
  size_t mer = 0;
  size_t multiply;
  size_t i;
  size_t end;

  // set our multiplier to be pow(4, mer) >> 2;
  multiply = kmer_pow_four(kmer_len) >> 2;

  end = *pos + kmer_len;

  // check if space or error, otherwise bump our multiplier and mer 
  for(i = *pos; i < end; ++i) {
    if(str[i] == SPACE) {
      // if our end has been expanded to the max string length
      // or if it is the first loop return the error length
      if(end == str_len || i == *pos) 
        return error_mer;

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

static size_t fread_save_n_bytes(char *buffer, FILE *fh, size_t save_size, size_t buffer_size, size_t len) {

  if(ftell(fh) == 0)  {
    fread(buffer, 1, buffer_size, fh);
    len = buffer_size - save_size;
  }
  else {

    size_t read_size = buffer_size - save_size;

    memcpy(buffer, buffer + len, save_size);

    len = fread(buffer + save_size, 1, read_size, fh);

  }

  return len;
}

unsigned long long *kmer_counts_from_file(FILE *fh, const unsigned int kmer) {
  char buffer[BUFFER_SIZE] = { 0 };
  bool header = false;
  bool started = false;

  size_t len = 0;

  const unsigned long width = kmer_pow_four(kmer);

  // malloc our return array
  unsigned long long *counts = calloc(width + 1, sizeof *counts);
  if(counts == NULL)  {
    fprintf(stderr, strerror(errno));
    exit(EXIT_FAILURE);
  }

  size_t save_size = kmer - 1;
  while((len = fread_save_n_bytes(buffer, fh, save_size, BUFFER_SIZE - 1, len)) != 0) {
    size_t i;
    char *ptr = buffer + save_size;
    size_t ptr_len = len;

    if(!started) {
      ptr_len = len + save_size;
      ptr = buffer;
    }

    header = translate_nucleotides_to_numbers(ptr, ptr_len, kmer_alpha, header);

    for(i = 0; i < (len - kmer + 1); i++) {
      size_t mer = calculate_mer(buffer, len, &i, kmer, width);
      counts[mer]++;
    }

    started = true;
  }

  return counts;
}
