// Kmer functions
void convert_kmer_to_num(char *str, unsigned long length);
unsigned long num_to_index(const char *str, const int kmer, long error_pos);
char *index_to_kmer(unsigned long long index, long kmer);
unsigned long long * get_kmer_counts_from_file(FILE *fh, int kmer);

// Utility functions
char *strnstrip(const char *s, char *dest, int c, int len);
unsigned long long pow_four(unsigned long long x);

extern const unsigned char alpha[256];
