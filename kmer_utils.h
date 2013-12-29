// Kmer functions
void kmer_convert_kmer_to_num(char *str, unsigned long length);
unsigned long kmer_num_to_index(const char *str, const int kmer, long error_pos);
char *kmer_index_to_kmer(unsigned long long index, long kmer);
unsigned long long * kmer_counts_from_file(FILE *fh, int kmer);

// Utility functions
unsigned long long kmer_pow_four(unsigned long long x);

extern const unsigned char kmer_alpha[256];
