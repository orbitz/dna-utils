// Kmer functions
void convert_kmer_to_num(char *str, const unsigned long length);
unsigned long num_to_index(const char *str, const int kmer, const long error_pos);
char *index_to_kmer(unsigned long long index, long kmer);
unsigned long long * get_kmer_counts_from_file(const char *fn, const int kmer);

// Utility functions
char *strnstrip(const char *s, char *dest, int c, int len);
inline unsigned long long pow_four(unsigned long long x);

// Variables
const unsigned char alpha[256]; 
