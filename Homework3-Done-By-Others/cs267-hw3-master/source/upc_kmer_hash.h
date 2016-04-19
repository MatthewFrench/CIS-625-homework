#ifndef UPC_KMER_HASH_H
#define UPC_KMER_HASH_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h> 
#include <string.h>
#include <upc.h>
#include "upc_contig_generation.h"

/* Auxiliary function for computing hash values */
int64_t hashseq(int64_t  hashtable_size, char *seq, int size)
{
   unsigned long hashval;
   hashval = 5381;
   for(int i = 0; i < size; i++) {
      hashval = seq[i] +  (hashval << 5) + hashval;
   }
   
   return hashval % hashtable_size;
}

/* Returns the hash value of a kmer */
int64_t hashkmer(int64_t  hashtable_size, char *seq)
{
   return hashseq(hashtable_size, seq, KMER_PACKED_LENGTH);
}

/* Looks up a kmer in the hash table and returns a pointer to that entry */
char lookup_kmer_r_ext(shared kmer_t* memory_heap, shared int64_t* next_index, 
                       shared int64_t* hash_table, int64_t tablesize, 
                       const unsigned char *kmer)
{
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
   
   // start traverse
   int64_t pos = hash_table[hashval];
   kmer_t tmp_kmer; // temporary storage
   for (; pos != -1; ) {
      // getch the global kmer to local memory
      upc_memget(&tmp_kmer, &memory_heap[pos], sizeof(kmer_t));
      
      // compare
      if ( memcmp(packedKmer, tmp_kmer.kmer, KMER_PACKED_LENGTH * sizeof(char)) == 0 ) {
         return tmp_kmer.r_ext;
      }
      
      pos = next_index[pos];
   }
   
   // ERROR! Should never reach here!!
   return 0;
}

/* Adds a kmer and its extensions in the hash table (note that a memory heap should be preallocated. ) */
int add_kmer(shared int64_t* next_index, int64_t k, 
             shared int64_t* hash_table, int tablesize, 
             shared kmer_t * memory_heap, 
             const unsigned char *kmer, char left_ext, char right_ext)
{
   /* Pack a k-mer sequence appropriately */
   char packedKmer[KMER_PACKED_LENGTH];
   packSequence(kmer, (unsigned char*) packedKmer, KMER_LENGTH);
   int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
   
   kmer_t tmp_kmer;// temp storage
   
   /* Add the contents to the appropriate kmer struct in the heap */
   memcpy(tmp_kmer.kmer, packedKmer, KMER_PACKED_LENGTH * sizeof(char));
   tmp_kmer.l_ext = left_ext;
   tmp_kmer.r_ext = right_ext;
   upc_memput(&memory_heap[k], &tmp_kmer, sizeof(kmer_t));
   
   /* insert k to the end of entry_list at hashtable[hashval] */
   int64_t ptr;
   ptr = bupc_atomicI64_cswap_strict((shared void*)(&hash_table[hashval]), -1, k);
   while(ptr != -1) {
     ptr = bupc_atomicI64_cswap_strict((shared void*)(&next_index[ptr]), -1, k);
   }
   
   return 0;
}

#endif // KMER_HASH_H
