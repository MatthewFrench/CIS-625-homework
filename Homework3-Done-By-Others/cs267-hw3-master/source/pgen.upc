#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <upc.h>
#include <upc_io.h>

#include "packingDNAseq.h"
#include "upc_kmer_hash.h"


//#define SINGLE_OUTPUT_FILE


int main(int argc, char *argv[]){

	/** Declarations **/
	double inputTime=0.0, constrTime=0.0, traversalTime=0.0;

	/** Read input **/
	upc_barrier;
	inputTime -= gettime();
	///////////////////////////////////////////
	// Your code for input file reading here //
	///////////////////////////////////////////
	
	
	char * input_file_name = argv[1];
	int64_t nKmers = getNumKmersInUFX(input_file_name);
	int64_t avg_n = nKmers / THREADS;
	int64_t rem_n = nKmers % THREADS;
	int64_t lines_to_read = avg_n + (MYTHREAD < rem_n);
	int64_t lines_to_ignore = avg_n * MYTHREAD;
	if(MYTHREAD <= rem_n) lines_to_ignore += MYTHREAD;
	else lines_to_ignore += rem_n;
	int64_t chars_to_read = lines_to_read * LINE_SIZE;
	int64_t chars_to_ignore = lines_to_ignore * LINE_SIZE;
	unsigned char* buffer = 
	  (unsigned char*) malloc((chars_to_read+5) * sizeof(unsigned char)); // local buffer
	
	upc_file_t *input_file;
	input_file = upc_all_fopen(input_file_name, UPC_RDONLY | UPC_INDIVIDUAL_FP, 0, NULL);
	upc_all_fseek(input_file, chars_to_ignore*sizeof(unsigned char), UPC_SEEK_SET);
	int64_t cur_chars_read = upc_all_fread_local(input_file, buffer, sizeof(unsigned char), chars_to_read,
	                                             UPC_IN_ALLSYNC | UPC_OUT_ALLSYNC);
	// All the data are stored in  * buffer *
	
	// close the open file
	upc_barrier;
	upc_all_fclose(input_file);
	
	
	//printf("Reading Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);
	
	
	upc_barrier;
	inputTime += gettime();

	/** Graph construction **/
	constrTime -= gettime();
	///////////////////////////////////////////
	// Your code for graph construction here //
	///////////////////////////////////////////
	
	//int ok=1;
	
	/* Initialize lookup table that will be used for the DNA packing routines */
   init_LookupTable();
	
	
	int64_t n = lines_to_read; // total kmers processed by the current thread
	// generate hash table and heap
	shared int64_t* next_index = (shared int64_t*)upc_all_alloc(nKmers, sizeof(int64_t));
  shared kmer_t* memory_heap = (shared kmer_t*)upc_all_alloc(nKmers, sizeof(kmer_t)); 
  
  // initialize hash_table
  int64_t tablesize = nKmers * LOAD_FACTOR;
  shared int64_t* hash_table = (shared int64_t*)upc_all_alloc(tablesize, sizeof(int64_t));
  int64_t i;
  upc_forall(i=0; i<tablesize; i++; &hash_table[i]) // initial hash table
    hash_table[i] = -1;
  upc_forall(i=0; i<nKmers; i++; &next_index[i]) // initial next index
    next_index[i] = -1;
  
  upc_barrier; // need to synchronize before we actually start!
  
  int64_t k = lines_to_ignore; // global kmer index
  int64_t ptr = 0;
  
  
//  i = 0; // test purpose
  
  
	while (ptr < cur_chars_read) {
	
    char left_ext = (char) buffer[ptr+KMER_LENGTH+1];
    char right_ext = (char) buffer[ptr+KMER_LENGTH+2];

    /* Add k-mer to hash table */
    add_kmer(next_index, k, hash_table, tablesize, memory_heap, &buffer[ptr], left_ext, right_ext);
    
    /*
    //// test whether it is inserted to the hash_table
    {
      char packedKmer[KMER_PACKED_LENGTH];
      packSequence(&buffer[ptr], (unsigned char*) packedKmer, KMER_LENGTH);
      int64_t hashval = hashkmer(tablesize, (char*) packedKmer);
      int64_t p = hash_table[hashval];
      while(p != k && p != -1) {
        p = next_index[p];
      }
      if(p != k) {
        printf("<ADD Failure!> Cannot find the kmer index [%d] in the hash_table on Thread#%d!\n", k, MYTHREAD);
      }
      
      if(ok==1) {
      
      kmer_t tmp_kmer;
      upc_memget(& tmp_kmer, &memory_heap[k], sizeof(kmer_t));
      
      
      char unpackedKmer[KMER_LENGTH+1];
      unpackSequence((unsigned char*) tmp_kmer.kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
      if(memcmp(unpackedKmer, &buffer[ptr], KMER_LENGTH * sizeof(char)) != 0) {
        printf("<Fetch ERROR> the fetched kmer cannot be correctly recovered on Thread#%d!\n", MYTHREAD);
      }

      }
    }
    
    */
    /* Move to the next k-mer in the input buffer */
    ptr += LINE_SIZE;
    k ++; 
    
    //i++;
  }
	
	
	//if(n != i) printf("<ERROR> on Thread #%d!\n >>> Expected Read %d kmers\n >>> Actual Read %d kmers\n", MYTHREAD, n, i);
	
	upc_barrier;
	constrTime += gettime();

  //printf("Building Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);

	/** Graph traversal **/
	traversalTime -= gettime();
	////////////////////////////////////////////////////////////
	// Your code for graph traversal and output printing here //
	// Save your output to "pgen.out"                         //
	////////////////////////////////////////////////////////////
	
	
	//ok=1;
	
	
	
	char cur_contig[MAXIMUM_CONTIG_SIZE+1], unpackedKmer[KMER_LENGTH+1]; 
	
	char output_file_name[50];
	sprintf(output_file_name, "pgen%d.out", MYTHREAD);
	FILE *out_file = fopen(output_file_name, "w"); // asynchronized & independent output
	
	i = 0; ptr = 0;
	for (; i < n; i ++, ptr += LINE_SIZE) {
	  char left_ext = (char) buffer[ptr+KMER_LENGTH+1];
	  if (left_ext != 'F') continue;
	  
	  //if(MYTHREAD == 0) {
	  //  printf("start find sequence, unpacking ... ");
	 // }
	  
	  // Start to find a sequence
//    unpackSequence((unsigned char*) cur_local_kmer->kmer,  (unsigned char*) unpackedKmer, KMER_LENGTH);
    
    //if(MYTHREAD == 0) {
	  //  printf("OK!\n");
	  //}
    
    
	  //memcpy(cur_contig, unpackedKmer, KMER_LENGTH * sizeof(char));
	  memcpy(cur_contig, &buffer[ptr], KMER_LENGTH * sizeof(char));
	  
	  int64_t posInContig = KMER_LENGTH;
    char right_ext = (char) buffer[ptr+KMER_LENGTH+2];
    
    /*
    if(ok==1) {
    // testing
    right_ext = lookup_kmer_r_ext(memory_heap, next_index, hash_table, tablesize,
                                 (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
    if(right_ext == 0) {
      printf("<ERROR at BEGIN> ca-nnot find the first kmer! on Thread#%d!\n", MYTHREAD);
      ok=0;
      continue;
    }
    }
    */
    
    // Keep traversing the graph
    while(right_ext != 'F' && right_ext != 0) {
      cur_contig[posInContig ++] = right_ext;
      
      // hash table lookup
      right_ext = lookup_kmer_r_ext(memory_heap, next_index, hash_table, tablesize,
                                   (const unsigned char *) &cur_contig[posInContig-KMER_LENGTH]);
    }
    
    if(right_ext == 0) {
      printf("<ERROR> cannot find the next kmer! on Thread#%d!\n", MYTHREAD);
      continue;
    }
    
    // print the contig
    cur_contig[posInContig] = '\0';
    fprintf(out_file, "%s\n", cur_contig);
	}
	
	// close the output file
	fclose(out_file);
	
	
	// clean the allocated memory
	upc_barrier;
	
	
  //printf("Traversal Finished on Thread#%d of %d threads.\n", MYTHREAD, THREADS);
	
	
	if(MYTHREAD == 0) {
	  upc_free(memory_heap);
	  upc_free(next_index);
	  upc_free(hash_table);
	}
	
	
	upc_barrier;
	traversalTime += gettime();

// hack implementation to produce a single output file for correctness checking
#ifdef SINGLE_OUTPUT_FILE

  if(MYTHREAD == 0) {
    out_file = fopen("pgen.out", "w");
    
    for(int t = 0; t < THREADS; ++ t) {
      char str[50];
      sprintf(str, "pgen%d.out", t);
      FILE*in_file = fopen(str, "r");
      
      while(fscanf(in_file, "%s", cur_contig) == 1)
        fprintf(out_file, "%s\n", cur_contig);
       
      fclose(in_file);
    }
    
    fclose(out_file);
  }

#endif


	/** Print timing and output info **/
	/***** DO NOT CHANGE THIS PART ****/
	if(MYTHREAD==0){
		printf("%s: Input set: %s\n", argv[0], argv[1]);
		printf("Number of UPC threads: %d\n", THREADS);
		printf("Input reading time: %f seconds\n", inputTime);
		printf("Graph construction time: %f seconds\n", constrTime);
		printf("Graph traversal time: %f seconds\n", traversalTime);
	}
	return 0;
}
