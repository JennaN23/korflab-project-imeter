# Project IMEter
## Program Descriptions
### at_get_seq.py
Filters data using a minimum IME score threshold. Prints in FASTA format.
### at_processing.py
  Prints the ID, entropy, and log2 expression score, as well as the categorized intron beginning, length of the intron, and IME scores for *A. thaliana*.
### avg_int_length.py
  Prints maximum expression, median intron length, average intron length, nucleotide frequencies, and "bins" showing the distribution of log2 expression scores.
### dm_get_100_seq.py
  Filters data using a minimum IME threshold. Prints in FASTA format.
### dm_processing.py
  Prints the ID and log2 expression score, as well as the categorized intron beginning, length of the intron, and IME scores for *D. melanogaster*.
### cluster.py
  Clusters a dataset and presents it as a graph.
### get_seq.py
   Filters data using minimum length, maximum length, maximum distance, and minimum expression as thresholds. Prints in FASTA format.
### ime_master_trainer.py
  Prints each kmer and its IME score. Includes arguments for k-mer size, donor length, acceptor length, proximal cutoff, distal cutoff, maximum intron length, and minimum expression.
### ind_tissue_xp.py
  Prints the average expression for each of the 11 *A. thaliana* tissues in at_ime_tissues.txt.gz, as well as how many introns are found in all tissues, a restricted number of tissues, or are tissue specific.
