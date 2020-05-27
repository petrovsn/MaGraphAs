# MaGraphAs
MaGraphAs. Aligner against genome graphs.

usage:

-g - path to gfa graph
-f - path to fasta-file with reads
-t - number of threads
-c - alignment treshhold
-o - output prefix

Output is a file with graph-based SNVs. 
-First column is a coordinates in the graph, 
-Second one is the alternative sequence, 
-Thierd one is the coverage of the SNV. 
