# MaGraphAs
MaGraphAs. Aligner against genome graphs.


usage:

./MaGraphAs -g file.gfa -f file.fasta -t 100 -c 0.8 -o output

-g - path to gfa graph
-f - path to fasta-file with reads
-t - number of threads
-c - alignment treshhold[default = 0.8]
-o - output prefix

Output is a file with graph-based SNVs. 
-First column is a coordinates in the graph, 
-Second one is the alternative sequence, 
-Thierd one is the coverage of the SNV. 
