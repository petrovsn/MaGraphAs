# MaGraphAs
MaGraphAs. Aligner against genome graphs.

MaGraphAs is an aligner and SNV caller on genome graphs. It takes an assembly or variance genome graph in .gfa format[http://gfa-spec.github.io/GFA-spec/GFA1.html]  and returns a list of founded SNVs. 

installation(requires C++11 and stl library):

bash compile.sh

usage:

./MaGraphAs -g file.gfa -f file.fasta -t 100 -c 0.8 -o output

-g - path to gfa graph

-f - path to fasta-file with reads

-t - number of threads

-c - alignment treshhold[default = 0.8]

-o - output prefix

-l - hash length. Default value is 13. Highly NOT RECOMMENDED to set this value lower than 6-7, it causes to decreasing of perfomance. 

Output is a tabulated file with graph-based SNVs. 

-First column is a coordinates in the graph, 

-Second one is the alternative sequence, 

-Thierd one is the coverage of the SNV. 

/sample folder contains artificial graph and set of reads to test the program. 


Tech support email: petrovsnwm@gmail.com
