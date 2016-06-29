# SparseAssembler

Exploiting sparseness in de novo genome assembly

http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-S6-S1

Basic command for the software:

./SparseAssembler g 10 k 51 LD 0 GS 200000000 NodeCovTh 1 EdgeCovTh 0 f frag_1.fastq f frag_2.fastq f frag_3.fastq  &

For memory usages, results and comparisons:

https://sites.google.com/site/sparseassembler/home/results-comparisons


Parameters:

k: kmer size, support 15~127.

g: number of skipped intermediate k-mers, support 1-25.

f: single end FASTA/FASTQ reads. Multiple inputs shall be independently imported with this parameter.

GS: genome size estimation in bp (used for memory pre-allocation), suggest a large value if possible.(e.g. ~ 3x genome size)

NodeCovTh: coverage threshold for spurious k-mers, support 0-16. (default 1)

EdgeCovTh: coverage threshold for spurious links, support 0-16. (default 0)

LD: load a saved k-mer graph.

PathCovTh: coverage threshold for spurious paths in the breadth-first search, support 0-100.

