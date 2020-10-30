[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## TSA pairHMM: Template switch alignment pair hidden Markov model


This repository contains code for the pairHMMs presented in the manuscript/preprint:

_Short-range template switching in great ape genomes explored using a pair hidden Markov model_
by Conor R. Walker, Aylwyn Scally, Nicola De Maio, and Nick Goldman <br/>
DOI: [insert DOI]

Code in this reposistory is modified from the four-point aligner (FPA) available here: https://github.com/ariloytynoja/fpa. Compilation and
usage is the same as for the FPA.

---

### Compilation

```sh
git clone https://github.com/conorwalker/tsa_pairhmm.git
cd tsa_pairhmm
make
```

### Typical usage

**Input**

A FASTA file (e.g. `aln.fa`) which contains a pairwise alignment of two nucleotide sequences, in which the
assumed descendant sequence is the first record and the assumed ancestral sequence is the second record, for example:

```
>descendant
ACGTAA-TG...
>ancestor
ACGT-ACTG...
```

**Scan for template switch events**

Mutation clusters within `aln.fa` are defined as >=2 nucleotide (nt) differences within a 10nt sliding window.
Whenever a mutation cluster is found, a local region around this mutation cluster is re-aligned under the
unidirectional and template switch pairHMMs, following the procedure described in the Methods and Supplementary
algorithms sections of [insert DOI]. 

To identify and re-align all mutation clusters within `aln.fa` under both pairHMMs, and create a CSV file `aln_scanned.csv` which contains a header line followed by one line per re-aligned mutation cluster. we use:

```sh
./tsa_pairhmm --scan --pair aln.fa > aln_scanned.csv
```

**Output**

For each mutation cluster, the associated line in `aln_scanned.csv` contains the following information:

```
chrom                   Chromosome
clus_start_chrom        Chromosome position for the start of the mutation cluster
                        (based on information in the input FASTA headers)
clus_start_align        Pairwise alignment position for the start of the mutation cluster
clust_start1            Starting position of the cluster in the descendant sequence
clust_end1              End position of the cluster in the descendant sequence
sp1_qry                 Descendant switch point 1 position
sp1_ref                 Ancestral switch point 1 position
sp2_ref                 Ancestral switch point 2 position
sp3_ref                 Ancestral switch point 3 position
sp4_ref                 Ancestral switch point 4 position
iden_up                 Identity of the region upstream of the re-aligned mutation cluster
ident_rep               2->3 fragment identity
ident_down              Identity of the region downstream of the re-aligned mutation cluster
ident_inv               TSA pairHMM identity
ident_fwd               Unidirectional pairHMM identity
ident_epo               Input alignment identity
masked                  Boolean: is the 2->3 fragment from a masked region? (0 = not masked, 1 = masked)
sum_ins                 Sum of insertions in the re-aligned region
sum_del                 Sum of deletions in the re-aligned region
sum_mis                 Sum of mismatches in the re-aligned region
sum_nuc                 Number of unique nucleotides in the 2->3 fragment
CpG                     CpG content of re-aligned region
clus_ins                Insertions in the mutation cluster
clus_del                Deletions in the mutation cluster
clus_mis                Mismatches in the mutation cluster
fwd_score               Unidirectional pairHMM log-probability score (in last i,j)
ts_score_local          TSA pairHMM log-probability, i.e. max(M3(•,m), I3(•,m), D3(•,m)) 
ts_score_global         Global TSA pairHMM log-probability i.e. max(M3(n,m), I3(n,m), D3(n,m))
                        (not used for any calculations)
ts_ref_seq_len          Ancestral sequence (x) length in the TSA pairHMM alignment
ts_qry_seq_len          Descendant sequence (y) length in the TSA pairHMM alignment
frag_L1_size            Length of L->1
frag_23_size            Length of 2->3
frag_4R_size            Length of 4->R
```
