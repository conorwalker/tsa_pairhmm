[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## TSA pairHMM: Template switch alignment pair hidden Markov model


This repository contains code for the pairHMMs presented in the manuscript/preprint:

_Short-range template switching in great ape genomes explored using a pair hidden Markov model_ </br>
by Conor R. Walker, Aylwyn Scally, Nicola De Maio, and Nick Goldman </br>
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

A FASTA file (e.g. `aln.fa` included here) which contains a pairwise alignment of two nucleotide sequences, in which the
assumed descendant sequence is the first record and the assumed ancestral sequence is the second record, for example:

```
>descendant
ATCTGGAA...
>ancestor
ATCTGGAA...
```

<br/>

**Scan for template switch events**

Mutation clusters within `aln.fa` are defined as >=2 nucleotide (nt) differences within a 10nt sliding window.
Whenever a mutation cluster is found, a local region around this mutation cluster is re-aligned under the
unidirectional and template switch pairHMMs, following the procedure described in the Methods and Supplementary
algorithms sections of [insert DOI]. 

To identify and re-align all mutation clusters within `aln.fa` under both pairHMMs, we use:

```sh
./tsa_pairhmm --scan --pair aln.fa > aln_scanned.csv
```
<br/>


**Output**

The above command creates a CSV file `aln_scanned.csv` which contains a header line followed by one line per re-aligned mutation cluster as follows:

```sh
head aln_scanned.csv

chrom,clus_start_chrom,clus_start_align,clust_start1,clust_end1,sp1_qry,sp1_ref,sp2_ref,sp3_ref,sp4_ref,iden_up,ident_rep,ident_down,ident_inv,ident_fwd,ident_epo,masked,sum_ins,sum_del,sum_mis,sum_nuc,CpG,clus_ins,clus_del,clus_mis,fwd_score,ts_score_local,ts_score_global,ts_ref_seq_len,ts_qry_seq_len,frag_L1_size,frag_23_size,frag_4R_size
10,123199418,601,601,603,602,602,740,740,605,1,1,1,1,0.988,0.976,0,0,1,0,1,0,0,1,1,-11,-15.4,-27.7,282,81,40,1,39
10,123201804,2988,2987,2994,2960,2961,3126,3126,2963,1,1,0.95,0.96,0.94,0.94,0,0,0,1,1,0,0,0,2,-18.7,-27.4,-39.4,282,82,9,1,40
10,123201987,3171,3170,3173,3171,3172,3302,3301,3175,1,1,1,1,0.975,0.976,3,0,0,2,2,0,0,0,2,-12.7,-15.5,-27.8,282,82,40,2,39
10,123202034,3218,3217,3221,3218,3219,3336,3334,3223,1,1,1,1,0.976,0.976,3,0,0,2,2,0,0,0,2,-12.7,-15.6,-28,282,82,40,3,39
10,123202123,3307,3306,3314,3307,3308,3452,3452,3310,1,1,0.975,0.987,0.974,0.974,3,0,0,1,1,0,0,0,2,-12.7,-21.4,-33.4,282,82,36,1,40
10,123202787,3971,3970,3972,3966,3967,4025,4025,3970,1,1,1,1,0.987,0.962,3,0,1,0,1,0,1,2,0,-11,-15.4,-27.7,282,81,36,1,40
10,123202852,4038,4035,4041,4041,4043,4137,4123,4146,1,0.933,1,0.987,0.473,0.511,3,0,87,-1,2,0,0,86,1,-25.6,-22.9,-35.2,371,84,40,15,23
10,123202871,4143,4054,4056,4041,4003,4137,4123,4146,1,0.933,1,0.982,0.671,0.947,3,0,27,-1,2,0,0,1,1,-22.3,-22.2,-34.5,282,55,1,15,39
10,123203007,4280,4190,4193,4186,4275,4413,4413,4275,1,1,1,1,0.974,0.974,3,2,0,0,1,0,2,0,0,-13.8,-15.4,-27.7,280,82,36,1,40
```

For each mutation cluster, the associated line in `aln_scanned.csv` contains the following information:

```
chrom                   Chromosome
clus_start_chrom        Chromosome position for the start of the mutation cluster
                        (based on information in the input FASTA headers)
clus_start_align        Pairwise alignment position for the start of the mutation cluster
clust_start1            Starting position of the cluster in the descendant sequence
clust_end1              End position of the cluster in the descendant sequence
sp1_qry                 Descendant switch point ① position
sp1_ref                 Ancestral switch point ① position
sp2_ref                 Ancestral switch point ② position
sp3_ref                 Ancestral switch point ③ position
sp4_ref                 Ancestral switch point ④ position
iden_up                 Identity of the region upstream of the re-aligned mutation cluster
ident_rep               ②→③ fragment identity
ident_down              Identity of the region downstream of the re-aligned mutation cluster
ident_inv               TSA pairHMM identity
ident_fwd               Unidirectional pairHMM identity
ident_epo               Input alignment identity
masked                  Boolean: is the ②→③ fragment from a masked region?
                        (0 = not masked, 1 = masked)
sum_ins                 Sum of insertions in the re-aligned region
sum_del                 Sum of deletions in the re-aligned region
sum_mis                 Sum of mismatches in the re-aligned region
sum_nuc                 Number of unique nucleotides in the ②→③ fragment
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
frag_L1_size            Length of Ⓛ→①
frag_23_size            Length of ②→③
frag_4R_size            Length of ④→Ⓡ
```

<br/>

**Assessing candidate events for significance and alignment quality**

Many of the mutation clusters re-aligned and output to `aln_scanned.csv` will not have arisen due to a template
switch event. Using our approach outlined in [insert DOI], we filter each candidate event line based on:

1. Statistical significance
2. A threshold on the per-base template switch alignment quality 
3. No masked sequence in the ②→③ fragment
4. All four unique nucleotides within the ②→③ fragment
5. A threshold on the maximum number of deletions in the mutation cluster

An example script `filter_csv.py` is included here to filter based on this criteria:

```sh
python filter_csv.py aln_scanned.csv
```

This will generate one file `aln_events.csv`, which contains any candidate events from `aln_scanned.csv` that pass our thresholds (one event in this case):
```sh
cat aln_events.csv

10,123237102,40311,38285,38298,38286,30379,30396,30385,30389,1,1,1,1,0.79,0.902,0,12,9,0,4,2,3,0,6,-38.8,-16.8,-29.1,289,92,40,12,39
```
<br/>


**Visualsing template switch alignments**

The event output contained within the above CSV files can be visualised using:
```sh
./tsa_pairhmm --pair aln.fa --print-file aln_events.csv

chr10:123237104-123237115

Template switch process:
F1: L CCAAACTCCATTTTTACTGACCATGCTAACACACACCCAAA 1
F3:                                                   4 GAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC R
RF:   CCAAACTCCATTTTTACTGACCATGCTAACACACACCCAAAAAATCCAAAGAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC
RR:   GGTTTGAGGTAAAAATGACTGGTACGATTGTGTGTGGGTTTTTTAGGTTTCTCGTATTAAGAACATGCAATAACGGTGTATGTCCTACAG
F2:                                               3 GTTTCTCGTATT 2


Unidirectional alignment (log-probability: -38.9)
 CAAACTCCATTTTTACTGACCATGCTAACACACAC---------CCAAAttatgctctttgGAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC
 CAAACTCCATTTTTACTGACCATGCTAACACACACCCAAAAAATCCAAA------------GAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC

Template switch alignment (log-probability: -16.8)
 CAAACTCCATTTTTACTGACCATGCTAACACACACCCAAA|TTATGCTCTTTG|GAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC
 CAAACTCCATTTTTACTGACCATGCTAACACACACCCAAA|TTATGCTCTTTG|GAGCATAATTCTTGTACGTTATTGCCACATACAGGATGTC
```
