This subdirectory contains standalone code for performing unidirectional pairHMM alignment, returning only the log-probability of the resulting alignments. This code was used for sampling genome-wide unidirectional pairHMM alignment log-probabilities.

### Compilation

```sh
make
```


### Usage

**Input**

A FASTA file (e.g. `tst_seqs.fa` included here) which contains a pair of aligned nucleotide sequences.

<br/>

**Calculate the unidirectional pairHMM alignment log-probability of the input pair of aligned sequences**

To calculate the log-probability of the whole alignment:
```sh
./unidirectional --scan --pair tst_seqs.fa
-56.9838
```

To calculate the log-probability of the whole alignment divided by the alignment length (the "per-base log-probability" in S8 Fig):
```sh
./unidirectional --scan --pair tst_seqs.fa
-0.518035
```
