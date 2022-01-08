This subdirectory contains standalone code for performing unidirectional pairHMM alignment, returning only the log-probability of the resulting alignments. This code was used for sampling genome-wide unidirectional pairHMM alignment log-probabilities.

### Compilation

```sh
make
```


### Usage

**Input**

A FASTA file (e.g. `tst_seqs.fa` included here) which contains a pair of nucleotide sequences.

<br/>

**Calculate the unidirectional pairHMM alignment log-probability of the input pair of sequences**

```sh
./unidirectional --scan --pair tst_seqs.fa
```

<br/>

**Output**

Return only the log-probability of the alignment.

```sh
-70.0504
```
