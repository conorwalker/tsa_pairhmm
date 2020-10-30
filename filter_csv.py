#!/usr/bin/env python3

"""
Copyright (C) 2020 EMBL - European Bioinformatics Institute
Contact: goldman@ebi.ac.uk, cwalker@ebi.ac.uk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License,
or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""

from sys import argv


def main():
    """
    Filter a TSA pairHMM output file.
    """

    # threshold for statistical significance
    logprob_thresh = 9

    # per base log probability threshold, determined using random
    # samples of Human/Chimp alignments
    per_base_logprob_thresh = -0.148574

    # input CSV
    scan_fi = argv[1]

    # output CSV
    filter_fi = "aln_events.csv"

    # store lines which meet criteria
    filtered_lines = []

    with open(scan_fi, "r") as scnfi:
        for line in scnfi:
            pairhmm_output = line.strip().split(",")
            # ignore header/footer lines
            if pairhmm_output[0] == "chrom" or pairhmm_output[0] == "scan finished":
                continue

            # unidirectional pairHMM log-probability
            uni_logprob = float(pairhmm_output[25])

            # TSA pairHMM log-probability
            tsa_logprob = float(pairhmm_output[26])

            # log of the ratio between the pairHMM probabilities
            log_probability_ratio = tsa_logprob - uni_logprob

            # length of the TSA pairHMM alignment
            tsa_align_len = float(pairhmm_output[29])
            # 14.78259 = (M1->M2) + (M2->M3)
            ts_per_base_logprob = (tsa_logprob + 14.78259) / tsa_align_len

            # sum of unique nucleotides in the 2->3 fragment
            sum_nuc = int(pairhmm_output[20])

            # mask boolean
            masked = int(pairhmm_output[16])

            # number of deletions in the focal mutation cluster
            clus_del = int(pairhmm_output[23])

            # assess event line
            if log_probability_ratio >= logprob_thresh and \
               ts_per_base_logprob >= per_base_logprob_thresh and \
               masked == 0 and \
               sum_nuc == 4 and \
               clus_del <= 50:
                filtered_lines.append(line)

    with open(filter_fi, "w") as out_fi:
        for line in filtered_lines:
            out_fi.write(line)


if __name__ == "__main__":
    main()
