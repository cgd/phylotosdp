==USAGE==

usage: phylotosdp.py [-h] --min-count-thresh MIN_COUNT_THRESH
                     [--all-strains ALL_STRAINS]
                     [--subset-strains SUBSET_STRAINS] --phylo-intervals
                     PHYLO_INTERVALS [PHYLO_INTERVALS ...] --unique-sdp-out
                     UNIQUE_SDP_OUT --sdp-interval-map-out
                     SDP_INTERVAL_MAP_OUT

calculates SDPs for all possible phylogeny tree splits

optional arguments:
  -h, --help            show this help message and exit
  --min-count-thresh MIN_COUNT_THRESH
                        ignore a tree split pattern where the "minor allele"
                        contains fewer than the given number of samples
  --all-strains ALL_STRAINS
                        file where each line contains a strain name. The order
                        of these strains corresponds to the 0-based indices
                        used in the phylogeny trees. This option should be
                        used in conjunction with the subset strains option
  --subset-strains SUBSET_STRAINS
                        file where each line contains a strain name. The
                        ordering of strains in this file will correspond to
                        the ordering of bits used in the output files
  --phylo-intervals PHYLO_INTERVALS [PHYLO_INTERVALS ...]
                        the input phylogeny interval files
  --unique-sdp-out UNIQUE_SDP_OUT
                        the unique SDP output file
  --sdp-interval-map-out SDP_INTERVAL_MAP_OUT
                        the SDP interval map output file

== PHYLOGENY INPUT FILES ==

The only requred input file is one or more phylogeny interval files. These
files will be formatted like:

middle_start,middle_end,start,end,op_start,op_end,sdp_count,tree
3125547,3160277,3125547,3261140,3125547,3262864,7,((0  1  2  3  4) (5  6  7) 0)
3262865,3301232,3180931,3337115,3160278,3337222,9,((0  1  2  3) (4  5  6  7) 0)

The data contained in all but the last column can be interval data or anything
else. It is simply copied over to the output in the SDP interval map output
file.

The phylogeny tree will be in the final column. Trees are represented by a
recursive structure of parentheses which will contain either subgroups or a list
of numbers separated by spaces. All tree splits are expected to be binary. Also
note that tree splits can be followed by a number (0 in the example given
above). This number is ignored for the purposes of constructing SDPs.

These files are similar to newick but not exactly the same. In particular the
zeros which follow groups of strain indices are ignored for the purposes of
constructing a phylogeny tree. Also all tree branches are expected to be binary.

== STRAIN NAME INPUT FILES ==

If you want to subset the output you will need to provide a file containing all
strain names represented in the phylogeny trees and another file that contains
the subset of strains that you would like to retain for the purpose of
generating SDPs. Both of these files have a simple format where there is
no header, and each strain gets its own line.

== EXAMPLE INVOCATION ==

The following shows how you might invoke this script when you are interested in
subsetting the strains:

./phylotosdp.py \
    --all-strains strains.txt --subset-strains strain-subset.txt \
    --phylo-intervals chr*maxk.csv --min-count-thresh 4 \
    --unique-sdp-out unique-sdps.txt --sdp-interval-map-out sdp-map.csv

