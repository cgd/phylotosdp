#!/usr/bin/env python -O

import argparse
import csv
import numbers
import os
import re

def is_num_list(xs):
    for x in xs:
        if not isinstance(x, numbers.Number):
            return False
    return True

def countLeaves(treeNode):
    if isinstance(treeNode, list):
        count = 0
        for subNode in treeNode:
            count += countLeaves(subNode)
        
        return count
    else:
        return 1

def all_tree_sdps(tree, min_count, subset_strain_indices):
    sdps = set()
    
    leaf_count = countLeaves(tree)
    def maybe_add_sdp(indexes):
        sdp = [False] * leaf_count
        for i in indexes:
            sdp[i] = True
        
        if subset_strain_indices is not None:
            sdp = [sdp[i] for i in subset_strain_indices]
        
        num_trues = 0
        for x in sdp:
            if x:
                num_trues += 1
        
        if num_trues * 2 > len(sdp) or (num_trues * 2 == len(sdp) and not sdp[0]):
            sdp = [not x for x in sdp]
            num_trues = len(sdp) - num_trues
        
        if num_trues >= min_count:
            sdps.add(tuple(sdp))
    
    def all_sdps_recursive(node):
        if isinstance(node, list):
            maybe_add_sdp(all_leaves(node))
            map(all_sdps_recursive, node)
    
    map(all_sdps_recursive, tree)
    
    return sdps

def all_leaves(tree):
    leaves = []
    for n in tree:
        if isinstance(n, list):
            leaves += all_leaves(n)
        else:
            leaves.append(n)
    return leaves

def branchLength(tree_string):
    if len(tree_string) == 0:
        raise Exception('tree_string is empty')
    
    paren_depth = 0
    i = 0
    for i in range(0, len(tree_string)):
        c = tree_string[i]
        if c == '(':
            paren_depth += 1
        elif c == ')':
            paren_depth -= 1
            if paren_depth == 0:
                break
            elif paren_depth < 0:
                raise Exception('unmatched closing paren')
        elif c.isspace():
            if paren_depth == 0:
                break
    
    if paren_depth != 0:
        raise Exception('unbalanced parens count: %i' % paren_depth)
    
    return i + 1

# TODO parsing algo is something like O(num_chars^2)
def parse_phylo_tree(tree_string):
    def parse_phylo_tree_recursive(tree_string):
        tree_string = tree_string.strip()
        
        if tree_string[0] == '(' and tree_string[len(tree_string) - 1] == ')':
            return [parse_phylo_tree_recursive(tree_string[1 : len(tree_string) - 1])]
        else:
            branches = []
            while len(tree_string) != 0:
                brLen = branchLength(tree_string)
                br = tree_string[0:brLen].strip()
                branches.append(br)
                
                tree_string = tree_string[brLen:len(tree_string)].strip()
            
            if len(branches) == 1:
                return int(branches[0])
            else:
                return [parse_phylo_tree_recursive(b) for b in branches]
    
    def clean_phylo_tree(tree):
        if(isinstance(tree, list)):
            # remove the ending zeros
            if len(tree) >= 2 \
            and all([isinstance(n, list) for n in tree[0 : len(tree) - 1]]) \
            and not isinstance(tree[len(tree) - 1], list):
                tree = tree[0 : len(tree) - 1]
            
            # remove branches that do nothing
            if len(tree) == 1 and isinstance(tree[0], list):
                return clean_phylo_tree(tree[0])
            else:
                return map(clean_phylo_tree, tree)
        else:
            return tree
    
    return clean_phylo_tree(parse_phylo_tree_recursive(tree_string))
    
def validate_tree(tree, leaf_count):
    if countLeaves(tree) != leaf_count:
        raise Exception('bad leaf count %i' % countLeaves(tree))
    
    indexSeen = [False] * leaf_count
    def validate_recursive(node):
        if is_num_list(node):
            for i in node:
                if indexSeen[i]:
                    raise Exception('error: observed leaf value twice')
                else:
                    indexSeen[i] = True
        elif len(node) != 2:
            raise Exception('expected all branches to be binary splits but found %i splits' % len(node))
        else:
            map(validate_recursive, node)
    validate_recursive(tree)
    
    if not all(indexSeen):
        raise Exception('not all leaves were seen: ' + str(indexSeen))

def sdp_to_str(sdp):
    sdp_str = ''
    for b in sdp:
        if b:
            sdp_str += '1'
        else:
            sdp_str += '0'
    return sdp_str

def get_subset_strain_indices(all_strains_file, subset_strains_file):
    # like the R match function. Returns a list of index matches for the first
    # argument in the second. If no match is found None is used for that item.
    def match(xs, ys):
        def maybe_index(x):
            try:
                return ys.index(x)
            except ValueError:
                return None
        
        return [maybe_index(x) for x in xs]
    
    all_strains = open(all_strains_file, 'r').readlines()
    subset_strains = open(subset_strains_file, 'r').readlines()
    
    matches = match(subset_strains, all_strains)
    non_matches = [subset_strains[i] for i in range(0, len(subset_strains)) if matches[i] is None]
    
    if len(non_matches) >= 1:
        raise Exception('Failed to find the following subset strains: ' + ', '.join(non_matches))
    else:
        return matches

# main entry point for script
def main():

    # parse command line arguments
    parser = argparse.ArgumentParser(description = 'calculates SDPs for all possible phylogeny tree splits')
    
    parser.add_argument(
        '--min-count-thresh',
        required = True,
        dest = 'min_count_thresh',
        help =
            'ignore a tree split pattern where the "minor allele" contains ' +
            'fewer than the given number of samples')
    parser.add_argument(
        '--all-strains',
        required = False,
        dest = 'all_strains',
        help =
            'file where each line contains a strain name. The order of ' +
            'these strains corresponds to the 0-based indices used in the ' +
            'phylogeny trees. This option should be used in conjunction with ' +
            'the subset strains option')
    parser.add_argument(
        '--subset-strains',
        required = False,
        dest = 'subset_strains',
        help =
            'file where each line contains a strain name. The ordering of strains ' +
            'in this file will correspond to the ordering of bits used in the output files')
    parser.add_argument(
        '--phylo-intervals',
        nargs = '+',
        required = True,
        dest = 'phylo_intervals',
        help = 'the input phylogeny interval files')
    parser.add_argument(
        '--unique-sdp-out',
        required = True,
        dest = 'unique_sdp_out',
        help = 'the unique SDP output file')
    parser.add_argument(
        '--sdp-interval-map-out',
        required = True,
        dest = 'sdp_interval_map_out',
        help = 'the SDP interval map output file')
    parser.add_argument(
        '--chr-capture-regex',
        required = False,
        dest = 'chr_capture_regex',
        help =
            'this argument should be a regular expression that will be matched ' +
            'against the name of each input phylogeny interval file. The first group in the match ' +
            'will be used as the chromosome string in the interval map output. ' +
            'It will be considered an error if any of the input phylogeny interval ' +
            'file names does not match this expression. An example of a regex that you ' +
            'could use for files named like "chr19maxk.csv" would be "chr(.+)maxk.csv"')
    
    args = parser.parse_args()
    
    subset_strain_indices = None
    if (args.all_strains is not None) and (args.subset_strains is not None):
        subset_strain_indices = get_subset_strain_indices(args.all_strains, args.subset_strains)
    
    chr_capture_regex = None
    if args.chr_capture_regex is not None:
        chr_capture_regex = re.compile(args.chr_capture_regex)
        if chr_capture_regex.groups != 1:
            raise Exception('--chr-capture-regex option was expected to have exactly one capture group')
    
    min_count = int(args.min_count_thresh)
    sdp_dict = dict()
    val_header = None
    for in_file in args.phylo_intervals:
        chr_name = None
        if chr_capture_regex is not None:
            match = chr_capture_regex.match(os.path.basename(in_file))
            if match:
                chr_name = match.group(1)
            else:
                raise Exception('%s does not match %s' % (os.path.basename(in_file), args.chr_capture_regex))
        
        with open(in_file, 'r') as open_in:
            csv_reader = csv.reader(open_in)
            val_header = csv_reader.next()
            val_header.pop()
            leaf_count = -1
            for row in csv_reader:
                tree = parse_phylo_tree(row[len(row) - 1])
                if leaf_count == -1:
                    leaf_count = countLeaves(tree)
                validate_tree(tree, leaf_count)
                
                for sdp in all_tree_sdps(tree, min_count, subset_strain_indices):
                    vals = []
                    if sdp in sdp_dict:
                        vals = sdp_dict[sdp]
                    else:
                        sdp_dict[sdp] = vals
                    data_rows = row[0 : len(row) - 1]
                    if chr_capture_regex is not None:
                        data_rows = [chr_name] + data_rows
                    vals.append(data_rows)
    
    with open(args.unique_sdp_out, 'wb') as unique_sdp_out:
        with open(args.sdp_interval_map_out, 'wb') as sdp_interval_map_out:
            unique_sdp_writer = csv.writer(unique_sdp_out)
            row = ['sdp']
            unique_sdp_writer.writerow(row)
            
            sdp_map_writer = csv.writer(sdp_interval_map_out)
            if chr_capture_regex is not None:
                row.append('chr')
            row += val_header
            sdp_map_writer.writerow(row)
            
            for sdp, intervals in sdp_dict.iteritems():
                sdp_str = sdp_to_str(sdp)
                unique_sdp_writer.writerow([sdp_str])
                for interval in intervals:
                    sdp_map_writer.writerow([sdp_str] + interval)

if __name__ == '__main__':
    main()
