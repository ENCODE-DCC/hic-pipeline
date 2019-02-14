#!/usr/bin/env python

import argparse
import json
import os


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'filepath',
        type=str,
        help='filepath to the file to convert to json'
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--library-complexity',
        help='A library_complexity.txt file',
        action='store_const',
        const=jsonify_library_complexity,
        dest='jsonify'
    )
    group.add_argument(
        '--library-stats',
        help='A statistics file produced during dedup',
        action='store_const',
        const=jsonify_library_stats,
        dest='jsonify'
    )
    group.add_argument(
        '--alignment-stats',
        help='An alignment_stats.txt file generated in fragment task',
        action='store_const',
        const=jsonify_alignment_stats,
        dest='jsonify'
    )
    args = parser.parse_args()
    return args


def parse_to_dict(filepath):
    """
    Convert a text file of key: value pairs to a dictionary, preprocessing the
    field names such that "Key Name" becomes "key_name" in the resulting
    dictionary
    """
    keys = []
    values = []
    with open(filepath) as f:
        data = f.readlines()
    for line in data:
        split_line = line.split(': ')
        key = '_'.join(split_line[0].strip().split(' '))
        keys.append(key.lower())
        values.append(split_line[1].strip())
    output = dict(zip(keys, values))
    return output


def jsonify_library_complexity(filepath):
    """
    An example libary_complexity.txt file looks like the following:
        Unique Reads: 12,101 
        PCR Duplicates: 120 
        Optical Duplicates: 0 
        Library Complexity Estimate: 618,223
    """
    jsonified = parse_to_dict(filepath)
    return jsonified


def jsonify_library_stats(filepath):
    """
    See https://github.com/aidenlab/juicer/blob/master/CPU/common/statistics.pl
    for more details on the meaning of the fields. An example stats.txt file
    looks like the following:
        Intra-fragment Reads: 6,969(57.59%)
        Hi-C Contacts: 5,132(42.41%)
        Ligation Motif Present: 3 (0.02%)
        3' Bias (Long Range): 65% - 35%
        Pair Type %(L-I-O-R): 25% - 23% - 27% - 25%
        Inter-chromosomal: 6 (0.05%)
        Intra-chromosomal: 5,126 (42.36%)
        Short Range (<20Kb): 4,537 (37.49%)
        Long Range (>20Kb): 589 (4.87%)
    """
    jsonified = parse_to_dict(filepath)
    value_sep = '-'
    for k, v in jsonified.items():
        if '%' in v and value_sep not in v:
            split_value = v.split('(')
            values = {
                'total': split_value[0].strip(),
                'percent': split_value[1].rstrip('%)')
            }
        else:
            split_value = v.split(value_sep)
            split_value = [i.strip() for i in split_value]
            if 'bias' in k:
                fields =['long_range', 'short_range']
            elif 'pair' in k:
                fields =['left', 'inner', 'outer', 'right']
            values = dict(zip(fields, split_value))
        jsonified[k] = values
    return jsonified


def jsonify_alignment_stats(filepath):
    """
    An example alignment_stats.txt file looks like the following:
        Sequenced Read Pairs:  332888
         Normal Paired: 297469 (89.36%)
         Chimeric Paired: 16 (0.00%)
         Collisions: 1 (0.00%)
         Low MAPQ Collisions: 2 (0.00%)
         Unmapped: 11303 (3.40%)
         MAPQ 0: 24097 (7.24%)
         Ligation Motif Present: 145 (0.04%)
        Alignable (Normal+Chimeric Paired): 297485 (89.36%)
    """
    jsonified = parse_to_dict(filepath)
    for k, v in jsonified.items():
        if '(' in v:
            split_value = v.split('(')
            values = {
                'total': split_value[0].strip(),
                'percent': split_value[1].rstrip('%)')
            }
        else:
            values = v.strip()
        jsonified[k] = values
    return jsonified


def main():
    args = parse_arguments()
    filepath = args.filepath
    jsonified = args.jsonify(filepath)
    filename = os.path.splitext(filepath)[0] + '.json'
    with open(filename, 'w') as f:
        f.write(json.dumps(jsonified, indent=4, sort_keys=True))


if __name__== '__main__':
    main()
