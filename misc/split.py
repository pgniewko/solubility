#! /usr/bin/env python


import sys

split_string = sys.argv[2]

with open(sys.argv[1], 'r') as fin:
    for line in fin:
        pairs = line.rstrip('\n').split(f'"{split_string}"')
        for entry in pairs:
            if len(entry) > 0:
                print(entry.replace(' ', ''))
