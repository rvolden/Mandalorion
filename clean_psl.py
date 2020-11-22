#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import numpy as np
import sys

def parse_contigs(psl_file, clean_psl_file):
    '''
    Adjusts some of the psl values
    '''

    out = open(clean_psl_file, 'w')
    for line in open(psl_file):
        a = line.strip().split('\t')
        start = int(a[15])
        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]
        readstarts = a[19].split(',')[:-1]

        blockstarts_clean = []
        readstarts_clean = []
        blocksizes_clean = []

        size_gap = []

        for x in range(0, len(blocksizes), 1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            readstart = int(readstarts[x])
            blockend = blockstart + blocksize

            size_gap.append(blocksize)

            try:
                next_blockstart = int(blockstarts[x+1])
                next_blocksize = int(blocksizes[x+1])
                next_readstart = int(readstarts[x+1])
                next_blockend = next_blockstart + next_blocksize

                gap = next_blockstart - blockend
                size_gap.append(gap)
            except:
                pass

        new_size_gap = []
        block = 0
        for index in range(0, len(size_gap), 1):
            if index%2 == 0:
                block += size_gap[index]
            if index%2 == 1:
                if size_gap[index] < 20:
                    block += size_gap[index]
                else:
                    new_size_gap.append(block)
                    new_size_gap.append(size_gap[index])
                    block = 0

        new_size_gap.append(block)

        current_position = start
        current_readposition = int(a[11])
        for index in range(0, len(new_size_gap), 1):
            inter = new_size_gap[index]
            if index%2 == 0:
                blockstarts_clean.append(str(current_position))
                blocksizes_clean.append(str(inter))
                readstarts_clean.append(str(current_readposition))

                current_position += inter
                current_readposition += inter

            if index%2 == 1:
                current_position += inter

        a[17] = str(len(blockstarts_clean))

        blockstarts_clean.append('')
        blocksizes_clean.append('')
        readstarts_clean.append('')

        a[18] = (',').join(blocksizes_clean)
        a[19] = (',').join(readstarts_clean)
        a[20] = (',').join(blockstarts_clean)

        new_line = ('\t').join(a)
        out.write(new_line + '\n')

def main():
    psl_file = sys.argv[1]
    clean_psl_file = sys.argv[2]
    parse_contigs(psl_file, clean_psl_file)

main()
