#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import os
import numpy as np
import argparse
import multiprocessing as mp


def revComp(seq):
    '''Returns the reverse complement of a seq'''
    bases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-'}
    return ''.join([bases[x] for x in list(seq)])[::-1]


def argParser():
    '''Parses command line arguments'''
    parser = argparse.ArgumentParser(
        description='Makes consensus sequences from R2C2 reads.',
        add_help=True, prefix_chars='-')
    parser.add_argument(
        '--path', '-p', type=str, action='store', default=os.getcwd(),
        help='Directory where all the files are/where they will end up.\
              Defaults to your current directory.')
    parser.add_argument('--subsample', '-s', type=int, action='store')
    parser.add_argument('--numThreads', '-n', type=int, action='store')
    parser.add_argument(
        '--config', '-c', type=str, action='store', default='',
        help='If you want to use a config file to specify paths to\
              programs, specify them here. Use for poa, racon, water,\
              blat, and minimap2 if they are not in your path.')
    return vars(parser.parse_args())


def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, racon, consensus, blat, and emtrey
    possible = set(['minimap2', 'consensus', 'racon',
                    'blat', 'emtrey', 'poa', 'medaka'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
    # check for missing programs
    # if missing, default to path
    for missing in possible - inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs


args = argParser()
path = args['path']
temp_folder = path + '/mp'
subsample = args['subsample']
numThreads = args['numThreads']

if args['config']:
    progs = configReader(args['config'])
    minimap2 = progs['minimap2']
    racon = progs['racon']
    consensus = progs['consensus']
    emtrey = progs['emtrey']
    blat = progs['blat']
    poa = progs['poa']
    medaka = progs['medaka']
else:
    minimap2 = 'minimap2'
    racon = 'racon'
    blat = 'blat'
    poa = 'abpoa'
    medaka = 'medaka_consensus'
    consensus = 'consensus.py'

consensus = 'python3 ' + consensus


def read_fastq_file(seq_file):
    '''
    Takes a FASTQ file and returns a list of tuples
    In each tuple:
        name : str, read ID
        seed : int, first occurrence of the splint
        seq : str, sequence
        qual : str, quality line
        average_quals : float, average quality of that line
        seq_length : int, length of the sequence
    '''
    read_list1 = []
    read_list2 = []

    length = 0
    for line in open(seq_file):
        length += 1
    lineNum = 0
    seq_file_open = open(seq_file, 'r')
    previous = ''
    burn = False
    while lineNum < length:
        name_root = seq_file_open.readline().strip()[1:].split('_')
        name = name_root[0]
        seq = seq_file_open.readline().strip()
        _ = seq_file_open.readline().strip()
        qual = seq_file_open.readline().strip()

        if previous != name:
            burn = False
            previous = name

        number = int(name_root[1])
        if number == 0:
            burn = True

        seq_length = len(seq)
        if burn:
            read_list2.append((name, seq, qual, seq_length))
        else:
            read_list1.append((name, seq, qual, seq_length))

        lineNum += 4
    return read_list1, read_list2


def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict, lastHead = {}, ''
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            if readDict:
                readDict[lastHead] = ''.join(readDict[lastHead])
            readDict[line[1:]] = []
            lastHead = line[1:]
        else:
            readDict[lastHead].append(line.upper())
    if readDict:
        readDict[lastHead] = ''.join(readDict[lastHead])
    return readDict


def determine_consensus(name, fasta, fastq, counter):
    '''Aligns and returns the consensus'''
    corrected_consensus = ''
    repeats = '0'
    fastq_reads_full, fastq_reads_partial = read_fastq_file(fastq)
    fasta_reads = []
    fasta_read_dict = read_fasta(fasta)
    for read, seq in fasta_read_dict.items():
        fasta_reads.append((read, seq))
    repeats = str(len(fasta_reads))

    out_Fq = temp_folder + '/' + counter + '_subsampled.fastq'
    out_F = temp_folder + '/' + counter + '_subsampled.fasta'
    Poa_Dir = temp_folder + '/' + counter + '_poa_direction.sam'
    PIR = temp_folder + '/' + counter + '_poa.fasta'
    combined_consensus_file = open(temp_folder + '/' + counter + '.fasta', 'w')
    out = open(out_Fq, 'w')
    outa = open(out_F, 'w')

    fastq_reads = fastq_reads_full + fastq_reads_partial
    if len(fastq_reads) > 0:
        if len(fastq_reads_full) < subsample:
            subsample_fastq_reads = fastq_reads
        else:
            indeces = np.random.choice(
                np.arange(0, len(fastq_reads_full)),
                min(len(fastq_reads_full), subsample), replace=False)
            subsample_fastq_reads = []
            for index in indeces:
                subsample_fastq_reads.append(fastq_reads_full[index])

        subread_counter = 0
        for read in subsample_fastq_reads:
            subread_counter += 1
            out.write('@' + read[0] + '_' + str(subread_counter) + '\n'
                      + read[1] + '\n+\n' + read[2] + '\n')
        out.close()

        indeces = np.random.choice(np.arange(0, len(fasta_reads)),
                                   min(len(fasta_reads), 20), replace=False)
        subsample_fasta_reads = []
        for index in indeces:
            subsample_fasta_reads.append(fasta_reads[index])
        for read, seq in subsample_fasta_reads:
            outa.write('>' + read + '\n' + seq + '\n')
        outa.close()

        out_First = temp_folder + '/' + counter + '_first.fasta'
        out_first = open(out_First, 'w')
        out_first.write('>' + subsample_fasta_reads[0][0]
                        + '\n' + subsample_fasta_reads[0][1] + '\n')
        out_first.close()
        os.system('%s --secondary=no -ax map-ont \
                  %s %s > %s 2> %s_minimap2_messages.txt'
                  % (minimap2, out_First, out_F,
                     Poa_Dir, temp_folder + '/' + counter))
        direction_dict = {}
        for line in open(Poa_Dir):
            if line[0] != '@':
                a = line.strip().split('\t')
                if len(a) > 1:
                    name1 = a[0]
                    direction = a[1]
                    direction_dict[name1] = direction
        outa = open(out_F, 'w')
        for read, seq in subsample_fasta_reads:
            if read in direction_dict:
                direction = direction_dict[read]
                if direction in ['0', '16']:
                    if direction == '0':
                        outa.write('>' + read + '\n' + seq + '\n')
                    else:
                        outa.write('>' + read + '\n' + revComp(seq) + '\n')
        outa.close()

        poa_cons = temp_folder + '/' + counter + '_consensus.fasta'
        final = temp_folder + '/' + counter + '_corrected_consensus.fasta'
        overlap = temp_folder + '/' + counter + '_overlaps.sam'

        os.system('%s %s > %s ' % (poa, out_F, PIR))
        reads = read_fasta(PIR)
        for read in reads:
            if 'Consensus' in read:
                out_cons_file = open(poa_cons, 'w')
                out_cons_file.write('>Consensus\n'
                                    + reads[read].replace('-', '') + '\n')
                out_cons_file.close()

        final = poa_cons
        for i in np.arange(1, 2):
            try:
                if i == 1:
                    input_cons = poa_cons
                    output_cons = poa_cons.replace('.fasta',
                                                   '_' + str(i) + '.fasta')
                else:
                    input_cons = poa_cons.replace('.fasta',
                                                  '_' + str(i - 1) + '.fasta')
                    output_cons = poa_cons.replace('.fasta',
                                                   '_' + str(i) + '.fasta')

                os.system('%s --secondary=no -ax map-ont \
                          %s %s > %s 2> %s_minimap2_messages.txt'
                          % (minimap2, input_cons, out_Fq,
                             overlap, temp_folder + '/' + counter))

                os.system('%s -q 5 -t 1 --no-trimming\
                          %s %s %s >%s  2>%s_racon_messages.txt'
                          % (racon, out_Fq, overlap, input_cons,
                             output_cons, temp_folder + '/' + counter))

                final = output_cons
            except:
                pass
        os.system('mkdir ' + temp_folder + '/' + counter)
        os.system('%s -f -i %s -d %s -o %s > %s_medaka_messages.txt 2>&1'
                  % (medaka, out_Fq, final,
                     temp_folder + '/' + counter, temp_folder + '/' + counter))
        final = temp_folder + '/' + counter + '/consensus.fasta'
        reads = read_fasta(final)
        if int(counter) % 100 == 0:
            print('\tfinished consensus %s' % (counter), end='\r')
        for read in reads:
            corrected_consensus = reads[read]
        combined_consensus_file = open(path + '/mp/' + counter + '.fasta', 'w')
        combined_consensus_file.write('>' + name + '_' + repeats + '\n'
                                      + corrected_consensus + '\n')
        combined_consensus_file.close()


def main():
    pool = mp.Pool(processes=numThreads)
    print('\tremoving files from previous run')
    os.system('rm -r ' + path + '/mp')
    print('\tdone removing files')
    os.system('mkdir ' + path + '/mp')
    if os.path.exists(path + '/Isoform_Consensi.fasta'):
        os.system('rm ' + path + '/Isoform_Consensi.fasta')
    os.system('touch {0}'.format(path + '/Isoform_Consensi.fasta'))
    counter = 0
    counter_list = []
    for line in open(path + '/isoform_list'):
        counter += 1
        fasta = line.split('\t')[0]
        fastq = line.split('\t')[1]
        name = line.split('\t')[2].strip()
        pool.apply_async(determine_consensus, [name, fasta, fastq, str(counter)])
        # determine_consensus(name, fasta, fastq, str(counter))
        counter_list.append(str(counter))
    print('\tmaking consensus sequences of', counter, 'isoforms')

    pool.close()
    pool.join()
    combined_consensus_file = open(path + '/Isoform_Consensi.fasta', 'w')

    for counter in counter_list:
        if os.path.exists(path + '/mp/' + counter + '.fasta'):
            for lines in open(path + '/mp/' + counter + '.fasta'):
                combined_consensus_file.write(lines)
    combined_consensus_file.close()


main()
