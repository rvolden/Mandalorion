import sys
import editdistance
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--reads_fasta_file', type=str)
parser.add_argument('-g', '--genome_fasta_file', type=str)
parser.add_argument('-f', '--hla_fasta_file', type=str)
parser.add_argument('-c', '--config_file', type=str)

args = parser.parse_args()
fasta_file = args.reads_fasta_file
genome_file = args.genome_fasta_file
hla_file = args.hla_fasta_file
config_file= args.config_file

def configReader(configIn):
    '''Parses the config file.'''
    progs = {}
    for line in open(configIn):
        if line.startswith('#') or not line.rstrip().split():
            continue
        line = line.rstrip().split('\t')
        progs[line[0]] = line[1]
    # should have minimap, poa, racon, water, consensus
    # check for extra programs that shouldn't be there
    possible = set(['poa', 'minimap2', 'gonk', 'consensus', 'racon', 'blat', 'emtrey', 'psl2pslx'])
    inConfig = set()
    for key in progs.keys():
        inConfig.add(key)
        # if key not in possible:
        #     raise Exception('Check config file')
    # check for missing programs
    # if missing, default to path
    for missing in possible-inConfig:
        if missing == 'consensus':
            path = 'consensus.py'
        else:
            path = missing
        progs[missing] = path
        sys.stderr.write('Using ' + str(missing)
                         + ' from your path, not the config file.\n')
    return progs

progs = configReader(config_file)
minimap2 = progs['minimap2']
emtrey = progs['emtrey']
psl2pslx = progs['psl2pslx']

def read_fasta(inFile):
    '''Reads in FASTA files, returns a dict of header:sequence'''
    readDict = {}
    for line in open(inFile):
        line = line.rstrip()
        if not line:
            continue
        if line.startswith('>'):
            readDict[line[1:]] = ''
            lastHead = line[1:]
        else:
            readDict[lastHead] += line.upper()
    return readDict

def find_gene_match(fasta_file, gene_list, genome_sequence):
    sam_file = fasta_file + '.sam'
    psl_file = fasta_file + '.psl'
    print('aligning isoforms to genome')
    os.system('%s --secondary=no -ax splice -t %s %s %s > %s ' %(minimap2, '1', genome_sequence, fasta_file, sam_file))
    print('converting sam to psl')
    os.system('%s -i %s > %s ' %(emtrey, sam_file, psl_file))
    print('Finding isoforms of HLA-genes')
    match_dict = {}
    for gene in gene_list:
        gene_name = gene[0]
        gene_chromosome = gene[1]
        gene_start = gene[2]
        gene_end = gene[3]
        for line in open(psl_file):
            a = line.strip().split('\t')
            name = a[9]
            chromosome = a[13]
            start = int(a[15])
            end = int(a[16])

            if chromosome == gene_chromosome:
                if gene_start < start < gene_end:
                    match_dict[name] = gene_name
    return match_dict

def collect_hla_names(hla_file):
    hla_seqs = read_fasta(hla_file)
    read_dict = {}
    for read in hla_seqs:
        name = read.split()[0]
        cat = read.split()[1]
        read_dict[name] = cat
    return read_dict

def find_allele_match(fasta_file, match_dict, hla_file):
    sam_file = fasta_file + '.hla.sam'
    psl_file = fasta_file + '.hla.psl'
    pslx_file = fasta_file + '.hla.pslx'
    print('aligning isoforms to hla_sequences')
    os.system('%s -ax map-ont -N 100 -t %s %s %s > %s ' %(minimap2, '1', hla_file, fasta_file, sam_file))
    print('converting sam to psl')
    os.system('%s -i %s > %s ' %(emtrey, sam_file, psl_file))
    print('converting psl to pslx')
    os.system('python3 %s -p %s -r %s -g %s -x %s ' %(psl2pslx, psl_file,
                                                      fasta_file, hla_file,
                                                      pslx_file))
    hla_dict = {}
    for line in open(pslx_file):
        a = line.strip().split('\t')

        name = a[9]
        hla_name = a[13]
        align_start = a[15]
        align_end = a[16]
        length = a[14]
        if align_start == '0' and align_end == length:
            cat = read_dict[hla_name]
            if name in match_dict:
                if name not in hla_dict:
                    hla_dict[name] = []
                hla_dict[name].append((int(a[1]), int(a[3])+int(a[5])+int(a[7]),
                                       hla_name, cat, match_dict[name], a[14]))
    for HLA_gene in ['A', 'B', 'C', 'DRA', 'DRB1', 'DPA1', 'DPB1', 'DQA1', 'DQB1']:
        for name, list1 in hla_dict.items():
            list1 = sorted(list1, key=lambda x :(x[0],x[1]))
            hla_name = list1[0][3]

            if hla_name.split('*')[0] == HLA_gene:
                print('Best Match', list1[0][3], 'Mismatches:', list1[0][0],
                      'Indels:', list1[0][1], 'Alignment length:', list1[0][5])

gene_list = []

gene_list.append(('HLA-A', 'chr6', 29942207, 29946087))
gene_list.append(('HLA-B', 'chr6', 31353043, 31358016))
gene_list.append(('HLA-C', 'chr6', 31268477, 31272311))
gene_list.append(('HLA-DRA1', 'chr6', 32439285, 32445414))
gene_list.append(('HLA-DRB1', 'chr6', 32578224, 32590373))
gene_list.append(('HLA-DRB5', 'chr6', 32516982, 32530918))
gene_list.append(('HLA-DPA1', 'chr6', 33064002, 33074116))
gene_list.append(('HLA-DPB1', 'chr6', 33075089, 33087875))
gene_list.append(('HLA-DQA1', 'chr6', 32637024, 32643690))
gene_list.append(('HLA-DQB1', 'chr6', 32658674, 32667200))

read_dict = collect_hla_names(hla_file)
match_dict = find_gene_match(fasta_file, gene_list, genome_file)
find_allele_match(fasta_file, match_dict, hla_file)
