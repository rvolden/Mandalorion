import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-v', '--vcf_file', type=str)
parser.add_argument('-r', '--pslx_file_for_phasing', type=str)
parser.add_argument('-s', '--pslx_files_to_be_sorted', type=str)
parser.add_argument('-p', '--output_path', type=str)

args = parser.parse_args()
vcf_file = args.vcf_file
read_file = args.pslx_file_for_phasing
sort_files = args.pslx_files_to_be_sorted
outfile_path = args.output_path

def read_vcf(vcf_file):
    snps = {}
    for line in open(vcf_file):
        if line[0] == '#':
            continue
        a = line.strip().split('\t')
        chromosome, position = a[0], a[1]
        ref, alt = a[3], a[4]

        status = a[9].split(':')[0]
        alt_split = alt.split(',')
        if len(alt_split) > 1:
            ref = alt = alt_split[0]

        first, second = status.split('/')[0], status.split('/')[1]
        if first != second:
            if chromosome not in snps:
                snps[chromosome] = {}
            if len(ref) == 1 and len(alt) == 1 and ref != alt:
                snps[chromosome][int(position)] = ((ref, first), (alt, second))
    return(snps)

def connecting_snps(comb_dict):
    for snp1 in comb_dict:
        add_set = set()
        for snp2 in comb_dict[snp1]:
            comb_dict[snp2].add(snp1)
            for snp3 in comb_dict[snp1]:
                comb_dict[snp2].add(snp3)
            for snp3 in comb_dict[snp2]:
                add_set.add(snp3)
        for snp3 in add_set:
            comb_dict[snp1].add(snp3)
    return comb_dict

def parse_reads(snps, read_file):
    comb_dict, con_dict, snp_dict = {}, {}, {}
    count, total = 0, 0
    reads, snp_count = {}, {}
    for line in open(read_file):
        count += 1
        if count % 100000 == 0:
            print('parsing read #', count)
        a = line.strip().split('\t')
        chromosome = a[13]

        start, end = int(a[15]), int(a[16])
        total += 1
        name = a[9]
        reads[name] = line
        offset_start = int(a[11])
        offset_end = int(a[10]) - int(a[12])

        blocksizes = np.array(a[18].split(',')[:-1], dtype=int)
        blockstarts = np.array(a[20].split(',')[:-1], dtype=int)

        start = min(blockstarts)
        end = max(blockstarts + blocksizes)

        readSequences = a[21].split(',')[:-1]
        genomeSequences = a[22].split(',')[:-1]
        readList, snpList = [], []

        for block in range(0, len(blockstarts), 1):
            genomeStart = blockstarts[block]
            readSequence = readSequences[block]
            genomeSequence = genomeSequences[block]
            for pos in range(0, len(genomeSequence), 1):
                readBase = readSequence[pos]
                genomeBase = genomeSequence[pos]
                if genomeStart + pos + 1 in snps[chromosome]:
                    ref, alt = snps[chromosome][genomeStart + pos + 1]
                    if readBase == ref[0] or readBase == alt[0]:
                        readList.append((chromosome, genomeStart+pos+1))
                        snpList.append((chromosome, genomeStart+pos+1, readBase))

                        if (chromosome, genomeStart+pos+1, readBase) not in snp_count:
                            snp_count[(chromosome, genomeStart+pos+1, readBase)] = 0
                        snp_count[(chromosome, genomeStart+pos+1, readBase)] += 1

        snp_dict[name] = snpList

        for snp1 in readList:
            if snp1 not in con_dict:
                con_dict[snp1] = {}
            for snp2 in readList:
                if snp1 != snp2:
                    if snp2 not in con_dict[snp1]:
                        con_dict[snp1][snp2] = 0
                    con_dict[snp1][snp2] += 1

    read_minimum = (1/1000000) * total

    for snp1, snps2 in con_dict.items():
        if snp1 not in comb_dict:
            comb_dict[snp1] = set()
        comb_counts = []
        for snp2 in snps2:
            comb_count = con_dict[snp1][snp2]
            comb_counts.append(comb_count)
        minimum = np.average(comb_counts) / 5
        for snp2 in snps2:
            comb_count = con_dict[snp1][snp2]
            if comb_count > minimum:
                 comb_dict[snp1].add(snp2)

    comb_dict = connecting_snps(comb_dict)
    comb_dict = connecting_snps(comb_dict)
    comb_dict = connecting_snps(comb_dict)

    used = {}
    groups = set()
    for snp1, group in comb_dict.items():
         groups.add(tuple(sorted(list(group))))

    connections = {}
    for name, snpList in snp_dict.items():
        for snp1 in snpList:
            if snp1 not in connections:
                connections[snp1] = {}
            for snp2 in snpList:
                if snp1 != snp2:
                    if snp2 not in connections[snp1]:
                         connections[snp1][snp2] = 0
                    connections[snp1][snp2] += 1
    return groups, connections, snp_dict, reads, snp_count, read_minimum

def filter_group(group, snps, snp_count):
    maximum, group_scores, group_filtered = 0, [], []
    for i in list(range(0, len(group), 1)):
        testing_snp = group[i]
        chromosome = testing_snp[0]
        position = testing_snp[1]
        options = snps[chromosome][position]

        for option in snps[chromosome][position]:
            if (chromosome, position, option[0]) not in snp_count:
                snp_count[(chromosome, position, option[0])] = 0
        score = np.abs(snp_count[(chromosome, position, options[0][0])] \
                       + snp_count[(chromosome, position, options[1][0])]) \
                - np.abs(snp_count[(chromosome, position, options[0][0])] \
                         - snp_count[(chromosome, position, options[1][0])])
        if score > maximum:
            maximum = score
        group_scores.append((testing_snp, score))
    for testing_snp, score in group_scores:
        if score / maximum > 0.2:
            group_filtered.append((testing_snp, score))

    group_filtered = sorted(group_filtered, key=lambda x: x[1], reverse=True)
    return group_filtered

def determine_haplo_snps(
        group_filtered, snps, haplo_snps, connections,
        starting_index, snp_count, read_minimum):

    for i in list(range(0, len(group_filtered), 1)):
        if i != starting_index:
            snp = group_filtered[i][0]
            chromosome, position = snp[0], snp[1]
            options = snps[chromosome][position]
            first, second = 0, 0

            for haplo_snp in haplo_snps:
                if haplo_snp in connections:
                    if (chromosome, position, options[0][0]) in connections[haplo_snp]:
                        first += connections[haplo_snp][(chromosome, position, options[0][0])]
                    if (chromosome, position, options[1][0]) in connections[haplo_snp]:
                        second += connections[haplo_snp][(chromosome, position, options[1][0])]

            if first/snp_count[(chromosome, position, options[0][0])] \
                    > ((second/snp_count[(chromosome, position, options[1][0])])*3):
                if first > read_minimum:
                    haplo_snps.add((chromosome, position, options[0][0]))
            elif second/snp_count[(chromosome, position, options[1][0])]
                    > ((first/snp_count[(chromosome, position, options[0][0])])*3):
                if second > read_minimum:
                    haplo_snps.add((chromosome, position, options[1][0]))
    return haplo_snps

def determine_haplo_group(
        group_filtered, snps, haplo_snps, connections,
        haplo_group, snp_count, hap, show, read_minimum):

    for i in list(range(0, len(group_filtered), 1)):
        snp = group_filtered[i][0]
        chromosome, position = snp[0], snp[1]
        options = snps[chromosome][position]
        first, second = 0, 0

        for haplo_snp in haplo_snps:
            if haplo_snp in connections:
                if (chromosome, position, options[0][0]) in connections[haplo_snp]:
                    first += connections[haplo_snp][(chromosome, position, options[0][0])]
                if (chromosome, position, options[1][0]) in connections[haplo_snp]:
                    second += connections[haplo_snp][(chromosome, position, options[1][0])]
        if show:
            print(chromosome, position, options, first, second,
                  snp_count[(chromosome, position, options[0][0])],
                  snp_count[(chromosome, position, options[1][0])])
        taken = False
        if first/snp_count[(chromosome, position, options[0][0])]
                > ((second/snp_count[(chromosome, position, options[1][0])])*3):
            if snp_count[(chromosome, position, options[0][0])] > read_minimum \
                    and first > (snp_count[(chromosome, position, options[0][0])]*len(haplo_snps))/5:
                haplo_group[(chromosome, position, options[0][0])] = ('f', hap, first, second, i)
                if show:
                    print('Taking first')
                    taken = True
        elif second/snp_count[(chromosome, position, options[1][0])]
                > ((first/snp_count[(chromosome, position, options[0][0])])*3):
            if snp_count[(chromosome, position, options[1][0])] > read_minimum \
                    and second > (snp_count[(chromosome, position, options[1][0])]*len(haplo_snps))/5:
                haplo_group[(chromosome, position, options[1][0])] = ('s', hap, first, second, i)
                if show:
                    print('Taking second')
                    taken = True
        if show and not taken:
                print('Not taking it')
    return haplo_group

def create_haplotypes(groups, connections, snps, snp_count, read_minimum):
    haplo_group = {}
    bed = open(outfile_path + '/snp.bed', 'w')

    print('read_minimum', read_minimum)
    for group in groups:
        group_filtered = filter_group(group, snps, snp_count)
        if len(group_filtered) > 0:
            starting_index, count = 0, 0
            starting_snp = group_filtered[starting_index][0]
            starting_chromosome = starting_snp[0]
            starting_position = starting_snp[1]
            starting_options = snps[starting_chromosome][starting_position]
            maternal = starting_options[0][0]
            paternal = starting_options[1][0]

            haplo_snps_maternal = set()
            haplo_snps_paternal = set()

            full_snp_maternal = (starting_chromosome, starting_position, maternal)
            full_snp_paternal = (starting_chromosome, starting_position, paternal)

            haplo_snps_maternal.add(full_snp_maternal)
            haplo_snps_paternal.add(full_snp_paternal)
            haplo_snps_maternal = determine_haplo_snps(group_filtered, snps,
                                                       haplo_snps_maternal,
                                                       connections, starting_index,
                                                       snp_count, read_minimum)
            haplo_snps_paternal = determine_haplo_snps(group_filtered, snps,
                                                       haplo_snps_paternal,
                                                       connections, starting_index,
                                                       snp_count, read_minimum)
            show = False
            haplo_group = determine_haplo_group(group_filtered, snps,
                                                haplo_snps_maternal,
                                                connections, haplo_group,
                                                snp_count, 0, show, read_minimum)
            haplo_group = determine_haplo_group(group_filtered, snps,
                                                haplo_snps_paternal,
                                                connections, haplo_group,
                                                snp_count, 1, show, read_minimum)

    for element, value in haplo_group.items():
        bed.write(element[0] + '\t' + str(element[1]) + '\t' \
                  + str(element[1]) + '\t+\t' + str(value[1]) + '\n')
    return haplo_group

def sort_reads(haplo_group, snp_dict, reads, sort_file):
    maternal = open(sort_file + '.allele1', 'w')
    paternal = open(sort_file + '.allele2', 'w')
    undetermined = open(sort_file + '.undetermined', 'w')

    file_dict = {}
    file_dict['0'] = maternal
    file_dict['1'] = paternal
    file_dict['U'] = undetermined

    for name, variants in snp_dict.items():
        var_list = []
        for var in variants:
            if var in haplo_group:
                var_list.append(haplo_group[var][1])
        length = len(set(var_list))
        if length == 1:
            haplotype = str(var_list[0])
        else:
            haplotype = 'U'
        if haplotype in ['0', '1']:
#        if haplotype in ['0', '1', 'U']: ## comment out line above and remove # at being of this line if you want to write out undetermined reads
            file_dict[haplotype].write(reads[name])

def main():
    print('reading snps')
    snps = read_vcf(vcf_file)
    print('parsing reads')
    groups, connections, snp_dict, reads, snp_count, read_minimum = parse_reads(snps, read_file)
    print('phasing snps')
    haplo_group = create_haplotypes(groups, connections, snps, snp_count, read_minimum)
    print('sorting reads')
    for sort_file in sort_files.split(','):
        groups, connections, snp_dict, reads, snp_count, read_minimum = parse_reads(snps, sort_file)
        sort_reads(haplo_group, snp_dict, reads, sort_file)

main()
