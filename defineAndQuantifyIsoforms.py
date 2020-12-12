#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import os
import sys
import numpy as np

path = sys.argv[2]
infile = sys.argv[1]
upstream_buffer = int(sys.argv[4])
downstream_buffer = int(sys.argv[3])
subreads = sys.argv[5]
fasta_files = sys.argv[6]
minimum_read_count = int(sys.argv[7])


def find_peaks(starts, ends):
    start_peaks, end_peaks = {}, {}

    start_count = {}
    end_count = {}
    for position in starts:
        if position not in start_count:
            start_count[position] = 0
        start_count[position] += 1

    for position in ends:
        if position not in end_count:
            end_count[position] = 0
        end_count[position] += 1

    for position in sorted(starts):
        if position not in start_peaks:
            window_count = 0
            for i in range(0, 10, 1):
                window_position = position + i
                if window_position in start_count:
                    window_count += start_count[window_position]
            if window_count >= minimum_read_count:
                for shift in np.arange(-upstream_buffer, downstream_buffer):
                    start_peaks[position + shift] = position

    for position in sorted(ends, reverse=True):
        if position not in end_peaks:
            window_count = 0
            for i in range(0, 10, 1):
                window_position = position - i
                if window_position in end_count:
                    window_count += end_count[window_position]
            if window_count >= minimum_read_count:
                for shift in np.arange(-downstream_buffer, upstream_buffer):
                    end_peaks[position + shift] = position

    return start_peaks, end_peaks


def collect_splice_events(path):
    splice_dict = {}
    for line in open(path + '/SS.bed'):
        a = line.strip().split('\t')
        chromosome = a[0]
        if not splice_dict.get(chromosome):
            splice_dict[chromosome] = {}
        splice_left = int(a[1])
        splice_right = int(a[2])
        for base in np.arange(splice_left, splice_right + 1):
            splice_dict[chromosome][base] = a[3].split('_')[0]
    return splice_dict


def sort_reads_into_splice_junctions(splice_dict, fasta_files, infile):
    start_end_dict, start_end_dict_mono, readDict = {}, {}, {}
    for fasta_file in fasta_files.split(','):
        tempSeqs, headers, sequences = [], [], []
        for line in open(fasta_file):
            line = line.rstrip()
            if not line:
                continue
            if line.startswith('>'):
                headers.append(line.split()[0][1:])
            # covers the case where the file ends while reading sequences
            if line.startswith('>'):
                sequences.append(''.join(tempSeqs).upper())
                tempSeqs = []
            else:
                tempSeqs.append(line)
        sequences.append(''.join(tempSeqs).upper())
        sequences = sequences[1:]
        for i in range(len(headers)):
            readDict[('-').join(headers[i].split('_')[0].split('-')[:5])] = [headers[i], sequences[i]]
    read_dict = readDict

    for line in open(infile):
        a = line.strip().split('\t')
        read_chromosome, read_direction = a[13], a[8]
        name = a[9].split('_')[0]
        name = ('-').join(name.split('-')[:5])
        show = False
        read_direction = '+'  # ignores read direction
        start, end = int(a[15]), int(a[16])
        if read_direction == '+':
            left_extra, right_extra = int(a[11]), int(a[10]) - int(a[12])
        if read_direction == '-':
            right_extra, left_extra = int(a[11]), int(a[10]) - int(a[12])

        failed = False
        identity = read_chromosome + '_'

        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]

        for x in range(0, len(blocksizes) - 1):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            left_splice = blockstart + blocksize
            right_splice = int(blockstarts[x + 1])
            if right_splice - left_splice > 50:
                if not splice_dict.get(read_chromosome):
                    failed = True
                    break
                else:
                    if not splice_dict[read_chromosome].get(left_splice) or not splice_dict[read_chromosome].get(right_splice):
                        failed = True
                        break
                left_splice_site = splice_dict[read_chromosome][left_splice]
                right_splice_site = splice_dict[read_chromosome][right_splice]
                identity += str(left_splice_site) + '-' + str(right_splice_site) + '~'
        if show:
            print(failed)
        if not failed:
            if identity.split('_')[1] != '':
                if not start_end_dict.get(identity):
                    start_end_dict[identity] = []
                start_end_dict[identity].append((start, end,
                                                 '>' + read_dict[name][0] + '\n'
                                                 + read_dict[name][1] + '\n',
                                                 left_extra,
                                                 right_extra,
                                                 read_direction))
            else:
                if not start_end_dict_mono.get(identity):
                    start_end_dict_mono[identity] = []

                start_end_dict_mono[identity].append((start, end,
                                                      '>' + read_dict[name][0] + '\n'
                                                      + read_dict[name][1] + '\n',
                                                      left_extra,
                                                      right_extra,
                                                      read_direction))

    return start_end_dict, start_end_dict_mono


def group_mono_exon_transcripts(start_end_dict, start_end_dict_mono):
    for identity in start_end_dict_mono:
        previous_end = 0
        iso_counter = 0
        new_identity = identity + 'M' + str(iso_counter)
        positions = start_end_dict_mono[identity]
        for position in sorted(positions):
            start = position[0]
            end = position[1]
            if start > previous_end:
                iso_counter += 1
                new_identity = identity + 'M' + str(iso_counter)
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = max(end, previous_end)
            else:
                if new_identity not in start_end_dict:
                    start_end_dict[new_identity] = []
                start_end_dict[new_identity].append(position)
                previous_end = end

    return start_end_dict


def define_start_end_sites(start_end_dict, start_end_dict_mono, individual_path):
    left_extras, right_extras = {}, {}
    file_set = set()
    isoform_counter, isoform_dict, subread_pointer = 0, {}, {}

    start_end_dict = group_mono_exon_transcripts(start_end_dict, start_end_dict_mono)

    number_of_isoforms = len(start_end_dict)
    counter = 0
    for identity in start_end_dict:
        counter += 1
        print('\tprocessing reads assigned to spliceform ', counter,
              ' of ', number_of_isoforms, ' into isoforms', end='\r')
        starts, ends = [], []
        ends = []
        positions = start_end_dict[identity]
        length_of_positions = len(positions)
        indexes = np.random.choice(range(0, length_of_positions, 1),
                                   size=min(10000, length_of_positions),
                                   replace=False)
        subsampled_positions = []
        for index in indexes:
            subsampled_positions.append(positions[index])
        for position in subsampled_positions:
            starts.append(int(position[0]))
            ends.append(int(position[1]))

        start_dict, end_dict = find_peaks(starts, ends)
        matched_positions = []
        left_extras[identity], right_extras[identity] = {}, {}
        for start, end, read, left_extra, right_extra, read_direction in positions:
            if int(start) in start_dict and int(end) in end_dict:
                left = start_dict[int(start)]
                right = end_dict[int(end)]
                if not left_extras[identity].get((left, right)):
                    left_extras[identity][(left, right)] = []
                    right_extras[identity][(left, right)] = []

                left_extras[identity][(left, right)].append(int(left_extra))
                right_extras[identity][(left, right)].append(int(right_extra))
                matched_positions.append((left, right, read, read_direction))

        median_left = {}
        median_right = {}
        for combination, values in left_extras[identity].items():
            median = np.median(values)
            median_left[combination] = median
        for combination, values in right_extras[identity].items():
            median = np.median(values)
            median_right[combination] = median

        for left, right, read, read_direction in matched_positions:
            medianLeft = median_left[(left, right)]
            medianRight = median_right[(left, right)]
            new_identity = identity + '_' + read_direction + '_' + str(left) \
                + '_' + str(right) + '_' \
                + str(round(medianLeft, 2)) \
                + '_' + str(round(medianRight, 2))
            if not isoform_dict.get(new_identity):
                isoform_counter += 1
                isoform_dict[new_identity] = isoform_counter

            subfolder = str(int(isoform_dict[new_identity] / 4000))
            if subfolder not in os.listdir(individual_path + '/parsed_reads/'):
                os.makedirs(individual_path + '/parsed_reads/' + subfolder)
            filename = subfolder + '/Isoform' + str(isoform_dict[new_identity])
            out_reads_fasta = open(individual_path + '/parsed_reads/' + filename + '.fasta', 'a')
            out_reads_subreads = open(individual_path + '/parsed_reads/' + filename + '_subreads.fastq', 'w')
            out_reads_fasta.write(read)
            out_reads_fasta.close()
            out_reads_subreads.close()

            file_set.add(individual_path + '/parsed_reads/' + filename
                         + '.fasta' + '\t' + individual_path
                         + '/parsed_reads/' + filename + '_subreads.fastq'
                         + '\t' + new_identity + '\n')
            read = read.split()[0].split('_')[0][1:]

            subread_pointer[read] = individual_path + '/parsed_reads/' + filename + '_subreads.fastq'
            # for subread, sequence, qual in subread_list:
            #     out_reads_subreads.write(subread + '\n' + sequence
            #                              + '\n+\n' + qual + '\n')

    out = open(individual_path + 'isoform_list', 'a')
    for item in file_set:
        out.write(item)
    out.close()
    return subread_pointer


def read_subreads(seq_file, infile, subread_pointer):
    lineNum, lastPlus, root_name = 0, False, ''
    name, seq, qual = '', '', ''
    for line in open(seq_file):
        line = line.rstrip()
        if not line:
            continue

        # make an entry as a list and append the header to that list
        if lineNum % 4 == 0 and line[0] == '@':
            if lastPlus:  # chrom_reads needs to contain root_names
                filepath = subread_pointer.get(root_name)
                if filepath:
                    outfile = open(filepath, 'a')
                    outfile.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))
                    outfile.close()
            name = line[1:]
            root_name = ('-').join(name.split('_')[0].split('-')[:5])

        if lineNum % 4 == 1:
            seq = line
        if lineNum % 4 == 2:
            lastPlus = True
        if lineNum % 4 == 3 and lastPlus:
            qual = line
        lineNum += 1


def main():
    individual_path = path
    print('\terasing files from previous run')
    os.system('rm -r ' + individual_path + '/parsed_reads/')
    os.system('rm -r ' + individual_path + '/mp/')
    os.system('mkdir ' + individual_path + '/parsed_reads')
    listDir = os.listdir(individual_path + '/parsed_reads')
    for file in listDir:
        print(file)
        os.system('rm -r ' + individual_path + '/parsed_reads/' + file)
    out = open(individual_path + 'isoform_list', 'w')
    out.close()

    print('\treading splice sites')
    splice_dict = splice_dict = collect_splice_events(path)
    print('\tdefining spliceforms')
    start_end_dict, start_end_dict_mono = sort_reads_into_splice_junctions(splice_dict, fasta_files, infile)
    print('\tdefining isoforms')
    subread_pointer = define_start_end_sites(start_end_dict, start_end_dict_mono, individual_path)
    print('\treading subreads')
    for subread_file in subreads.split(','):
        print('    ' + subread_file)
        read_subreads(subread_file, infile, subread_pointer)


if __name__ == '__main__':
    main()
