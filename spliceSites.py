#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import numpy as np
import time

infile =  sys.argv[1]
out_path = sys.argv[2]
cutoff = float(sys.argv[3])
genome_file = sys.argv[4]
refine = sys.argv[5]
sam_file =  sys.argv[6]
splice_site_width = int(sys.argv[7])
minimum_read_count = int(sys.argv[8])
out = open(out_path + '/SS.bed', 'w')

def scan_for_best_bin(entry, distance_range, iterator_shift, density_dict,
                      peak_areas, chromosome, side):
    '''
    Stuff about the inputs, outputs, and how it works
    '''

    best_extra_list, peak_center, bases = 0, 0, 0
    coverage_area, best_directions = {}, {}
    best_directions['+']=0
    best_directions['-']=0
    for x in distance_range:
        extra_list_expression = 0
        directions = {}
        directions['+'] = 0
        directions['-'] = 0
        coverage_set = {}
        called = False
        for y in distance_range:
            if entry+x+y in peak_areas[chromosome][side]:
                called = True
        if not called:
            for y in distance_range:
                if entry+x+y in density_dict:
                    for item in density_dict[entry+x+y]:
                        extra_list_expression+=1
                        directions[item[4]] += 1
                        for covered_position in item[3]:
                            if covered_position not in coverage_set:
                                coverage_set[covered_position] = 0
                            coverage_set[covered_position] += 1

        if extra_list_expression > best_extra_list:
            best_extra_list = extra_list_expression
            peak_center = entry + x
            coverage_area = coverage_set
            best_directions = directions

    return best_extra_list, peak_center, coverage_area, best_directions

def determine_coverage(coverage_area, chromosome, reverse,
                       peak_center, histo_coverage):
    coverage = [0]
    coverage_area2 = []
    for covered_position,count in coverage_area.items():
        if count > 1:
            coverage_area2.append(covered_position)

    coverage_area = sorted(coverage_area2, reverse=reverse)
    counter = 0
    for base_f in coverage_area:
        count = False
        if not reverse:
            if base_f > peak_center:
                count = True
        elif reverse:
            if base_f < peak_center:
                count = True
        if count:
            if counter <= 3:
                counter += 1
                base_f = myround(base_f)
                if base_f in histo_coverage[chromosome]:
                    coverage.append(histo_coverage[chromosome][base_f])
            else:
                break
    coverage = max(coverage)
    return coverage, coverage_area

def myround(x, base=10):
    '''Rounds to the nearest base'''
    return int(base * round(float(x)/base))

def find_peaks(density_dict, out, peaks, reverse, cutoff,
               histo_coverage, side, peak_areas,chromosome):
    '''
    Insert docstring
    '''

    if not reverse:
        distance_range = range(-splice_site_width, splice_site_width)
        iterator_shift = 1
    if reverse:
        distance_range = range(splice_site_width, -splice_site_width, -1)
        iterator_shift =- 1

    entry_list = []
    for entry in density_dict:
        entry_list.append([entry, density_dict[entry]])
    for entry, density in sorted(entry_list,
                                 key=lambda x: len(x[1]),
                                 reverse=True):

        if len(density) >= minimum_read_count:
            if entry not in peak_areas[chromosome][side]:

                best_extra_list, peak_center, \
                coverage_area, best_directions \
                = scan_for_best_bin(entry, distance_range, iterator_shift,
                                    density_dict, peak_areas,
                                    chromosome, side)

                coverage, coverage_area \
                = determine_coverage(coverage_area, chromosome, reverse,
                                     peak_center, histo_coverage)


                if coverage > 0:
                    proportion = round(best_extra_list/coverage, 3)
                    if proportion > cutoff:
                        Left_to_Right = best_directions['+']
                        Right_to_Left = best_directions['-']
                        Type = '-'

                        if Left_to_Right < Right_to_Left:
                            if reverse:
                                Type = '3'
                            else:
                                Type = '5'

                        elif Left_to_Right > Right_to_Left:
                            if reverse:
                                Type = '5'
                            else:
                                Type = '3'

                        if Type != '-':
                            peaks += 1
                            out.write(chromosome + '\t'
                                      + str(peak_center-splice_site_width)
                                      + '\t' + str(peak_center+splice_site_width)
                                      + '\t' + str(Type) + side + str(peaks) + '_'
                                      + str(peak_center-splice_site_width) + '_'
                                      + str(peak_center+splice_site_width) + '_'
                                      + str(proportion) + '\t' + str(peaks) + '\n')
                            for base in range(peak_center - splice_site_width,
                                              peak_center + splice_site_width +1):
                                peak_areas[chromosome][side][base] = 1
        else:
            break

    return peaks, peak_areas

def collect_reads(infile, sam_file,direction_dict,target_chromosome):
    histo_left_bases, histo_right_bases, histo_coverage = {}, {}, {}
    histo_coverage[target_chromosome] = {}
    histo_left_bases[target_chromosome] = {}
    histo_right_bases[target_chromosome] = {}




    for line in open(infile):
        a = line.strip().split('\t')
        chromosome = a[13]
        if chromosome == target_chromosome:
            name = a[9].split('_')[0]
            direction = a[8]
            if name in direction_dict:
                direction = direction_dict[name]

            length = int(a[10])
            begin, span = int(a[15]), int(a[16])
            blocksizes = a[18].split(',')[:-1]
            blockstarts = a[20].split(',')[:-1]

            coverage_set = set()
            low_bounds, up_bounds = [], []
            aligned_bases = 0
            for x in range(0, len(blocksizes)):
                blockstart = int(blockstarts[x])
                blocksize = int(blocksizes[x])
                aligned_bases += blocksize
                blockend = blockstart + blocksize
                for y in range(0, blocksize, 10):
                    rounded = myround(blockstart + y)
                    coverage_set.add(rounded)
                for yy in range(y, blocksize):
                    rounded = myround(blockstart + yy)
                    coverage_set.add(rounded)
                if blockstart != begin:
                    up_bounds.append(blockstart)
                if blockend != span:
                    low_bounds.append(blockend)

            for rounded in coverage_set:
                if rounded not in histo_coverage[chromosome]:
                    histo_coverage[chromosome][rounded] = 0
                histo_coverage[chromosome][rounded] += 1

            if aligned_bases/length > 0.70:
                for low_bound in low_bounds:
                    if low_bound not in histo_left_bases[chromosome]:
                        histo_left_bases[chromosome][low_bound] = []
                    histo_left_bases[chromosome][low_bound].append([0, begin,
                                                                    span, coverage_set,
                                                                    direction])
                for up_bound in up_bounds:
                    if up_bound not in histo_right_bases[chromosome]:
                        histo_right_bases[chromosome][up_bound] = []
                    histo_right_bases[chromosome][up_bound].append([0, begin,
                                                                    span, coverage_set,
                                                                    direction])

    return histo_left_bases, histo_right_bases, histo_coverage

def parse_genome(input_file, left_bounds, right_bounds):
    chromosome_list = set()
    gene_dict = {}
    for line in open(input_file):
        a = line.strip().split('\t')
        if len(a) > 7:
             if a[2] == 'exon':
                 testKey = a[8].split('transcript_id "')[1].split('"')[0]
                 if not gene_dict.get(testKey):
                     gene_dict[testKey] = []
                 gene_dict[testKey].append((a[0], a[3], a[4], a[6]))

    read_list = []
    for transcript_id in gene_dict:
        transcript_data = gene_dict[transcript_id]

        chromosome = transcript_data[0][0]
        chromosome_list.add(chromosome)
        if chromosome not in right_bounds:
            left_bounds[chromosome], right_bounds[chromosome]= {}, {}
            left_bounds[chromosome]['5'], right_bounds[chromosome]['5'] = [], []
            left_bounds[chromosome]['3'], right_bounds[chromosome]['3'] = [], []

        start = sorted(transcript_data, key=lambda x: int(x[1]))[0][1]
        end = sorted(transcript_data, key=lambda x: int(x[2]), reverse=True)[0][2]

        for entry in transcript_data:
            if entry[1] != start:
                if entry[3] == '+':
                    right_bounds[chromosome]['3'].append(int(entry[1])-1)
                elif entry[3] == '-':
                    right_bounds[chromosome]['5'].append(int(entry[1])-1)
            if entry[2] != end:
                if entry[3] == '+':
                    left_bounds[chromosome]['5'].append(int(entry[2]))
                if entry[3] == '-':
                    left_bounds[chromosome]['3'].append(int(entry[2]))
    return chromosome_list, left_bounds, right_bounds

def make_genome_bins(bounds, side, peaks, chromosome, peak_areas):
    for type1 in ['5', '3']:
        covered = {}
        position_list = sorted(bounds[type1], key=int)
        for index1 in range(0, len(position_list)):
            if index1 not in covered:
                sub_list = []
                sub_list.append(position_list[index1])
                for index2 in range(index1, len(position_list)):
                    if position_list[index2] - max(sub_list) <= splice_site_width:
                        sub_list.append(position_list[index2])
                        covered[index2] = 1
                    else:
                        break
                single = False
                if len(sub_list) > 1:
                    splice_distances = []
                    for splice_pos in range(0, len(sub_list)-1):
                        splice_distances.append(int(sub_list[splice_pos+1])
                                                - int(sub_list[splice_pos]))
                    if min(splice_distances) > 3:
                        for x in range(0, len(sub_list), 1):
                            if x != 0:
                                start = int(sub_list[x]
                                            - ((sub_list[x] - sub_list[x-1])/2))
                            else:
                                start = int(sub_list[x]) - 1
                            if x != len(sub_list) - 1:
                                end = int(sub_list[x]
                                          + ((sub_list[x+1] - sub_list[x])/2))
                            else:
                                end = int(sub_list[x]) + 1

                            out.write(chromosome + '\t' + str(start) + '\t'
                                      + str(end) + '\t' + type1 + side
                                      + str(peaks) + '_' + str(start) + '_'
                                      + str(end) + '_A' + '\t' + str(peaks) + '\n')
                            for base in range(start, end+1):
                                peak_areas[chromosome][side][base] = 1
                            peaks += 1
                    else:
                         single = True
                else:
                    single = True
                if single:
                    start = min(sub_list) - splice_site_width
                    end = max(sub_list) + splice_site_width
                    out.write(chromosome + '\t' + str(start) + '\t' + str(end)
                              + '\t' + type1 + side + str(peaks) + '_'
                              + str(start) + '_' + str(end) + '_A'
                              + '\t' + str(peaks) + '\n')
                    for base in range(start, end+1):
                        peak_areas[chromosome][side][base] = 1
                    peaks += 1

    return peaks, peak_areas

def get_alignment_direction(sam_file):
    direction_dict = {}
    direction_conversion = {}
    direction_conversion[('0','+')] = '+'
    direction_conversion[('0','-')] = '-'
    direction_conversion[('16','+')] = '-'
    direction_conversion[('16','-')] = '+'

    for line in open(sam_file):
        if line[0] != '@':
            a = line.strip().split('\t')
            read_name = a[0].split('_')[0]
            read_direction = a[1]
            if read_direction in ['0', '16']:
                for entry in a[10:]:
                    if 'ts:A:' in entry:
                        direction = entry.split('ts:A:')[1]
                        direction_dict[read_name] = direction_conversion[(read_direction, direction)]
    return direction_dict

def collect_chromosomes(isoform_psl,chromosomes):
    for line in open(isoform_psl):
        a = line.strip().split('\t')
        chromosome = a[13]
        chromosomes.add(chromosome)
    return chromosomes






def main():
    left_bounds = {}
    right_bounds = {}
    print('    parsing annotated splice sites')
    chromosome_list,left_bounds, right_bounds = parse_genome(genome_file, left_bounds, right_bounds)

    Left_Peaks = 0
    Right_Peaks = 0

    chromosome_list = collect_chromosomes(infile,chromosome_list)
    direction_dict = get_alignment_direction(sam_file)

    peak_areas = {}
    print(sorted(list(chromosome_list)))
    for chromosome in sorted(list(chromosome_list)):
        print('    now processing ', chromosome)
        print('    collecting reads')
        histo_left_bases, histo_right_bases, \
        histo_coverage = collect_reads(infile,sam_file,direction_dict,chromosome)


        peak_areas[chromosome] = {}
        peak_areas[chromosome]['l'] = {}
        peak_areas[chromosome]['r'] = {}
        if chromosome not in left_bounds:
            left_bounds[chromosome]={}
            left_bounds[chromosome]['5'], left_bounds[chromosome]['3'] = [], []
        if chromosome not in right_bounds:
            right_bounds[chromosome] = {}
            right_bounds[chromosome]['5'], right_bounds[chromosome]['3'] = [], []
        if 'g' in refine:
            Left_Peaks_old = Left_Peaks
            Right_Peaks_old = Right_Peaks
            Left_Peaks, peak_areas = make_genome_bins(left_bounds[chromosome], 'l',
                                                      Left_Peaks, chromosome, peak_areas)
            Right_Peaks, peak_areas = make_genome_bins(right_bounds[chromosome], 'r',
                                                       Right_Peaks, chromosome, peak_areas)

            print('    parsed annotation-based splice-sites',
                  Left_Peaks - Left_Peaks_old,
                  Right_Peaks - Right_Peaks_old)
        Left_Peaks_old = Left_Peaks
        Right_Peaks_old = Right_Peaks
        Left_Peaks, peak_areas = find_peaks(histo_left_bases[chromosome],
                                            out, Left_Peaks, True, cutoff,
                                            histo_coverage, 'l',
                                            peak_areas, chromosome)
        Right_Peaks, peak_areas = find_peaks(histo_right_bases[chromosome],
                                             out, Right_Peaks, False, cutoff,
                                             histo_coverage, 'r',
                                             peak_areas, chromosome)


        print('    detected read-based splice-sites',
              Left_Peaks - Left_Peaks_old,
              Right_Peaks - Right_Peaks_old)



main()
