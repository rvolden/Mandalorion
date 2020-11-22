import sys
import os

path=sys.argv[1]


combined_consensus_file = open(path + '/Isoform_Consensi.fasta', 'w')

for file1 in os.listdir(path + '/mp/'):
    try:
        counter=int(file1.split('.fasta')[0])
        print(counter,file1)
        for lines in open(path + '/mp/' + file1):
            combined_consensus_file.write(lines)
    except:
        print(file1)
combined_consensus_file.close()
