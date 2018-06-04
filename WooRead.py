import argparse
import os


parser = argparse.ArgumentParser(description='WooRead - A tool to filt unequal reads from FastQ Screen. '
                                             'For research team of Dr. Liu in iSynBio, SIAT only. '
                                             'Usage: python WooRead.py '
                                             '-1 file/name1.tagged_filter.fastq -2 file/name2.tagged_filter.fastq '
                                             'Output: file/name_1.fastq file/name_2.fastq '
                                             'Copyleft: Jacky Woo from ZHK Research team, iSynBio, SIAT.')
parser.add_argument('-1', '--input_file_1', type=str, help='Input FASTQ file R1')
parser.add_argument('-2', '--input_file_2', type=str, help='Input FASTQ file R2')
parser.add_argument('-o', '--output_file', type=str, required=False,
                    help='Output FASTQ file prefix (Default: input name and file)')

args = parser.parse_args()
if args.output_file is not None:
    outfile_1 = args.output_file + '_1.fastq'
    outfile_2 = args.output_file + '_2.fastq'
else:
    outfile_1 = args.input_file_1.rsplit('.tagged_filter.fastq')[0] + '_1.fastq'
    outfile_2 = args.input_file_2.rsplit('.tagged_filter.fastq')[0] + '_2.fastq'
print 'Thank you for use WooRead! You set input %s, %s and output %s, %s' % (
    args.input_file_1, args.input_file_2, outfile_1, outfile_2)
main_dict = {}
outlist = []
seqfile = open(args.input_file_1)
counter = 0
count = 0
outer = ['', '', '']
out = ''
for line in seqfile:
    count += 1
    if count == 1:
        out = line.split(' ')[0]
        outer[0] = line.split('#')[0]
    elif count == 2:
        outer[1] = line.split('\n')[0]
    elif count == 4:
        outer[2] = line.split('\n')[0]
        counter += 1
        main_dict.update({out: [outer, None]})
        outlist += [(counter, out)]
        count = 0
        outer = ['', '', '']
seqfile.close()
seqfile = open(args.input_file_2)
count = 0
outer = ['', '', '']
out = ''
for line in seqfile:
    count += 1
    if count == 1 and line.split(' ')[0] in main_dict:
        out = line.split(' ')[0]
        outer[0] = line.split('#')[0]
    elif count == 2 and outer[0]:
        outer[1] = line.split('\n')[0]
    elif count == 4:
        if outer[0]:
            outer[2] = line.split('\n')[0]
            main_dict[out][1] = outer
        count = 0
        outer = ['', '', '']
seqfile.close()
outlist.sort()
outfh_1 = open(outfile_1, 'w')
outfh_2 = open(outfile_2, 'w')
for item in outlist:
    line = main_dict[item[1]]
    if line[1]:
        outfh_1.write(line[0][0] + '\n' + line[0][1] + '\n+\n' + line[0][2] + '\n')
        outfh_2.write(line[1][0] + '\n' + line[1][1] + '\n+\n' + line[1][2] + '\n')
outfh_1.close()
outfh_2.close()
print 'WooRead finished!'
