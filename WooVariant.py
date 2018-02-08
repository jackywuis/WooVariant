import pysam
import argparse


parser = argparse.ArgumentParser(description='WooVariants - A tool to find SNPs and indels in prokaryotic BAM file')
parser.add_argument('-r', '--reference', type=str, help='Reference FASTA file')
parser.add_argument('-i', '--input_file', type=str, help='Input bam file')
parser.add_argument('-o', '--output_file', type=str, required=False,
                    help='Output vcf file prefix (Default: input name and file)')
parser.add_argument('-s', '--is_min', type=int, default=10, required=False,
                    help='Minimum threshold for maximum number of reads supporting a variant (Default: 10)')
parser.add_argument('-d', '--depth', type=int, default=20, required=False,
                    help='Minimum threshold for raw read depth (Default: 20)')
parser.add_argument('-m', '--ma', type=int, default=5, required=False,
                    help='Minimum threshold of mutation abundance as percentage (Default: 5)')
parser.add_argument('-a', '--only_first', action='store_false', required=False,
                    help='Set to output only when majority of a position is variant')

args = parser.parse_args()
template = {}
infh = open(args.reference, 'r')
for line in infh:
    if '>' in line:
        node = line.split('>')[1].split('\n')[0].split(' ')[0]
        template.update({node: ''})
    else:
        template[node] += line.split('\n')[0]
infh.close()
if args.output_file:
    output_file = args.output_file + '_woo.vcf'
else:
    output_file = args.input_file.rsplit('.bam')[0] + '_woo.vcf'
print 'Thank you for use WooVariant! You set reference %s, input %s, output %s, minimum read number %s, '\ 
  'minmum read depth %s, minimum mutation abundance %s' % (
    args.reference, args.input_file, output_file, args.is_min, args.depth, args.ma)
outfh = open(output_file, 'w')
outfh.write('##fileformat=VCFv4.1\n##source=WooNTvariant\n##reference=file://' + args.reference +
            '\n##FILTER=<ID=IS,Number=1,Type=Integer,Description="Maximum number of reads supporting a mutation>"\n' +
            '##FILTER=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">\n' +
            '##FILTER=<ID=AF,Number=1,Type=Float,Description="Mutation abundance">\n' + 
            '##INFO=<ID=TYPE,Number=1,Type=String,Description="Variant type">\n' + 
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
samfile = pysam.AlignmentFile(args.input_file, 'rb')
for line in samfile.text.split('\n'):
    if 'SQ' in line:
        dele = [False, '', 0, 0]
        mx = 0
        for row in samfile.pileup(line.split('\t')[1].split(':')[1], 0, int(line.split('\t')[2].split(':')[1])):
            if row.n >= args.is_min:
                counter = {}
                for read in row.pileups:
                    if read.indel > 0:
                        position = read.alignment.query_sequence[read.query_position]
                        count = 0
                        while count < read.indel:
                            count += 1
                            position += read.alignment.query_sequence[read.query_position + count]
                    elif not read.query_position:
                        position = 'D'
                    else:
                        position = read.alignment.query_sequence[read.query_position]
                    if position not in counter:
                        counter.update({position: 0})
                    counter[position] += 1
                count = 0
                sort_list = []
                for item in counter:
                    if (float(counter[item]) * 100 / float(row.n)) > args.ma:
                        count += 1
                        sort_list += [(counter[item], item)]
                if not sort_list:
                    sort_list.sort(reverse=True)
                    output_nt = ''
                    if (sort_list[0][1] != template[line.split('\t')[1].split(':')[1]][row.pos]) and (sort_list[0][0] >= args.depth):
                        output_nt = sort_list[0][1]
                        mx = sort_list[0][0]
                    elif (count > 1) and (sort_list[1][0] >= args.depth) and args.only_first:
                        output_nt = sort_list[1][1]
                        mx = sort_list[1][0]
                    if output_nt != 'D' and output_nt != '':
                        outfh.write(line.split('\t')[1].split(':')[1] + '\t' + str(row.pos + 1) + '\t.\t' +
                                    template[line.split('\t')[1].split(':')[1]][
                                        row.pos] + '\t' + output_nt + '\t.\tPASS\tIS=' + str(mx) + ';DP=' +
                                    str(row.n) + ';AF=%.4f' % (float(mx) / float(row.n)))
                        if len(output_nt) > 1:
                            outfh.write(';TYPE=indel\n')
                        else:
                            outfh.write(';TYPE=SNP\n')
                    if output_nt == 'D':
                        if not dele[0]:
                            dele[2] = row.pos
                        dele[1] += last_nt
                        dele[0] = True
                        if mx > dele[3]:
                            dele[3] = mx
                    elif dele[0]:
                        dele[1] += last_nt
                        outfh.write(
                            line.split('\t')[1].split(':')[1] + '\t' + str(dele[2]) + '\t.\t' + dele[1] + '\t' +
                            dele[1][0] + '\t.\tPASS\tIS=' + str(dele[3]) + ';DP=' + str(row.n) +
                            ';AF=%.4f;TYPE=indel\n' % (float(dele[3]) / float(row.n)))
                        dele = [False, '', 0, 0]
            last_nt = template[line.split('\t')[1].split(':')[1]][row.pos]
samfile.close()
outfh.close()
