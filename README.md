# WooVariant
This is a small python program for SNV calling based on pysam package. Prokaryotic might be preferred.

A sorted BAM file and a reference sequence are required. Output format is VCF.

In the output VCF, you will get information including IS (Maximum number of reads supporting a mutation), DP (Raw reads depth), AF (Mutation abundance), and TYPE (Variant types, including SNP and indel).

The major usage is:

python woovariant.py -r reference.fasta -i query.sorted.bam

Enter -h or --help would print following messages:

-r --reference    Reference FASTA file

-i --input_file   Input bam file

-o --output_file  Output vcf file prefix (Default: input name and file)

-s --is_min       Minimum threshold for maximum number of reads supporting a variant (Default: 10)

-d --depth        Minimum threshold for raw read depth (Default: 20)

-m --ma           Minimum threshold of mutation abundance as percentage (Default: 5)

-a --only_first   Set to output only when majority of a position is variant
