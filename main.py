import argparse
import os
import sys
from FileOperator import FileOperator
from utils import assert_existed
from utils import assert_gtf_file
from utils import is_existed
from utils import print_line_split


def main():
    parser = argparse.ArgumentParser(
        description='This software calculates each gene coverage using sequencing data and the gtf file.')
    parser.add_argument('-g', '--gtf', 
                        help='The gtf file used for locating gene postions', required=True)
    parser.add_argument('-d', '--depth', 
                        help='The path of the file containing depth information for each site.', default='')
    parser.add_argument('-t', '--is_transcript', 
                        help='Calculate coverage based on transcript. Default is false, parse "True" if you want to change.', default='False')
    parser.add_argument('-f', '--filter_depth', 
                        help='The depth lower than this will be filtered. Default is 3', default=3)
    parser.add_argument('-s', '--sambamba_dir',
                        help='The path of the directory containing the sambamba exectutive files. Please specify if it is not in $PATH', default='')
    parser.add_argument('-m', '--samtools_dir',
                        help='The path of the directory containing the samtools, exectutive files. Please specify if it is not in $PATH', default='')
    parser.add_argument('-z', '--samtools_only', help='Using samtools to sort bam file. default is False, set "True" if you don\'t install sambamba', default='False')
    parser.add_argument('-x', '--threads', help='The threads used in this program when needed. Default is 8', default=8)
    
    parser.add_argument('-b', '--bam',
                        help='The path of the bam file, if give, the software will sort it and calculate the depth of each location. default is None.', default='')
    parser.add_argument('-i', '--sorted_bam', 
                        help='The sorted bam file, if give, the software will calculate the depth of each location. Default is None', default='')
    parser.add_argument('-o', '--out_dir', 
                        help='The output directory where the output files will be stored, default is the current directory.', default='')
    parser.add_argument('-p', '--prefix',
                        help='The prefix of output dir.', default='PAVGene')
    
    args = parser.parse_args()
    
    gtf_fpath = os.path.abspath(args.gtf)
    assert_gtf_file(gtf_fpath)
    
    depth_fpath = args.depth
    bam_fpath = args.bam
    sorted_bam_fpath = args.sorted_bam
    prefix = args.prefix
    threads = args.threads
    filter_depth = args.filter_depth
    sambamba_dir = args.sambamba_dir
    samtools_dir = args.samtools_dir
    out_dir = args.out_dir
    samtools_only = args.samtools_only.lower()
    is_transcript = args.is_transcript.lower()
    
    gene_or_trans = ('transcript' if 't' in is_transcript else 'gene')  
    
    fo = FileOperator(gtf_fpath=gtf_fpath,
                      bam_fpath=bam_fpath,
                      depth_fpath=depth_fpath,
                      sambamba_dir=sambamba_dir,
                      samtools_dir=samtools_dir,
                      out_dir=out_dir,
                      samtools_only=samtools_only,
                      sorted_bam_fpath=sorted_bam_fpath,
                      screen_depth=filter_depth,
                      gene_or_trans=gene_or_trans,
                      threads=threads,
                      pre=prefix)
    
    out_fpath = fo.gtf2bed()
    print_line_split('The bed file is {}'.format(out_fpath))
    out_pos_fpath = fo.bed2pos()
    print_line_split('The position file is {}'.format(out_pos_fpath))
    
    out_depth_fpath = fo.bam2depth()
    print_line_split('The depth file is {}'.format(out_depth_fpath))
    out_pos_depth_filter_fpath = fo.filter_depth()
    out_gene_cov_fpath = fo.calculate_id_cov()
    print_line_split('The filter depth file is {}'.format(
        out_pos_depth_filter_fpath))
    print_line_split('The coverage file is {}'.format(out_gene_cov_fpath))
    fo.generate_stats()
    fo.generate_final_rep()
    
    
if __name__ == '__main__':
    main()
