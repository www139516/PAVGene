import os
import sys
from FileProcessor import FileProcessor
from DataFrame import DataFrame
from utils import print_line_split
from sys_cmd import get_sorted_aligned_file_using_samtools
from sys_cmd import get_sorted_aligned_file_using_sambamba
from sys_cmd import get_pos_depth


class FileOperator:
    
    def __init__(self, 
                 gtf_fpath, 
                 sambamba_dir='', 
                 samtools_dir='',  
                 samtools_only='', 
                 bam_fpath='', 
                 sorted_bam_fpath='', 
                 depth_fpath='', 
                 screen_depth=3, 
                 gene_or_trans='gene', 
                 threads=8, 
                 out_dir='', 
                 pre='') -> None:
        
        print_line_split('Initializing FileOperator...', '#')
        
        print_line_split('Setting the depth for screening to {}'.format(screen_depth))
        self.screen_depth = screen_depth
        print_line_split('Setting the thread number to {}'.format(threads))
        self.threads = threads
        
        self.samtools_only = (True if samtools_only.startswith('t') else False)
        if samtools_only:
            print_line_split('Will use samtools for sorting and depth calculation.')
        else:
            print_line_split('Will use sambamba for sorting and samtools for depth calculation.')
        
        self.gtf_fpath = os.path.abspath(gtf_fpath)
        self.gtf_fname = os.path.basename(self.gtf_fpath)
        
        self.cur_dir = os.getcwd()
        self.pre = pre
        self.out_dir = ('{}_out'.format(self.pre) if not out_dir else out_dir) 
        
        if '/' in self.out_dir:
            self.out_dpath = os.path.abspath(self.out_dir)
        else:
            self.out_dpath = os.path.join(self.cur_dir, self.out_dir)
        if not os.path.exists(self.out_dpath):
            os.mkdir(self.out_dpath)
        print_line_split('The output files will be in the {}'.format(self.out_dpath))
        
        self.gene_or_trans = gene_or_trans
        print_line_split('Will calculate coverage based on per {}'.format(self.gene_or_trans))
        
        if bam_fpath:
        
            self.bam_fpath = os.path.abspath(bam_fpath)
            self.bam_fname = os.path.basename(self.bam_fpath)
            if not self.pre:
                self.pre = os.path.splitext(self.bam_fname)[0]
            self.sorted_bam_fname = self.pre + '.sorted.bam'
            # print(self.sorted_bam_fname)
            self.sorted_bam_fpath = os.path.join(self.out_dpath, self.sorted_bam_fname)
            
        if sorted_bam_fpath:
            if bam_fpath:
                print_line_split('unsorted bam file will not be used since you have given the sorted bam file.')
                
            self.bam_fpath = ''
            self.bam_fname = ''
            
            self.sorted_bam_fpath = os.path.abspath(sorted_bam_fpath)
            self.sorted_bam_fname = os.path.basename(self.sorted_bam_fpath)
            if not self.pre:
                self.pre = os.path.splitext(self.sorted_bam_fname)[0]
        
        if depth_fpath:
            if bam_fpath:
                print_line_split(
                    'Unsorted bam file will not be used since you have given the depth file.')
            if sorted_bam_fpath:
                print_line_split(
                    'Sorted bam file will not be used since you have given the depth file.')

            self.pos_depth_fpath = os.path.abspath(depth_fpath)
            self.pos_depth_fname = os.path.basename(self.pos_depth_fpath)
            
            self.bam_fpath = ''
            self.bam_fname = ''
            self.sorted_bam_fpath = ''
            self.sorted_bam_fname = ''
            
            if not self.pre:
                self.pre = os.path.splitext(self.pos_depth_fname)[0]
        else:
            self.pos_depth_fpath = os.path.join(
                self.out_dpath, self.pre + '.bed.depth')
            
        self.sambamba_dpath = ( os.path.abspath(sambamba_dir) if sambamba_dir else '' )
        self.samtools_dpath = ( os.path.abspath(samtools_dir) if samtools_dir else '' )
        
        
            
        self.gene_bed_fpath = os.path.join(
            self.out_dpath, os.path.splitext(self.gtf_fname)[0] + '.gene.bed') ## The file with bed format containing gene id
        self.trans_bed_fpath = os.path.join(
            self.out_dpath, os.path.splitext(self.gtf_fname)[0] + '.trans.bed')  # The file with bed format containing transcript id
        
        self.gene_pos_fpath = os.path.join(
            self.out_dpath, os.path.splitext(self.gtf_fname)[0] + '.gene.pos')  ## The file with collapsed gene postion
        self.trans_pos_fpath = os.path.join(
            self.out_dpath, os.path.splitext(self.gtf_fname)[0] + '.trans.pos') ## The file with collapsed transcript position
        
        self.merged_fpath = os.path.join(
            self.out_dpath, self.pre + '.merged.pos.depth.txt')  # The file merged pos and depth
                    
        self.pos_depth_filter_fpath = os.path.join(
            self.out_dpath, self.pre + '.bed.filtered.depth')
        
        self.gene_cov_fpath = os.path.join(
            self.out_dpath, self.pre + '.gene.cov')  # The file with gene coverage
        self.trans_cov_fpath = os.path.join(
            self.out_dpath, self.pre + '.trans.cov')  # The file with transcript coverage
        self.stats_fpath = os.path.join(self.out_dpath, self.pre + '.stats.txt') # The file with stats info
        self.res_fpath = os.path.join(
            self.out_dpath, self.pre + 'final.txt')
        
    def bam2depth(self) -> str:
        if self.bam_fpath:
            print_line_split('Start sorting bam...')
            self.__bam2sortedbam()
            self.__sortedbam2depth()
            print_line_split('Sorting complete.')
        elif self.sorted_bam_fpath:
            print_line_split('Start generating depth file...')
            self.__sortedbam2depth()
            print_line_split('Depth file is completed.')
        elif self.pos_depth_fpath:
            print_line_split('Depth file is given.')
        else:
            print_line_split('Data file is missing, please check.')
            sys.exit(1)
        return self.pos_depth_fpath
    
    def __bam2sortedbam(self):
        if self.samtools_only:
            get_sorted_aligned_file_using_samtools(self.bam_fpath, self.sorted_bam_fpath, self.samtools_dpath, self.threads)
        else:
            get_sorted_aligned_file_using_sambamba(self.bam_fpath, self.sorted_bam_fpath, self.sambamba_dpath, self.threads)
            
    
    def __sortedbam2depth(self):
        
        if 'gene' in self.gene_or_trans:
            get_pos_depth(self.sorted_bam_fpath, self.samtools_dpath, self.gene_bed_fpath, self.pos_depth_fpath)
        else:
            get_pos_depth(self.sorted_bam_fpath, self.samtools_dpath, self.trans_bed_fpath, self.pos_depth_fpath)
        
        
    def gtf2bed(self):
        out_fpath = ''
        if 'gene' in self.gene_or_trans:
            self.__gtf2bed_with_genes()
            out_fpath = self.gene_bed_fpath
        else:
            self.__gtf2bed_with_transcripts()
            out_fpath = self.trans_bed_fpath
        return out_fpath
    
    def __gtf2bed_with_genes( self ) -> str:
        
        print_line_split('Converting the gtf file to bed file with gene id...')
        
        fp = FileProcessor(self.gtf_fpath, 'gene')
        fp.gtf2bed( self.gene_bed_fpath )
        
        print_line_split('Finished converting to bed file with gene id.')     
   
    def __gtf2bed_with_transcripts(self) -> str:
        
        print_line_split('Coverting the gtf file to bed file with transcript id...')
        
        fp = FileProcessor(self.gtf_fpath, 'transcript')
        fp.gtf2bed(self.trans_bed_fpath)
        
    def bed2pos(self):
        
        out_fpath = ''
        if 'gene' in self.gene_or_trans:
            self.__bed2mapping_gene_file()
        else:
            self.__bed2mapping_trans_file()                
        
        if 'gene' in self.gene_or_trans:
            out_fpath = self.gene_pos_fpath
        else:
            out_fpath = self.trans_pos_fpath
        return out_fpath
        
    
    def __bed2mapping_gene_file(self):
        
        print_line_split('collapsing gene bed file to mapping file with each position...')
        FileProcessor.bed_to_mapping(self.gene_bed_fpath, self.gene_pos_fpath)
        return self.gene_pos_fpath
    
    def __bed2mapping_trans_file(self):
        
        print_line_split(
            'collapsing transcript bed file to mapping file with each position...')
        FileProcessor.bed_to_mapping(
            self.trans_bed_fpath, self.trans_pos_fpath)
        return self.trans_pos_fpath
    
    def filter_depth(self):
        
        print_line_split('Filter the depths that are too lower...')
        
        FileProcessor.depth_filter(self.pos_depth_fpath, self.pos_depth_filter_fpath, self.screen_depth)
        return self.pos_depth_filter_fpath
    
    
    def calculate_id_cov(self):
        
        print_line_split('Start calculating id coverage...')
        
        if self.gene_or_trans.lower() == 'gene':
            print_line_split('Calculating coverage for each gene...')    
            df_opt = DataFrame(self.gene_pos_fpath, self.pos_depth_filter_fpath)
            df_opt.merge_data(self.merged_fpath)
            df_opt.calculate_id_coverage(self.gene_cov_fpath)
            ret = self.gene_cov_fpath
        elif 'trans' in self.gene_or_trans.lower():
            print_line_split('Calculating coverage each transcript ...')
            df_opt = DataFrame(self.trans_pos_fpath, self.pos_depth_filter_fpath)
            df_opt.merge_data(self.merged_fpath)
            df_opt.calculate_id_coverage(self.trans_cov_fpath)
            ret = self.trans_cov_fpath
        else:
            print_line_split('Please enter the right request gene or transcript...', 'X')
            sys.exit(1)
        
        print_line_split('Finished coverage calculating. Please open {} to check the result files.'.format(self.out_dpath))
            
        return ret
    
    def generate_stats(self):
        if 'gene' in self.gene_or_trans:
            DataFrame.stat(self.gene_cov_fpath, self.stats_fpath)
        else:
            DataFrame.stat(self.trans_cov_fpath, self.stats_fpath)
            
    def generate_final_rep(self):
        if 'gene' in self.gene_or_trans:
            DataFrame.merge_id_pos_cov(self.gene_bed_fpath, self.gene_cov_fpath, self.res_fpath)
        else:
            DataFrame.merge_id_pos_cov(self.trans_bed_fpath, self.trans_cov_fpath, self.res_fpath)

if __name__ == '__main__':
    fo = FileOperator('./data/gtf.txt', depth_fpath='./data/pos.depth.txt')
    out_fpath = fo.__gtf2bed_with_genes()
    out_tran_fpath = fo.__gtf2bed_with_transcripts()
    out_gene_pos_fpath = fo.__bed2mapping_gene_file()
    out_trans_pos_fpath = fo.__bed2mapping_gene_file()
    out_pos_depth_filter_fpath = fo.filter_depth()
    out_gene_cov_fpath = fo.calculate_id_cov()
    print(out_fpath)
    print(out_tran_fpath)
    print(out_gene_pos_fpath)
    print(out_pos_depth_filter_fpath)
    print(out_gene_cov_fpath)
    
