import subprocess
import os
from utils import assert_existed, assert_samtools
from utils import assert_bam_file
from utils import assert_samtools
from utils import assert_sambamba
from utils import print_line_split

def get_sorted_aligned_file_using_sambamba(bam_fpath, out_sorted_bam, sambamba_dir='', threads=8):
    
    assert_sambamba(sambamba_dir)
    
    bam_fpath = os.path.abspath(bam_fpath)
    assert_bam_file(bam_fpath)
    
    if sambamba_dir:
        sambamba_dir = os.path.abspath(sambamba_dir)
        assert_existed(sambamba_dir)
    
    sambamba = os.path.join(sambamba_dir, 'sambamba')
    cmd = '{sambamba} sort -t {threads} -o {out_bam} -p {bam}'.format(sambamba=sambamba, 
                                                                      threads = threads,
                                                                      out_bam=out_sorted_bam,
                                                                      bam=bam_fpath)
    
    print_line_split('Sorting bam file using sambamba...')
    print(cmd)
    if subprocess.check_call(cmd, shell=True) != 0:
        raise SystemCommandError
    print_line_split('Finished sorting bam.')
    

def get_sorted_aligned_file_using_samtools(bam_fpath, out_bam_fpath, samtools_dir='', threads=8):
    
    bam_fpath = os.path.abspath(bam_fpath)
    
    assert_bam_file(bam_fpath)
    if samtools_dir:
        samtools_dir = os.path.abspath(samtools_dir)
        assert_existed(samtools_dir)
        
    samtools = os.path.join(samtools_dir, 'samtools')
    cmd = '{samtools} sort -o {out} -@ {threads} {in_bam}'.format(samtools=samtools,
                                                                  out=out_bam_fpath,
                                                                  threads = threads,
                                                                  in_bam=bam_fpath)
    
    print_line_split('Sorting bam file using samtools...')
    print(cmd)
    if subprocess.check_call(cmd, shell=True) != 0:
        raise SystemCommandError
    print_line_split('Finished sorting bam.')
    
 
def get_pos_depth(sorted_bam_fpath, samtools_dir, bed_fpath, out_depth_fpath):
    print_line_split('Checking files...')
    
    assert_existed(sorted_bam_fpath)
    assert_existed(bed_fpath)
    assert_samtools(samtools_dpath=samtools_dir)
    
    samtools = os.path.join(samtools_dir, 'samtools')

    cmd = '{samtools} depth -b {bed_fpath} {bam} >{out}'.format(samtools=samtools,
                                                         bed_fpath=bed_fpath,
                                                         bam=sorted_bam_fpath,
                                                         out=out_depth_fpath)
    if subprocess.check_call(cmd, shell=True) != 0:
        raise SystemCommandError
    
