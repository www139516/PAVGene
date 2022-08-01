import os
import sys


def assert_existed(fpath):
    
    fpath = os.path.abspath(fpath)
    print_line_split('Checking the file {}...'.format(fpath))

    if not os.path.exists(fpath):
        print_line_split('The file {path} is not exist, please check.'.format(path=fpath), 'X')
        sys.exit(1)


def assert_gtf_file(fpath):
    
    fpath = os.path.abspath(fpath)
    assert_existed(fpath)
    
    if not 'gtf' in os.path.basename(fpath).lower():
        print_line_split('The gtf file provided might has problems, please check if "gtf" is in the file name', 'X')
        sys.exit(1)
        
        
def assert_bam_file(fpath):
    fpath = os.path.abspath(fpath)
    assert_existed(fpath)
    
    if not 'bam' in os.path.basename(fpath).lower():
        print_line_split('The bam file provided might has problems, please check if "bam" is in the file name.', 'X')
        sys.exit(1)
        

def assert_sorted_bam_file(fpath):
    fpath = os.path.abspath(fpath)
    assert_existed(fpath)
    
    if not 'sort' in os.path.basename(fpath) or 'bam' in os.path.basename(fpath):
        print_line_split('The sorted bam file provided might has problems, please check if ".sorted.bam" is the file name.', 'X')
        sys.exit(1)
    
def assert_execute_file(dpath, fname):
    flag = False
    if dpath:
        dpath = os.path.abspath(dpath)
        fpath = os.path.join(dpath, fname)
        flag = os.path.exists(fpath)
    else:
        env_dpaths = os.getenv('PATH').split(':')
        for env_dpath in env_dpaths:
            fpath = os.path.join(env_dpath, fname)
            print(fpath)
            if os.path.exists(fpath):
                flag = True
                break
    return flag
        
def assert_samtools(samtools_dpath):
    print_line_split('Checking samtools...')
    if not assert_execute_file(samtools_dpath, 'samtools'):
        print_line_split('The samtools is not in the designated paths, please ensure it is installed correctly.')
        sys.exit(1)
    else:
        print_line_split('The samtools is ready.')

def assert_sambamba(sambamba_dpath):
    print_line_split('Checking sambamba...')
    if not assert_execute_file(sambamba_dpath, 'sambamba'):
        print_line_split(
            'The sambamba is not in the designated paths, please ensure it is installed correctly.')
        sys.exit(1)
    else:
        print_line_split('The sambamba is ready.')

def is_existed(fpath):
    fpath = os.path.abspath(fpath)
    return os.path.exists(fpath)


def print_line_split(out_line='', symb='#') -> None:
    
    print('\n')
    print('\n')
    print(out_line)
    symbs = symb * 80
    print(symbs)
    print('\n')


