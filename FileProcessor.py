import re
import sys
from utils import print_line_split

class FileProcessor:
    
    def __init__( self, fin, gene_or_trans ) -> None:
        """Used for processing files and convert input file to the required file format

        Args:
            fin (str): gtf file path
            gene_or_trans (str): gene or trans 
        """
        
        print_line_split('Initializing FileProcesser...')
        
        self.in_gtf = fin
        self.out_bed = None
        self.gene_or_trans = gene_or_trans
        
    def __get_bed_line_from_gtf( self, line ) -> str:
        """get info required for bed file from gtf file.

        Args:
            line (str): line read from gtf file
            gene_or_trans: transcript or gene
        Returns:
            str: trans_id or gene_id from gtf file
        """
        out_line_lst = line.split()
        id_content = ' '.join(out_line_lst[8:])
        
        # print(id_content)
        
        if 'trans' in self.gene_or_trans.lower():
            name_pat = re.compile(r'(?<=transcript_id ").*?(?=")')
            # id = re.findall(name_pat, out_line_lst[-2])[0]
        else:
            name_pat = re.compile(r'(?<=gene_id ").*?(?=")')
            # id = re.findall(name_pat, out_line_lst[-1])[0]
        id = re.findall(name_pat, id_content)[0]
        if (len(out_line_lst) < 1):
            return ''
        
        return out_line_lst[0] + ' ' + out_line_lst[3] + ' ' + out_line_lst[4] + ' ' + id + '\n'

    def gtf2bed( self, fout_bed ):
        self.out_bed = fout_bed
        with open( self.in_gtf, 'r' ) as fin:
            with open(self.out_bed, 'w') as fout:
                for line in fin:
                    out_line = ''
                    out_line = self.__get_bed_line_from_gtf(line)
                    fout.writelines(out_line)

    
    @staticmethod
    def bed_to_mapping(fin_bed, fout_mapping):
        id_pos = list() # storing genne/trans id for each position [chr, set(start, .., end), id]
        line_lst = list() 
        print_line_split('Reading ids from the first block of bed file, each "*" represents 1000 lines.')
        with open(fin_bed, 'r') as fin:
            with open( fout_mapping, 'w' ) as fo:
                
                count = 0
                cnt_star = 0
                print_line_split(
                    'Reading ids from the next block of bed file, each "*" represents 1000 lines.')
                for line in fin:
                    line_lst = line.split()
                    # print(line)
                    
                    if not line_lst[2] in id_pos:
                        
                        if id_pos:
                            for i in id_pos[1]:
                                fo.writelines("{chr}\t{pos}\t{id}\n".format(chr=id_pos[0], pos=i, id=id_pos[-1]))
                        id_pos = [line_lst[0], set(i for i in range(int(line_lst[1]), int(line_lst[2])+1)), line_lst[-1]]
                        # print(id_pos)
                        count += 1
                        if (count % 1000 == 0):
                            print('*', end=' ')
                            sys.stdout.flush()
                            cnt_star += 1
                            if (cnt_star % 20 == 0):
                                print('', end='\n')
                    else:
                        
                        id_pos[1] = id_pos[1] | set(i for i in range(int(line_lst[1]), int(line_lst[2])+1))
        print('\n')
        print_line_split('Finished reading id information.')
                        
    @staticmethod
    def depth_filter(cov_in_fpath, filtered_cov_fpath, depth_threshold=3):
        line_lst = list()
        count_total = 0
        count_fil = 0
        count_retain = 0
        print_line_split('Start filtering the coverage file, the depth lower than {} will be discarded'.format(depth_threshold))
        
        with open (cov_in_fpath, 'r') as fin:
            with open( filtered_cov_fpath, 'w') as fout:
                for line in fin:
                    line_lst = line.split()
                    count_total += 1
                    if int(line_lst[-1]) >= depth_threshold:
                        fout.writelines(line)
                        count_retain += 1
                    else:
                        count_fil += 1
                        
        print_line_split(
            "Finished filtering, {fil} lines were discarded, \nkept {retain} lines out of {total} lines. \nThe percentage is {rate}%".format(
                fil = count_fil, 
                retain = count_retain,
                total = count_total,
                rate = count_retain / count_total * 100
            ))
