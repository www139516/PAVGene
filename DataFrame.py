import pandas as pd 
from utils import print_line_split


class DataFrame():
    
    '''
    depth_id file format
    chr    pos                   id
0    1  44289  gene:Zm00001d027230
1    1  44290  gene:Zm00001d027230
2    1  44291  gene:Zm00001d027230
3    1  44292  gene:Zm00001d027230
4    1  44293  gene:Zm00001d027230

    filtered_depth_fpath
    chr    pos  depth
0    1  44352     11
1    1  44353     11
2    1  44355     10
3    1  44356     10
4    1  44357     12

    merged_df
    chr     pos     id      chr_pos depth   not_covered
1       44289   gene:Zm00001d027230     1_44289         True
1       44290   gene:Zm00001d027230     1_44290         True
1       44291   gene:Zm00001d027230     1_44291         True
1       44292   gene:Zm00001d027230     1_44292         True
1       44293   gene:Zm00001d027230     1_44293         True
1       44294   gene:Zm00001d027230     1_44294         True
1       44295   gene:Zm00001d027230     1_44295         True

    '''
    
    
    def __init__(self, depth_id_fpath, filtered_depth_fpath) -> None:
        
        print_line_split('Initializing the coverage calculation module...')
        
        self.depth_id_df = pd.read_csv(depth_id_fpath, sep='\t',names=['chr', 'pos', 'id'], usecols=[0, 1, 2], dtype={'chr': 'str'})
        # print( self.depth_id_df.head() )
        self.filtered_depth_df = pd.read_csv(filtered_depth_fpath, sep='\s+', names=[
                                             'chr', 'pos', 'depth'], usecols=[0, 1, 2], dtype={'chr': 'str'})
        # print(self.filtered_depth_df.head())
        self.merged_df = None
        self.df_coverage = None

    def merge_data(self, merged_fpath):
        
        print_line_split("Start merging information...")
        
        self.depth_id_df.chr = self.depth_id_df.chr.astype(str)
        self.depth_id_df.pos = self.depth_id_df.pos.astype(str)
        self.filtered_depth_df.chr = self.filtered_depth_df.chr.astype(str)
        self.filtered_depth_df.pos = self.filtered_depth_df.pos.astype(str)
        
        self.depth_id_df['chr_pos'] = self.depth_id_df.chr + '_' + self.depth_id_df.pos
        self.filtered_depth_df['chr_pos'] = self.filtered_depth_df['chr'] + '_' + self.filtered_depth_df['pos']
        self.filtered_depth_df.drop(['chr', 'pos'], axis=1, inplace=True)
        
        self.merged_df = pd.merge(self.depth_id_df, self.filtered_depth_df, on='chr_pos', how='left')
        self.merged_df['not_covered'] = self.merged_df.isnull().T.any()
        
        self.merged_df.to_csv(merged_fpath, sep='\t', index=None)
        
        print_line_split('Finished initialization.')
        

    def calculate_id_coverage(self, out_cov_fpath):
        
        print_line_split('Calculating coverage of each id...')
        
        self.merged_df.drop('depth', axis=1, inplace=True)
        print(self.merged_df.head())
        self.df_coverage = self.merged_df.groupby('id').not_covered.value_counts()
        print(self.df_coverage.head())
        
        id_cov_lst = list()
        line = ''
        
        with open(out_cov_fpath, 'w') as fo:
            for k_v in self.df_coverage.items():
                if not id_cov_lst or id_cov_lst[0][0][0] == k_v[0][0]:
                    id_cov_lst.append(k_v)
                elif len(id_cov_lst) == 1 and id_cov_lst[0][0][1] == True:
                    line = '{id}\t{cov}'.format(id=id_cov_lst[0][0][0], cov=100)
                    print(line)
                    id_cov_lst = list()
                    fo.writelines("{out_line}\n".format(out_line=line))

                elif len(id_cov_lst) == 1 and id_cov_lst[0][0][1] == False:
                    line = '{id}\t{cov}'.format(id=id_cov_lst[0][0][0], cov=0)
                    print(line)
                    id_cov_lst = list()
                    fo.writelines("{out_line}\n".format(out_line=line))

                if len(id_cov_lst) > 1:
                    line = '{id}\t{cov}'.format(
                        id=id_cov_lst[0][0][0], cov=id_cov_lst[0][1] / (id_cov_lst[0][1] + id_cov_lst[1][1]) * 100)
                    print(line)
                    id_cov_lst = list()
                    fo.writelines("{out_line}\n".format(out_line=line))


            if len(id_cov_lst) == 1 and id_cov_lst[0][0][1] == True:
                line = '{id}\t{cov}'.format(id=id_cov_lst[0][0][0], cov=100)
                print(line)
                id_cov_lst = list()
            elif len(id_cov_lst) == 1 and id_cov_lst[0][0][1] == False:
                line = '{id}\t{cov}'.format(id=id_cov_lst[0][0][0], cov=0)
                print(line)
                id_cov_lst = list()
            elif len(id_cov_lst) > 1:
                line = '{id}\t{cov}'.format(
                    id=id_cov_lst[0][0][0], cov=id_cov_lst[0][1] / (id_cov_lst[0][1] + id_cov_lst[1][1]) * 100)
                print(line)
                id_cov_lst = list()
            else:
                print("Finished calculating coverage.")
            
            
            fo.writelines("{out_line}\n".format(out_line = line))
            
            print_line_split('Finished calculating coverage.')
            
    @classmethod
    def stat(cls, cov_fath, out_fpath):
        print_line_split('Start statistical analysis...')
        df_cov = pd.read_csv(cov_fath, sep='\t', names=['id', 'coverage'])
        df_cov['cov_bin_counts'] = pd.cut(df_cov['coverage'], [0, 20, 40, 60, 80, 100], labels=['0 to 20', '21 to 40', '41 to 60',
                                                                                                '61 to 80', '81 to 100'])
        
        df_cov_counts_by_bin = df_cov.cov_bin_counts.value_counts()
        with open(out_fpath, 'w') as fo:
            print_line_split('Writing statistical results to the {}'.format(out_fpath))
            fo.writelines('Bins\tCounts\n')

            for k, v in df_cov_counts_by_bin.items():
                
                print('{key} >>> {value}'.format(
                    key=k,
                    value=v
                ))
                fo.writelines('{key}\t{value}\n'.format(
                    key=k,
                    value=v
                ))   
        
        print_line_split('Finished statistical analysis.')
    
    @classmethod
    def merge_id_pos_cov(cls, pos_fpath, cov_fpath, out_fpath):
        df_pos = pd.read_csv(pos_fpath, names=['chr', 'start', 'end', 'id'], sep='\t|\s+')
        df_cov = pd.read_csv(cov_fpath, names=['id', 'coverage'], sep='\t|\s+')
        df_merge = pd.merge(df_pos, df_cov, on='id')
        df_merge.to_csv(out_fpath, sep='\t', index=None)   

                
                
    
