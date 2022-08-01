# PAVGene
Calculating gene coverage using aligned file.
## Dependencies
python >= 3

pandas

samtools

sambamba (optional)

## Usage:
Get help

python /your/path/to/PAVGene/main.py --help

$ python /share//nas1/smp1/tifu/scripts/PAVGenV2/main.py --help

usage: main.py [-h] -g GTF [-d DEPTH] [-t IS_TRANSCRIPT] [-f FILTER_DEPTH]
               [-s SAMBAMBA_DIR] [-m SAMTOOLS_DIR] [-z SAMTOOLS_ONLY]
               [-x THREADS] [-b BAM] [-i SORTED_BAM] [-o OUT_DIR] [-p PREFIX]

This software calculates each gene coverage using sequencing data and the gtf
file.

optional arguments:

  -h, --help            show this help message and exit
  
  -g GTF, --gtf GTF     The gtf file used for locating gene postions
  
  -d DEPTH, --depth DEPTH
                        The path of the file containing depth information for
                        each site.
                        
  -t IS_TRANSCRIPT, --is_transcript IS_TRANSCRIPT
                        Calculate coverage based on transcript. Default is
                        false, parse "True" if you want to change.
                        
  -f FILTER_DEPTH, --filter_depth FILTER_DEPTH
                        The depth lower than this will be filtered. Default is
                        3
                        
  -s SAMBAMBA_DIR, --sambamba_dir SAMBAMBA_DIR
                        The path of the directory containing the sambamba
                        exectutive files. Please specify if it is not in $PATH
                        
  -m SAMTOOLS_DIR, --samtools_dir SAMTOOLS_DIR
                        The path of the directory containing the samtools,
                        exectutive files. Please specify if it is not in $PATH
                        
  -z SAMTOOLS_ONLY, --samtools_only SAMTOOLS_ONLY
                        Using samtools to sort bam file. default is False, set
                        "True" if you don't install sambamba
                        
  -x THREADS, --threads THREADS
                        The threads used in this program when needed. Default
                        is 8
                        
  -b BAM, --bam BAM     The path of the bam file, if give, the software will
                        sort it and calculate the depth of each location.
                        default is None.
                        
  -i SORTED_BAM, --sorted_bam SORTED_BAM
                        The sorted bam file, if give, the software will
                        calculate the depth of each location. Default is None
                        
  -o OUT_DIR, --out_dir OUT_DIR
                        The output directory where the output files will be
                        stored, default is the current directory.
                        
  -p PREFIX, --prefix PREFIX
                        The prefix of output dir.
                        

## use with bam
python /share/nas1/smp1/tifu/scripts/PAVGenV2/main.py -g path/to/gtf -b bam -p prefix
