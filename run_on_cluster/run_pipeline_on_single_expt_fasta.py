#!/home/lucas_pkuhpc/anaconda3/bin/python
import sys
import os
import argparse
import subprocess
import time
import pprint



# parse cmd line args
ap = argparse.ArgumentParser()
ap.add_argument("-bn", "--base_name", required=True, help="basename for all output files")
ap.add_argument("-g", "--genome_fasta_file", required=True, help="FASTA file with the genome / amplified sequence")
ap.add_argument("-r1", "--read1", required=True,  help="FASTQ file first read")
ap.add_argument("-r2", "--read2", required=False, help="FASTQ file second read" , default=' ')
ap.add_argument("-fmh", "--find_mh", required=False, help="path to find_mh binary" , default='/lustre2/lucas_pkuhpc/bin/find_mh')
ap.add_argument("-cs", "--catch_signatures", required=False, help="cmd for catch_signatures.awk" , default='/usr/bin/gawk -v SZ=25 -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/catch_signatures.awk')
ap.add_argument("-gs", "--generate_signatures", required=False, help="cmd for generate_signatures.awk" , default='/usr/bin/gawk -v SZ=25 -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/generate_signatures.awk')
ap.add_argument("-nc", "--normalize_count", required=False, help="cmd for normalize_count.awk" , default='/usr/bin/gawk -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/normalize_count.awk')
ap.add_argument("-sb", "--sigtobed", required=False, help="cmd for sigtobed.awk" , default='/usr/bin/gawk -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/sigtobed.awk')
ap.add_argument("-bwa", "--bwa", required=False, help="path to bwa binary" , default='/lustre2/lucas_pkuhpc/bin/bwa')
ap.add_argument("-t", "--nthreads", required=False, help="number of threads" , default=1 , type=int)
ap.add_argument("-f", "--force_overwrite_flag", help="run all cmds even if the output files already exist" , action="store_true" )
ap.add_argument("-rn", "--run_counts_normalization_flag", help="generates .norm.tsv from .counts.tsv. needed for milllion X amplicon seq, but not 10k coverage WGS" , action="store_true" )
ap.add_argument("-no2unix", "--no_run_dos2unix_on_fasta", help="don't run dos2unix to make sure the newlines are correct" , action="store_true" )
ap.add_argument("-q", "--mapq_min_quality", help="use samtoosl to remove reads with mapq < Q" , default=0 , required=False  , type=int )
ap.add_argument("-m", "--ram_per_thread", help="how many GB of RAM to use for each thread (int)" , default=5 , required=False  , type=int )
args = ap.parse_args()


# make sure we can find all the commands & files
#  also prints the commands used, for logging & debugging
die_with_error_flag = False
if not os.path.isfile( args.find_mh ):
	print( "find_mh not found, exiting\n   looked in:\t" ,  args.find_mh  , "\n" )
	die_with_error_flag = True
if not os.path.isfile( args.bwa ):
	print( "bwa not found, exiting\n   looked in:\t" ,  args.bwa  , "\n" )
	die_with_error_flag = True

if not os.path.isfile( args.read1 ):
	print( "read1 not found, exiting\n   looked in:\t" ,  args.read1  , "\n" )
	die_with_error_flag = True
#if not os.path.isfile( args.read2 ):
#	print( "read2 not found, exiting\n   looked in:\t" ,  args.read2  , "\n" )
#	die_with_error_flag = True
if not os.path.isfile( args.genome_fasta_file ):
	print( "genome_fasta_file not found, exiting\n   looked in:\t" ,  args.genome_fasta_file  , "\n" )
	die_with_error_flag = True

pprint.pprint( vars(args) ) 
print('\n')

if die_with_error_flag:
	sys.exit(1)

# file names to create
bamcramfile = args.base_name + '.bam'
probesdir = 'probes___' + os.path.basename(args.genome_fasta_file)
signatures_file = os.path.join( probesdir ,  os.path.basename(args.genome_fasta_file) + '.signatures' )

# ## ## ## start by indexing the genome and making the .cram file
print("Start time: ", time.asctime(time.localtime()))
start_time = time.perf_counter()


if not os.path.isfile( args.genome_fasta_file + '.pac'):
    print('--- indexing genome --- ')
    cmd = args.bwa + ' index ' + args.genome_fasta_file 
    print( cmd  )
    subprocess.run(  cmd  ,  shell=True , check=True)
    print(' ---- ---- ')

print(' ---- aligning reads to the genome ---- ')
mapq_filter = ' ' 
if (args.mapq_min_quality > 0):
    mapq_filter = ' -q ' + str(args.mapq_min_quality) + ' ' 

ram_arg = ' -m' + str(args.ram_per_thread) + 'G '
cmd = args.bwa + ' mem -Y ' + ' -t ' + str(args.nthreads) + ' ' + args.genome_fasta_file + ' ' + args.read1 + ' ' + args.read2 \
	+ ' | samtools view -@ 3 -F 2048 -Sb ' + mapq_filter \
	+ ' | samtools sort -O BAM  --reference ' + args.genome_fasta_file + ' -@ 4 ' + ram_arg + ' -o ' + bamcramfile
print( cmd  )
if os.path.isfile( bamcramfile):
	if (os.path.getsize(bamcramfile) > 1000):
		print( bamcramfile + " already exists. Skipping bwa mem\n")
else:
	subprocess.run(  cmd  ,  shell=True , check=True)
	subprocess.run(  'samtools index '  + bamcramfile ,  shell=True , check=True)
end_time_1 = time.perf_counter()
print(' ---- bwa mem took ' + str( end_time_1 - start_time)   + ' seconds.      ---- \n')

# ## ##  run find_mh & generate_signatures.awk
print(' ---- find_mh & extract signatures ------ ')

# don't do anything if signatures_files exists, and if force_overwrite is False
if ( os.path.isfile( signatures_file + '.bed') & (args.force_overwrite_flag == False) ):
	print( signatures_file + '.bed    exists. not running\n')
else:
    # by default, run dos2unix on the fasta file to fix problems
    if not (args.no_run_dos2unix_on_fasta):
        cmd = 'dos2unix ' + args.genome_fasta_file
        subprocess.run(  cmd  ,  shell=True , check=True)

    # (1) generate .mh file
    subprocess.run(  'mkdir -p ' + probesdir ,  shell=True , check=True)
    probefiles = os.path.join( probesdir , os.path.basename(args.genome_fasta_file) + '.mh')
    cmd = args.find_mh + ' ' + probefiles + ' < ' + args.genome_fasta_file 
    print( cmd  )
    subprocess.run(  cmd  ,  shell=True , check=True)
	
    # (2) generate signatures file
    list_of_mh_files = subprocess.run( 'ls ' + probefiles  +  '*' + '  | sort -t. -k1n,1 ' , shell=True , capture_output=True , text=True )
    list_of_mh_files = list_of_mh_files.stdout.replace( '\n' , ' ')
    cmd = args.generate_signatures + ' ' + args.genome_fasta_file + ' ' + list_of_mh_files + ' > ' + signatures_file
    print( cmd )
    subprocess.run(  cmd  ,  shell=True , check=True)
	
    # (3) generate .bed file for calculating depth (coverage) and normalized read counts
    cmd = args.sigtobed + ' ' + signatures_file  + ' > ' + signatures_file + '.bed'
    print (cmd ) 
    subprocess.run(  cmd  ,  shell=True , check=True)

end_time_2 = time.perf_counter()
print( ' ---------  find_mh & extract signatures took : ' +  str( end_time_2 - end_time_1) + ' seconds -------\n')


# ## ## catch_signatures.awk & calculate relative coverage (normalized read counts) : 
# (1) samtools view to count reads
if (os.path.isfile( args.base_name + '.sign.count.tsv' ) or os.path.isfile( args.base_name + '.sign.count.tsv.gz' )) :
    print( args.base_name + '.sign.count.tsv[.gz] found, skipping. ' ) 
else:
    cmd = 'samtools view -T ' + args.genome_fasta_file + ' ' + bamcramfile + ' | ' + args.catch_signatures + ' ' + signatures_file + ' - > ' + args.base_name + '.sign.count.tsv'
    print( cmd ) 
    subprocess.run(  cmd  ,  shell=True , check=True)
#extract the sites found & compress the original
    cmd2 = 'gzip --force --best ' + args.base_name + '.sign.count.tsv'
    subprocess.run(  cmd2  ,  shell=True , check=True)
    cmd = 'gunzip -c ' + args.base_name + '.sign.count.tsv.gz | grep -Pv \'\t0\t0$\' ' + args.base_name + '.sign.count.tsv | perl -ne \'chomp ; print "$_\t' + args.base_name + '\n" ; \' > ' + args.base_name + '.sites_found.txt'
    print( cmd ) 
    subprocess.run(  cmd  ,  shell=True , check=True)



if (args.run_counts_normalization_flag):
    coverage_file = args.base_name + '.cov'
    if os.path.isfile( coverage_file ) and (os.path.getsize(coverage_file) > 1000) :
        print( coverage_file + " already exists. Skipping samtools depth\n")
    else:
        # (2) samtools depth  to count total reads at each position w/MH
        cmd = 'samtools depth -b ' + signatures_file + '.bed' + ' -d 0 ' + bamcramfile + ' > ' + coverage_file
        print( cmd ) 
        subprocess.run(  cmd  ,  shell=True , check=True)

    # (3) normalize_count creates the final output file
    #gunzip -c AmpliconLib_E3_Library_4.sign.count.tsv.gz | /usr/bin/gawk -f /home/lucas_pkuhpc/Develop/MicroHomologyMediatedIndels/modules/normalize_count.awk AmpliconLib_E3_Library_4.cov  - 
    #cmd = args.normalize_count + ' ' + args.base_name + '.cov' +  '   ' + args.base_name + '.sign.count.tsv'  + ' > ' + args.base_name + '.sign.norm.tsv'
    if os.path.isfile( args.base_name + '.sign.norm.tsv'):
        if os.path.getsize(args.base_name + '.sign.norm.tsv') < 1000 :
            os.remove(args.base_name + '.sign.norm.tsv')
            print( '.sign.norm.tsv exists but is too small. deleting!')
    if not os.path.isfile( args.base_name + '.sign.norm.tsv') :
        cmd = 'gunzip -c ' + args.base_name + '.sign.count.tsv.gz | ' + args.normalize_count + ' ' + coverage_file +  ' - > ' + args.base_name + '.sign.norm.tsv'
        print( cmd )
        subprocess.run(  cmd  ,  shell=True , check=True)
        cmd2 = 'gzip --force --best ' + args.base_name + '.sign.norm.tsv'
        subprocess.run(  cmd2  ,  shell=True , check=True)
    else:
        print( args.base_name + '.sign.norm.tsv' + 'already exits, skipping')




print( ' --------- catch_signatures & normalized coverage took: ' +  str( time.perf_counter() - end_time_2) + ' seconds ------\n')
print( '  TOTAL time = ' +  str( time.perf_counter() - start_time) + ' seconds ------\n')


