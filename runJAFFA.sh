#! /bin/csh -f
# set echo

##  qsub -q <queue> -w e -N <job_name> \
##    -l h_vmem=<memory, e.g. 4G> -l h_rt=<time> -l s_rt=<time> \
##    -pe smp <num_slots> -o <outputlogfile> -e <errorlogfile> <pathtoScript> <arg1> <arg2>

#$ -q jrf.q
#$ -w e -N jaffa_fx
#$ -l h_vmem=4G
#$ -pe smp 8 -o smap.log -j y
#$ -V -cwd
## -t 1-16 -tc 4

#<usage>
# as array :
# qsub ./discasm.qsh ../path/to/fastq/ ./path/to/output  
# for single run
# SGE_TASK_ID=1 NSLOTS=8 ./runJAFFA.qsh ../path/to/fastq/ ./path/to/output
#</usage>

source /ifs/e63data/reis-filho/genomes/homo_sapiens/gVars.csh

/home/gularter/src/JAFFA/tools/bin/bpipe run -n 1 -r -p fastqInputFormat="%.*.fastq.gz" ~/src/JAFFA/JAFFA_assembly.groovy $1

