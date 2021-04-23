DSBS analyzer
==========

DSBS analyzer is a pipeline to analyzing Double Strand Bisulfite Sequencing data, which could simultaneously identify SNVs and evaluate DNA methylation levels in a single base resolution. In DSBS, bisulfite-converted Watson strand and reverse complement of bisulfite-converted Crick strand derived from the same double-strand DNA fragment were sequenced in read 1 and read 2, and aligned to the same position on reference genome. By simultaneous analyzing the sequence of read 1 and read 2, the sequence and DNA methylation state of DNA fragment could be deduced.

Schematic
---------

<img src="https://github.com/tianguolangzi/pic/blob/main/raw.png"  width="50%" height="50%" /><br/>


DSBS analyzer Usage
-------------------
### QuickStart

`python3 DSBS_analyzer.py --config DSBS.config  -a /data_path/*_1.fq.gz  -b /data_path/*_2.fq.gz   -o  /outdir  & `

`python3 DSBS_analyzer.py`

```
                                                
             m                                      m                           
             |                                      |                           
        AGCTACGT  bisulfite treat   -->  A  G {T} TACGT        AGTTACGT         
        ||||||||-------------------->|   |  |  |  |||||                         
        TCGATGCA        (C>T)       -->  T {T} G  ATGCA        AACTAGCT         
              |                                      |                        
              m                                      m                        
        C>T  1) bisulfite conversion ?                                          
             2) mutation ?                                                      
             3) sequencing errors ?                                             
    
usage: DSBS_analyzer.py [-h] [--config CONFIG] -a [FQ1 [FQ1 ...]] -b
                                  [FQ2 [FQ2 ...]] [-r REF] 
                                  [-k [KNOWNSITES [KNOWNSITES ...]]] [-d DEPTH] 
				  [-q CLEANQUALITY] [--gapsize GAPSIZE] [--fastquniq]  
				  [--remove_tmp] [--bissnp] [-t CPU] [--pp PP] -o
                                  OUTDIR [--qsub]

  -h, --help            show this help message and exit
  --config CONFIG       config file
  -a [FQ1 [FQ1 ...]], --fq1 [FQ1 [FQ1 ...]]
                        Fastq1(s), Space separation
  -b [FQ2 [FQ2 ...]], --fq2 [FQ2 [FQ2 ...]]
                        Fastq2(s), Space separation
  -r REF, --ref REF     refrence
  -k [KNOWNSITES [KNOWNSITES ...]], --knownsites [KNOWNSITES [KNOWNSITES ...]]
                        knownsites
  -d DEPTH, --depth DEPTH
                        the minimum depth , default is 4
  -q CLEANQUALITY, --cleanquality CLEANQUALITY
                        quality when clening the fastq
  --gapsize GAPSIZE     the size of insert or del, default is 3
  --fastquniq           rmdup fastq before cleanning fastq
  --remove_tmp          rm temp file
  --bissnp              bisSNP
  -t CPU, --cpu CPU     thread, default is 20
  --pp PP               1,execute ,0,skip. 1)qualimap, 2)clean, 3)align,
                        4)dealBam, 5)call mutation and methylation , default 11111
  -o OUTDIR, --outdir OUTDIR
                        outdir
  --qsub                use qsub
```



#### config

`cat  DSBS.config`

Here is the example of config file which configured the path of softwares and databases and  software operating environment. Y For those necessary softwares and databases, you are on your own. It is important to note that each software needs to have execute permission(`chmod a+x soft`).

```
#soft

java:/public/soft/jdk1.8.0_111/bin/java:soft
java1.6:/usr/bin/java:soft
python3:/public/soft/python3:soft

fastqc:/public/soft/fastqc:soft
qualimap:/public/soft/qualimap_v2.2.1/qualimap:soft

#clean,rmdup
fastuniq:/public/soft/FastUniq/fastuniq:soft
cutadapt:/public/soft/cutadapt:soft
filter_fq:/public/soft/filter_fq.py:soft
stdmap_awk:/public/soft/stdmap.awk:soft
diffmap_awk:/public/soft/diffmap.awk:soft

bsmap:/public/soft/bsmap2.7/bsmap:soft
bissnp:/public/soft/BisSNP/BisSNP-0.82.2.jar:soft

#qsub
job_sub_py_2:/public/soft/job_sub_py_2.py:soft

#index, sort, merge, split 
samtools:/public/soft/samtools-1.8/samtools:soft
picard:/public/soft/picard-2.17.3/picard.jar:soft


#call SNP and  methy
DSBS:/public/soft/DSBS.py:soft
#resource
refDir:/public/database/hg19/:resource
ref:/public/database/hg19/hg19.fa:resource
1000G_omni:/public/database/1000G_omni2.5.hg19.sites.vcf:resource
1000G_phase1:/public/database/1000G_phase1.snps.high_confidence.hg19.sites.vcf:resource
dbsnp:/public/database/dbsnp_138.hg19.vcf:resource
dbsnp_gz:/public/database/dbsnp_138.hg19.vcf.gz:resource
hapmap:/public//database/hapmap_3.3.hg19.sites.vcf:resource
```



Dependencies
------------

DSBS depends on 
* [python3+](https://www.python.org/)
* [termcolor](https://pypi.python.org/pypi/termcolor/1.1.0)
* [pysam-0.11.0](https://github.com/pysam-developers/pysam/archive/v0.11.0.tar.gz)
* [pyfasta](https://pypi.python.org/pypi/pyfasta/0.5.2)
* [tabix](https://)

## DSBS Usage

### QuickStart

`python DSBS.py  example.bam -g chr1.fa -d dbsnp138.vcf.gz --chr chr1 -o outdir -q` 


### Index 
One time only, it is necessary to index  a reference sequence and the dbsnp file。
The commands:

	bgzip dbsnp_138.hg19.vcf
	tabix -s 1 -b 2 dbsnp_138.hg19.vcf.gz (tabix -s 1 -b 2 -p vcf dbsnp_138.hg19.vcf.gz )
 	

`python DSBS.py`

```
usage: DSBS.py [-h] [-v] [--maxDistance MAXDISTANCE] [--maxLen MAXLEN]
               [--mixLen MIXLEN] [--mixReadQual MIXREADQUAL]
               [--mixReadQualN MIXREADQUALN] [--maxN MAXN]
               [--maxSeqErr MAXSEQERR] [--maxSnp MAXSNP] [--maxIndel MAXINDEL]
               [--maxBp MAXBP] [--maxMut MAXMUT] [--minQual MINQUAL]
               [--minVaf MINVAF] [--secAlign] [--debug] [-q] [-c COVERAGE]
               [-p CPU] [-o OUTDIR] -g GENOMEFILE -d DBSNP [--Chr CHR]
               bamFile
	       
positional arguments:
bam                           input BAM file
optional arguments:
-h, --help                    show this help message and exit
-v, --version                 show program's version number and exit
--maxDistance MAXDISTANCE     the maxmum distance of the start postions of the
                              paired reads in reference genome, default is 50
--maxLen MAXLEN               the maxmum length of a read, default is 200
--minLen MINLEN               the minimum length of a read, default is 50
--minReadQual MINREADQUAL     the minimum quality for assessing the read 
                              overall qualities, default is 20
--minReadQualN MINREADQUALN   the minimum ratio which quality is greater than
                              the minReadQual, default is 0.6
--maxN MAXN                   the maxmum number of N within a read, default is 5
--maxSeqErr MAXSEQERR         the maxmum number of sequencing errors within a paired
                              reads, default is 10
--maxSnp MAXSNP               the maxmum number of snp within a paired reads,
                              default is 5
--maxIndel MAXINDEL           the maxmum number of indel within a paired reads,
                              default is 1
--maxBp MAXBP                 a certain range for allowing the currence of certain
                              mutation ,default is 50
--maxMut MAXMUT               the maxmum number of allowing the currence of
                              mutation within a certain range ,default is 3
--minQual MINQUAL             the minimum quality for call Variation, default is 20
--minVaf MINVAF               the minimum value of allele, default is 0.5
--window WINDOW               ignoring the snps in a certain range of a indel, the
                              range is 5bp
--minCoverage MINCOVERAGE     coverage or minimum number of reads desired,
                              default is 8
--maxCoverage MAXCOVERAGE     coverage or maxmum number of reads desired,
                              default is 500
--strand                      only consider the correct comparison direction,   
                              read1=++ && read2=-- || read1=-+ && read2=+-
--secAlign                    consider non-optimal alignments
--CpG                         output CpG
--CHG                         output CHG
--CHH                         output CHH
--debug                       dubug mode
-q, --quite                   quiet mode
-p CPU, --cpu CPU             the number of working threads, default is 50
-o OUTDIR, --outDir OUTDIR    the outDir
-g GENOMEFILE, --genomeFile GENOMEFILE
                            the chromosome reference fasta
-d DBSNP, --dbsnp DBSNP       the dbsnp file
--Chr CHR                    chromosome
```

####  Single jobs
`
python3 DSBS.py   chr1.merge.sorted.rmdup.realigned.recal.bam  --maxBp 50  --minVaf 0.1  -q -p 40  -o outdir  -g  chr1.fa --Chr chr1 -d dbsnp_138.hg19.vcf.gz  --secAlign 
`

#### Batch jobs
use job_sub_py_2.py which based on PBS to  deliver the jobs

`
for i in {1..22} X Y M;
do 
python3 job_sub_py_2.py --jobname chr$i --cpu 2 --work "python3 /public/soft/DSBS.py 	/datapath/chr$i.merge.sorted.rmdup.realigned.recal.bam --maxBp 50 --minVaf 0.1  -q --cpu 40  -o /outdir -g   ~/public/database/hg19/chr$i.fa  --Chr chr$i -d 	/public/database/dbsnp_138.hg19.vcf.gz "
done
`



Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/tianguolangzi/DSBS/issues) for any of these, including example reports where possible.



