
![image](https://github.com/tianguolangzi/pic/blob/main/1.png)
|Software|version|website|
|---------|---------|-------|
|Fastuniq|1.1|https://homes.cs.washington.edu/~dcjones/fastq-tools|
|Fastqc|0.11.5|http://www.bioinformatics.babraham.ac.uk/projects/FastQC|
|Cutadapt|1.11|https://github.com/marcelm/cutadapt|
|Trim_galore|0.4.2|http://www.bioinformatics.babraham.ac.uk/projects/trim_galore|
|Bsmap|2.74|https://code.google.com/archive/p/bsmap/|
|Samtools|1.8|http://samtools.sourceforge.net/|
|Picard|2.9|http://broadinstitute.github.io/picard/|
|Bissnp|0.82.2|https://sourceforge.net/projects/bissnp/|
|Bamqc||https://github.com/s-andrews/BamQC|
|Qualimap|2.2.1|http://qualimap.bioinfo.cipf.es/|
|VarScan|2.3.9|http://varscan.sourceforge.net/|
|GATK|3.8|https://software.broadinstitute.org/gatk/|
|Freebayes|1.1.0-50|https://github.com/ekg/freebayes/releases|
|MethylDackel|0.3|https://github.com/dpryan79/MethylDackel/releases|
|Bcftools|0.1.19|http://www.htslib.org/doc/bcftools.html|
|BS-Snper||https://github.com/hellbelly/BS-Snper|
|BWA|0.7.7|https://github.com/lh3/bwa/releases|
|Tabix|0.2.5|http://www.htslib.org/doc/tabix.html|
|hg19.fa|hg19|ftp://ftp.broadinstitute.org/bundle/hg19/|
|dbsnp|138|ftp://ftp.broadinstitute.org/bundle/hg19/|
|Hapmap|3.3|ftp://ftp.broadinstitute.org/bundle/hg19/|
|Mills_and_1000G_gold_standard.indels|hg19|ftp://ftp.broadinstitute.org/bundle/hg19/|
|1000G_phase1.snps|hg19|ftp://ftp.broadinstitute.org/bundle/hg19/|
|1000G_phase1.indels|hg19|ftp://ftp.broadinstitute.org/bundle/hg19/|
|1000G_omni2.5|hg19|ftp://ftp.broadinstitute.org/bundle/hg19/|
|R|3.2.3|https://www.r-project.org/|
|Vcftools|0.1.13|http://vcftools.sourceforge.net/|
|Perl|5.20.3|https://www.perl.org/|
|Java|1.6&&1.8|https://java.com/download|
|Multiqc|0.9|https://github.com/ewels/MultiQC/releases|



## DSBS

DSBS is a special tool to call variants and methylations for the DSBS project.
DSBS is written in Python3.5 (tested with v3.4 / v3.5 ).

## Citation


### WGS vs DSBS
* Coverage in whole genome
* Coverage in non-repeat region
* Coverage in Exome region
* Snp Overlap between WGS and DSBS

### WGBS vs DSBS
* Coverage in all CpG sites of whole genome
* Coverage in all CpG sites of non-repeat region
* Coverage in all CpG sites of repeat region
* Methylation Pearson in WGBS and DSBS	

### BS-snper vs FreeBayes vs bcftools vs VarScan vs GATK vs DSBS
* [BS-snper](https://github.com/hellbelly/BS-Snper) 
* [FreeBayes](https://github.com/ekg/freebayes) 
* [bcftools](https://github.com/samtools/bcftools) 
* [VarScan](https://github.com/dkoboldt/varscan) 
* [GATK](https://github.com/broadinstitute/gatk) 
* [DSBS](https://github.com/tianguolangzi/DSBS)
	


## Dependencies
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
use [job_sub_py_2](https://github.com/tianguolangzi/ZK-Tools)  based on [PBS](https://github.com/pbspro/pbspro) to  deliver the jobs

`
for i in {1..22} X Y M;
do 
job_sub_py_2 --jobname chr$i --cpu 2 --work "python3 ~/zhangkun/bin/DSBS.py 	/public/home/zhangkun/work/DSBS/hg19/DSBS_NEW_merge/bam/chr$i.merge.sorted.rmdup.realigned.recal.bam --maxBp 50 --minVaf 0.1  -q --cpu 40  -o /public/home/zhangkun/work/DSBS/hg19/DSBS_NEW_merge/script6 -g   ~/public/database/hg19/chr$i.fa  --Chr chr$i -d 	dbsnp_138.hg19.vcf.gz "
done
`

## DSBS pipeline Usage

### QuickStart
`python3 DSBS_pipeline.V1.py -a fq1(s) -b fq2(s) --ref hg19.fa --config  DSBS.config -o outdir`

`python3 DSBS_pipeline.V1.py`

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
    
usage: DSBS_pipeline.V1.py [-h] [--config CONFIG] -a [FQ1 [FQ1 ...]] -b
                                  [FQ2 [FQ2 ...]] [-r REF] 
                                  [-k [KNOWNSITES [KNOWNSITES ...]]] [-d DEPTH] 
				  [-q CLEANQUALITY] [--gapsize GAPSIZE] [--fastquniq]  
				  [--remove_tmp] [--bissnp] [-t CPU] [--pp PP] -o
                                  OUTDIR [--qsub]
DSBS_pipeline.V1.py: error: the following arguments are required: -a/--fq1, -b/--fq2, -o/--outdir
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
                        4)dealBam, 5)call mutation and methylation 默认 11111
  -o OUTDIR, --outdir OUTDIR
                        outdir
  --qsub                use qsub
```

#### Run DSBS_pipeline.V1.py
`python3 DSBS_pipeline.V1.py --config DSBS.config  -a /data_path/*_1.fq.gz  -b /data_path/*_2.fq.gz  --fastquniq   --bissnp  -t 40  --pp 11111 -o  /outdir/ --qsub  & `

#### config

`cat  DSBS.config`

Here is the example of config file which configured the path of softwares and databases and  software operating environment. Y For those necessary softwares and databases, you are on your own. It is important to note that each software needs to have execute permission(`chmod a+x soft`).

```
#soft

java:/public/home/zhangkun/bin/java_bin/jdk1.8.0_111/bin/java:soft
java1.6:/usr/bin/java:soft
python3:/public/bin/python3:soft

fastqc:/public/home/public/bin/fastqc:soft
qualimap:/public/home/zhangkun/bin/qualimap_v2.2.1/qualimap:soft
#muiltqc:/public/software/Python3.5/bin/multiqc:soft

#clean,rmdup
fastuniq:/public/home/zhangkun/bin/bin/FastUniq/fastuniq:soft
cutadapt:/public/bin/cutadapt:soft
filter_fq:/public/home/zhangkun/bin/filter_fq.py:soft
stdmap_awk:/public/home/zhangkun/work/DSBS/bin/stdmap.awk:soft
diffmap_awk:/public/home/zhangkun/work/DSBS/bin/diffmap.awk:soft

bsmap:/public/home/zhangkun/work/DSBS/bin/bsmap2.7/bsmap:soft
bissnp:/public/home/zhangkun/work/DSBS/bin/BisSNP/BisSNP-0.82.2.jar:soft

#qsub
job_sub_py_2:/public/home/zhangkun/bin/python3_bin/job_sub_py_2.py:soft

#index, sort, merge, split 
samtools:/public/home/zhangkun/bin/samtools:soft
picard:/public/home/zhangkun/bin/picard.jar:soft


#call SNP and  methy
DSBS:/public/home/zhangkun/bin/DSBS.py:soft
#resource
refDir:/public/home/zhangkun/database/hg19/:resource
ref:/public/home/zhangkun/database/hg19/hg19.fa:resource
1000G_omni:/public/home/zhangkun/database/1000G_omni2.5.hg19.sites.vcf:resource
1000G_phase1:/public/home/zhangkun/database/1000G_phase1.snps.high_confidence.hg19.sites.vcf:resource
dbsnp:/public/home/zhangkun/database/dbsnp_138.hg19.vcf:resource
dbsnp_gz:/public/home/zhangkun/database/dbsnp_138.hg19.vcf.gz:resource
hapmap:/public/home/zhangkun/database/hapmap_3.3.hg19.sites.vcf:resource
```

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/tianguolangzi/DSBS/issues) for any of these, including example reports where possible.

If any question, please :
[@Zhang Kun](https://github.com/tianguolangzi) (tianguolangzi@yahoo.com)



## Contributors
* Project lead and main author: [@Liang Jialong](https://github.com/lll)
* Project lead and main author: [@Zhang Kun](https://github.com/tianguolangzi)

Thanks to my friends who give me a hand with project:
* [@Li Xianfeng](https://github.com/xflicsu)
* [@Teng Huajing]()


