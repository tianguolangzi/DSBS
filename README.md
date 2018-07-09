
## DSBS pipeline depends on 

|Software|version|website|
|---------|---------|-------|
|Fastuniq|1.1|https://homes.cs.washington.edu/~dcjones/fastq-tools|
|Fastqc|0.11.5|http://www.bioinformatics.babraham.ac.uk/projects/FastQC|
|Cutadapt|1.11|https://github.com/marcelm/cutadapt|
|Trim_galore|0.4.2|http://www.bioinformatics.babraham.ac.uk/projects/trim_galore|
|Bsmap|2.74|https://code.google.com/archive/p/bsmap/|
|Samtools|0.1.19|http://samtools.sourceforge.net/|
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
|Dbsnp|138|ftp://ftp.broadinstitute.org/bundle/hg19/|
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

## License

## Accuracy


### WGS vs DSBS
* Coverage in whole genome
	[此处要放DSBS和WGS图]
* Coverage in non-repeat region
* Coverage in Exome region
* Snp Overlap between WGS and DSBS

### WGBS vs DSBS
* Coverage in all CpG sites of whole genome
* Coverage in all CpG sites of non-repeat region
* Coverage in all CpG sites of repeat region

	[此处要放甲基化图]
* Methylation Pearson in WGBS and DSBS	

### BS-snper vs FreeBayes vs bcftools vs VarScan vs GATK vs DSBS
* [BS-snper](https://github.com/hellbelly/BS-Snper) 
* [FreeBayes](https://github.com/ekg/freebayes) 
* [bcftools](https://github.com/samtools/bcftools) 
* [VarScan](https://github.com/dkoboldt/varscan) 
* [GATK](https://github.com/broadinstitute/gatk) 
* [DSBS](https://github.com/tianguolangzi/DSBS)
	
[软件对比图]	


## Dependencies
DSBS depends on 
* [python3+](https://www.python.org/)
* [termcolor](https://pypi.python.org/pypi/termcolor/1.1.0)
* [pysam](https://pypi.python.org/pypi/pysam)
* [pyfasta](https://pypi.python.org/pypi/pyfasta/0.5.2)
* [tabix](https://)

## DSBS Usage

### QuickStart

`python DSBS.py  example.bam -g chr1.fa -d dbsnp138.vcf.gz --chr chr1 -o outdir -q` 


### Index 
One time only, it is necessary to index  a reference sequence and the dbsnp file。
The commands:

	bgzip dbsnp_138.hg19.vcf
	tabix -s 1 -b 2 dbsnp_138.hg19.vcf.gz

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
```

#### positional arguments:
 * `bamFile`               the input BAM file

#### optional arguments:
 * `-h, --help`             show this help message and exit
 * `-v, --version`          show program's version number and exit
 * `--maxDistance MAXDISTANCE`  配对read间比对到参考基因组上起始位置间的差,默认50
 * `--maxLen MAXLEN`       read最长长度,默认200
 * `--mixLen MIXLEN`       read最短长度,默认50
 * `--mixReadQual MIXREADQUAL` 筛选read整体质量值时最小质量值,默认20
 * `--mixReadQualN MIXREADQUALN`  筛选read整体质量值时,大于最小质量值比例,默认0.6
 * `--maxN MAXN`           一条read最多出现N的个数,默认4
 * `--maxSeqErr MAXSEQERR` 配对read最多允许出现测序错误个数,默认10
 * `--maxSnp MAXSNP`       配对read最多允许出现Snp的个数,默认5
 * `--maxIndel MAXINDEL`   配对read最多允许出现Indel的个数,默认1
 * `--maxBp MAXBP`         筛选一定范围多个突变时,限定的范围,默认是50bp
 * `--maxMut MAXMUT`       一定范围内,最多允许出现突变的个数,默认3个
 * `--minQual MINQUAL`     支持SNP的最小质量值,默认20
 * `--minVaf MINVAF`       等位基因频率最小值,默认0.1
 * `--secAlign`            找突变的时候考虑非最优比对，默认否
 * `--debug`               debug输出模式,调试用的
 * `-q, --quite`           安静模式,默认开启
 * `-c COVERAGE, --coverage COVERAGE`  coverage or minimum number of reads desired
 * `-p CPU, --cpu CPU`     子进程,默认50
 * `-o OUTDIR, --outdir OUTDIR` outdir
 * `-g GENOMEFILE, --genomeFile GENOMEFILE` input FASTA file
 * `-d DBSNP, --dbsnp DBSNP`   dbsnp
 * `--Chr CHR`            染色体号

####  Single jobs

`
python3 DSBS.py   chr1.merge.sorted.rmdup.realigned.recal.bam  --maxBp 50  --minVaf 0.1  -q -p 40  -o outdir  -g  chr1.fa --Chr chr1 -d dbsnp_138.hg19.vcf.gz  --secAlign 
`

#### Batch jobs
use [job_sub_py_2](https://github.com/tianguolangzi/ZK-Tools)  based on [PBS](https://github.com/pbspro/pbspro) to  deliver the jobs

`
for i in {1..22} X Y M;
do 
job_sub_py_2 --jobname chr$i --cpu 2 --work "python3 ~/zhangkun/bin/python3_bin/DSBS.py 	/public/home/jcli/zhangkun/work/DSBS/hg19/DSBS_NEW_merge/bam/chr$i.merge.sorted.rmdup.realigned.recal.bam --maxBp 50 --minVaf 0.1  -q --cpu 40  -o /public/home/jcli/zhangkun/work/DSBS/hg19/DSBS_NEW_merge/script6 -g   ~/public/database/hg19/chr$i.fa  --Chr chr$i -d 	dbsnp_138.hg19.vcf.gz "
done
`

## DSBS pipeline Usage

### QuickStart
`python3 DSBS_pipeline.py -a fq1(s) -b fq2(s) --ref hg19.fa --config  DSBS_pipeline_config -o outdir`

`python3 ~/zhangkun/work/DSBS/bin/DSBS_pipeline.2017.6.14.py`

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
    
usage: DSBS_pipeline.2017.6.14.py [-h] [--config CONFIG] -a [FQ1 [FQ1 ...]] -b
                                  [FQ2 [FQ2 ...]] [-r REF] [-g GTF]
                                  [-k [KNOWNSITES [KNOWNSITES ...]]]
                                  [-d DEPTH] [-p PMTSIZE] [-w WINSIZE]
                                  [-q CLEANQUALITY] [--gapsize GAPSIZE]
                                  [--fastquniq] [--revise_sam] [--remove_tmp]
                                  [--bs2real] [--bissnp] [-t CPU] [--pp PP] -o
                                  OUTDIR [--qsub]
DSBS_pipeline.2017.6.14.py: error: the following arguments are required: -a/--fq1, -b/--fq2, -o/--outdir
```

#### optional arguments:
   * `-h, --help`            show this help message and exit
  * `--config CONFIG`       config file 
  * `-a [FQ1 [FQ1 ...]], --fq1 [FQ1 [FQ1 ...]]` Fastq1(s), 空格分隔
  * `-b [FQ2 [FQ2 ...]], --fq2 [FQ2 [FQ2 ...]]` Fastq2(s), 空格分隔
  * `-r REF, --ref REF`     参考基因组
  * `-g GTF, --gtf GTF`     参考基因
  * `-k [KNOWNSITES [KNOWNSITES ...]], --knownsites [KNOWNSITES [KNOWNSITES ...]]` 已知数据库
  * `-d DEPTH, --depth DEPTH`     最小深度,默认4
  * `-p PMTSIZE, --pmtsize PMTSIZE` 启动子区间大小,默认2000
  * `-w WINSIZE, --winsize WINSIZE`  窗口大小,默认200000
  * `-q CLEANQUALITY, --cleanquality CLEANQUALITY`  clean时的质量值控制
  * `--gapsize GAPSIZE`     插入缺失大小,默认3
  * `--fastquniq`           clean前rmdup一遍,默认False
  * `--revise_sam`          修正可能是由于测序错误导致的多测或少测一个碱基的情况
  * `--remove_tmp`          删除临时文件
  * `--bs2real`             未配对的转换重新bwa对比一遍,默认False
  * `--bissnp`              局部重新比对
  * `-t CPU, --cpu CPU`     进程, 默认 40
  * `--pp PP`               1,执行,0,跳过. 1)质控, 2)清洗, 3)比对, 4)处理, 5)突变和甲基化 , 6)统计作图 ,7)交集, 默认 1111111
  * `-o OUTDIR, --outdir OUTDIR` 输出目录     
  * `--qsub`                如果有多对fqs,是否使用qsub来多节点比对

#### Run DSBS_pipeline.2017.06.14.py
`python3 ~/zhangkun/work/DSBS/bin/DSBS_pipeline.2017.6.14.py --config ~/zhangkun/work/DSBS/bin/anzhen.config  -a /public/home/jcli/zhangkun/work/DSBS/DSBS_wangwen/*_1.fq.gz  -b /public/home/jcli/zhangkun/work/DSBS/DSBS_wangwen/*_2.fq.gz  --fastquniq   --bissnp  -t 40  --pp 1111100 -o  /public/home/jcli/zhangkun/work/DSBS/hg19/DSBS_20180223_wangwen/ --qsub  & `

#### config

`[jcli@loginI1 ~]$ cat  ~/zhangkun/work/DSBS/bin/DSBS.anzhen.config`

```
#这是安贞服务器上DSBS运行环境配置
#soft

#运行环境
java:/public/home/jcli/zhangkun/bin/java_bin/jdk1.8.0_111/bin/java:soft
java1.6:/usr/bin/java:soft
python3:/public/home/jcli/public/bin/python3:soft

#质控
fastqc:/public/home/jcli/public/bin/fastqc:soft
qualimap:/public/home/jcli/zhangkun/bin/bin/qualimap_v2.2.1/qualimap:soft
#muiltqc:/public/home/jcli/public/software/Python3.5/bin/multiqc:soft

#clean,rmdup
fastuniq:/public/home/jcli/zhangkun/bin/bin/FastUniq/fastuniq:soft
cutadapt:/public/home/jcli/public/bin/cutadapt:soft
filter_fq:/public/home/jcli/zhangkun/bin/python3_bin/filter_fq.py:soft
SAMdeduplicate:/public/home/jcli/zhangkun/work/DSBS/bin/SAMdeduplicate.pl:soft
stdmap_awk:/public/home/jcli/zhangkun/work/DSBS/bin/stdmap.awk:soft
diffmap_awk:/public/home/jcli/zhangkun/work/DSBS/bin/diffmap.awk:soft

Bs2Real:/public/home/jcli/zhangkun/work/DSBS/bin/Bs2Real.py:soft
bwa:/public/home/jcli/public/bin/bwa:soft
bsmap:/public/home/jcli/zhangkun/work/DSBS/bin/bsmap2.7/bsmap:soft
bissnp:/public/home/jcli/zhangkun/work/DSBS/bin/BisSNP/BisSNP-0.82.2.jar:soft
#bwameth:/public/home/jcli/public/software/Python3.5/bin/bwameth.py:soft

#qsub
job_sub_py_2:/public/home/jcli/zhangkun/bin/python3_bin/job_sub_py_2.py:soft

#index, sort, merge, split 
samtools:/public/home/jcli/public/bin/samtools:soft
picard:/public/home/jcli/zhangkun/bin/java_bin/picard.jar:soft
#picard:/histor/sun/wujinyu/zhangkun/DSBS/bin/picard2.9/picard.jar:soft


#call SNP and  methy
#DSBS:/public/home/jcli/zhangkun/work/DSBS/bin/DSBS.py:soft
DSBS:/public/home/jcli/zhangkun/bin/python3_bin/DSBS.py:soft
#resource
refDir:/public/home/jcli/public/database/hg19/:resource
ref:/public/home/jcli/public/database/hg19/hg19.fa:resource
gtf:/public/home/jcli/zhangkun/data_base/UCSC_gtf/hg19.gene.gtf.gz:resource
1000G_omni:/public/home/jcli/zhangkun/work/DSBS/bin/database/1000G_omni2.5.hg19.sites.vcf:resource
1000G_phase1:/public/home/jcli/zhangkun/work/DSBS/bin/database/1000G_phase1.snps.high_confidence.hg19.sites.vcf:resource
dbsnp:/public/home/jcli/zhangkun/work/DSBS/bin/database/dbsnp_138.hg19.vcf:resource
dbsnp_gz:/public/home/jcli/zhangkun/work/DSBS/bin/database/dbsnp_138.hg19.vcf.gz:resource
hapmap:/public/home/jcli/zhangkun/work/DSBS/bin/database/hapmap_3.3.hg19.sites.vcf:resource
```

## Update 
* V1.1 (2017-09-10)
  * Add the dbsnp id in outfile
* V1.2 (2018-01-29)
  * Add the methylation level in CpG/CHG/CHH in outfile
* V1.3 (2018-03-26)
  * Add the hemimethylation level in CpG/CHG/CHH in outfile
* V1.3 (2018-07-09)
  * Add the CNV information in outfile
## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/tianguolangzi/DSBS/issues) for any of these, including example reports where possible.

If any question, please :
[@Zhang Kun](https://github.com/tianguolangzi) (tianguolangzi@yahoo.com; tianguolangzi@gmail.com)



## Contributors
* Project lead and main author: [@Liang Jialong](https://github.com/lll)
* Project lead and main author: [@Zhang Kun](https://github.com/tianguolangzi)

Thanks to my friends who give me a hand with project:
* [@Shi Xiaohui]()
* [@Li Xianfeng](https://github.com/xflicsu)
* [@Teng Huajing]()


