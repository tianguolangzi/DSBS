## DSBS

DSBS is a special tool to call variants and methylations for the DSBS project.
DSBS is written in Python3.5 (tested with v3.4 / v3.5 ).

## Accuracy

### WGS vs DSBS

	[此处要放DSBS和WGS图]
	
### WGBS vs DSBS

	[此处要放甲基化图]
	
### bissnip vs DSBS
	
	[软件对比图]	


## Dependencies
DSBS depends on 
* [python3+](https://www.python.org/)
* [termcolor](https://pypi.python.org/pypi/termcolor/1.1.0)
* [pysam](https://pypi.python.org/pypi/pysam)
* [pyfasta](https://pypi.python.org/pypi/pyfasta/0.5.2)
* [tabix](https://)

## Usage

### QuickStart

`python DSBS.py  example.bam -g chr1.fa -d dbsnp138.vcf.gz --chr chr1 -o outdir -q` 


### Index 
One time only, it is necessary to index  a reference sequence and the dbsnp file。
The commands:

	bgzip dbsnp138.vcf 
	tabix -s 1 -b 2 dbsnp138.vcf.gz 

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
 * `bamFile`               the input BAM file

#### optional arguments:
 * `-h, --help`             show this help message and exit
 * `-v, --version`          show program's version number and exit
 * `--maxDistance MAXDISTANCE`  配对read间比对到参考基因组上起始位置间的差,默认50
 * `--maxLen MAXLEN`       read最长长度,默认200
 * `--mixLen MIXLEN`       read最短长度,默认50
 * `--mixReadQual MIXREADQUAL` 筛选read整体质量值时最小质量值,默认20
 * `--mixReadQualN MIXREADQUALN` 筛选read整体质量值时,大于最小质量值比例,默认0.6
 * `--maxN MAXN`           一条read最多出现N的个数,默认4
 * `--maxSeqErr MAXSEQERR` 配对read最多允许出现测序错误个数,默认10
 * `--maxSnp MAXSNP`       配对read最多允许出现Snp的个数,默认5
 * `--maxIndel MAXINDEL`   配对read最多允许出现Indel的个数,默认1
 * `--maxBp MAXBP`         筛选一定范围多个突变时,限定的范围,默认是50bp
 * `--maxMut MAXMUT`       一定范围内,最多允许出现突变的个数,默认3个
 * `--minQual MINQUAL`     支持SNP的最小质量值,默认20
 * `--minVaf MINVAF`       等位基因频率最小值,默认0.1
 * `--secAlign`            找突变的时候考虑非最优比对，默认否
 * `--debug`               dubug输出模式,调试用的
 * `-q, --quite`           安静模式,默认开启
 * `-c COVERAGE, --coverage COVERAGE`  coverage or minimum number of reads desired
 * `-p CPU, --cpu CPU`     子进程,默认50
 * `-o OUTDIR, --outdir OUTDIR` outdir
 * `-g GENOMEFILE, --genomeFile GENOMEFILE` input FASTA file
 * `-d DBSNP, --dbsnp DBSNP`   dbsnp
 * `--Chr CHR`            chr1.fa

## Update 
* V1.1 (2017-09-10)
  * Add the dbsnp id in outfile 
 
## Contributions & Support

Contributions and suggestions for new features are welcome, as are bug reports! Please create a new [issue](https://github.com/tianguolangzi/DSBS/issues) for any of these, including example reports where possible.

If any question, please :
[@Zhang Kun](https://github.com/tianguolangzi) (tianguolangzi@yahoo.com; tianguolangzi@gmail.com)



## Contributors
* Project lead and main author: [@Liang Jialong](https://github.com/lll)
* Project lead and main author: [@Zhang Kun](https://github.com/tianguolangzi)
Thanks to my friends who give me a hand with project:
* [@Shi Xiaohui]()
* [@Li Xianfeng]()
* [@Teng Huajing]()
