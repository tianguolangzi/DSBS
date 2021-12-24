#!/public/home/jcli/public/bin/python3
# -*- coding: utf-8 -*-

#strand: 
#   ++: forward strand of Watson of reference (BSW)   
#   +-: reverse strand of Watson of reference (BSWC)
#   -+: forward strand of Crick of reference (BSC)  
#   --: reverse strand of Crick of reference (BSCC) 


from termcolor import colored
import os,sys
import re
import string
import argparse
from pyfasta import Fasta
from collections import defaultdict
import pysam
from  multiprocessing import Pool,Manager
import time

__version__ = "2.0"
__date__='20200825'
__author__='Kun Zhang'
__email__='tianguolangzi@yahoo.com'


def checkDir(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except :
            sys.stderr.write(colored("Error: failed to mkdir {} !\n".format(path),'red'))
            sys.exit(1)

def checkFile(File):
    if not os.path.exists(File):
        sys.stderr.write(colored('Error: {} is not existed! \n'.format(File),'red'))
        return 0
    else:
        return 1

def checkContinueBase(seq):
    seq=seq.upper()
    length=len(seq)
    return max([seq.count(i)/length for i in 'ATGC'])


def faseq(fa, chrom, start, end, cache=[None]):
    '''
    this is called by pileup which is ordered by chrom
    so we can speed things up by reading in a chrom at
    a time into memory
    '''
    if cache[0] is None or cache[0][0] != chrom:
        seq = "".join(x.strip() for i, x in enumerate(nopen("|samtools faidx %s %s" % (fa, chrom))) if i >0).upper()
        cache[0] = (chrom, seq)
    chrom, seq = cache[0]
    return seq[start - 1: end]
    
    


def getTime():
    nowTime=time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())
    return nowTime

def getDbsnp(dbsnp):
    #2017-09-10
    if not checkFile(dbsnp+".tbi"):
        if  os.system('tabix -s 1 -b 2 {dbsnp}'.format(dbsnp=dbsnp)):
            raise Exception("Please prepare a tabix file for dbsnp by the command 'tabix -s 1 -b 2 {dbsnp}' ".format(dbsnp=dbsnp))
    dbsnp_tb=pysam.TabixFile(dbsnp)
    return dbsnp_tb


def vcf2annovar(start0,end0,ref0,alt0):
    # AAC>A ==> AC>-
    minL = min(len(ref0),len(alt0))
    maxL = max(len(ref0),len(alt0))
    for i0 in range(minL):
        if ref0[0]==alt0[0]:
            ref0 = ref0[1:]
            alt0 = alt0[1:]
            start0 += 1
        else:
            break
    # GC>TC ==> G>T
    minL = min(len(ref0),len(alt0))
    maxL = max(len(ref0),len(alt0))
    for i0 in range(minL):
        if ref0[-1]==alt0[-1]:
            ref0 = ref0[:-1]
            alt0 = alt0[:-1]
    if ref0=="": 
        ref0 = "-"
        start0 -= 1
    if alt0=="": alt0 = "-"
    end0 = start0+len(ref0)-1
    return(start0,end0,ref0,alt0)


def pairCheckSnp(Mut_stand,baseDict,refNuc,IsSNP,Ismethy,minVaf):
    #
    ''''
    2016-11-20 
    +/-
    read1 read2 real_+  real_-
      A     A     A      T
      T     T     T      A
      C     C     mC     G
      G     G     G      mC
      #T/C  C/T   C      G
      T     C     C      G
      #G/A  A/G   G      C
      G     A     G      C
    
    2017-04-13
    ref->  alt   judgement condition
     A  ->  G       GG 1 or GA 2
     A  ->  C       CC 1 or TC 1
     A  ->  T       TT 1
     G  ->  A       AA 2
     G  ->  C       CC 1 or TC 1
     G  ->  T       TT 1 
     C  ->  A       AA 1
     C  ->  G       GG 1 or GA 1
     C  ->  T       TT 2
     T  ->  A       AA 1
     T  ->  G       GG 1 or GA 1
     T  ->  C       CC 1 or TC 2
    '''
    
    
    pairedCov=sum(baseDict.values())
    if pairedCov ==0:
        return ''
    
    if not IsSNP:
        if not Ismethy:return ''
        #no mutation
        if   refNuc =='G':
            if baseDict[('G','mC')]>=Mut_stand['G'][('G','mC')]['N']:
                methy=baseDict[('G','mC')]/pairedCov
                return ['-','-','-','-mC','%0.2f'%methy]
            else:
                return ''
        elif refNuc =='C':
            if baseDict[('mC','G')]>=Mut_stand['C'][('mC','G')]['N']:
                methy=baseDict[('mC','G')]/pairedCov
                return ['-','-','-','+mC','%0.2f'%methy]
            else:
                return ''
    else:
        mut=[]
        mut_VAF=[]
        if not Ismethy:
            #mutation but no methy
            for k in baseDict:
                if 'mC' in k:continue
                if k[0] == refNuc:continue
                if baseDict[k]>=Mut_stand[refNuc][k]['N']:
                    vaf=baseDict[k]/pairedCov
                    if vaf > minVaf:  
                        mut.append(k[0])
                        mut_VAF.append('%0.2f'%(vaf))
            if len(mut) ==0:return ''
            if len(mut)==1:
                if float(mut_VAF[0]) == 1:
                    return    [mut[0],mut_VAF[0],'homo','-','-']
                else:
                    return    [mut[0],mut_VAF[0],'het','-','-']
            else:
                return [",".join(mut),",".join(mut_VAF),"het","-","-"]
        else:

            #mutation and methy
            #mut=[]
            #mut_VAF=[]
            met=[]
            met_VAF=[]
            if refNuc in  ('A','T'):
                _C_flag=0
                _G_flag=0
                #Methylation due to mutation
                if baseDict[('T','A')] >=Mut_stand[refNuc][('T','A')]['N']:
                    #A>T
                    if refNuc!='T':
                        vaf=baseDict[('T','A')]/pairedCov
                        if  vaf > minVaf:
                            mut.append('T')
                            mut_VAF.append('%0.2f'%(vaf))
                
                if baseDict[('A','T')] >=Mut_stand[refNuc][('A','T')]['N']:
                    #T>A
                    if refNuc!='A':
                        vaf=baseDict[('A','T')]/pairedCov
                        if  vaf > minVaf:
                            mut.append('A')
                            mut_VAF.append('%0.2f'%(vaf))
                
                if baseDict[('G','mC')]>=Mut_stand[refNuc][('G','mC')]['N']:
                    met.append('-mC')
                    met_VAF.append('%0.2f'%(baseDict[('G','mC')]/(baseDict[('G','mC')]+baseDict[('G','C')])))
                    _G_flag=1
                if baseDict[('mC','G')]>=Mut_stand[refNuc][('mC','G')]['N']:
                    met.append('+mC')
                    met_VAF.append('%0.2f'%(baseDict[('mC','G')]/(baseDict[('mC','G')]+baseDict[('C','G')])))
                    _C_flag=1
                if baseDict[('C','G')]>=Mut_stand[refNuc][('C','G')]['N'] or _C_flag:
                    vaf=(baseDict[('mC','G')]+baseDict[('C','G')])/pairedCov
                    if vaf >minVaf:
                        mut.append('C')
                        mut_VAF.append('%0.2f'%(vaf))
                if baseDict[('G','C')]>=Mut_stand[refNuc][('G','C')]['N'] or _G_flag:
                    vaf=(baseDict[('G','mC')]+baseDict[('G','C')])/pairedCov
                    if vaf >minVaf:
                        mut.append('G')
                        mut_VAF.append('%0.2f'%(vaf))

                if len(mut)==0 :
                    if len(met)==0:
                        return ''
                    else:
                        return ['-','-',"-" ,','.join(met),','.join(met_VAF)]
                if len(mut)==1:
                    if float(mut_VAF[0]) ==1:
                        if len(met) ==0:
                            return [mut[0],mut_VAF[0],"homo",'-','-']
                        else:
                            return [mut[0],mut_VAF[0],"homo",','.join(met),','.join(met_VAF)]
                    else:
                        if len(met) ==0:
                            return [mut[0],mut_VAF[0],"het" ,'-','-']
                        else:
                            return [mut[0],mut_VAF[0],"het" ,','.join(met),','.join(met_VAF)]
                else:   
                    if len(met) ==0:
                        return [",".join(mut), ",".join(mut_VAF), "het", "-", "-"]
                    else:
                        return [",".join(mut), ",".join(mut_VAF), "het", ",".join(met), ",".join(met_VAF)]
            else:
                #G,C
                G_C_flag=0
                if refNuc=='G': 
                    if baseDict[('G','mC')]>=Mut_stand['G'][('G','mC')]['N']:
                        met.append('-mC')
                        met_VAF.append('%0.2f'%(baseDict[('G','mC')]/(baseDict[('G','mC')]+baseDict[('G','C')])))
                    if baseDict[('mC','G')]>=Mut_stand['G'][('mC','G')]['N']:
                        met.append('+mC')
                        met_VAF.append('%0.2f'%(baseDict[('mC','G')]/(baseDict[('mC','G')]+baseDict[('C','G')])))
                        G_C_flag=1
                    #mut.append('C')
                    if baseDict[('C','G')]>=Mut_stand['G'][('C','G')]['N'] or G_C_flag:
                        vaf=(baseDict[('mC','G')]+baseDict[('C','G')])/pairedCov
                        if vaf >minVaf:
                            mut.append('C')
                            mut_VAF.append('%0.2f'%(vaf))
                    #else:
                    #    mut_VAF.append('%0.2f'%(baseDict[('mC','G')]/pairedCov))
                    
                    if baseDict[('A','T')]>=Mut_stand['G'][('A','T')]['N']:
                        vaf=baseDict[('A','T')]/pairedCov
                        if vaf >minVaf:
                            mut.append('A')
                            mut_VAF.append('%0.2f'%(vaf))
                    if baseDict[('T','A')]>=Mut_stand['G'][('T','A')]['N']:
                        vaf=baseDict[('T','A')]/pairedCov
                        if vaf >minVaf:
                            mut.append('T')
                            mut_VAF.append('%0.2f'%(vaf))
                    #return ",".join(mut)+"\t"+",".join(map(str,mut_VAF))+"\t"+",".join(met)+"\t"+",".join(met_VAF)
                else:
                    C_G_flag=0
                    if baseDict[('mC','G')]>=Mut_stand['C'][('mC','G')]['N']:
                        met.append('+mC')
                        met_VAF.append('%0.2f'%(baseDict[('mC','G')]/(baseDict[('mC','G')]+baseDict[('C','G')])))
                    if baseDict[('G','mC')]>=Mut_stand['C'][('G','mC')]['N']:
                        met.append('-mC')
                        met_VAF.append('%0.2f'%(baseDict[('G','mC')]/(baseDict[('G','mC')]+baseDict[('G','C')])))
                        C_G_flag=1
                        #mut.append('G')

                    #2017-08-23
                    if baseDict[('G','C')]>=Mut_stand['C'][('G','C')]['N'] or C_G_flag:
                        vaf=(baseDict[('G','mC')]+baseDict[('G','C')])/pairedCov
                        if vaf >minVaf:
                            mut.append('G')
                            mut_VAF.append('%0.2f'%(vaf))
                    #else:
                    #    mut_VAF.append('%0.2f'%(baseDict[('G','mC')]/pairedCov))
                    
                    if baseDict[('A','T')]>=Mut_stand['C'][('A','T')]['N']:
                        vaf=baseDict[('A','T')]/pairedCov
                        if vaf>minVaf:
                            mut.append('A')
                            mut_VAF.append('%0.2f'%(vaf))
                    if baseDict[('T','A')]>=Mut_stand['C'][('T','A')]['N']:
                        vaf=baseDict[('T','A')]/pairedCov
                        if vaf >minVaf:
                            mut.append('T')
                            mut_VAF.append('%0.2f'%(vaf))
                
                if len(mut)==0:
                    if  len(met)==0: 
                        return '' 
                    else:
                        return ['-','-','-',",".join(met),",".join(met_VAF)]
                if len(mut) ==1:
                    if float(mut_VAF[0])==1:
                        if len(met)==0:
                            return [",".join(mut),mut_VAF[0],'homo',"-","-"]
                        else:
                            return [",".join(mut),mut_VAF[0],'homo',",".join(met),",".join(met_VAF)]
                    else:
                        if len(met)==0:
                            return [",".join(mut),mut_VAF[0],'het',"-","-"]
                        else:
                            return [",".join(mut),mut_VAF[0],'het',",".join(met),",".join(met_VAF)]
                else:
                    if len(met) == 0:
                        return [",".join(mut),",".join(mut_VAF),'het',"-","-"]
                    else:
                        return [",".join(mut),",".join(mut_VAF),'het',",".join(met),",".join(met_VAF)]

                    
                    
def sub_snpOutput(ID,START,END,BAM, GENOME,dbsnp, CHR, snpList,methyList,indelList,cpgList,chgList,chhList,cpgAppearList,cpgDisappearList,hemiMethyList,Mut_stand,maxDistance,maxLen,minLen,minReadQual,minReadQualN,maxN,maxSeqErr,maxSnp,maxIndel,maxBp,maxMut,minQual,minVaf,secAlign,strand,quite,debug):

    #global methyList,snpList
    global dbsnp_tb 
    if not quite:sys.stdout.write(colored('INOFR: {}\t{}\t{}\t{} is processesing\n'.format(CHR,ID,START,END),'green'))
    startT=time.time()
    openFile = Fasta(GENOME)
    
    SeqErrPoss={}
    
    allDict = defaultdict(lambda: {'gene_type': {('A','T'): 0, ('T','A'): 0, ('mC','G'): 0, ('G','mC'): 0, ('C','G'): 0,( 'G','C'): 0}, 'ref': '', "C_type": '-', 'Hemimethy': 0, 'all_Dep': 0, 'pair_dep': 0, 'mut_dep': 0, 'methy': 0})
    indelDict={k0:{k:{} for k in [CHR] } for k0 in ['insert','del']}

    #reflactNeuCount={('A','T'):0,('T','A'):0,('mC','G'):0,('G','mC'):0,('C','G'):0,('G','C'):0}
    reflactNeu={
    '+':{
         ('A','A'): ('A','T'),
         ('T','T'): ('T','A'),
         ('C','C'): ('mC','G'),
         ('G','G'): ('G','mC'),
         ('T','C'): ('C','G'),
         ('G','A'): ('G','C')},
    '-':{
         ('A','A'): ('A','T'),
         ('T','T'): ('T','A'),
         ('C','C'): ('mC','G'),
         ('G','G'): ('G','mC'),
         ('C','T'): ('C','G'),
         ('A','G'): ('G','C')}}

    #2017-09-10: dbsnp ID dict
    dbsnp_tb=getDbsnp(dbsnp)
    db_infors=defaultdict(lambda:".")
    for db_infor in dbsnp_tb.fetch(CHR, START, END):
        db_chr,db_pos,db_id,db_ref,db_alt=db_infor.strip().split()[:5]
        if len(db_ref)==len(db_alt):
            start = end = db_pos
            ref,alt=db_ref,db_alt
        else:
            start,end,ref,alt = vcf2annovar(int(db_pos),int(db_pos),db_ref,db_alt)
        db_infors[db_chr+":"+str(start)+":"+str(end)+":"+ref+":"+alt]=db_id
    if debug:print(ID,'db_infors is ok')


    C_type=defaultdict(lambda:'-')
    #sequence0=faseq(GENOME, CHR, START, END+3)  if END%100000==0  else faseq(GENOME, CHR, START, END)
    sequence0 = openFile[CHR][START:END+3].upper() if END%100000==0  else openFile[CHR][START:END].upper()
    sequence=sequence0[:-3]
    
    for i,s in enumerate(sequence):
        pos=START+i
        #CG:["CG"]
        #CHG:["CCG","CAG","CTG"]
        #CHH:["CCC","CCA","CCT","CAC","CAA","CAT","CTC","CTA","CTT"]
        if s!="C":continue
        triNeu=sequence0[i:i+3]
        if triNeu[1]=="G":C_type[pos]='CpG'    #  CG.append((pos,pos+1))
        elif triNeu[2]=="G":C_type[pos]='CHG' # ((pos,pos+2))
        else:C_type[pos]='CHH'                #  CHH.append((pos,pos+2))


    cigar_vlaue={0:'M',1:'I',2:'D',3:'S',4:'H'}
    samfile = pysam.AlignmentFile(BAM, "rb")
    tmp_dict={}

    #store duplicate 
    duplicate={'+':[],'-':[]}
    AA=0
    if debug:print('{} test{}'.format(ID,1))
    for read in samfile.fetch(CHR, START, END):
        #flush tmp dict
        site_tmp_dict=defaultdict(lambda:{'gene_type':{
        ('A','T'):0,('T','A'):0,('mC','G'):0,
        ('G','mC'):0,('C','G'):0,('G','C'):0},
        'ref':'','C_type':'-','Hemimethy':0,
        'seq_err':0,'all_Dep':0,'pair_dep':0,
        'mut_dep':0,'methy':0,'qual':(20,20)})
        
        ref_start=read.reference_start
        next_ref_start=read.next_reference_start
        
        
        # filter 1
        #check f optical or PCR duplicate
        if read.is_duplicate:
            if debug:
                print(read.qname,"f optical or PCR duplicate")
            continue

        #filter 2 
        #check is secondary
        if (not secAlign) and read.is_secondary:
            if debug:
                print(read.qname,"non-optimal")
            continue
        
        #filter 3
        #check paired reads  overlap

        if abs(next_ref_start-ref_start) >maxDistance: 
            if debug:
                print(read.qname," therer is no overlap between paired reads")
            continue
            #
            #to make False positive smaller, maxDistance should be smaller
            
        #filter 4
        # a read which contains more than a certain percentage of low-quality bases should be discarded
        qualities=read.query_qualities.tolist()    
        if len(list(filter(lambda x:x<minReadQual ,qualities))) / len(qualities) >minReadQualN:
            if debug:
                print(read.name,"too much low-quality bases")
            continue   
        
        #filter 5
        #a read which cigar contains S or H should be discarded
        cigarstring=read.cigarstring
        if 'S' in cigarstring or 'H' in cigarstring :
            #if debug:
            #    print(read.qname,'S or H in cigar')
            continue  
        
        #filter 6
        #a read read can have at most one insert or del
        if cigarstring.count('I') >maxIndel or cigarstring.count('D') >maxIndel:
            #if debug:
            #    print(read.qname,'too much indel')
            continue
        
        #filter 7
        #a read which contains too much N  should be discarded
        realseq=read.query_sequence.upper()  
        if realseq.count('N') >maxN:
            #if debug:
            #    print(read.qname,'too much N')
            continue
        
        #filter 8
        #a read that is too long or too short must be discarded
        readlen=read.rlen
        if readlen <= minLen or readlen>=maxLen :            
            #if debug:
            #    print(read.qname,'too long or short')
            continue
        
        
        #complete sequence and quality value through cigar
        cigar=read.cigar #[(0, 55), (2, 2), (0, 89)]
        cigarlsit=[]
        if len(cigar)>1:
            if cigar[1][0]=='2': 
                # del
                # quality 30
                # base    *
                for c in cigar:
                    for v,l in zip(c):
                        cigarlsit.extend([cigar_vlaue[v]]*l)
                qualities=qualities[:cigar[0][1]]+[30]*cigar[1][1]+qualities[-cigar[2][1]:]
            elif cigar[1][0]=='1': 
                #insert
                for c in cigar:
                    for v,l in zip(c):
                        cigarlsit.extend([cigar_vlaue[v]]*l)
                qualities=qualities[:cigar[0][1]]+qualities[-cigar[2][1]:]
                        
        else:
            cigarlsit=[cigar_vlaue[cigar[0][0]]]*cigar[0][1]


        chr=samfile.getrname(read.tid)

        readName=read.qname

        #print(read.query_alignment_sequence,read.get_reference_sequence())
        #refseq=read.get_reference_sequence()     #the refseq  displays del, hides insert
        refseq=read.query_alignment_sequence
        Stand = '-' if read.is_reverse else '+'
        ref_end=read.reference_end
        isRead1=read.is_read1
        mappingstrand=read.get_tags()[-1][-1]  #(bsmap ZS stand)
        #ST-E00144:480:H2KYYCCXY:4:2115:31913:53961 144 146 144 902187 902333 902189 255 [(0, 55), (2, 2), (0, 89)] 0 144
        #qstat qend  o based
        #[(0, 12), (1, 1), (0, 131)]    1 means I
        if readName not in tmp_dict:
            tmp_dict[readName]={
            'ref_start':ref_start,'ref_end':ref_end,
            'chr':chr,'realseq':realseq,'refseq':refseq,
            'qualities':qualities,'Stand':Stand,'readlen':readlen,
            'cigar':cigar,'cigarstring':cigarstring,
            'cigarlsit':cigarlsit,'mappingstrand':mappingstrand}
        else:
            #filter 9 
            # paired reads which have different stand should be discard
            if Stand != tmp_dict[readName]['Stand']:
                #if debug:
                #    print(readName,'paired read diff stand',Stand,tmp_dict[readName]['Stand'])
                del tmp_dict[readName]
                continue
            
            #filter 10
            if strand:
                Map_ERR_Stand=1
                #2017-10-23  #
                #When read1 is aligned to the W chain processed by BS, read2 should be aligned ro the complementary chain of C chain processed by BS
                #When read1 is aligned to the C chain processed by BS, read2 should be aligned ro the complementary chain of W chain processed by BS
                if isRead1:
                    if mappingstrand =="++" and  tmp_dict[readName]['mappingstrand'] =='--':
                        Map_ERR_Stand=0
                        
                    elif mappingstrand =="-+" and  tmp_dict[readName]['mappingstrand'] =='+-':
                        Map_ERR_Stand=0
                    #if debug:
                    #    print(readName+"\t"+mappingstrand+"\t"+tmp_dict[readName]['mappingstrand'])
                else:
                    if mappingstrand =="--" and  tmp_dict[readName]['mappingstrand'] =='++':
                        Map_ERR_Stand=0
                    elif mappingstrand =="+-" and  tmp_dict[readName]['mappingstrand'] =='-+':
                        Map_ERR_Stand=0
                    #if debug:
                    #    print(readName+"\t"+tmp_dict[readName]['mappingstrand']+"\t"+mappingstrand)
                if Map_ERR_Stand:
                    #if debug:
                    #    print(readName,"wrong aligned stand")
                    continue
                
            
            #filter 11
            # duplicate  again 
            #add a  duplicate function
            
            if (ref_start,ref_end) in duplicate[Stand]  or (tmp_dict[readName]['ref_start'],tmp_dict[readName]['ref_start']) in duplicate[Stand]  :
                del tmp_dict[readName]
                continue
            else:
                duplicate[Stand].append((ref_start,ref_end))
                duplicate[Stand].append((tmp_dict[readName]['ref_start'],tmp_dict[readName]['ref_start']))
                
            
            # get overlap
            #1 based [ )
            overlap_s = max(ref_start, tmp_dict[readName]['ref_start'])
            overlap_e = min(ref_end,   tmp_dict[readName]['ref_end']  )
            human_ref_seq = openFile[chr][overlap_s:overlap_e].upper()
            
            indel=0
            #if the paired read is inserted at the same time, it is judged that the inserted position is different.
            if 'I' in cigarstring and 'I' in tmp_dict[readName]['cigarstring']:
                
                pos=ref_start+cigar[0][1]
                if pos == tmp_dict[readName]['ref_start']+tmp_dict[readName]['cigar'][0][1] and cigar[1][1] ==tmp_dict[readName]['cigar'][1][1]:
                    if pos not in indelDict['insert'][chr]:indelDict['insert'][chr][pos]={'allcoverage':0,'paircoverage':0,'mutcoverage':0,'insetcon':{}}
                    if isRead1:
                        con1=realseq[cigar[0][1]:cigar[0][1]+cigar[1][1]]
                        con2=tmp_dict[readName]['realseq'][tmp_dict[readName]['cigar'][0][1]:tmp_dict[readName]['cigar'][0][1]+tmp_dict[readName]['cigar'][1][1]]
                    else:
                        con2=realseq[cigar[0][1]:cigar[0][1]+cigar[1][1]]
                        con1=tmp_dict[readName]['realseq'][tmp_dict[readName]['cigar'][0][1]:tmp_dict[readName]['cigar'][0][1]+tmp_dict[readName]['cigar'][1][1]]
                    insetcon=''
                    for r1,r2 in zip(con1,con2):
                        #print(Stand,r1,r2)
                        if (r1,r2) not in reflactNeu[Stand]:
                            if debug :
                                sys.stderr.write(colored('Warning: find a inser ERROR in read {}, {}, {}, {}, {}\n'.format(readName,pos,Stand,r1,r2),"yellow"))
                        else:
                            insetcon+=reflactNeu[Stand][(r1,r2)][0].replace('m','')
                    
                    indelDict['insert'][chr][pos]['allcoverage'] +=2
                    indelDict['insert'][chr][pos]['paircoverage']+=2
                    indelDict['insert'][chr][pos]['mutcoverage'] +=2
                    if insetcon  not in indelDict['insert'][chr][pos]['insetcon']:
                        indelDict['insert'][chr][pos]['insetcon'][insetcon]=2
                    else:
                        indelDict['insert'][chr][pos]['insetcon'][insetcon]+=2
                    indel=1
            elif 'D' in cigarstring and 'D' in tmp_dict[readName]['cigarstring']:
                pos=ref_start+cigar[0][1]
                if  pos == tmp_dict[readName]['ref_start']+tmp_dict[readName]['cigar'][0][1] and cigar[1][1] ==tmp_dict[readName]['cigar'][1][1]:
                    if pos not in indelDict['del'][chr]:indelDict['del'][chr][pos]={'allcoverage':0,'paircoverage':0,'mutcoverage':0,'delcon':{}}
                    delcon=ref_seq= openFile[chr][pos:pos+cigar[1][1]]
                    indelDict['del'][chr][pos]['allcoverage'] +=2
                    indelDict['del'][chr][pos]['paircoverage']+=2
                    indelDict['del'][chr][pos]['mutcoverage'] +=2
                    if delcon  not in indelDict['del'][chr][pos]['delcon']:
                        indelDict['del'][chr][pos]['delcon'][delcon]=2
                    else:
                        indelDict['del'][chr][pos]['delcon'][delcon]+=2
                    indel=1

                    
            # get new read1 and read2  from overlap
            poss=range(overlap_s,overlap_e)
            if isRead1:
                read1=refseq[overlap_s-ref_start:readlen-(ref_end-overlap_e)]  #[ref_start,ref_end] [0,readlen]
                read2=tmp_dict[readName]['refseq'][overlap_s-tmp_dict[readName]['ref_start']:tmp_dict[readName]['readlen']-(tmp_dict[readName]['ref_end']-overlap_e)]      
                #[tmp_dict[readName]['ref_start'],tmp_dict[readName]['ref_end']] 
                #[0,                              tmp_dict[readName]['readlen']]
                
                read1Qua=qualities[overlap_s-ref_start:readlen-(ref_end-overlap_e)]
                read2Qua=tmp_dict[readName]['qualities'][overlap_s-tmp_dict[readName]['ref_start']:tmp_dict[readName]['readlen']-(tmp_dict[readName]['ref_end']-overlap_e)]  
            else:
                read2=refseq[overlap_s-ref_start:readlen-(ref_end-overlap_e)]   #[ref_start,ref_end] [0,readlen]
                read1=tmp_dict[readName]['refseq'][overlap_s-tmp_dict[readName]['ref_start']:tmp_dict[readName]['readlen']-(tmp_dict[readName]['ref_end']-overlap_e)]
                
                read2Qua=qualities[overlap_s-ref_start:readlen-(ref_end-overlap_e)]
                read1Qua=tmp_dict[readName]['qualities'][overlap_s-tmp_dict[readName]['ref_start']:tmp_dict[readName]['readlen']-(tmp_dict[readName]['ref_end']-overlap_e)]
            

            SNP=0
            methy=0
            seqERROR=0
            SNPsite=[]
            tmpLen=len(read1)
            
            
            #check the quality of 10 bases in the head and 10 bases in the tail
            Start_filter_qual= 0 if len(list(filter(lambda x:x<minReadQual,read1Qua[:10])))/10 >minReadQualN or len(list(filter(lambda x:x<minReadQual,read2Qua[:10])))/10 >minReadQualN else 1
            End_filter_qual =0 if len(list(filter(lambda x:x<minReadQual,read1Qua[-10:])))/10 >minReadQualN or len(list(filter(lambda x:x<minReadQual,read2Qua[-10:])))/10 >minReadQualN else 1
            
            tmpDelPos=[]
            for n,r1,r2,r,p,q1,q2 in zip (range(len(read1)),read1,read2,human_ref_seq,poss,read1Qua,read2Qua):
                #n index , r1  read1 base, r2 read2 base, r ref base, p pos, q1  read1 quality, q2 read2 quality
                if r1 == '-' or  r2 == '-':
                    tmpDelPos.append(p)
                    continue
                    
                if (r1,r2) not in reflactNeu[Stand]:
                    seqERROR+=1
                    if p not in SeqErrPoss:SeqErrPoss[p]=''
                    site_tmp_dict[p]['seq_err']=1
                    if debug:
                        sys.stderr.write(colored('Warning: find a seqERROR in read {}, {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(readName,ref_start,n+1,Stand,r1,r2,q1,q2,read1,read2),"yellow"))
                    continue
                    
                realNeu=reflactNeu[Stand][(r1,r2)]
                

                site_tmp_dict[p]['gene_type'][realNeu]+=2
                site_tmp_dict[p]['ref']=r
                

                if p in indelDict['del'][chr]:
                    indelDict['del'][chr][p]['paircoverage']+=2
                if p in indelDict['insert'][chr]:
                    indelDict['insert'][chr][p]['paircoverage']+=2
                site_tmp_dict[p]['pair_dep']+=2
                

                if  r not in realNeu[0].replace('m',''):
                    SNP+=1
                    SNPsite.append(p)
                    #print(readName,p,r,realNeu)
                    #site_tmp_dict[p]['SNP']=1
                    site_tmp_dict[p]['qual']=(q1,q2)
                    site_tmp_dict[p]['mut_dep']+=2

                    #if n<=9:
                    #    #n 10 
                    #    if  Start_filter_qual:
                    #        SNP+=1
                    #        SNPsite.append(p)
                    #        site_tmp_dict[p]['qual']=(q1,q2)
                    #        site_tmp_dict[p]['mut_dep']+=2
                    #    else:
                    #        continue
                    #elif  (tmpLen-n)<10:
                    #    if End_filter_qual:
                    #        SNP+=1
                    #        SNPsite.append(p)
                    #        site_tmp_dict[p]['qual']=(q1,q2)
                    #        site_tmp_dict[p]['mut_dep']+=2
                    #    else:
                    #        continue
                    #elif checkContinueBase(human_ref_seq[n-7:n]) >=0.8:
                    #    #print('')
                    #    #print(list(filter(lambda x:x<minReadQual,read1Qua[n-5:n+3])),list(filter(lambda x:x<minReadQual,read2Qua[n-5:n+3])))
                    #    if len(list(filter(lambda x:x<minReadQual,read1Qua[n-5:n+3])))/8 >minReadQualN or len(list(filter(lambda x:x<minReadQual,read2Qua[n-5:n+3])))/8 >minReadQualN:
                    #        continue
                    #    else:
                    #        SNP+=1
                    #        SNPsite.append(p)
                    #        site_tmp_dict[p]['qual']=(q1,q2)
                    #        site_tmp_dict[p]['mut_dep']+=2
                    #else:
                    #    SNP+=1
                    #    SNPsite.append(p)
                    #    site_tmp_dict[p]['qual']=(q1,q2)
                    #    site_tmp_dict[p]['mut_dep']+=2

                if 'mC' in realNeu:
                    site_tmp_dict[p]['methy']+=2

                if  r == "G" and C_type[p-1]=="CpG":
                    if site_tmp_dict[p]['mut_dep']>0:continue  #  Prevent G mutate
                    if site_tmp_dict[p-1]['methy'] > 0 and 'mC' not in realNeu:
                        #C is methy but G not methy 
                        site_tmp_dict[p-1]['Hemimethy']+=2
                        
                    elif 'mC' in realNeu and  site_tmp_dict[p-1]['methy']== 0:
                        if site_tmp_dict[p-1]['seq_err']:                        
                            #Prevent C from having sequencing errors, causing C not to be methylated and G to have
                            continue
                        else:
                            
                            if n==0 or p-1 in tmpDelPos or site_tmp_dict[p-1]['mut_dep']>0:
                               #the first position of paired read shoul  not be G
                               #the c postion of CG should not be del(-)
                               #the C of CG shoul not mutated
                               #防止CG的C出现del(-)
                                continue
                            else:
                                #C is not methy but G is methy
                                site_tmp_dict[p]['Hemimethy']+=2

            #Filter 12
            #the read which contains much seqERRORs  should be discarded
            if seqERROR >= maxSeqErr:
                continue

            #Filter 13
            #a read which contains much mutations should be discarded
            if SNP >= maxSnp:
                continue

            #Filter 14
            #the number of mutations in a certain range should be smaller than a value
            if len(SNPsite) >=maxMut:
                tmppos=[]
                for i in range(len(SNPsite) -maxMut+1):
                    if SNPsite[i+maxMut-1]-SNPsite[i]  <=maxBp:
                        tmppos.extend(SNPsite[i:i+maxMut])
                tmppos=set(tmppos)
                if len(tmppos)>0:
                    for p in tmppos:
                        del site_tmp_dict[p]
                    if debug: sys.stdout.write(colored('Warning :the number of snps in {}bp range in read {} is greater than {}, {}\t{}\t{}\n'.format(maxBp,readName,maxMut),'yellow'))

            ## Filter 15
            #the quality of snp should be greater than a certain value
            for p in site_tmp_dict:
                if site_tmp_dict[p]['ref']=="":
                    continue
                allDict[p]['mut_dep']   += site_tmp_dict[p]['mut_dep']
                allDict[p]['pair_dep']  += site_tmp_dict[p]['pair_dep']
                allDict[p]['methy']     += site_tmp_dict[p]['methy']
                allDict[p]['Hemimethy'] += site_tmp_dict[p]['Hemimethy']
                allDict[p]['ref']        = site_tmp_dict[p]['ref']
                #-10*log(0.025,10) ~ 16.02059991327962
                #-10*log(0.05,10) ~ 13.01029995663981
                if site_tmp_dict[p]['qual'][0] >=minQual and site_tmp_dict[p]['qual'][1] >=minQual:
                    for gt in site_tmp_dict[p]['gene_type']:
                        allDict[p]['gene_type'][gt]+=site_tmp_dict[p]['gene_type'][gt]
                    #allDict[p]['ref']=site_tmp_dict[p]['ref']
                else:
                    if  debug:
                        sys.stderr.write(colored('Warning : the quality {} of the {}st postion in read {} is less than the quality {} of supporting call snp , discard\n'.format(site_tmp_dict[p]['qual'], p+1, readName, minQual),'yellow'))
                           
    SNPCon=''
    INDELCon=''
    MethyCon=''
    CpGCon='' 
    CHHCon=''
    CHGCon=''
    CpG_appear=''
    CpG_disappear=''
    HemimethyCon=''
    SeqErr_methy=0
    SeqErr_Snp=0
    SeqErr_hemimethy=0
    
    for pos in sorted(allDict.keys()):
        if pos<START or pos>END:continue
        if allDict[pos]['ref']=='N':continue 
        con=pairCheckSnp(Mut_stand,allDict[pos]['gene_type'],allDict[pos]['ref'],allDict[pos]['mut_dep'],allDict[pos]['methy'],minVaf)

        if allDict[pos]['ref'] == 'C':
            if C_type[pos]=="CpG":
                #the depth of C in CpG is 0
                if allDict[pos]['pair_dep'] ==0: 
                    if allDict[pos+1]['pair_dep'] !=0:
                        # #the depth of G in CpG is gt 0
                        tmp_methy=allDict[pos+1]['methy']/allDict[pos+1]['pair_dep']
                        tmp_dep=allDict[pos+1]['pair_dep']
                        CpGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G",'%0.2f'%(tmp_methy),int(tmp_dep))
                    else:
                        #the depth of the CpG is 0
                        CpGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G","NA","NA")
                else:
                    if allDict[pos+1]['pair_dep'] == 0:
                        tmp_methy=allDict[pos]['methy']/allDict[pos]['pair_dep']
                        tmp_dep=allDict[pos]['pair_dep']
                        CpGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G",'%0.2f'%(tmp_methy),int(tmp_dep))
                    else:
                        #the depth of C and G in CpG are both gt 0
                        tmp_methy=(allDict[pos]['methy']/allDict[pos]['pair_dep']+allDict[pos+1]['methy']/allDict[pos+1]['pair_dep'])/2
                        tmp_dep=(allDict[pos+1]['pair_dep']+allDict[pos]['pair_dep'])/2
                        CpGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G",'%0.2f'%(tmp_methy),int(tmp_dep))
                        if allDict[pos]['Hemimethy']>0 :
                            if allDict[pos+1]['Hemimethy']>0:
                                #C and G  both  Hemimethy 
                                #tmp_Hemimethy1 = allDict[pos]['Hemimethy']/allDict[pos]['pair_dep']
                                #tmp_Hemimethy2 = allDict[pos+1]['Hemimethy']/allDict[pos+1]['pair_dep']
                                tmp_Hemimethy=(allDict[pos]['Hemimethy']+allDict[pos+1]['Hemimethy'])/allDict[pos]['pair_dep']
                                HemimethyCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G","Hemimethy_C_G",'%0.2f'%(tmp_Hemimethy),allDict[pos]['pair_dep'])
                                if pos in SeqErrPoss:SeqErr_hemimethy+=1
                            else:
                                #C Hemimethy
                                tmp_Hemimethy=allDict[pos]['Hemimethy']/allDict[pos]['pair_dep']
                                HemimethyCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G","Hemimethy_C",'%0.2f'%(tmp_Hemimethy),allDict[pos]['pair_dep'])
                                if pos in SeqErrPoss:SeqErr_hemimethy+=1
   
                        elif allDict[pos+1]['Hemimethy']>0:
                            #G Hemimethy
                            tmp_Hemimethy=allDict[pos+1]['Hemimethy']/allDict[pos+1]['pair_dep']
                            HemimethyCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G","Hemimethy_G",'%0.2f'%(tmp_Hemimethy),allDict[pos+1]['pair_dep'])
                            if pos in SeqErrPoss:SeqErr_hemimethy+=1
                        else:
                            #there is no hemimethy in the CpG
                            HemimethyCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","G","non_Hemimethy",'0',int((allDict[pos]['pair_dep']+allDict[pos+1]['pair_dep'])/2))
            elif C_type[pos]=='CHG':
                if allDict[pos]['pair_dep'] !=0:
                    tmp_methy=allDict[pos]['methy']/allDict[pos]['pair_dep']
                    CHGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","CHG",'%0.2f'%(tmp_methy),allDict[pos]['pair_dep'])
                else:
                    CHGCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","CHG","NA","NA")
            else:
                #C_type[pos]=='CHH':
                if allDict[pos]['pair_dep'] !=0:
                    tmp_methy=allDict[pos]['methy']/allDict[pos]['pair_dep']
                    CHHCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","CHH",'%0.2f'%(tmp_methy),allDict[pos+1]['pair_dep'])
                else:
                    CHHCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+2,"C","CHH","NA","NA")
        if con:
            if con[0]=='-':
                MethyCon+='{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]))
                if pos in SeqErrPoss:SeqErr_methy+=1
            else:
            
                #20201002
                #filter 16
                #homo  See the true and false positive distribution graph of chromosome X for the reason
                #het  See the true and false positive distribution graph of ATGC for the reason
                if con[2] =="homo":
                    if sum([allDict[pos]['gene_type'][gt] for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]) <=4:
                        continue
                elif con[2] =="het":
                    if allDict[pos]['ref']=="A":
                        if sum([allDict[pos]['gene_type'][gt] for gt in [('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]) <=4:
                            continue
                    elif allDict[pos]['ref']=="G":
                        if sum([allDict[pos]['gene_type'][gt] for gt in [('A','T'),('T','A'),('mC','G'),('C','G')]]) <=4:
                            continue
                    elif allDict[pos]['ref']=="C":
                        if sum([allDict[pos]['gene_type'][gt] for gt in [('A','T'),('T','A'),('G','mC'),('G','C')]]) <=4:
                            continue
                    elif allDict[pos]['ref']=="T":
                        if sum([allDict[pos]['gene_type'][gt] for gt in [('A','T'),('mC','G'),('G','mC'),('C','G'),('G','C')]]) <=4:
                            continue
                db_key=chr+":"+str(pos+1)+":"+str(pos+1)+":"+allDict[pos]['ref']+":"+con[0]   
                db_id=db_infors[db_key]
                if  debug:
                    debugCout=[0]*8
                    debugCon="+A:+T:+C:+G:-A:-T:-C:-G"
                    debugCout[0] =allDict[pos]['gene_type'][('A','T')]    #+A
                    debugCout[5] =allDict[pos]['gene_type'][('A','T')]    #-T
                    debugCout[1] =allDict[pos]['gene_type'][('T','A')]    #+T
                    debugCout[4] =allDict[pos]['gene_type'][('T','A')]    #-A
                    debugCout[2]+=allDict[pos]['gene_type'][('mC','G')]  #+C
                    debugCout[2]+=allDict[pos]['gene_type'][('C','G')]   #+C
                    debugCout[3]+=allDict[pos]['gene_type'][('G','mC')]  #+G
                    debugCout[3]+=allDict[pos]['gene_type'][('G','C')]   #+G
                    debugCout[6]+=allDict[pos]['gene_type'][('G','C')]   #-C
                    debugCout[6]+=allDict[pos]['gene_type'][('G','mC')]  #-C
                    debugCout[7]+=allDict[pos]['gene_type'][('C','G')]   #-G
                    debugCout[7]+=allDict[pos]['gene_type'][('mC','G')]  #-G
                    debugCout=":".join([str(iii) for iii in debugCout])

                    SNPCon+='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,db_id,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]),debugCon,debugCout)
                else:
                    if pos in SeqErrPoss:SeqErr_Snp+=1
                    SNPCon+='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,db_id,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]))

                    if C_type[pos]=="CpG" or C_type[pos-1]=="CpG" :
                        CpG_disappear+='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,db_id,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]),"".join([allDict[pp]['ref'] for pp in range(pos-2,pos+3)]))     #2018-01-31

                    if allDict[pos-1]['ref']=="C" and "G" in con[0] or allDict[pos+1]['ref']=="G" and "C" in con[0]:
                        CpG_appear+=   '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,db_id,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]),"".join([allDict[pp]['ref'] for pp in range(pos-2,pos+3)]))            #2018-01-29

                    #old CpG disapper,new CpG apper  TCGGT > TCCGT

                if con[3] != '-':
                    MethyCon+='{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,allDict[pos]['ref'],"\t".join(con),'A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]))
                    if pos in SeqErrPoss:SeqErr_methy+=1
        else:
            if debug:
                debugCout=[0]*8
                debugCon="+A:+T:+C:+G:-A:-T:-C:-G"
                debugCout[0]= allDict[pos]['gene_type'][('A','T')]    #+A
                debugCout[5]= allDict[pos]['gene_type'][('A','T')]    #-T
                debugCout[1]= allDict[pos]['gene_type'][('T','A')]    #+T
                debugCout[4]= allDict[pos]['gene_type'][('T','A')]    #-A
                debugCout[2]+=allDict[pos]['gene_type'][('mC','G')]  #+C
                debugCout[2]+=allDict[pos]['gene_type'][('C','G')]   #+C
                debugCout[3]+=allDict[pos]['gene_type'][('G','mC')]  #+G
                debugCout[3]+=allDict[pos]['gene_type'][('G','C')]   #+G
                debugCout[6]+=allDict[pos]['gene_type'][('G','C')]   #-C
                debugCout[6]+=allDict[pos]['gene_type'][('G','mC')]  #-C
                debugCout[7]+=allDict[pos]['gene_type'][('C','G')]   #-G
                debugCout[7]+=allDict[pos]['gene_type'][('mC','G')]  #-G
                debugCout=":".join([str(iii) for iii in debugCout])
                SNPCon+='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chr,pos+1,".",allDict[pos]['ref'],'-\t-\t-\t-\t-','A_T:T_A:mC_G:G_mC:C_G:G_C',":".join([str(allDict[pos]['gene_type'][gt]) for gt in [('A','T'),('T','A'),('mC','G'),('G','mC'),('C','G'),('G','C')]]),debugCon,debugCout)



    indel='insert'
    for chr in indelDict[indel]:
        if len(indelDict[indel][chr]) ==0 :continue
        for pos in sorted(indelDict[indel][chr].keys()):
            if pos<START or pos>END:continue
            for insetcon in indelDict[indel][chr][pos]['insetcon']:
                db_key=chr+":"+str(pos+1)+":"+str(pos+1)+":"+'-'+":"+insetcon
                db_id=db_infors[db_key]
                INDELCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,db_id,'-',insetcon,indelDict[indel][chr][pos]['mutcoverage'],indelDict[indel][chr][pos]['paircoverage'],indelDict[indel][chr][pos]['allcoverage'])
    indel='del'
    for chr in indelDict[indel]:
        if len(indelDict[indel][chr]) ==0 :continue
        for pos in sorted(indelDict[indel][chr].keys()):
            if pos<START or pos>END:continue 
            for delcon in indelDict[indel][chr][pos]['delcon']:
                db_key=chr+":"+str(pos+1)+":"+str(pos+1+len(delcon)-1)+":"+delcon+":"+"-"
                db_id=db_infors[db_key]
                #INDELCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,pos+1+len(delcon)-1,delcon,'-',indelDict[indel][chr][pos]['mutcoverage'],indelDict[indel][chr][pos]['paircoverage'],indelDict[indel][chr][pos]['allcoverage'])
                INDELCon+="{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,pos+1,db_id,delcon,'-',indelDict[indel][chr][pos]['mutcoverage'],indelDict[indel][chr][pos]['paircoverage'],indelDict[indel][chr][pos]['allcoverage'])


    if not quite:
        sys.stdout.write(colored('INFOR: thread {},  snpcons                            {}.\n'.format(ID,SNPCon.count("\n")),       'green'))
        sys.stdout.write(colored('INFOR: thread {},  indelcons                          {}.\n'.format(ID,INDELCon.count("\n")),     'green'))
        sys.stdout.write(colored('INFOR: thread {},  methycons                          {}.\n'.format(ID,MethyCon.count("\n")),     'green'))
        sys.stdout.write(colored('INFOR: thread {},  CpGCon                             {}.\n'.format(ID,CpGCon.count("\n")),       'green'))
        sys.stdout.write(colored('INFOR: thread {},  CHGCon                             {}.\n'.format(ID,CHGCon.count("\n")),       'green'))
        sys.stdout.write(colored('INFOR: thread {},  CHHCon                             {}.\n'.format(ID,CHHCon.count("\n")),       'green'))
        sys.stdout.write(colored('INFOR: thread {},  CpG_appear                         {}.\n'.format(ID,CpG_appear.count("\n")),   'green'))
        sys.stdout.write(colored('INFOR: thread {},  CpG_disappear                      {}.\n'.format(ID,CpG_disappear.count("\n")),'green'))
        sys.stdout.write(colored('INFOR: thread {},  HemimethyCon                       {}.\n'.format(ID,HemimethyCon.count("\n")), 'green'))
        sys.stdout.write(colored('INFOR: thread {},  seqERROR                           {}.\n'.format(ID,len(SeqErrPoss)),          'green'))
        sys.stdout.write(colored('INFOR: thread {},  the seqERROR number in Hemimethy   {}.\n'.format(ID,SeqErr_hemimethy),         'green'))
        sys.stdout.write(colored('INFOR: thread {},  the seqERROR number in variation   {}.\n'.format(ID,SeqErr_Snp),               'green'))
        sys.stdout.write(colored('INFOR: thread {},  the seqERROR number in methyaltion {}.\n'.format(ID,SeqErr_methy),             'green'))


    snpList[ID]          = SNPCon
    methyList[ID]        = MethyCon
    indelList[ID]        = INDELCon
    cpgList[ID]          = CpGCon
    chgList[ID]          = CHGCon
    chhList[ID]          = CHHCon
    cpgAppearList[ID]   = CpG_appear
    cpgDisappearList[ID]= CpG_disappear
    hemiMethyList[ID]    = HemimethyCon

    if not quite:sys.stdout.write(colored('INFOR: {} {} was finished! Used time are {}S!\n'.format(CHR,ID,int(time.time()-startT)), 'green'))




def snpOutput(ref, bam, genomeFile, dbsnp, CHR, outDir, cpu=50, maxDistance=50, maxLen=200, minLen=50, minReadQual=20, minReadQualN=0.6, maxN=5, maxSeqErr=10, maxSnp=5, maxIndel=1, maxBp=50, maxMut=3, minQual=20, minVaf=0.01, secAlign=False, strand=False, quite=False, debug=False):

    global methyList,snpList
    #hg19 chromosome  length
    if  ref == "hg19":
        CHRs={'chr1':249250621,'chr2':243199373,'chr3':198022430,'chr4':191154276,'chr5':180915260,'chr6':171115067,'chr7':159138663,'chr8':146364022,'chr9':141213431,'chr10':135534747,'chr11':135006516,'chr12':133851895,'chr13':115169878,'chr14':107349540,'chr15':102531392,'chr16':90354753,'chr17':81195210,'chr18':78077248,'chr19':59128983,'chr20':63025520,'chr21':48129895,'chr22':51304566,'chrX':155270560,'chrY':59373566,'chrM':11850}    
    elif ref == "hg38":
        CHRs={"chr1":248956422, "chr2":242193529, "chr3":198295559, "chr4":190214555, "chr5":181538259, "chr6":170805979, "chr7":159345973, "chr8":145138636, "chr9":138394717, "chr10":133797422, "chr11":135086622, "chr12":133275309, "chr13":114364328, "chr14":107043718, "chr15":101991189, "chr16":90338345, "chr17":83257441, "chr18":80373285, "chr19":58617616, "chr20":64444167, "chr21":46709983, "chr22":50818468, "chrX":156040895, "chrY":57227415, "chrM":16569}



    checkFile(bam)
    checkFile(genomeFile)
    checkFile(dbsnp)

    root = os.path.splitext(os.path.basename(bam))[0]
    if outDir =='':
        outDir=os.path.dirname(bam)
        outDir=os.path.dirname(outDir)+"/snpmeth" if outDir !='' else './snpmeth'
    checkDir(outDir)


    #getDbsnp(dbsnp)

    Mut_stand={
    'A':{
        ('A','T'): {'N':2,'con':''            },
        ('G','mC'):{'N':2,'con':'A\t>\tG\t-mC'},
        ('G','C'): {'N':4,'con':'A\t>\tG\t-'  },
        ('mC','G'):{'N':2,'con':'A\t>\tC\t+mC'},
        ('C','G'): {'N':2,'con':'A\t>\tC\t-'  },
        ('T','A'): {'N':2,'con':'A\t>\tT\t-'  }
        },
    'G':{
        ('G','C'): {'N':2,'con':''            },
        ('A','T'): {'N':4,'con':'G\t>\tA\t-'  },
        ('mC','G'):{'N':2,'con':'G\t>\tC\t+mC'},
        ('C','G'): {'N':2,'con':'G\t>\tC\t-'  },
        ('T','A'): {'N':2,'con':'G\t>\tT\t-'  },
        ('G','mC'):{'N':4,'con':'G\t-\t-\t-mC'}
        },
    'C':{
        ('C','G'): {'N':2,'con':''            },
        ('A','T'): {'N':2,'con':'C\t>\tA\t-'  },
        ('G','mC'):{'N':2,'con':'C\t>\tG\t-mC'},
        ('G','C'): {'N':2,'con':'C\t>\tG\t-'  },
        ('T','A'): {'N':4,'con':'C\t>\tT\t-'  },
        ('mC','G'):{'N':4,'con':'C\t-\t-\t+mC'}
    },
    'T':{
        ('T','A'): {'N':2,'con':''            },
        ('A','T'): {'N':2,'con':'T\t>\tA\t-'  },
        ('G','mC'):{'N':2,'con':'T\t>\tG\t-mC'},
        ('G','C'): {'N':2,'con':'T\t>\tG\t-'  },
        ('mC','G'):{'N':2,'con':'T\t>\tC\t+mC'},
        ('C','G'): {'N':4,'con':'T\t>\tC\t-'  }
        }
    }



    #bam = pysam.Samfile(bamfile, 'rb')

    sys.stdout.write("\n")
    sys.stdout.write(colored('INFOR:  version:                     {}\n'.format(__version__ ), 'green'))
    sys.stdout.write(colored('INFOR:  date:                        {}\n'.format(__date__    ), 'green'))
    sys.stdout.write(colored('INFOR:  author:                      {}\n'.format(__author__  ), 'green'))
    sys.stdout.write(colored('INFOR:  email:                       {}\n'.format(__email__   ), 'green'))
    sys.stdout.write(colored("INFOR:  bam:                         {}\n".format(bam         ), 'green'))
    sys.stdout.write(colored("INFOR:  reference genome:            {}\n".format(genomeFile  ), 'green'))
    sys.stdout.write(colored("INFOR:  dbsnp:                       {}\n".format(dbsnp       ), 'green'))
    sys.stdout.write(colored("INFOR:  chromosome:                  {}\n".format(CHR         ), 'green'))
    sys.stdout.write(colored("INFOR:  outDir:                      {}\n".format(outDir      ), 'green'))
    sys.stdout.write(colored("INFOR:  maxDistance:                 {}\n".format(maxDistance ), 'green'))
    sys.stdout.write(colored("INFOR:  maxLen:                      {}\n".format(maxLen      ), 'green'))
    sys.stdout.write(colored("INFOR:  minLen:                      {}\n".format(minLen      ), 'green'))
    sys.stdout.write(colored("INFOR:  minReadQual:                 {}\n".format(minReadQual ), 'green'))
    sys.stdout.write(colored("INFOR:  minReadQualN:                {}\n".format(minReadQualN), 'green'))
    sys.stdout.write(colored("INFOR:  maxN:                        {}\n".format(maxN        ), 'green'))
    sys.stdout.write(colored("INFOR:  maxSeqErr:                   {}\n".format(maxSeqErr   ), 'green'))
    sys.stdout.write(colored("INFOR:  maxSnp:                      {}\n".format(maxSnp      ), 'green'))
    sys.stdout.write(colored("INFOR:  maxIndel:                    {}\n".format(maxIndel    ), 'green'))
    sys.stdout.write(colored("INFOR:  maxBp:                       {}\n".format(maxBp       ), 'green'))
    sys.stdout.write(colored("INFOR:  maxMut:                      {}\n".format(maxMut      ), 'green'))
    sys.stdout.write(colored("INFOR:  minQual:                     {}\n".format(minQual     ), 'green'))
    sys.stdout.write(colored("INFOR:  minVaf:                      {}\n".format(minVaf      ), 'green'))
    sys.stdout.write(colored("INFOR:  secAlign:                    {}\n".format(secAlign    ), 'green'))
    sys.stdout.write(colored("INFOR:  strand:                      {}\n".format(strand      ), 'green'))
    sys.stdout.write(colored("INFOR:  quite:                       {}\n".format(quite       ), 'green'))
    sys.stdout.write(colored("INFOR:  debug:                       {}\n".format(debug       ), 'green'))
    sys.stdout.write(colored("INFOR:  starting time:               {}\n".format(getTime()  ), 'green'))
    sys.stdout.write(colored("INFOR:  starting to deal chromosome  {}\n".format(CHR         ), 'green'))


    #10000bp as a region
    nu                 = int(CHRs[CHR]/100000)
    IDs                = range(nu+1)
    manager            = Manager()
    snpList            = manager.list(['-']*(nu+1))
    methyList          = manager.list(['-']*(nu+1))
    indelList          = manager.list(['-']*(nu+1))
    cpgList            = manager.list(['-']*(nu+1))
    chgList            = manager.list(['-']*(nu+1))
    chhList            = manager.list(['-']*(nu+1))
    cpgAppearList      = manager.list(['-']*(nu+1))
    cpgDisappearList   = manager.list(['-']*(nu+1))
    hemiMethyList      = manager.list(['-']*(nu+1))

    pool=Pool(processes=cpu)
    for ID in IDs:
        #0 base
        START=ID*100000
        if ID != nu:
            END=(ID+1)*100000
        else:
            END=CHRs[CHR]
        pool.apply_async(sub_snpOutput, (ID,START,END,bam,genomeFile,dbsnp,CHR,snpList,methyList,indelList,cpgList,chgList,chhList,cpgAppearList,cpgDisappearList,hemiMethyList,Mut_stand, maxDistance,maxLen,minLen,minReadQual,minReadQualN,maxN,maxSeqErr,maxSnp,maxIndel,maxBp,maxMut,minQual,minVaf,secAlign,strand,quite,debug))
    pool.close()
    pool.join()



    snpFile          = open('{}/{}.Snp.txt'.format(outDir,root),            'w')
    indelFile        = open('{}/{}.Indel.txt'.format(outDir,root),          'w')
    MethFile         = open('{}/{}.Methylation.txt'.format(outDir,root),    'w')
    hemiMethyFile    = open('{}/{}.Hemimethylation.txt'.format(outDir,root),'w')
    cpgAppearFile    = open('{}/{}.CpG_appear.txt'.format(outDir,root),     'w')
    cpgDisappearFile = open('{}/{}.CpG_disappear.txt'.format(outDir,root),  'w')
    cpgFile          = open('{}/{}.cpgFile.txt'.format(outDir,root),        'w')
    chgFile          = open('{}/{}.chgFile.txt'.format(outDir,root),        'w')
    chhFile          = open('{}/{}.chhFile.txt'.format(outDir,root),        'w')


    for i in range(nu+1):
        #if i !=3:continue
        if snpList[i] =='-' :
            if i ==nu:
                sys.stderr.write(colored('Error: snpList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: snpList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
            #continue
        else:
            snpFile.write(snpList[i])

        if indelList[i] =='-' :
            if i ==nu:
                sys.stderr.write(colored('Error: indelList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: indelList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
            #continue
        else:
            indelFile.write(indelList[i])
       
        if  methyList[i] == '-':
            if i ==nu:
                sys.stderr.write(colored('Error: methyList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: methyList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            MethFile.write(methyList[i])
        if cpgList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: cpgList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: cpgList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            cpgFile.write(cpgList[i])
        if chgList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: chgList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: chgList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            chgFile.write(chgList[i])
        if chhList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            chhFile.write(chhList[i])

        if cpgAppearList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            cpgAppearFile.write(cpgAppearList[i])

        if cpgDisappearList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: chhList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            cpgDisappearFile.write(cpgDisappearList[i])
            
        if hemiMethyList[i] == '-':
            if i==nu:
                sys.stderr.write(colored('Error: hemiMethyList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,CHRs[CHR]),'red'))
            else:
                sys.stderr.write(colored('Error: hemiMethyList, {} {} {}-{} is wrong!\n'.format(CHR,i,i*100000,(i+1)*100000),'red'))
        else:
            hemiMethyFile.write(hemiMethyList[i])
            

    indelFile.close()
    snpFile.close()
    MethFile.close()
    cpgFile.close()
    chgFile.close()
    chhFile.close()
    cpgAppearFile.close()
    cpgDisappearFile.close()
    hemiMethyFile.close()


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v','--version', action='version', version='DSBS.py {}'.format(__version__))
    parser.add_argument('bam',  help='input BAM file')
    parser.add_argument('--maxDistance',      type=int,   default=50,   help='the maxmum distance of the start postions of the paired reads in reference genome, default is 50')
    parser.add_argument('--maxLen',           type=int,   default=200,  help='the maxmum length of a read, default is 200')
    parser.add_argument('--minLen',           type=int,   default=50 ,  help='the minimum length of a read, default is 50')
    parser.add_argument('--minReadQual',      type=int,   default=20 ,  help='the minimum quality for assessing the read overall qualities, default is 20')
    parser.add_argument('--minReadQualN',     type=float, default=0.6,  help='the minimum ratio which quality is greater than the minReadQual, default is 0.6')
    parser.add_argument('--maxN',             type=int,   default=5 ,   help='the maxmum number of N within a read, default is 5')
    parser.add_argument('--maxSeqErr',        type=int,   default=10,   help='the maxmum number of sequencing errors within a paired reads, default is 10')
    parser.add_argument('--maxSnp',           type=int,   default=5 ,   help='the maxmum number of snp within a paired reads, default is 5')
    parser.add_argument('--maxIndel',         type=int,   default=1 ,   help='the maxmum number of indel within a paired reads, default is 1')
    parser.add_argument('--maxBp',            type=int,   default=50 ,  help='a certain range for allowing the currence of certain mutation ,default is 50')
    parser.add_argument('--maxMut',           type=int,   default=3 ,   help='the maxmum number of allowing the currence of mutation within a certain range ,default is 3')
    parser.add_argument('--minQual',          type=int,   default=20 ,  help='the minimum quality for call Variation, default is 20')
    parser.add_argument('--minVaf',           type=float, default=0.1,  help='The minimum value of allele, default is 0.1')
    parser.add_argument('--window',           type=int,   default=5 ,   help='ignoring the snps in a certain range of a indel, the range is 5bp')
    parser.add_argument('--minCoverage',      type=int,   default=8,    help='coverage or minimum number of reads desired')
    parser.add_argument('--maxCoverage',      type=int,   default=500,  help='coverage or maxmum number of reads desired')
    parser.add_argument('--strand',      action='store_true',default=False , help='only consider the correct comparison direction, read1==++ && read2==-- || read1==-+ && read2==+-')
    parser.add_argument('--secAlign',    action='store_true',default=False , help='consider non-optimal alignments')
    parser.add_argument('--CpG',         action='store_true',default=False , help='output CpG')
    parser.add_argument('--CHG',         action='store_true',default=False , help='output CHG')
    parser.add_argument('--CHH',         action='store_true',default=False , help='output CHH')
    parser.add_argument('--debug',       action='store_true',default=False , help='dubug mode')
    parser.add_argument('-q','--quite',  action='store_true',default=False , help='quiet mode')
    parser.add_argument('-p', '--cpu',    type=int, default=50, help='the number of working threads, default is 50')
    parser.add_argument('-o', '--outDir', default='', help='the outDir')
    parser.add_argument('-g', '--genomeFile', required=True, help='the chromosome reference fasta')
    parser.add_argument('-d', '--dbsnp',      required=True, help='the dbsnp file')
    parser.add_argument('--Chr',  default='chr1', help='chromosome')
    parser.add_argument('--ref',  default='hg19',choices=['hg19','hg38'], help='genome reference')
    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    if args.bam=='':                                parser.error("Missing bam file")
    if len(args.genomeFile) == 0:                   parser.error("Missing reference file, use -g or --genomeFile option.")
    if len(args.dbsnp) == 0:                        parser.error("Missing dbsnp file, use -d or --dbsnp option.")
    if args.cpu <0 :                                parser.error("Invalid -p value,             must >= 1")
    if args.maxDistance<0 :                         parser.error("Invalid --maxDistance value,  must >= 1")
    if args.maxLen <0 or args.maxLen >300 :         parser.error("Invalid --maxLen value,       must >= 1 and <300")
    if args.minLen<0 or args.minLen >300:           parser.error("Invalid --minLen value,       must >= 1 and <300")
    if args.minReadQual <0 :                        parser.error("Invalid --minReadQual value,  must >= 0")
    if 1<=args.minReadQualN or args.minReadQualN<0: parser.error("Invalid --minReadQualN value, must >= 0 and <1")
    if args.maxN <0 or args.maxN >50:               parser.error("Invalid --maxN value,         must >= 0 and <50")
    if args.maxSeqErr<0:                            parser.error("Invalid --maxSeqErr value,    must >= 0 and <50")
    if args.maxSnp<0   :                            parser.error("Invalid --maxSnp    value,    must >= 0 and <50")
    if args.maxIndel<0 :                            parser.error("Invalid --maxIndel  value,    must >= 0 and <50")
    if args.maxBp<0    :                            parser.error("Invalid --maxBp     value,    must >= 0 and <50")
    if args.maxMut<0   :                            parser.error("Invalid --maxMut    value,    must >= 0 and <50")
    if args.minQual<0  :                            parser.error("Invalid --minQual   value,    must >= 0 ")
    if args.minVaf <0 or args.minVaf >1:            parser.error("Invalid --minVaf    value,    must >= 0 and <1")
    snpOutput(args.ref, args.bam, args.genomeFile,args.dbsnp, args.Chr,args.outDir,args.cpu,args.maxDistance,args.maxLen,args.minLen,args.minReadQual,args.minReadQualN,args.maxN,args.maxSeqErr,args.maxSnp,args.maxIndel,args.maxBp,args.maxMut,args.minQual,args.minVaf,args.secAlign,args.strand,args.quite,args.debug)
    
if __name__ == '__main__':
    main()
