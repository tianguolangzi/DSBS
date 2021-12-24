#!/public/software/python3
#coding=utf-8


from termcolor import colored
import sys, os, gzip, time
from os import path as op
from collections import defaultdict
import subprocess
import argparse
from multiprocessing import Pool
def purpose():
    sys.stdout.write(colored('''                                                
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
    \n''','green'))
    




def checkDir(path):
    if not op.exists(path):
        try:
            os.makedirs(path)
        except :
            sys.stderr.write(colored("ERROR: failed to mkdir {} !\n".format(path),'red'))
            sys.exit(1)

def checkFile(File):
    if not op.exists(File):
        sys.stderr.write(colored('ERROR: {} is not existed! \n'.format(File),'red'))
        sys.exit(1)
    else:
        return 1


def writeLog(File,Con):
    with open(File,'w') as A1:
        A1.write(Con)

def systemRun(command,quite=0):
    if not quite:sys.stdout.write(colored('command: {}\n'.format(command),'green'))
    if os.system(command):
        sys.stderr.write(colored('EEROR: {} \n'.format(command),'red'))
        sys.exit(1)
    else:
        pass

def systemRun2(logs,commands,pp=1):

    print(logs,commands)
    popens={}
    tmppopen=[]
    i=0
    Flag=0
    Nu=len(commands)
    if pp >len(commands):
        pp=len(commands)
    while Nu > 0:
        if i != len(commands) :
            for file in logs:
                if len(popens)<pp:
                    if file not in tmppopen:
                        print(commands[file])
                        sys.stdout.write(colored('command: {}\n'.format(commands[file]),'green'))
                        popens[file]=subprocess.Popen(commands[file], shell=True)
                        tmppopen.append(file)
                        i+=1
        time.sleep(5)
        for file  in tmppopen:
            if file not in popens:continue
            if popens[file].poll() ==0:
                if popens[file].returncode != 0:
                    sys.stderr.write("ERROR: "+ commands[file] +" \n")
                    Flag=1
                else:
                    sys.stdout.write(colored("{} was finished\n".format(file),'green'))
                    writeLog(logs[file],'it was finished\n')
                del  popens[file]
                Nu=Nu-1
                if Nu==0:break
            elif popens[file].poll() ==None:
                continue
            elif popens[file].poll() >0:
                sys.stderr.write(colored("ERROR: "+ commands[file] +" \n",'red'))
                Flag=1
                Nu=Nu-1
                if Nu==0:break
    if Flag:sys.exit(3)



def qsubRun(logs,commands,cpu=1,check_time=3600):
    tmplogs,jobID,tmpcommands={},{},{}
    for file in logs:
        tmplogs[file]=logs[file]
        sys.stdout.write(colored('command: {}\n'.format(commands[file]),'green'))
        tmpcommand=""" {python3} {job_sub_py_2}  --cpu {cpu} --work  " {command}  &&  echo it was finished >  {tmplog}" """.format(python3=configDict['python3'], job_sub_py_2=configDict['job_sub_py_2'], cpu=cpu,command=commands[file], tmplog=logs[file])
        tmpcommands[file]=tmpcommand
        jobID[file]=os.popen(tmpcommand).read().strip().split(".")[0]   #'836836.admin1\n'
        sys.stdout.write(colored("INOFR: file {} jobID is {}\n".format(file,jobID[file]),"green"))
        time.sleep(2)
    start=time.time()
    while 1:
        length=len(tmplogs)
        if  length== 0 :break
        time.sleep(1200) 
        if time.time() - start < check_time: 
            jobs_state=os.popen("qstat").read().strip().split("\n")
            for file in tmplogs:
                if not op.exists(tmplogs[file]):
                    if os.popen("qstat|grep {ID}".format(ID=jobID[file])).read().strip().split("\n")[0]=="":
                        sys.stdout.write(colored('command: {}\n'.format(commands[file]),'green'))
                        jobID[file]=os.popen(tmpcommands[file]).read().strip().split(".")[0]
                        sys.stdout.write(colored("INOFR: Because of some err, the file {} jobID is {}\n".format(file,jobID[file]),"green"))
                        time.sleep(2)
        else:
            sys.stderr.write(colored('ERROR: the followed tasks are not finished in time.\n','red'))
            for file in tmplogs:
                if not op.exists(tmplogs[file]):
                    sys.stderr.write(colored('ERROR: {}\n'.format(commands[file]),'red'))
            sys.exit(1)
        tmplogs_1={}
        for file in tmplogs:
            if not op.exists(tmplogs[file]):
                tmplogs_1[file]=tmplogs[file]
        tmplogs=tmplogs_1


def subStr(s1, s2):   
    m=[[0 for i in range(len(s2)+1)]  for j in range(len(s1)+1)]
    mmax=0
    p=0
    for i in range(len(s1)):  
        for j in range(len(s2)):  
            if s1[i]==s2[j]:  
                m[i+1][j+1]=m[i][j]+1  
                if m[i+1][j+1]>mmax:  
                    mmax=m[i+1][j+1]  
                    p=i+1  

    ss=''
    if p!=0:
        ss=s1[p-mmax:p]
        while 1:
            if ss[-1]  in  ("-","_","."):
                ss=ss[:-1]
            else:
                break
    return ss


def fastQC(QcDir,FqList):
    checkDir(QcDir)
    popen={}
    tmppopen=[]
    command={}
    tmpFqList=[fq for fq in FqList  if not op.exists(outdir+'/log/'+op.basename(fq)+".fastQC.log") ]
    FqList=tmpFqList
    if len(FqList) != 0:
        for  fq in FqList:
            command[fq]='{fastqc} -o {QcDir} -t 10  {fq} '.format(fastqc=configDict['fastqc'],QcDir=QcDir,fq=fq)
            popen[fq]=subprocess.Popen(command[fq], shell=True)
            tmppopen.append(fq)
        err=[]
        time.sleep(120)
        while len(popen) > 0:
            time.sleep(60)
            for fq in tmppopen:
                if fq not in popen:continue
                if popen[fq].poll() ==0:
                    if popen[fq].returncode != 0:
                        err.append(fq)
                    else:
                        writeLog(outdir+'/log/'+op.basename(fq)+".fastQC.log",'it was finished\n')
                    del popen[fq]
                else:
                    popen[fq].wait()
        if len(err) >0:
            sys.stderr.write(colored("ERROR: "+ "   ".join([command[fq] for fq in err]) +" \n",'red'))
            sys.exit(1)

def rmdupFq(fq1,fq2,cleandir):
    #2017-4-11
    #add fastqUniq 
    popen={}
    command={}
    ##log
    tmplogs=[]
    tmpfq1=[]
    tmpfq2=[]
    for a,b in zip(fq1,fq2):
        subfq=subStr(a,b)
        if  not op.exists(outdir+'/log/'+op.basename(subfq)+".rmdupFq.log"):
            tmpfq1.append(a)
            tmpfq2.append(b)
    fq1,fq2=tmpfq1,tmpfq2
    
    for a,b in zip(fq1,fq2):
        subfq=op.basename(subStr(a,b))
        tmplog=outdir+'/log/'+subfq+".rmdupFq.log"
        if  op.exists(tmplog):continue
        sys.stdout.write(colored("INFOR: start rmdup fastq:\n",'green'))
        sys.stdout.write(colored(a + "\t" + b +"\n",'blue'))
        if '.gz' in a : 
            systemRun('gzip -cd {} > {}/rawdata/{}'.format(a,outdir,op.basename(a[:-3])))
            #a=a[:-3]
            a='{}/rawdata/{}'.format(outdir,op.basename(a[:-3]))
        if '.gz' in b : 
            systemRun('gzip -cd {} > {}/rawdata/{}'.format(b,outdir,op.basename(b[:-3])) )
            #b=b[:-3]
            b='{}/rawdata/{}'.format(outdir,op.basename(b[:-3]))

        A2=open(cleandir+'/{}.txt'.format(subfq),'w')
        A2.write(a+"\n"+b+"\n")
        A2.close()
        fq1rmdup=cleandir+"/"+subfq +".1.rmdup.fq"
        fq2rmdup=cleandir+"/"+subfq +".2.rmdup.fq"
        if qsub:
            tmplogs.append(tmplog)
            command='{fastuniq}  -i {input_list} -t q -o {fq1rmdup} -p {fq2rmdup} -c 0'.format(fastuniq=configDict['fastuniq'], input_list=cleandir+"/{}.txt".format(subfq), fq1rmdup=fq1rmdup, fq2rmdup=fq2rmdup)
            command=" {python3} {job_sub_py_2}  --work  ' {command}  &&  echo it was finished >  {tmplog}' ".format(python3=configDict['python3'], job_sub_py_2=configDict['job_sub_py_2'], command=command, tmplog=tmplog)
            systemRun(command)
        else:
            command[subfq]='{fastuniq}  -i {input_list} -t q -o {fq1rmdup} -p {fq2rmdup} -c 0'.format(fastuniq=configDict['fastuniq'], input_list=cleandir+"/{}.txt".format(subfq), fq1rmdup=fq1rmdup, fq2rmdup=fq2rmdup)
            popen[subfq]=subprocess.Popen(command[subfq], shell=True)
        #systemRun(command)
    if qsub:
        start=time.time()
        while 1:
            length=len(tmplogs)
            if  length== 0 :break
            time.sleep(60)
            if time.time() - start > 3600*4: 
                sys.stderr.write(colored('ERROR: the followed files which should be rmdup by fastqUniq are not finished in time.\n','red'))
                for log in tmplogs:
                    if not op.exists(log):
                        sys.stderr.write(colored('ERROR: {}\n'.format(op.basename(log)),'red'))
                sys.exit(1)
            tmplogs_1=[]
            for i in range(length):
                if not op.exists(tmplogs[i]):
                    tmplogs_1.append(tmplogs[i])
            tmplogs=tmplogs_1
    else:
        if len(fq1) != 0:
            for p in popen:
                popen[p].wait()
            err=[]
            for p in popen:
                if popen[p].returncode != 0:
                    err.append(p)
                else:
                    writeLog(outdir+'/log/'+p+".rmdupFq.log",'it was finished\n')
            if len(err) >0:
                sys.stderr.write(colored("ERROR: "+ "   ".join([command[p] +"\n"  for p in err]) +" \n",'red'))
                sys.exit(1)
    


#step 1 clean fastq()
def cleanFQ(fq1,fq2,cleanDir):
    global fastquniq,CFq1s,CFq2s,cleanquality
    CFq1s=[]
    CFq2s=[]
    if fastquniq: 
        rmdupFq(fq1,fq2,cleanDir) 
        tmp1,tmp2=[],[]
        for a,b in zip(fq1,fq2):
            subfq=op.basename(subStr(a,b))
            tmp1.append(cleanDir+"/"+subfq +".1.rmdup.fq")
            tmp2.append(cleanDir+"/"+subfq +".2.rmdup.fq")
        fq1,fq2=tmp1,tmp2
    popen={}
    command={}
    ##log
    tmpfq1=[]
    tmpfq2=[]
    #if 
    for a,b in zip(fq1,fq2):
        subfq=subStr(a,b)
        fq1Clean=cleanDir+"/"+op.basename(subStr(a,b))+".1.clean.fq"
        fq2Clean=cleanDir+"/"+op.basename(subStr(a,b))+".2.clean.fq"
        CFq1s.append(fq1Clean)
        CFq2s.append(fq2Clean)
        if  not op.exists(outdir+'/log/'+op.basename(subfq)+".clean.log"):
            tmpfq1.append(a)
            tmpfq2.append(b)
    fq1,fq2=tmpfq1,tmpfq2
    
    for a,b in zip(fq1,fq2):
        subfq=op.basename(subStr(a,b))
        sys.stdout.write(colored("INFOR: start clean fastq:\n",'green'))
        sys.stdout.write(colored(a + "\t" + b +"\n",'blue'))
        fq1Clean=cleanDir+"/"+subfq +".1.clean.fq"
        fq2Clean=cleanDir+"/"+subfq +".2.clean.fq"
        command[subfq]='{cutadapt} -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT   -G CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   -a TGCCAGGTGGTAAGTGAAGTTATTTGGTGT  -A ACACCAAATAACTTCACTTACCACCTGGCA  -a GATCGGAAGAGCACACGTCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGACGTGTGCTCTTCCGATC   -q {cleanquality},{cleanquality} --minimum-length 50   -o {fq1c} -p {fq2c}  {Fq1}   {Fq2}'.format(cutadapt=configDict['cutadapt'], cleanquality=cleanquality,fq1c=fq1Clean, fq2c=fq2Clean, Fq1=a ,Fq2=b )  
        print(command[subfq])

        popen[subfq]=subprocess.Popen(command[subfq], shell=True)
        #systemRun(command)
    if len(fq1) != 0:
        for p in popen:
            popen[p].wait()
        err=[]
        for p in popen:
            if popen[p].returncode != 0:
                err.append(p)
            else:
                writeLog(outdir+'/log/'+p+".clean.log",'it was finished\n')
        if len(err) >0:
            sys.stderr.write(colored("ERROR: "+ "   ".join([command[p] for p in err]) +" \n",'red'))
            sys.exit(1)


def filterFq(Fq1,Fq2):
    import filter_fq
    subfq=op.basename(subStr(Fq1,Fq2))
    if  not op.exists(outdir+'/log/'+subfq+".filterFq.log"): 
        if not remove_tmp:             
            signal =filter_fq.start(40,Fq1,Fq2,False)
        else:
            signal =filter_fq.start(40,Fq1,Fq2,True)
        if not signal:
            writeLog(outdir+'/log/'+subfq+".filterFq.log",'it was finished\n')
        #command='{python3} {filter_fq}  {Fq1} {Fq2} -t 40 '.format(python3=configDict['python3'], filter_fq=configDict['filter_fq'], Fq1=Fq1, Fq2=Fq2)
        #systemRun(command)


def bsmap(Fq1,Fq2,AlignDir):
    tmplogs,commands={},{}
    for a,b in zip(Fq1,Fq2):
        subfq=subStr(a,b)
        #print('subfq',subfq)
        subfq=os.path.basename(subfq)
        cleanq1=outdir+"/clean/"+subfq +".1.clean.fq"
        cleanq2=outdir+"/clean/"+subfq +".2.clean.fq"
        subfq=op.basename(subStr(cleanq1,cleanq2))
        bam='{}/{}.bsmap.bam'.format(AlignDir,subfq)
        filterFq(cleanq1,cleanq2)
        if  not op.exists(outdir+'/log/'+subfq+".bsmap.log"): 
            command= """{bsmap}  -a  {FqClean1}  -b  {FqClean2}  -d  {Ref}    -m 10 -x 500  -S 1 -z 33  -s 16  -g 3 -n 1 -q 0 -f 5 -p 20 -u -r 1 -v 0.08  |awk ' $6 !~/-/ && $6 !~ /[DI]0M$/'| {samtools} view -b -S - -o  {Bam}""".format(bsmap=configDict['bsmap'], FqClean1=cleanq1+".passed", FqClean2=cleanq2+".passed", Ref=ref, Bam=bam, samtools=configDict['samtools'])
            tmplog=outdir+'/log/'+subfq+".bsmap.log"
            tmplogs[bam]=tmplog
            commands[bam]=command
    if len(commands) >=1:
        systemRun2(tmplogs,commands,2)
    

def pickUpSam(bams):
    #2017--6-15
    tmplogs={}
    commands={}
    for bam in bams:
        log=outdir+"/log/"+op.basename(bam)[:-3]+"stdmap.log"
        if not op.exists(log):
            command='''{samtools} view -h -F  12 {bam}  | awk -f {stdmap_awk}  |{samtools} view -bS - > {stdBAM}'''.format(samtools=configDict['samtools'],stdmap_awk=configDict['stdmap_awk'],bam=bam,stdBAM=bam[:-3]+'std.map.bam')
            tmplogs[bam+'std']=log
            commands[bam+'std']=command
        # log=outdir+"/log/"+op.basename(bam)[:-3]+"diffmap.log"
        # if not op.exists(log):
            # command='''{samtools} view -F 12 {bam}  | awk -f {diffmap_awk} >{diffmap}'''.format(samtools=configDict['samtools'] ,diffmap_awk=configDict['diffmap_awk'],bam=bam,diffmap=bam[:-3]+"diffReadsName" )
            # tmplogs[bam+'err']=log
            # commands[bam+'err']=command
        # log=outdir+"/log/"+op.basename(bam)[:-3]+"unmap.log"
        # if not op.exists(log):
            # command='''{samtools} view  -f 4  {bam}  | cut -f 1 >{unmap1} && samtools view  -f 8 {bam}  | cut -f 1 >{unmap2} && cat {unmap1} {unmap2} >{unmap}'''.format(samtools=configDict['samtools'],bam=bam,unmap1=bam[:-3]+'UnMap1',unmap2=bam[:-3]+"UnMap2",unmap=bam[:-3]+"unmapReadName")
            # tmplogs[bam+"unmap"]=log
            # commands[bam+"unmap"]=command
    if qsub:
        qsubRun(tmplogs,commands,3,3600*4)
    else:
        if len(commands) >= 1:
            systemRun2(tmplogs,commands,4)

def bisSNP(inbams):
    '''
    Bis-SNP can only work under Sun JDK 1.6 currently.
    3. Perl: Bis-SNP perl scripts require Perl v 5.8.8 or later.
    '''
    #Add read group tag to BAM file

    tmplogs={}
    commands={}
    for inbam in inbams:
        withRGBam=inbam[:-3]+"withRG.bam"
        tmplog=outdir+'/log/'+op.basename(withRGBam)[:-3]+"log"
        if  not op.exists(tmplog): 
            command='''{java} -Xmx4g -jar {picard} AddOrReplaceReadGroups  I={inbam} O={withRGBam} ID=ceshi LB=lib1 PL=illumina PU=unit1 SM=S103 VALIDATION_STRINGENCY=SILENT'''.format(java=configDict['java'],picard=configDict['picard'],inbam=inbam,withRGBam=withRGBam)
            #samtools index S103.withRG.bam
            tmplogs[withRGBam]=tmplog
            commands[withRGBam]=command
    if len(commands)>0:
        if qsub:
            qsubRun(tmplogs,commands,5,3*3600)
        else:
            systemRun2(tmplogs,commands,8)
        
    tmplogs={}
    commands={}
    for inbam in inbams:
        withRGBam=inbam[:-3]+"withRG.bam"
        tmplog=outdir+'/log/'+op.basename(withRGBam)[:-3]+"index.log"
        if  not op.exists(tmplog):
            command='''{samtools} index {withRGBam}'''.format(withRGBam=withRGBam,samtools=configDict['samtools'])
            tmplogs[withRGBam]=tmplog
            commands[withRGBam]=command
    if len(commands) >0:
        systemRun2(tmplogs,commands,8)
    
    #Find indel region
    tmplogs={}
    commands={}
    for inbam in inbams:
        intervals=inbam[:-3]+'indel_target.intervals'
        withRGBam=inbam[:-3]+"withRG.bam"
        tmplog=outdir+'/log/'+op.basename(intervals)+'.log'
        if  not op.exists(tmplog): 
            command= """{java} -Xmx2g -jar {bissnp} -R {ref} -I {withRGBam} -T BisulfiteRealignerTargetCreator {known}  -o {intervals} -nt 10""".format(java=configDict['java1.6'],bissnp=configDict['bissnp'],ref=ref,known="-known\t" + "\t-known\t".join(knownsites),withRGBam=withRGBam,intervals=intervals)
            tmplogs[intervals]=tmplog
            commands[intervals]=command
    if len(commands)>0:
        if qsub:
            qsubRun(tmplogs,commands,5,10*3600)
        else:
            systemRun2(tmplogs,commands,4)

    tmplogs={}
    commands={}         
    for inbam in inbams:
        withRGBam=inbam[:-3]+"withRG.bam"
        intervals=inbam[:-3]+'indel_target.intervals'
        realignedBam=withRGBam[:-3]+"realigned.bam"
        tmplog=outdir+'/log/'+op.basename(realignedBam)[:-3]+'log'
        if  not op.exists(tmplog): 
            command="""{java} -Xmx10g -jar  {bissnp} -R {ref} -I {withRGBam} -T BisulfiteIndelRealigner -targetIntervals {intervals} {known} -cigar -o {realignedBam}""".format(java=configDict['java1.6'],bissnp=configDict['bissnp'],ref=ref,withRGBam=withRGBam,intervals=intervals ,known="-known\t" + "\t-known\t".join(knownsites),realignedBam=realignedBam)
            ##may add  --maxReadsInMemory 1000000 
            tmplogs[realignedBam]=tmplog
            commands[realignedBam]=command
    if len(commands)>0:
        if qsub:
            qsubRun(tmplogs,commands,5,24*3600)
        else:
            systemRun2(tmplogs,commands,4)
        
        
    #Base quality recalibration
    #Count Covariant
    tmplogs={}
    commands={}         
    for inbam in inbams:
        withRGBam=inbam[:-3]+"withRG.bam"   
        recalFile=withRGBam[:-3]+"recalFile.csv"
        realignedBam=withRGBam[:-3]+"realigned.bam"
        tmplog=outdir+'/log/'+op.basename(recalFile)[:-3]+'log'
        if  not op.exists(tmplog): 
            command="""{java} -Xmx2g -jar {bissnp} -R {ref} -I {realignedBam} -T BisulfiteCountCovariates {knownSites}  -cov ReadGroupCovariate -cov QualityScoreCovariate  -cov CycleCovariate  -recalFile {recalFile} -nt 10""".format(java=configDict['java1.6'],bissnp=configDict['bissnp'],realignedBam=realignedBam,ref=ref,knownSites="-knownSites\t" + "\t-knownSites\t".join(knownsites),recalFile=recalFile)
            tmplogs[recalFile]=tmplog
            commands[recalFile]=command
    if len(commands)>0:
        if qsub:
            qsubRun(tmplogs,commands,5,24*3600)
        else:
            systemRun2(tmplogs,commands,4)        

    tmplogs={}
    commands={}         
    for inbam in inbams:
        recalBam=inbam[:-3]+"realigned.recal.bam"
        #Write recalibrated base quality score into BAM file
        tmplog=outdir+'/log/'+op.basename(recalBam)[:-3]+'log'
        withRGBam=inbam[:-3]+"withRG.bam"   
        recalFile=withRGBam[:-3]+"recalFile.csv"
        realignedBam=withRGBam[:-3]+"realigned.bam"
        if  not op.exists(tmplog): 
            command="""{java} -Xmx10g -jar  {bissnp} -R  {ref}  -I  {realignedBam} -o {recalBam} -T BisulfiteTableRecalibration -recalFile {recalFile} -maxQ 40""".format(java=configDict['java1.6'],bissnp=configDict['bissnp'],known="-known\t" + "\t-known\t".join(knownsites),ref=ref,realignedBam=realignedBam,recalBam=recalBam,recalFile=recalFile)
            tmplogs[recalBam]=tmplog
            commands[recalBam]=command
    if len(commands)>0:
        if qsub:
            qsubRun(tmplogs,commands,5,24*3600)
        else:
            systemRun2(tmplogs,commands,4)   


    

def rmdupBam(Bams,rmdupBAMs):
    logs,commands={},{}
    for Bam,rmdupBAM in zip(Bams,rmdupBAMs):
        log=outdir+'/log/'+op.basename(rmdupBAM)+".PicardRmdupBAM.log"
        if  not op.exists(log):
            logs[Bam]=log
            commands[Bam]='{java} -jar {picard} MarkDuplicates I={sortedBAM} O={rmdupBAM} M={outdir}/marked_dup_metrics.txt  REMOVE_DUPLICATES=true'.format(java=configDict['java'], picard=configDict['picard'], sortedBAM=Bam, rmdupBAM=rmdupBAM, outdir=op.dirname(rmdupBAM))
    if qsub:
        qsubRun(logs,commands,5,3600*2)
    else:
        systemRun2(logs,commands,4)
    

def rmdupBam2(Bams,rmdupBAMs):
    logs,commands={},{}
    for Bam,rmdupBAM in zip(Bams,rmdupBAMs):
        log=outdir+'/log/'+op.basename(rmdupBAM)+".SamtoolsRmdupBAM.log"
        if  not op.exists(log):
            logs[Bam]=log
            commands[Bam]='{samtools} rmdup {Bam} {rmdupBAM}'.format(samtools=configDict['samtools'], Bam=Bam, rmdupBAM=rmdupBAM)
    if qsub:
        qsubRun(logs,commands,5,3600*2)
    else:
        systemRun2(logs,commands,4)


def cleanSam(Bams,cleanBams):
    logs,commands={},{}
    for Bam,cleanBam in zip(Bams,cleanBams):
        log=outdir+'/log/'+op.basename(Bam)+".cleanSam.log"
        if  not op.exists(log):
            logs[Bam]=log
            commands[Bam] = '{java} -jar {picard} CleanSam  I={Bam}  O={cleanBam} '.format(java=configDict['java'], picard=configDict['picard'], Bam=Bam, cleanBam=cleanBam)
    if qsub:
        qsubRun(logs,commands,5,3600*2)
    else:
        systemRun2(logs,commands,8)



def mergeBam(Bams,mergeBam):
    if  not op.exists(outdir+'/log/'+op.basename(mergeBam)+".mergeBam.log"): 
        command = '{samtools} merge -@ 30 {mergebam} {bamlist} '.format(samtools=configDict['samtools'], mergebam=mergeBam, bamlist=" ".join(Bams) )
        systemRun(command)
        writeLog(outdir+'/log/'+op.basename(mergeBam)+".mergeBam.log",'it was finished\n')

def sortedBam(Bam,sortedBAM):
    if not  op.exists(outdir+'/log/'+op.basename(sortedBAM)+".sortedBAM.log"): 
        command = ' {samtools} sort  -@ 30 {bam}  -o {sortedBam}'.format(samtools=configDict['samtools'], bam=Bam,sortedBam=sortedBAM)
        systemRun(command)
        writeLog(outdir+'/log/'+op.basename(sortedBAM)+".sortedBAM.log",'it was finished\n')

#step 3 devied the bam file to two types file()
def splitBam(Bam):
    '''
    Then devided to 23 or 24 files by chromsome
        1-22,X,Y
    '''
    if not  op.exists(outdir+'/log/'+op.basename(Bam)+".index.log"): 
        command="{samtools} index {BAM} ".format(samtools=configDict['samtools'], BAM=Bam)
        systemRun(command)
        writeLog(outdir+'/log/'+op.basename(Bam)+".index.log",'it was finished\n')
    popen={}
    command={}
    if  not op.exists(outdir+'/log/'+op.basename(Bam)+".splitBam.log"): 
        CHR= ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
        for chr in CHR:
            checkDir(outdir+'/chr_BAM/'+chr)
            checkDir(outdir+'/nomal/'+chr)
            command[chr]='{samtools} view -b  {BAM} {chr} > {path}/{chr}/{chr}.merge.sorted.rmdup.bam  && {samtools} index {path}/{chr}/{chr}.merge.sorted.rmdup.bam '.format(samtools=configDict['samtools'], BAM=Bam, chr=chr, path=outdir+"/chr_BAM" if not "/nomal/Bs.sorted.bam" in Bam else outdir+"/nomal")
            popen[chr]=subprocess.Popen(command[chr],shell=True)
        for chr in CHR:popen[chr].wait()
        for chr in CHR:
            if popen[chr].returncode != 0:  
                sys.stderr.write(colored("ERROR: "+ command[chr] +" \n",'red'))
                sys.exit(3)
        writeLog(outdir+'/log/'+op.basename(Bam)+".splitBam.log",'it was finished\n')

    #File=out_dir+'/'+sample+'/'+'chr_BAM/'+sample+'_split_BAM.log'
    #Con='step3  finished\n'
    #writeLog(File,Con)

def bamQC(Bam):
    if not  op.exists(outdir+'/log/'+op.basename(Bam)+".bamQC.log"): 
        command="{qualimap} bamqc -bam {Bam} --java-mem-size=10G ".format(qualimap=configDict['qualimap'], Bam=Bam)
        systemRun(command)
        writeLog(outdir+'/log/'+op.basename(Bam)+".bamQC.log",'it was finished\n')
    
def dealBam(Bams):


    pickUpSam(Bams)
    
    for sam in Bams:
        bam=sam[:-3]+"std.map.bam"
        sortedBam(bam,sam[:-3]+"sorted.bam")
        
    #cleanSam
    #picard rmdup
    tmpBams,tmpRmDupBams=[],[]
    for sam in Bams:
        tmpBams.append(sam[:-3]+"sorted.bam")
        tmpRmDupBams.append(sam[:-3]+"sorted.cleanSam.bam")
    cleanSam(tmpBams,tmpRmDupBams)


    #picard rmdup
    tmpBams,tmpRmDupBams=[],[]
    for sam in Bams:
        tmpBams.append(sam[:-3]+"sorted.cleanSam.bam")
        tmpRmDupBams.append(sam[:-3]+"sorted.Prmdup.bam")
    rmdupBam(tmpBams,tmpRmDupBams)
    
    #samtools rmdup
    tmpBams,tmpRmDupBams=[],[]
    for sam in Bams:
        tmpBams.append(sam[:-3]+"sorted.Prmdup.bam")
        tmpRmDupBams.append(sam[:-3]+"sorted.Prmdup.Srmdup.bam")
    rmdupBam2(tmpBams,tmpRmDupBams)
    
    if len(Bams) >1:
        bam=[sam[:-3]+"sorted.Prmdup.Srmdup.bam"  for  sam in Bams]
        subbam=subStr(bam[0],bam[1])
        if len(bam) >2:
            for b in bam:
                subbam=subStr(subbam,b)
            mergebam=[]
            #if subbam[-1]  in  ("-","_","."):subbam=subbam[:-1]
        if subbam=='':subbam=time.strftime("%Y_%m_%d_%H",time.localtime())
        mergeBam(bam,subbam+".sorted.rmdup.merge.bam")
        sortedBam(subbam+".sorted.rmdup.merge.bam",subbam+".sorted.rmdup.merge.sorted.bam")
        bamQC(subbam+".sorted.rmdup.merge.sorted.bam")
        # if bs2real:
            # bwa(subbam+".sorted.rmdup.merge.sorted.bam") #
            # splitBam(outdir+"/nomal/Bs.sorted.bam")
        if not bissnp:
            splitBam(subbam+".sorted.rmdup.merge.sorted.bam")
        else:

            #bisSNP(subbam+".sorted.rmdup.merge.sorted.bam")
            #splitBam(subbam+".sorted.rmdup.merge.sorted.realigned.recal.bam")
                        
            splitBam(subbam+".sorted.rmdup.merge.sorted.bam")
            chr_bams=['{path}/chr_BAM/{chr}/{chr}.merge.sorted.rmdup.bam'.format(path=outdir,chr=chr) for chr in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY'] ]
            bisSNP(chr_bams)
    else:
        bamQC(Bams[0][:-3]+"sorted.Prmdup.Srmdup.bam")
        # if bs2real:
            # bwa(Bams[0][:-3]+"sorted.Prmdup.Srmdup.bam")
            # splitBam(outdir+"/nomal/Bs.sorted.bam")
        if not bissnp:
            splitBam(Bams[0][:-3]+"sorted.Prmdup.Srmdup.bam")
        else:
            #splitBam([Sam[0][:-3]+"sorted.Prmdup.Srmdup.bam"])
            bisSNP([Bams[0][:-3]+"sorted.Prmdup.Srmdup.bam"])
            splitBam(Bams[0][:-3]+"sorted.Prmdup.Srmdup.realigned.recal.bam")

#step 4 get mutation file()
def getSnpAndMeth(outdir):

    global qsub
    CHR= ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']
    snpCommand=methyCommand=''
    tmplogs=[]

    # tmpChrBamPath=['chr_BAM','nomal'] if  bs2real else ['chr_BAM']
    tmpChrBamPath= ['chr_BAM']
    for Path in tmpChrBamPath:
        tmpCHR=[ chr for chr in CHR if not op.exists(outdir+'/log/'+chr+'.'+Path+".snpmeth.log")]
        #CHR=tmpCHR
        snpCommand  ='cat   ' if not snpCommand   else snpCommand   + '  &&  cat  '
        methyCommand='cat   ' if not methyCommand else methyCommand + '  &&  cat  '
        for chr in CHR:
            if Path!='nomal':
                if bissnp and len(fq1)>1:
                    bamFile='{}/{}/{}.merge.sorted.rmdup.realigned.recal.bam'.format(outdir+"/"+Path,chr,chr)
                else:
                    bamFile='{}/{}/{}.merge.sorted.rmdup.bam'.format(outdir+"/"+Path,chr,chr)
            else:
                bamFile='{}/{}/{}.merge.sorted.rmdup.bam'.format(outdir+"/"+Path,chr,chr)
            genomeFile='{refDir}/{chr}.fa'.format(refDir=configDict['refDir'], chr=chr)
            command='''{python3} {DSBS}  {bam} --ref {ref} --Chr {chr} --maxBp 50 --minVaf 0.15 -q -g {genomeFile} -o {outdir}/{Path}/ --cpu {cpu} -d {dbsnp} '''.format(python3=configDict['python3'], DSBS=configDict['DSBS'], bam=bamFile, ref='hg38' if 'hg38' in configDict['ref'] else 'hg19' ,genomeFile=genomeFile, chr=chr, outdir=outdir,Path=Path, cpu=cpu, dbsnp=configDict['dbsnp_gz'])
            tmplog=outdir+'/log/'+chr+'.'+Path+".snpmeth.log"
            tmplogs.append(tmplog)
            snpCommand   +=outdir +'/' +Path + '/' + op.splitext(os.path.basename(bamFile))[0]+'.Snp.txt   ' 
            methyCommand +=outdir +'/' +Path + '/' + op.splitext(os.path.basename(bamFile))[0]+'.Methylation.txt   ' 
            if chr in tmpCHR:
                #systemRun(command)
                #writeLog(tmplog,'it was finished\n')
                if not qsub:
                    systemRun(command)
                    writeLog(tmplog,'it was finished\n')
                else:
                    #if  'anzhen' in config:
                    command1=" {python3} {job_sub_py_2}  --work  '{command}  &&  echo it was finished >  {tmplog}' ".format(python3=configDict['python3'], job_sub_py_2=configDict['job_sub_py_2'],  command=command, tmplog =tmplog)
                    systemRun(command1)
                    time.sleep(10)
        snpCommand  +="> {}/{}/{}.Snp.vcf  ".format(outdir,Path,Path)
        methyCommand+='>{}/{}/{}.Methylation.txt  '.format(outdir,Path,Path)
    if qsub:
        start=time.time()
        while 1:
            length=len(tmplogs)
            if  length== 0 :break
            time.sleep(60)
            if time.time() - start > 3600*8: 
                sys.stderr('ERROR: the followed files called by DSBS are not finished in time.\n')
                for log in tmplogs:
                    if not op.exists(log):
                        sys.stderr.write('ERROR: {}\n'.format(op.basename(log)))
                sys.exit(1)
            tmplogs_1=[]
            for i in range(length):
                if not op.exists(tmplogs[i]):
                    tmplogs_1.append(tmplogs[i])
            tmplogs=tmplogs_1
    
    if snpCommand:
        systemRun(snpCommand)
        systemRun(methyCommand)
    for Path in tmpChrBamPath:
        systemRun("cp {outdir}/{Path}/{Path}.Snp.vcf  {outdir}/snpmeth/".format(outdir=outdir,Path=Path))
        systemRun("cp {outdir}/{Path}/{Path}.Methylation.txt  {outdir}/snpmeth/".format(outdir=outdir,Path=Path))
    if len(tmpChrBamPath) ==1:
        systemRun("ln -s {outdir}/snpmeth/{path}.Snp.vcf  {outdir}/snpmeth/SNP.vcf ".format(outdir=outdir,path=tmpChrBamPath[0]))
        systemRun("ln -s {outdir}/snpmeth/{path}.Methylation.txt  {outdir}/snpmeth/METHY.txt ".format(outdir=outdir,path=tmpChrBamPath[0]))
    else:
        #if  op.exists('{outdir}/{Path}/{Path}.Snp.vcf'.format(outdir=outdir,Path=Path)):return
        SNP={chr:defaultdict(lambda:{'con':'','dep':[0,0,0]}) for chr in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']}
        METHY={chr:defaultdict(lambda:{'con':'','dep':[0,0,0]}) for chr in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']}
        for Path in tmpChrBamPath:
            with open ('{outdir}/{Path}/{Path}.Snp.vcf'.format(outdir=outdir,Path=Path),'r')  as A1:
                for line in A1.readlines():
                    temp=line.strip().split()
                    SNP[temp[0]][temp[1]]['con']="\t".join(temp[2:6])
                    for i in range(3):SNP[temp[0]][temp[1]]['dep'][i]+=int(temp[i+6])
            with open ('{outdir}/{Path}/{Path}.Methylation.txt'.format(outdir=outdir,Path=Path),'r')  as A2:
                for line in A2.readlines():
                    temp=line.strip().split()
                    METHY[temp[0]][temp[1]]['con']="\t".join(temp[2:7])
                    for i in range(3):METHY[temp[0]][temp[1]]['dep'][i]+=int(temp[i+7])
        with open('{outdir}/snpmeth/SNP.vcf'.format(outdir=outdir),'w') as A1:
            for chr in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']:
                for k in sorted(SNP[chr].keys(),key=int):
                    A1.write('{}\t{}\t{}\t{}\n'.format(chr,k,SNP[chr][k]['con'],"\t".join(list(map(str,SNP[chr][k]['dep'])))))
        with open('{outdir}/snpmeth/METHY.txt'.format(outdir=outdir),'w') as A1:
            for chr in ['chr'+str(i) for i in range(1,23)] + ['chrX', 'chrY']:
                for k in sorted(METHY[chr].keys(),key=int):
                    A1.write('{}\t{}\t{}\t{}\n'.format(chr,k,METHY[chr][k]['con'],"\t".join(list(map(str,METHY[chr][k]['dep'])))))
                    




    

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config'     ,   default='DSBS.config' ,help='config file')
    parser.add_argument('-a','--fq1',      nargs='*',   required='true', help='Fastq1(s), Space separation')
    parser.add_argument('-b','--fq2',      nargs='*',   required='true', help='Fastq2(s), Space separation')
    parser.add_argument('-r','--ref',      type=str,    default='',      help='refrence')
    parser.add_argument('-k','--knownsites', nargs='*'  ,default=['1000G_omni','1000G_phase1','dbsnp','hapmap'],     help='knownsites')
    parser.add_argument('-d', '--depth',   type=int,    default=4,       help='the minimum depth , default is 4')
    parser.add_argument('-q','--cleanquality',  type=int,    default=30, help='quality when clening the fastq')
    parser.add_argument('--gapsize',       type=int,    default=3,       help='the size of insert or del, default is 3')
    parser.add_argument('--fastquniq',     action='store_true',default=False, help='rmdup fastq before cleanning fastq') 
    parser.add_argument('--remove_tmp', action='store_true',default='False',help='rm temp file')
    parser.add_argument('--bissnp',      action='store_true',default=False, help='bisSNP')
    parser.add_argument('-t','--cpu',      type=int,    default=20,      help='thread, default is 20')
    parser.add_argument('--pp',            default='11111',            help='1,execute ,0,skip. 1)qualimap, 2)clean, 3)align, 4)dealBam, 5)call mutation and methylation, default is 11111')
    parser.add_argument('-o',"--outdir",   required='true',              help='outdir')
    parser.add_argument('--qsub'    ,      action='store_true',          help='use qsub ')
    return parser
 
def main():
    purpose()
    parser = get_parser()
    args = parser.parse_args()
    global config,fq1,fq2,ref,adapter,cpu,pp,outdir,qsub,depth,fastquniq,configDict,bissnp,knownsites,gapsize,cleanquality,remove_tmp
    
    config      =  args.config
    fq1         =  args.fq1
    fq2         =  args.fq2
    ref         =  args.ref
    cpu         =  args.cpu
    pp          =  args.pp
    outdir      =  args.outdir
    qsub        =  args.qsub
    depth       =  args.depth  
    fastquniq   =  args.fastquniq
    bissnp      =  args.bissnp
    knownsites  =  args.knownsites
    gapsize     =  args.gapsize
    cleanquality=  args.cleanquality
    remove_tmp  =  args.remove_tmp
    

    configDict={}
    resDict={}
    softDict={}
    if not op.exists(config):
        sys.stderr.write(colored("ERROR: config {} is not existed!\n".format(config),'red'))
    else:
        with open(config,'r') as A1:
            for line in A1.readlines():
                #print(line.strip())
                if line[0] == "#" :   continue
                if line.strip() =="": continue
                temp=line.strip().split(":")
                configDict[temp[0]]=temp[1]
                if temp[2]=='soft':
                    softDict[temp[0]]=temp[1]
                elif temp[2]=='resource':
                    resDict[temp[0]]=temp[1]

    for kk in resDict:configDict[kk]=resDict[kk]
    sys.stdout.write(colored("StartTime: "+ time.strftime("%Y-%m-%d %H:%M:%S",time.localtime())+"\n",'green'))
    sys.stdout.write(colored(sys.argv[0]+"   :version 1.0 "+"\n",'green'))
    #purpose()
    if len (fq1) != len(fq2): sys.stderr.write(colored("ERROR: Fq should be paired! \n\t{}\n\t{}\n".format(" ".join(fq1)," ".join(fq2)),'red'))

    #checkFile
    sys.stdout.write(colored("fq:\n",'green'))
    for a,b in zip(fq1,fq2):
        if checkFile(a) and checkFile(b):
            sys.stdout.write(colored("\t{}\n\t{}\n".format(a,b),'blue'))
    
    ref=ref or resDict['ref']
    if checkFile(ref):
        sys.stdout.write(colored("Reference:\n\t{}\n".format(ref),'blue'))
        configDict['ref']=ref
        configDict['refDir']=op.dirname(ref)

    
    print('resDict',resDict)
    tmpKnownsite=[]
    for ks in knownsites:
        
        if ks in resDict:
            if checkFile(resDict[ks]):
                configDict[ks]=resDict[ks]
                tmpKnownsite.append(resDict[ks])
                sys.stdout.write(colored("knownSite:\n\t{}\n".format(resDict[ks]),'blue'))
        else:
            if checkFile(ks):
                sys.stdout.write(colored("knownSite:\n\t{}\n".format(ks),'blue'))
            tmpKnownsite.append(ks)
    knownsites=tmpKnownsite

    checkDir(outdir)
    for d in ['rawdata','fastqc','clean','align','chr_BAM','snpmeth','log']:checkDir(outdir+"/"+d)

    #step
    step=[fastQC,cleanFQ,bsmap,dealBam,getSnpAndMeth]
    for i,j in enumerate(str(pp)):
        if j == '1':
            if i==0:
                step[i](outdir+"/fastqc",fq1+fq2)
            elif i==1:
                step[i](fq1,fq2,outdir+"/clean")
            elif i==2:
                step[i](fq1,fq2,outdir+"/align/")
            elif i==3:
                bams=[]
                for a,b in zip(fq1,fq2):
                    subfq=subStr(a,b)
                    subfq=os.path.basename(subfq)
                    bams.append(outdir+"/align/"+subfq+".bsmap.bam")
                step[i](bams)
            elif i==4:
                step[i](outdir)

if __name__ == '__main__':
    main()
