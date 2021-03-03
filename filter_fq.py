#-*-coding=utf-8 -*-#
import os,sys,gzip
import time
import argparse
from collections import defaultdict
from multiprocessing import Manager,Pool
__version__=1.0


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action = 'version', version = 'filter_fq.py {}'.format(__version__))
    parser.add_argument('fq1',  help = '' )
    parser.add_argument('fq2',  help = '' )
    parser.add_argument('-t', '--cpu', type = int, default=20, help = '')
    return parser

def readfile(File):
    A1=gzip.open(File) if File[-3:] == '.gz' else open(File)
    lines=A1.readlines()
    A1.close()
    return lines

def writefile(j,fq1Name,fq2Name):

    F1=readfile(fq1Name)
    F2=readfile(fq2Name)
    A1=open("{}.passed".format(fq1Name),"w") 
    A2=open("{}.passed".format(fq2Name),"w") 
    A3=open("{}.nopassed".format(fq1Name),"w")
    A4=open("{}.nopassed".format(fq2Name),"w")
    print("start",j,len(F1),len(F2))
    Con1=[]
    Con2=[]
    i=0
    while i<len(F1)//4:
        if F1[i*4+1] == "\n" or F2[i*4+1] == "\n":
            A3.writelines(F1[i*4:i*4+4])
            A4.writelines(F2[i*4:i*4+4])
            i+=1
            continue
        Con1.extend(F1[i*4:i*4+4])
        Con2.extend(F2[i*4:i*4+4])
        i+=1
    print("end",j,len(Con1),len(Con2))
    A1.writelines(Con1)
    A2.writelines(Con2)
    A1.close()
    A2.close()
    A3.close()
    A4.close()


    
def start(cpu,fq1,fq2,AA):

    pool=Pool(processes=cpu)
    #length=len(fq1Read)
    #Num=int(length//4000000)+1
    print(cpu,fq1,fq2)
    C1='cat '
    C2='cat '
    N1='cat '
    N2='cat '
    Rm="rm  "
    j=0
    if not os.system("split -l 2000000 {} -d -a 3 {}.tmp && split -l 2000000 {} -d -a 3 {}.tmp ".format(fq1,fq1,fq2,fq2)):
        path=os.path.dirname(fq1)
        if path=='':path=os.getcwd()
        #print(os.listdir(path))
        for i in os.listdir(path):
            if os.path.basename(fq1) + ".tmp" in i:
                j+=1
        for n in range(j):
            nu=str(n)
            nu="0"*(3-len(nu))+nu
            fq1Name=fq1+".tmp" +nu
            fq2Name=fq2+".tmp" +nu
            #print(fq1Name,fq2Name)
            C1+= fq1Name+".passed" +"  " 
            C2+= fq2Name+".passed" +"  " 
            N1+= fq1Name+".nopassed" +"  " 
            N2+= fq2Name+".nopassed" +"  " 
            Rm+=fq1Name + " " + fq1Name+".passed" +"  "  +fq2Name+ " "+ fq2Name+".passed" +"  " + fq1Name+".nopassed" +"  " +fq2Name+".nopassed" +"  " 

            pool.apply_async(writefile, (n,fq1Name,fq2Name))
    pool.close()
    pool.join()
    command1 =  C2 +" > " +fq2 +".passed"+ " && " +N1+" > "+fq1 +".nopassed"+" && " +N2 +" > "+ fq2+".nopassed"+ " && " + C1 +" > "+fq1 +".passed"  
    print(command1)
    if not os.system(command1):
        command2=Rm
        if not os.system(command2):
            print("merger is OVER")
            return 0 

def work():
    parser = get_parser()
    args = parser.parse_args()
    cpu=args.cpu
    fq1=args.fq1
    fq2=args.fq2
    start(cpu,fq1,fq2)
if __name__=='__main__':
    work()



