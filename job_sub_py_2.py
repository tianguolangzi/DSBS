#!/public/home/jcli/public/bin/python3
# -*- coding: utf-8 -*-
#author: Zhang Kun
#vesion: 2.0
#date:2016-8-20
#email:tianguolangzi@yahoo.com
#too yang too naive
import sys,getopt,os,gzip,time,re
def usage():
    print('''\033[01;33m
        usage: python %s -option  <argement> 
                    --node      deaufalt is 1
                                    1: which is chosed by the program
                                    node_name:the one of all nodes
                    --cpu       deaufalt is 1
                    --jobname   deaufalt is /public/home/jcli/temp/jobs/%s/%s
                    --work      is necessary
                    --mem       deaufalt is none
                    --config    deaufalt is /public/home/jcli/zhangkun/bin/python3_bin/node_config
                    --auto      deaufalt is Y
                                    Y :node which is chosed by the PBS
                                    N :node which is chosed by the program
                    -h/--help'''%(sys.argv[0],time.strftime('%Y-%m-%d-%H'),time.strftime('%H-%M-%S')))
    print('''\033[01;32m
eg:
   python-->  job_sub_py --node fat6 --cpu 40 --jobname /public/home/jcli/temp/jobs/test --work "python3 /public/home/jcli/temp/python3/timer.py"
   shell--->  job_sub_py --node fat6 --cpu 40 --jobname /public/home/jcli/temp/jobs/test --work "for i in {1..10}; do echo \$i ;done"
   python-->  job_sub_py --cpu 40 --auto N ---mem 40G --work "python3 /public/home/jcli/temp/python3/timer.py"
        ''')

def option():
    global node,cpu,jobname,work,nodes,states,allmems,tasks,mem,config,auto
    if len(sys.argv) < 2:
        usage()
        sys.exit(2)
    try :
        opts ,args = getopt.getopt(sys.argv[1:],"h",["help","node=","cpu=","work=","jobname=","mem=","config=","auto="])
        #opts is (('-i','infile'),('-o','outflie'),('-c','no'),('-h',''))
    except getopt.GetoptError:
        print("\033[01;31m"+ "get option error!")
        usage()
        sys.exit(2)
    for opt ,val in opts:
        if opt in ("-h","--help"):
            usage()
            sys.exit(1)
        else:
            if opt in ('--node',):
                node = val
            if opt in ("--cpu",):
                cpu = val
            if opt in ("--work",):
                work = val
            if opt in ("--jobname",):
                jobname = val
            if opt in ("--mem",):
                mem = val
            if opt in ("--config",):
                config = val
            if opt in ("--auto",):
                auto = val
    if "work" not in globals():
        usage()
        sys.exit(2)

    if "jobname" not in globals():
        path = "/public/home/jcli/temp/jobs/%s"%time.strftime('%Y-%m-%d-%H')
        if not os.path.exists(path):
            os.makedirs(path)
        jobname = path + "/" + time.strftime('%H-%M-%S')
        ID=0
        while 1:
            if not os.path.exists(jobname) :
                break
            ID+=1
            jobname=jobname+"_"+str(ID)
    else:
        if not os.path.dirname(jobname):
            jobname = os.getcwd() + "/" + jobname

    if "cpu" not in globals():
        cpu = "1"
    else:
        if not cpu.isdigit():
            cpu = re.match("[0-9]+",cpu)
            if not cpu:
                cpu = "1"
            else:
                cpu = cpu.group()
        if int(cpu) >48:
            cpu = "48"

    if "config" not in globals():
        if os.path.exists("/public/home/zhangkun/bin/python3_bin/node_config"):
            config = "/public/home/zhangkun/bin/python3_bin/node_config"
    node_infor = get_node_infor()
    if "auto" not in globals():
        auto = "Y"
    if "node" not in globals():
        node = "1"
    if node != "1":
        if node not in node_infor:
            print("\033[01;31m"+"%s is not existed"%node)
            node = get_node(node_infor,0)
            print ("\033[01;31m"+"%s is chosed for %s"%(node,jobname))
        else:
            node=get_node(node_infor,[node])
            if node == "1": node = get_node(node_infor,0)
    else:
        if auto in ("N","n","NO","no","No","nO") :
            node = get_node(node_infor,0)
        else:
            node = "1"


def get_node(node_infor,nodes):
    if not nodes:
        if "mem" in globals() and int(mem[:-1]) >= 50:
            nodes=[key for key in node_infor.keys() if node_infor[key]["priority"]=="high"] +[key for key in node_infor.keys() if node_infor[key]["priority"]=="low"]
        else:
            nodes=[key for key in node_infor.keys() if node_infor[key]["priority"]=="low"] +[key for key in node_infor.keys() if node_infor[key]["priority"]=="high"]

    for temp_node in nodes:
        print(node_infor[temp_node])
        if node_infor[temp_node]["able_state"] == "disable":
            continue
        if node_infor[temp_node]["state"] != "free":
            continue
        if int(node_infor[temp_node]["ncpus"]) - int(node_infor[temp_node]["tasks"]) < int(cpu):
            continue
        if "mem" in globals():
            if abs((int(node_infor[temp_node]["allmem"])-int(node_infor[temp_node]["resi"]))/1024 - int(mem[:-1])) <= 10:
                return temp_node
        else:
            return temp_node
    else:
        return "1"

def get_node_infor():
    #nodes=os.popen("pestat |grep -v state |cut -c 5-15 ").read().split("\n")
    #nodes=[n.strip() for n in nodes if n!=""]
    #states=os.popen("pestat |grep -v state |cut -c 18-22 ").read().split("\n")
    #states=[s.strip() for s in states if s!=""]
    #node state phymem ncpus allmem resi tasks
    information = os.popen("pestat |grep -v state |awk '{print $1,$2,$4,$5,$6,$7,$9}'").read().strip().split("\n")
    node_infor = {}
    for infor in information:
        temp_infor = infor.split()
        if temp_infor[0] not in node_infor: 
            node_infor[temp_infor[0]] = {"state":"","phymem":"","ncpus":"","allmem":"","resi":"","tasks":"","max_cpu":"0","max_mem":"0","able_state":"able","priority":"low"}
        for i,j in zip(["state","phymem","ncpus","allmem","resi","tasks"],temp_infor[1:]): 
            node_infor[temp_infor[0]][i] = j
    #print(config)
    if config:
        with open(config,"r") as A1:
            #fat6   48  516 able/disable    high/low
            temp = [line.strip().split() for line in A1.readlines()[1:]]
            for infor  in temp:
                if infor[0] in node_infor: 
                    for i,j in zip(["max_cpu","max_mem","able_state","priority"],infor[1:]): node_infor[infor[0]][i]=j
    return node_infor

def creat_job_sub_sh():
    with open(jobname + ".sh","w") as A1:
        if "mem" in globals():
            A1.write("""#!/bin/bash
#PBS -N %s
#PBS -o %s.log
#PBS -e %s.err
#PBS -q new
#PBS -l nodes=%s:ppn=%s,mem=%s
%s\n"""%(os.path.basename(jobname),jobname,jobname,node,cpu,mem,work))
        else:
            A1.write("""#!/bin/bash
#PBS -N %s
#PBS -o %s.log
#PBS -e %s.err
#PBS -q new
#PBS -l nodes=%s:ppn=%s
%s\n"""%(os.path.basename(jobname),jobname,jobname,node,cpu,work))



def INT():
    option()
    creat_job_sub_sh()
    os.system("qsub %s.sh"%jobname) 
INT()

