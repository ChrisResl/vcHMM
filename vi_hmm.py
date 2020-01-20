import argparse
from Bio import SeqIO

import numpy as np
import pandas as pd
import itertools as it
import math

from matplotlib import pyplot as plt
#import numpy as np
#from numpy import linalg as la
from scipy import stats
#import math
import warnings
warnings.filterwarnings("ignore")



def draw_lines(val_param, name_param, num_softwares):
    colors = ["b-", "y-", "r-", "k-", "p-"]
    markers=['o','v',"s","p", "h"]
    softwares=["Freebayes", "Samtools", "Varscan", "vcHMM", "viHMM"]
    fig = plt.figure(figsize=(10, 10))
    for i in range(num_softwares):
        covs = [7, 14, 21]

        plt.plot(covs, val_param[3 *i: 3*i+3], colors[i], label=softwares[i], marker=markers[i])
    plt.title("Comparison of "+name_param)
    plt.xlabel("Coverages")
    plt.ylabel(name_param)
    plt.legend()
    #plt.show()
    plt.savefig("/Users/chin-hua-huang/Desktop/"+name_param)

def read_stats_3(stats_file):
    f = open(stats_file, 'r')
    collect = []
    names=[]
    for x in f:
        spl = x.split("  ")
        # print(spl)
        tmp = [float(spl[1]), float(spl[2]), float(spl[3])]
        names.append(spl[0])
        collect.append(tmp)

    return collect, names


def avg_perform(arr):
    sen=[]
    prec=[]
    f1=[]

    for i in range(len(arr)):
        sen.append(arr[i][0])
        prec.append(arr[i][1])
        f1.append(arr[i][2])

    sen=np.mean(np.asarray(sen))
    prec=np.mean(np.asarray(prec))
    f1=np.mean(np.asarray(f1))
    avg=[sen, prec, f1]
    return avg


def read_stats(stats_file):
    f = open(stats_file, 'r')
    collect=[]

    for x in f:
        spl=x.split("  ")

        tmp=[float(spl[1]), float(spl[2]), float(spl[3])]

        collect.append(tmp)



    return collect


def read_file(file_path):
    ''''''

    collect = []
    f=open(file_path, 'r')

    for x in f:
        if x[0:6]=='simref':
            tmp=x.split("\t") #tmp: the splitted line
            if len(tmp) == 4:
                tmp2= [tmp[0], tmp[1], tmp[2], tmp[3][0:-1]] #tmp2: the splitted line without \n
                collect.append(tmp2)
                print("-+-+-")
                print(tmp2)
            if len(tmp) == 5: #If the file has 5 columns
                if tmp[4][0:-1] == "+" or tmp[4][0:-1] == "-":
                    tmp3=[tmp[0], tmp[1], tmp[2], tmp[3]]
                    print("-truth-")
                else:

                    tmp3= [tmp[0], tmp[1], tmp[3], tmp[4][0:-1]]

                if tmp3[2] == "-":
                    tmp3[2] = ""
                if tmp3[3] == "-":
                    tmp3[3] = ""
                print("---")
                print(tmp3)
                collect.append(tmp3)
    return collect


def eval_result(gt, vcf):
    ''''''
    '''gt: information of ground truth'''
    '''vcf: information of vcf file'''
    tp = 0.  # True Positive: If the position in the ground truth is also predicted in the vcf
    fp = 0.  # False Positive: If the position is present in the vcf file but not in ground truth
    fn = 0.  # False Negative: If the position is not present in the vcf file but in ground truth

    '''Number of items that have the matching locations'''
    v_locs = []
    g_locs = []
    v_locs_snp = []
    v_locs_indel = []
    g_locs_snp = []
    g_locs_indel=[]
    for v in vcf:
        v_locs.append(v)
        if len(v[2]) == len(v[3]):
            v_locs_snp.append(int(v[1]))
        else:
            v_locs_indel.append(int(v[1]))

    for g in gt:
        g_locs.append(g)
        if len(g[2]) == len(g[3]):
            g_locs_snp.append(int(g[1]))
        else:
            g_locs_indel.append(int(g[1]))

    for vs in v_locs_snp:
        if vs in g_locs_snp:
            tp += 1.0
        else:
            fp += 1.0

    for vi in v_locs_indel:
        if vi in g_locs_indel:
            tp += 1.0
        else:
            fp += 1.0

    for g in g_locs:
        if g not in v_locs:
            fn += 1.0

    #print(tp)
    #print(fp)
    #print(fn)

    '''sen: Sensitivity'''
    '''prec: Precision'''
    '''f1: F1 Score'''
    sen = tp / (tp + fn)
    prec = tp / (tp + fp)
    f1 = 2. * tp / (2. * tp + fp + fn)

    return sen, prec, f1




def write_file(results, output_path, software_name):
    file=open(output_path, 'a')

    file.write("{}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(software_name, results[0], results[1], results[2]))
    file.close()
    print("Complete")



def parsing():
    parser = argparse.ArgumentParser()
    '''Positional argument: '''
    parser.add_argument('i', type=str, nargs=2 ,help='First position is for the ground truth \n the second are vcf-files.')


    '''Optional argument: '''
    parser.add_argument('-stats_in', type=str, help='<statistics.txt> as input.') #2
    parser.add_argument('-stats_out', type=str, help='Average performance of a tool as input.') #2
    parser.add_argument('-o', type=str, help='Output path')
    args = parser.parse_args()
    return args




def main():
    ''''''

    args = parsing()



    '''--------------First stage: analyse the vcf file--------------------'''

    software_name="simref"
    #For every vcf file, analyse the it(The for-loop is driven by a bash script(eval_vcf.sh))
    gt=read_file(args.i[0])
    vcf=read_file(args.i[1])

    results=eval_result(gt, vcf)

    print("Write the confusion table into the file...")
    write_file(results, args.o, software_name)


    print("--End--")

    '''--------------Second stage: average the statistics for every software in every coverage--------------------'''

    '''f = open(args.stats_out, 'a')
    file_name = args.stats_in
    spl=file_name.split("/")
    software_cov=spl[-1][0:-10]
    


    arr = read_stats(args.stats_in)
    avg=avg_perform(arr) #The average of every parameter
        

    f.write("{}   {:.4f}   {:.4f}   {:.4f}\n".format(software_cov,avg[0], avg[1], avg[2]))
    f.close()'''
    '''-------------Third stage: Plot the analysis on the canvas------------------------------------------------'''
    '''file_name = args.stats_in
    values, names=read_stats_3(file_name) #values: the table of real number; names: the list of software names
    n=len(names)//3
    values = np.asarray(values)

    sens=values[:,0]
    precs=values[:,1]
    f1s=values[:,2]

    draw_lines(sens, "Sensitivities", n)
    draw_lines(precs, "Precisions", n)
    draw_lines(f1s, "F1-scores", n)'''


    '''----Post processing----'''



main()