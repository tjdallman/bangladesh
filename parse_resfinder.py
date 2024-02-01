 #imports
import sys, re, math
from os import system
from datetime import datetime
from numpy  import *
from math import sin, cos, sqrt, atan2, radians
from pyproj import Proj, transform
import pandas as pd
from matplotlib import pyplot
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import glob


def usage():

    sys.exit()

def get_opts():
    if len(sys.argv) != 6:
        usage()
    else:
        return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

def read_genes(list_file):

    inhandle = open(list_file, 'r')
    yield_c = {}
    name = ""
    count = ""
    for each_line in inhandle.readlines():
        if each_line.startswith('>'):
            name = each_line[1:].strip()
            #print (name)
        else:
            count = len(each_line)
            #print (count)
            yield_c[name] = count

    return yield_c

def read_counts(count_file_dir):
    read_counts = {}
    for name in glob.glob(count_file_dir+"/*.txt"):
        temp = name.split("/")
        temp2 = temp[4].split(".")
        temp3 = temp2[0].split("_")
        samp = temp3[0]

        inhandle = open(name,'r')
        for line in inhandle.readlines():
            temp = line.split("\t")
            if int(temp[2]) > 0:
                try:
                    read_counts[samp][temp[0]] = temp[2]
                except:
                    read_counts[samp] = {}
                    read_counts[samp][temp[0]] = temp[2]
    return read_counts                

def read_res_db(res_db):
    class_map = {}
    class_names = []
    for name in glob.glob(res_db+"/*.fsa"):   
        temp = name.split("/")
        temp2 = temp[4].split(".")
        class_name = temp2[0]
        class_names.append(class_name)
        inhandle = open(name,'r')
        for each_line in inhandle.readlines():
            if each_line.startswith('>'):
                name = each_line[1:].strip()
                class_map[name] = class_name

    return class_map, class_names


def read_clusters(clust_file):
    inhandle = open(clust_file, 'r')
    clusters = {}
    for each_line in inhandle.readlines():
        if each_line.startswith('>'):
            cluster_name = each_line[1:].strip()
        else:
            temp = each_line.split()
            samp_name = temp[2][1:-3]
            try:
                clusters[cluster_name].append(samp_name)
            except:
                clusters[cluster_name] = []
                clusters[cluster_name].append(samp_name)
    return clusters            

def convert_to_fpkm(read_counts, gene_c, yield_c):
    fpkm_counts = {}
    for sample in read_counts:
        fpkm_counts[sample] = {}
        for hit in read_counts[sample]:
            fpkm = (float(read_counts[sample][hit]) * 1000000000) / (int(yield_c[sample]) * int(gene_c[hit]))
            fpkm_counts[sample][hit] = fpkm
    return fpkm_counts        

def read_yield(list_file):

    inhandle = open(list_file, 'r')
    yield_c = {}
    for each_line in inhandle.readlines():
        temp = (each_line.strip()).split('\t')
        yield_c[temp[0]] = temp[1]
    return yield_c


def aggregate_to_clusters(fpkm_counts, cluster_map):
    cluster_counts={}
    for cluster in cluster_map:
        for sample in fpkm_counts:
            try:
                cluster_counts[sample][cluster] = 0
            except:
                cluster_counts[sample] = {}
                cluster_counts[sample][cluster] = 0               
            for hit in fpkm_counts[sample]:
                if hit in cluster_map[cluster]:
                    cluster_counts[sample][cluster] = cluster_counts[sample][cluster] + float(fpkm_counts[sample][hit])                    
    return cluster_counts

def aggregate_to_class(fpkm_counts, class_map, classes):
    cluster_counts={}
    for class_name in classes:
        for sample in fpkm_counts:
            try:
                cluster_counts[sample][class_name] = 0
            except:
                cluster_counts[sample] = {}
                cluster_counts[sample][class_name] = 0               
            for hit in fpkm_counts[sample]:
                if class_map[hit] == class_name:
                    cluster_counts[sample][class_name] = cluster_counts[sample][class_name] + float(fpkm_counts[sample][hit])                    
    return cluster_counts    

def print_fpkm_table(cluster_counts):
    df = pd.DataFrame(cluster_counts)
    df = df.loc[~(df==0.0).all(axis=1)]
    #df.to_csv(r'resfinder_fkpm.tsv', sep='\t')
    df.to_csv(r'resfinder_fpkm_prot.tsv', sep='\t')


def rename_clusters(cluster_map, class_map):
    new_cluster_map = {}
    for cluster in cluster_map:
        if class_map[cluster_map[cluster][0]] != 'disinfectant':
            name = cluster+"|"+cluster_map[cluster][0]+"|"+class_map[cluster_map[cluster][0]]
            new_cluster_map[name] = cluster_map[cluster]

    return new_cluster_map
###### MAIN ###
 
gene_file, count_file_dir, clust_file, res_db, yields = get_opts()

gene_c = read_genes(gene_file)

read_counts = read_counts(count_file_dir)

class_map, classes = read_res_db(res_db)

cluster_map = read_clusters(clust_file)

cluster_map = rename_clusters(cluster_map, class_map)

yield_c = read_yield(yields)

fpkm_counts = convert_to_fpkm(read_counts, gene_c, yield_c)

#fpkm_counts = read_counts

cluster_counts = aggregate_to_clusters(fpkm_counts, cluster_map)

#cluster_counts = aggregate_to_class(fpkm_counts, class_map, classes)

print_fpkm_table(cluster_counts)






