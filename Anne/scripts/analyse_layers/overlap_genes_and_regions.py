#!/usr/bin/env python3

#Imports
import sys
import gzip

sys.path.append('D:/Hanze_Groningen/STAGE/00git/non-coding-somatic-mutations-in-cancer/Anne/scripts/')
from config import get_config


def overlap(c1, b1, e1, c2, b2, e2):
    """
	Get overlap 
    :param c1: chromosome of the gene
    :param b1:  start position of the gene
    :param e1: stop positio of the gene
	:param c2: chromosome of the region
    :param b2:  start position of the region
    :param e2:  stop position of the region      
    :return:    
    """
    # return 0 if there is no overlap
    # return 1 if there is overlap
    if c1 != c2:
        return 0
    if e2 < b1 or b2 > e1:
        return 0
    if e2 > b1 and e2 < e1:
        return 1
    if b2 > b1 and e2 < e1:
        return 1
    if b2 < e1 and e2 > e1:
        return 1
    return 0


def get_info_genes(config):
    """
	Gets the information about the genes
	:param config: Dictionary with as keys the name of the paths and as value the paths   
    :return: objs: a list with chromosome, start, stop and name of the gene
	         genedata: dictionary with as key chromosome and as value a list (chrobjs)
	"""
    genedata = {}
    objs = []
    fh = gzip.open(config['allgenes'], 'rt')
    fh.readline()
    for line in fh:
        elems = line.strip().split("\t")
        # Get information out of file
        chr = elems[0]
        sta = int(elems[1])
        sto = int(elems[2])
        nam = elems[3]
        obj = [chr, sta, sto, nam, 0, 0, 0]
        chrobjs = genedata.get(chr)
        if chrobjs is None:
            chrobjs = []
            genedata[chr] = chrobjs
        chrobjs.append(obj)
        objs.append(obj)
    fh.close()
    return objs, genedata


def get_info_regions(genedata, config):
    """
	Gets the information about the regions (chromosome divided into 2000 bp pieces)
	:param genedata: dictionary with as key chromosome and as value a list (chrobjs)
	:param config: Dictionary with as keys the name of the paths and as value the paths  
    :return: df_whole: Dataframe with all the information out of the database
	"""
    fh = gzip.open(config['Region_2000'], 'rt')
    fh.readline()
    lnctr = 0
    for line in fh:
        elems = line.strip().split("\t")
        # Get information out of file
        chr = elems[0]
        sta = round(float(elems[1]))
        sto = int(elems[2])
        ct1 = float(elems[3])
        ct2 = float(elems[4])
        chrobjs = genedata.get(chr)
        if chrobjs is not None:
            if ct1 > 0 or ct2 > 0:
                for gene in chrobjs:
                    # Search for overlap
                    if overlap(gene[0], gene[1], gene[2], chr, sta, sto) == 1:
                        # Adds the amount of breast counts
                        gene[4] += ct1
                        # Adds the amount of nonbreast counts
                        gene[5] += ct2
                        # Keeps track of how many regions a gene consists of
                        gene[6] += 1
        lnctr += 1
        if lnctr % 10000 == 0:
            print("{} lines parsed".format(lnctr), end='\r')
    fh.close()


def write_output(objs, config):
    """
	Writes the output
	:param objs: a list with chromosome, start, stop and name of the gene
	:param config: Dictionary with as keys the name of the paths and as value the paths  
    :return: df_whole: Dataframe with all the information out of the database
	"""
    fh = open(config['output_overlap'], 'w')
    for gene in objs:
        chr = gene[0]
        sta = gene[1]
        sto = gene[2]
        nam = gene[3]
        ct1 = gene[4]
        ct2 = gene[5]
        ct3 = gene[6]
        outln = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr, sta, sto, nam, ct1, ct2, ct3)
        fh.write(outln)
    fh.close()


def main():
    # Call config
    config = get_config('Anne')
    objs, genedata = get_info_genes(config)
    get_info_regions(genedata, config)
    write_output(objs, config)


if __name__ == '__main__':
    main()
