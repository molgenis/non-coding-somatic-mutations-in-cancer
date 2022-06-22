import sys
import gzip

def overlap(c1,b1,e1,c2,b2,e2):
	"""

    :param : 
    :param :  
    :param :        
    :return:    
    """
	if c1 != c2:
		return 0
	if e2 < b1 or b2 > e1:
		return 0
	if e2 > b1 and e2 < e1:
		return 1
	if b2>b1 and e2<e1:
		return 1
	if b2<e1 and e2>e1:
		return 1
	return 0

genedata = {}
objs = []
fh = gzip.open('allgenes.tsv.gz','rt')
fh.readline()
for line in fh:
	elems = line.strip().split("\t")
	chr = elems[0]
	sta = int(elems[1])
	sto = int(elems[2])
	nam = elems[3]
	obj = [chr,sta,sto,nam,0,0,0]
	chrobjs = genedata.get(chr)
	if chrobjs is None:
		chrobjs = []
		genedata[chr] = chrobjs
	chrobjs.append(obj)
	objs.append(obj)
fh.close()



fh = gzip.open('Region_2000_ALL_both_0_TESTS.tsv.gz','rt')
fh.readline()
lnctr = 0
for line in fh:
	elems = line.strip().split("\t")
	chr = elems[0]
	sta = round(float(elems[1]))
	sto = int(elems[2])
	ct1 = float(elems[3])
	ct2 = float(elems[4])
	chrobjs = genedata.get(chr)
	if chrobjs is not None:
		if ct1 > 0 or ct2 > 0:
			for gene in chrobjs:
				if overlap(gene[0],gene[1],gene[2],chr,sta,sto) == 1:
					gene[4] += ct1
					gene[5] += ct2
					gene[6] += 1
	lnctr += 1
	if lnctr % 10000 == 0:
		print("{} lines parsed".format(lnctr), end='\r')
fh.close()

fh = open('output.txt','w')
for gene in objs:
	chr = gene[0]
	sta = gene[1]
	sto = gene[2]
	nam = gene[3]
	ct1 = gene[4]
	ct2 = gene[5]
	ct3 = gene[6]
	outln = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr,sta,sto,nam,ct1,ct2,ct3)
	fh.write(outln)
fh.close()
