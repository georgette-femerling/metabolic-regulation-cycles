#!/usr/bin/env python
# coding: utf-8

###################################################
##  get_cycles_output2                           ##
## Author: Georgette Femerling                   ##
## Version: v2                                   ##
## Date: 7-09-2020                               ##
## Description: Takes output report from         ##
##              get_cycles.pl from regulon db    ##
##              and from new network, parses and ##
##              compares them to generate a      ##
##              comparison table                 ##
## Input: [1] Start TF                           ##
##        [2] Report from new Network            ##
##        [3] Report from RegulonDB Network      ##
##        [4] Output Directory                   ##
## Output: One table enlisting in columns        ##
##       [1] Dominant MF(genes in MF)            ##
##       [2] Multifun Term(genes in MF)          ##
##       [3] Direct Regulation known             ##
##       [4] Direct Regulation unknown           ##
##       [5] Indirect Regulation(TF/gene) known  ##
##       [6] Indirect Regulation(TF/gene) unknown##
##  *** All known/unknown interactions are       ##
##      relative to the regulonDB report         ##
###################################################

# Example run as:
# python bin/get_cycles_output2.v2.py Nac Cycles/Nac_Cycles_Report_o1_v3.txt Cycles/Nac_Cycles_Report_Regulon_o1_v3.txt Cycles/

# Import Modules
from collections import defaultdict 
from datetime import date
import sys
import os
import re

start_tf = sys.argv[1]
report_new = sys.argv[2]
report_regulon = sys.argv[3]
output_dir = sys.argv[4]

new = defaultdict(lambda:defaultdict(lambda: [set(), set()])) #dictionary of dictionary of list of two lists
with open(report_new) as file2:
    for line in file2:
        if re.match("\n",line):
            continue
        elif re.match("#", line):
            if re.match("# TF name:", line):
                tf = line.split("\t")[1].strip("\n")
            elif re.match("# Multifun ID:", line):
                MFid = line.split("\t")[1].strip("\n")
            elif re.match("# Multifun:", line):
                MF = line.split("\t")[1].strip("\n")
                dom = (MFid,MF)
        else:
            TF = line.split("\t")[0]
            gene = line.split("\t")[1]
            MF = line.split("\t")[2].strip("\n")
            if TF == start_tf:
                new[dom][MF][0].add(gene)
            else:
                new[dom][MF][1].add((TF,gene))

regulon = defaultdict(lambda:defaultdict(lambda: [set(), set()])) #dictionary of dictionary of list of two lists
with open(report_regulon) as file1:
    for line in file1:
        if re.match("\n",line):
            continue
        elif re.match("#", line):
            if re.match("# TF name:", line):
                tf = line.split("\t")[1].strip("\n")
            elif re.match("# Multifun ID:", line):
                MFid = line.split("\t")[1].strip("\n")
            elif re.match("# Multifun:", line):
                MF = line.split("\t")[1].strip("\n")
                dom = (MFid,MF)
        else:
            TF = line.split("\t")[0]
            gene = line.split("\t")[1]
            MF = line.split("\t")[2].strip("\n")
            if TF == start_tf:
                regulon[dom][MF][0].add(gene)
            else:
                regulon[dom][MF][1].add((TF,gene))


out = open(str(output_dir+start_tf+"_cycles_report_o2.txt"),"w")
out.write("# From get_cycles_output2.py \n")
out.write("# Date: "+str(date.today())+"\n")
out.write("# Genes marked with * participate in a reaction with a DeltaG less than -58.5\n")
out.write("# Multifun Terms marked with * are new relative to RegulonDB network\n")
out.write("#\n")
out.write("# Input TF: Nac\n")
out.write("#\n")
out.write("Dominant MF(genes in MF)\tMultifun Term(genes in MF)\tDirect Regulation known\tDirect Regulation unknown\tIndirect Regulation(TF/gene) known\tIndirect Regulation(TF/gene) unknown\n")
for dom in new.keys():
    if dom in regulon.keys():
        for mf in new[dom]:
            if mf in regulon[dom]:
                knowdir = regulon[dom][mf][0] & new[dom][mf][0]
                knowdir.update(regulon[dom][mf][0].difference(new[dom][mf][0]))
                unknowndir = new[dom][mf][0].difference(regulon[dom][mf][0])
                knownindir = regulon[dom][mf][1] & new[dom][mf][1]
                knownindir.update(regulon[dom][mf][1].difference(new[dom][mf][1]))
                unknownindir = new[dom][mf][1].difference(regulon[dom][mf][1])
                out.write(":".join(dom)+"\t"+mf+"\t"+",".join(knowdir)+"\t"+",".join(unknowndir)+"\t"+",".join("/".join(indir) for indir in knownindir)+"\t"+",".join("/".join(indir) for indir in unknownindir)+"\n")
            else:
                out.write(":".join(dom)+"\t"+"*"+mf+"\t"+"\t"+",".join(new[dom][mf][0])+"\t"+"\t"+",".join("/".join(indir) for indir in new[dom][mf][1])+"\n")
    else:
        for mf in new[dom]:
            if mf in regulon.values():
                out.write(":".join(dom)+"\t"+mf+"\t"+"\t"+",".join(new[dom][mf][0])+"\t"+"\t"+",".join("/".join(indir) for indir in new[dom][mf][1])+"\n")
            else:
                out.write(":".join(dom)+"\t"+"*"+mf+"\t"+"\t"+",".join(new[dom][mf][0])+"\t"+"\t"+",".join("/".join(indir) for indir in new[dom][mf][1])+"\n")
out.close()