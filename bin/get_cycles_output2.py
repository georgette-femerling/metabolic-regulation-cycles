from collections import defaultdict 
from datetime import date
import re

start_tf = "Nac"

direct = defaultdict(list)
indirect = defaultdict(list)

with open("Nac_Cycles_Report_o1_v3") as file2:
    for line in file2:
        if re.match("#", line):
            

        else:
            TF = line.split("\t")[0]
            gene = line.split("\t")[1]
            MF = line.split("\t")[2].strip("\n")
            if TF == start_tf:
                direct[MF].append(gene)
            else:
                indirect[MF].append((TF,gene))


out = open("Nac_Cycles_Report3_PotInt.txt","w")
out.write("# From get_cycles_output2.py \n")
out.write("# Date: "+str(date.today())+"\n")
out.write("# Genes marked with * participate in a reaction with a DeltaG less than -58.5\n")
out.write("#\n")
out.write("# Input TF: Nac\n")
out.write("# Start Multifun: BC-1.8.3: nitrogen metabolism\n")
out.write("#\n")
out.write("Multifun Term(genes in MF)\tDirect Regulation\tIndirect Regulation(TF/gene)\n")
for mf in set(direct.keys()+indirect.keys()):
    if mf in direct.keys() and mf not in indirect.keys():
         out.write(mf+"\t"+",".join(direct[mf])+"\tNA\n")
    elif mf not in direct.keys() and mf in indirect.keys():
        out.write(mf+"\tNA\t"+",".join( "/".join(indir) for indir in indirect[mf])+"\n")
    elif mf in direct.keys() and mf in indirect.keys():
        out.write(mf+"\t"+",".join(direct[mf])+"\t"+",".join( "/".join(indir) for indir in indirect[mf])+"\n")
out.close()