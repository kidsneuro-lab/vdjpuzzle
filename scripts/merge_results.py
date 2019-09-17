chains=['TRA','TRB','TRG','TRD','IGH','IGK','IGL']

import os
import sys


hasIGH=os.path.isfile(sys.argv[1] + "/IGH.tsv") 
hasIGK=os.path.isfile(sys.argv[1] + "/IGK.tsv") 
hasIGL=os.path.isfile(sys.argv[1] + "/IGL.tsv") 
hasTRA=os.path.isfile(sys.argv[1] + "/TRA.tsv") 
hasTRB=os.path.isfile(sys.argv[1] + "/TRB.tsv") 
hasTRG=os.path.isfile(sys.argv[1] + "/TRG.tsv") 
hasTRD=os.path.isfile(sys.argv[1] + "/TRD.tsv")

files=[]
if (hasIGH):
    files.append(("IGH","summary_corrected/IGH.tsv"))
if (hasIGK):
    files.append(("IGK","summary_corrected/IGK.tsv"))
if (hasIGL):
    files.append(("IGL","summary_corrected/IGL.tsv"))
if (hasTRA):
    files.append(("TRA","summary_corrected/TRA.tsv"))
if (hasTRB):
    files.append(("TRB","summary_corrected/TRB.tsv"))
if (hasTRG):
    files.append(("TRG","summary_corrected/TRG.tsv"))
if (hasTRD):
    files.append(("TRD","summary_corrected/TRD.tsv"))


Nreads={}
#read number of reads in a cell
with open("nreads.csv", 'r') as seqCell:
    for line in seqCell:
        tokens=line.split('\t')
        Nreads[tokens[0]]=float(tokens[1].replace("\n",""))/4
    


for f in files:
    with open(f[1], 'r') as receptorFile:
        linecount=0
        extraCSV=open(f[1]+".extra", 'w')
        extraCSV.write("Isotype/Constant\tMembrane/Secreted\tExpression\n")
        for line in receptorFile:
            linecount+=1
            if linecount>1:
                cellID=line.split('\t')[0]
                seqID=line.split('\t')[1]
                #seqID=seqID.split(' ')[0] #remove everything after space
                #seqID=seqID[1::] #remove > symbol
                cellID=cellID[4::]
                #find constant
                blastfile="blast/"+f[0]+"_"+cellID+".out"
                hmmfile="blast/IGH_"+cellID+"membrane.out"
                constant=""
                max_identity=0
                mem_sec="NA"
                if f[0]=="IGH":
                    mem_sec="Not found"
                    hasIGH=os.path.isfile(hmmfile) 
                    if hasIGH:
                        with open(hmmfile, 'r') as hmmOutput:
                            for line3 in hmmOutput:
                                if seqID in line3:
                                    if 'membrane' in line3:
                                        mem_sec="Membrane" 
                                    else:
                                        mem_sec="Secreted"
                if os.path.isfile(blastfile):
                    with open(blastfile, 'r') as blastOutput:
                        for line2 in blastOutput:
                            if not line2[0]=="#":
                                token=line2.split('\t')
                                identity=float(token[2])
                                if (token[0]==seqID) and (max_identity<identity):
                                    max_identity=identity
                                    constant=token[1].split('|')[1]
                #find expression
                expressionfile="kallisto/"+f[0]+"_"+cellID+"/"+f[0]+"_out/abundance.tsv"
                tpm="0"
                if os.path.isfile(expressionfile):
                    with open(expressionfile, 'r') as expressionOutput:
                        for line2 in expressionOutput:
                            if not line2[0]=="#":
                                token=line2.split('\t')
                                if (token[0]==seqID) and (tpm=="0"):
                                    rpm=float(token[3].replace("\n",""))/(Nreads[cellID]/1000000) 
                                    rpkm=rpm/(float(token[2].replace("\n","")))
                extraCSV.write(constant+"\t"+mem_sec+"\t"+str(rpkm)+"\n")
        extraCSV.flush()
        extraCSV.close()
    cmd="paste "+f[1]+" "+f[1]+".extra > "+f[1]+".final"
    os.system(cmd)
    os.system("mv "+f[1]+".final "+f[1])
os.system("rm summary_corrected/*.extra")
