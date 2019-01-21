chains=['TRA','TRB','TRG','TRD','IGH','IGK','IGL']

import os
import sys


hasIGH=os.path.isfile("summary_corrected/IGH.csv") 
hasIGK=os.path.isfile("summary_corrected/IGK.csv") 
hasIGL=os.path.isfile("summary_corrected/IGL.csv") 
hasTRA=os.path.isfile("summary_corrected/TRA.csv") 
hasTRB=os.path.isfile("summary_corrected/TRB.csv") 
hasTRG=os.path.isfile("summary_corrected/TRG.csv") 
hasTRD=os.path.isfile("summary_corrected/TRD.csv") 

files=[]
if (hasIGH):
    files.append(("IGH","summary_corrected/IGH.csv"))
if (hasIGK):
    files.append(("IGK","summary_corrected/IGK.csv"))
if (hasIGL):
    files.append(("IGL","summary_corrected/IGL.csv"))
if (hasTRA):
    files.append(("TRA","summary_corrected/TRA.csv"))
if (hasTRB):
    files.append(("TRB","summary_corrected/TRB.csv"))
if (hasTRG):
    files.append(("TRG","summary_corrected/TRG.csv"))
if (hasTRD):
    files.append(("TRD","summary_corrected/TRD.csv"))


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
        extraCSV.write("\tIsotype/Constant\tMembrane/Secreted\tExpression\n")
        for line in receptorFile:
            linecount+=1
            if linecount>1:
                cellID=line.split('\t')[0]
                seqID=line.split('\t')[1]
                seqID=seqID.split(' ')[0] #remove everything after space
                seqID=seqID[1::] #remove > symbol
                cellID=cellID[4::]
                #find constant
                blastfile="blast/"+f[0]+"_"+cellID+".out"
                constant=""
                max_identity=0
                mem_sec="NA"
                if f[0]=="IGH":
                    mem_sec="Secreted"
                if os.path.isfile(blastfile):
                    with open(blastfile, 'r') as blastOutput:
                        for line2 in blastOutput:
                            if not line2[0]=="#":
                                token=line2.split('\t')
                                identity=float(token[2])
                                if (token[0]==seqID) and (max_identity<identity):
                                    constant=token[1].split('|')[1]
                                    if ("_M" in token[1].split('|')[0]) and (f[0]=="IGH"): #grepl("_M", token[1].split('|')[0]):
                                        mem_sec="Membrane-bound"
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
                extraCSV.write("\t"+constant+"\t"+mem_sec+"\t"+str(rpkm)+"\n")
        extraCSV.flush()
        extraCSV.close()
    cmd="paste "+f[1]+" "+f[1]+".extra > "+f[1]+".final"
    os.system(cmd)
    os.system("mv "+f[1]+".final "+f[1])
os.system("rm summary_corrected/*.extra")





