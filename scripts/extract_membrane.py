import sys

with open(sys.argv[2], 'r') as blast:
    #for each contigs
    #extract qend (maybe it is better to do q.end-10)
    readID=""
    qEnd={}
    MaxIdentity={}
    for line in blast.readlines():
        if line.startswith("# Query"):
            readID=line.split("# Query: ")[1]
            readID=readID.replace("\n","").replace("\r","")
            qEnd[readID]=""
            MaxIdentity[readID]=0
        elif line.startswith("#"):
            #do nothing
            pass
        else:
            tokens=line.split("\t")
            identity=float(tokens[2])
            if identity>MaxIdentity[readID]:
                MaxIdentity[readID]=identity
                qEnd[readID]=int(tokens[7])

with open(sys.argv[1], 'r') as fasta:
    membrane=""
    for line in fasta.readlines():
        if line[0]==">":
            if len(membrane)>20:
                print readID
                print membrane
            readID=line.replace("\n","").replace("\r","")
            membraneStart=qEnd[readID[1::]]
            counting=0
            membrane=""
            isMembrane=False
        else:
            line2=line.replace("\n","").replace("\r","")
            if isMembrane or ((len(line2)+counting)<membraneStart):
                counting+=len(line2)
            else:
                membrane+=line2[(membraneStart-counting)::]
    
    if len(membrane)>20:
        print readID
        print membrane