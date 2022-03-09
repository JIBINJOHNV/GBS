
import pandas as pd
import glob
import argparse
import subprocess as sp
import os
import subprocess

print("Please index your reference genome using bowtie; eg bowtie2-build -f genome.fasta genome")


parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
parser.add_argument('-FastaPath','--FastaPath', help="Fasta file path should end with/: eg /data/reference/", required=True)
parser.add_argument('-FastFile','--FastFile', help="Reference fasta file", required=True)
parser.add_argument('-R1Site','--R1Site', help="Restriction Enzyme 1 site: eg G^AATTC", required=True)
parser.add_argument('-R2Site','--R2Site', help="Restriction Enzyme 2 site: eg T^TAA", required=True)
parser.add_argument('-R1Enzyme','--R1Enzyme', help="First Restriction Enzyme name: eg EcoRI", required=True)
parser.add_argument('-R2Enzyme','--R2Enzyme', help="Second Restriction Enzyme name: eg MspI", required=True)
parser.add_argument('-LSize','--LSize', help="First Restriction Enzyme name: eg 90", required=True)
parser.add_argument('-USize','--USize', help="Second Restriction Enzyme name: eg 130", required=True)
parser.add_argument('-MQuality','--MQuality', help="Minimum mapping quality: eg 20", required=True)
parser.add_argument('-OutName','--OutName', help="Name of output file/folder: eg Mouse", required=True)


args=parser.parse_args()

R1Enzyme=args.R1Enzyme
R2Enzyme=args.R2Enzyme
R1Site=args.R1Site
R2Site=args.R2Site
LSize=args.LSize
USize=args.USize
FastFile=args.FastFile
MQuality=args.MQuality
OutFile=args.OutName
FastaPath=args.FastaPath


LSize=int(LSize)
USize=int(USize)
MQuality=int(MQuality)

###########
Ans=input("Is reference genome indexed: ")
if Ans.upper()=="YES":
    pass
else:
    sys.exit()

############



fragmatic="/data/Software_Resources/fragmatic/" #Reference: https://github.com/tkchafin/fragmatic/blob/master/fragmatic.pl


#R1Enzyme="EcoRI"
#R2Enzyme="MsPI"
#R1Site="G^AATTC"
#R2Site="T^TAA"
#LSize=90
#USize=130
#FastFile="genome.fas"
#OutFile="TEstOut"
#MQuality=20

#Split the genome fasta file by chromosome/scafolde
os.system("faidx -x "+ FastaPath+FastFile)

#dOUBLE DIGESTION 
for fasta in glob.glob("*"+FastFile.split(".")[-1]):
    if fasta==FastFile:
        pass
    else:
        cmd="perl "+fragmatic+"fragmatic.pl -i " + fasta + " -r " + '"'+R1Site + " " + R2Site +'"'+" -f -o " + OutFile+fasta
        os.system(cmd)


#mERGING CHROMOSOME WISE DOUBLE DIGESTED FASTA FILE
F1OutFile=OutFile+"_"+R1Enzyme+"_"+R2Enzyme+"_DoubleDigested.fa"
F2OutFile=OutFile+"_"+R2Enzyme+"_"+R1Enzyme+"_DoubleDigested.fa"

os.system("cat *"+R1Site.replace("^", '')+"-"+R2Site.replace("^", '')+"*" + "> "+F1OutFile)
os.system("cat *"+R2Site.replace("^", '')+"-"+R1Site.replace("^", '')+"*" + "> "+F2OutFile)

os.system("cat "+ F1OutFile +" " +F2OutFile +" > temp") 

os.system('awk \'BEGIN{RS=">"}{if(NR>1)print ">"$1"_"(NR-1)"\\n"$2}\' temp > DoubleDigested.fasta' )
os.system("rm temp")


#Remove intermediate files
os.system("rm *"+R1Site.replace("^", '')+"* "+ "*"+R2Site.replace("^", '')+"* " +"; rm *Missing*" +"; rm *tsv*" )


#Extract by size
DoubleDigeste="DoubleDigested.fasta"
cmd="awk -v LSize="+str(LSize) + " -v USize=" + str(USize) +" -F_ " +"'$3> "+str(LSize) +" && $3< " + str(USize)+" {print $0}' " + DoubleDigeste

FirstPrefix=DoubleDigeste[:-6]+str(LSize)+"_"+str(USize)
os.system(cmd +"> "+FirstPrefix+".txt")


os.system("sed -i 's/>//g' " + FirstPrefix+".txt" )
os.system("seqtk subseq "+DoubleDigeste + " "+ FirstPrefix+".txt" + " > "+FirstPrefix+".fasta")
os.system("fastx_collapser < "+FirstPrefix+".fasta > "+FirstPrefix+"_collapsed.fasta")




FirsName='DoubleDigested_' + R1Enzyme + "_"+R1Site.split("^")[1] +"_" +R2Enzyme+"_"+R2Site.split("^")[0]+"_"+str(LSize)+"_"+str(USize)+".unique.fasta"
Seconame='DoubleDigested_' + R2Enzyme + "_" +R2Site.split("^")[1] +"_"+R1Enzyme+"_"+R1Site.split("^")[0]+"_"+str(LSize)+"_"+str(USize)+".unique.fasta"


#Extracting fasta sequence by RE site order; like EcoR-MseI and MseI-EcoRI
cmd="cat "+FirstPrefix+"_collapsed.fasta | grep "+str('^')+R1Site.split("^")[1] +" | grep "+R2Site.split("^")[0]+"$ |" + 'awk \'{print \">\"NR\"\\n"$0}\' > ' +FirsName
os.system(cmd)
cmd="cat "+FirstPrefix+"_collapsed.fasta | grep "+str('^')+R2Site.split("^")[1] +" | grep "+R1Site.split("^")[0]+"$ |" + 'awk \'{print \">\"NR\"\\n"$0}\' > ' +Seconame
os.system(cmd)



#Align the digested products to the genome
cmd="bowtie2 -p 8 -x " + FastFile.split(".")[0] +" -f " + FirsName  + " -S "+ FirsName[:-6]+".sam 2> " +FirsName[:-6]+"_stats.txt"
os.system(cmd)
cmd="bowtie2 -p 8 -x " + FastFile.split(".")[0] +" -f " + Seconame  + " -S "+ Seconame[:-6]+".sam 2> " +Seconame[:-6]+"_stats.txt"
os.system(cmd)


#Extract the chromosome coordinate of double digested fragments
cmd='grep -v "@"  '+ FirsName[:-6]+".sam " + \
"| awk -v var="+str(MQuality) +" -F\"\t\" '$5>var" +"{print $3,$4,$4+length($10),$4+length($10)-$4,$2}' | sort -g -k1,2  > " + FirsName[:-6]+".bed "
os.system(cmd)

cmd='grep -v "@"  '+ Seconame[:-6]+".sam " + \
"| awk -v var="+str(MQuality) +" -F\"\t\" '$5>var" +"{print $3,$4,$4+length($10),$4+length($10)-$4,$2}' | sort -g -k1,2  > " + Seconame[:-6]+".bed "
os.system(cmd)



#Genome coverage calculation
Gsize = sp.getoutput('grep -v ">" '+FastFile +" | wc -c ")


Fsize=sp.getoutput("awk '{s+=$4}END{print s}' "+ FirsName[:-6]+".bed ")
FTotalfragment=sp.getoutput("cat " +FirsName[:-6]+".bed | wc -l")
FgenomeCovered=int(Fsize)/int(Gsize)*100

os.system("awk '$5==0{print $0}' "+ FirsName[:-6]+".bed >" + FirsName[:-6]+"_NoRevComp.bed" )
F_NR_size=sp.getoutput("awk '{s+=$4}END{print s}' "+ FirsName[:-6]+"_NoRevComp.bed")
F_NR_Totalfragment=sp.getoutput("cat " +FirsName[:-6]+"_NoRevComp.bed | wc -l")
F_NR_genomeCovered=int(F_NR_size)/int(Gsize)*100



Ssize=sp.getoutput("awk '{s+=$4}END{print s}' "+ Seconame[:-6]+".bed ")
STotalfragment=sp.getoutput("cat " +Seconame[:-6]+".bed | wc -l")
SgenomeCovered=int(Ssize)/int(Gsize)*100


os.system("awk '$5==0{print $0}' "+ Seconame[:-6]+".bed >" + Seconame[:-6]+"_NoRevComp.bed" )
S_NR_size=sp.getoutput("awk '{s+=$4}END{print s}' "+ Seconame[:-6]+"_NoRevComp.bed")
S_NR_Totalfragment=sp.getoutput("cat " +Seconame[:-6]+"_NoRevComp.bed | wc -l")
S_NR_genomeCovered=int(S_NR_size)/int(Gsize)*100





Tsize=int(Fsize)+int(Ssize)
TTotalfragment=int(FTotalfragment)+int(STotalfragment)
TgenomeCovered=int(Tsize)/int(Gsize)*100



T_NR_size=int(F_NR_size)+int(S_NR_size)
T_NR_Totalfragment=int(F_NR_Totalfragment)+int(S_NR_Totalfragment)
T_NR_genomeCovered=int(T_NR_size)/int(Gsize)*100



##Printing output

for file in ["Enzyme1-2.csv",'Enzyme2-1.csv',"Enzyme1-2-and-2-1.csv"]:
    if os.path.exists(file):
        pass
    else:
        os.system("echo GenomeSize,FragmentSize,FragmentCount,Genomecovered,EnzymeCombination-Size >"+file)


os.system("echo "+str(Gsize)+","+str(Fsize) +","+str(FTotalfragment)+","+ str(FgenomeCovered)+","+FirsName[15:-13] +">>Enzyme1-2.csv")
os.system("echo "+str(Gsize)+","+str(Ssize) +","+str(STotalfragment)+","+ str(SgenomeCovered)+","+Seconame[15:-13] +">>Enzyme2-1.csv")
os.system("echo "+str(Gsize)+","+str(Tsize) +","+str(TTotalfragment)+","+ str(TgenomeCovered)+","+Seconame[15:-19]+FirsName[15:-13]  +">> Enzyme1-2-and-2-1.csv")


if os.path.exists("dd-RAD.csv"):
    pass
else:
    os.system("echo GenomeSize,E1-2Fragment,E1-2FragementCount,E1-2GenomeCovered,E1-E2Names-Size,\
    E2-1Fragment,E2-1FragementCount,E2-1GenomeCovered,E2-E1Names-Size, \
    E12and21Fragment,E12and21FragementCount,E12and21GenomeCovered,E1E2andE2E1Combination-Size >"+"dd-RAD.csv")


os.system("echo "+str(Gsize)+","+str(F_NR_size)+","+str(F_NR_Totalfragment)+",\
"+str(F_NR_genomeCovered)+","+str(FirsName[15:-13])+","+str(S_NR_size)+","+str(S_NR_Totalfragment)+",\
"+str(S_NR_genomeCovered)+","+str(Seconame[15:-13])+","+str(T_NR_size)+","+str(T_NR_Totalfragment)+",\
"+str(T_NR_genomeCovered)+","+Seconame[15:-19]+FirsName[15:-13] +">>dd-RAD.csv")


###
os.system("mkdir -p "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))
os.system(" mv *bed "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))
os.system(" mv *_stats.txt "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))



os.system("samtools view -u "+FirsName[:-6]+".sam | samtools sort > "+FirsName[:-6]+".bam ")
os.system("samtools view -u "+Seconame[:-6]+".sam | samtools sort > "+Seconame[:-6]+".bam ")
os.system("samtools index " + FirsName[:-6]+".bam " + FirsName[:-6]+".bai ")
os.system("samtools index " + Seconame[:-6]+".bam " + Seconame[:-6]+".bai ")



os.system(" mv *bam "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))
os.system(" mv *bai "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))

os.system(" mv DoubleDigested.fasta "+R1Enzyme+"_"+R2Enzyme+"_"+str(LSize)+"_"+str(USize))

os.system("rm *DoubleDigested*") 


for fasta in glob.glob("*"+FastFile.split(".")[-1]):
    if fasta==FastFile:
        pass
    else:
        os.system("rm " + fasta)
