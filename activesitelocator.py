#!/usr/bin/python

from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import reverse_complement, translate

import sys
import csv
import os
# Finds genomic DNA coordinates from protein coordinates.
# inputs: genbank file of protein, AA seq of protein, Active site of protein #, mutation to be made (Amino acid to amino acid)
# output: genbank file of protein, Start location of changed codon, current codon, end location, codon to change to
#
# kinase d->m
# kinase may need to have HRD seq to be active. Search for this
#



# function to retrieve genbank data from NCBI website
def genbank_retriever(accession): 
	Entrez.email = "gjordanray@gmail.com"     # tells NCBI who you are
	handle = Entrez.efetch(db="nucleotide", id = accession, rettype = "gb", retmode="xml") # rettype can be 'gb' or 'fasta'
	return handle

#function to separate a large genbank file into smaller 1 entry ones
def genbank_breakup(genbank_file):
	for rec in SeqIO.parse(genbank_file, "genbank"):
		SeqIO.write([rec], open(rec.name + ".gb", "w"), "genbank")

#function to read a genbank file and return the DNAseq, and all of the CDS featues
def genbank_read(genbank_file):
	gfile= SeqIO.read(genbank_file, "genbank")
	seq = gfile.seq
	#seq_str = str(file.seq)
	CDS= list()
	for feature in gfile.features:
		newfeature = feature.type.strip()
		if newfeature == "CDS":
			CDS.append(feature)
	return seq, CDS

# function to parse a long genbank file of multiple entries
def genbank_parse(genbank_file):
	names = list()
	seqs = list()
	CDSs = list()
	for record in SeqIO.parse(genbank_file, "genbank"):
		seq= record.seq
		CDS = list()
		protid = record.id
		protname = str(record.name)[:record.name.find('_')]
		name = str(protid)
		for feature in record.features:
			newfeature = feature.type.strip()
			if newfeature == "CDS":
				CDS.append(feature)
		seqs.append(str(seq))
		CDSs.append(CDS)
		names.append(name)
	return names, seqs, CDSs

def csv_parse(csv_file):
	names = list()
	ProteinSeqs = list()
	ActiveSites = list()
	InitialAminoAcids = list()
	FinalAminoAcids = list()
	counter =0;
	with open(csv_file, 'rU') as csvfile:
		csvreader = csv.reader(csvfile, delimiter = ',')
		for row in csvreader:
			counter +=1
			if counter ==1:
				continue
			names.append(row[0]),
			ProteinSeqs.append(row[1]),
			ActiveSites.append(row[2]),
			InitialAminoAcids.append(row[3]),
			FinalAminoAcids.append(row[4])
	return names, ProteinSeqs, ActiveSites, InitialAminoAcids, FinalAminoAcids


# Removes * at end of peptide if its there
def adjustAminoAcid(aminoAcidSeq):
	if aminoAcidSeq[-1]== '*':
		aminoAcidSeq = aminoAcidSeq[:-1]
	return aminoAcidSeq	


# return correct CDS of multiple given in genbankfile
def cdsFilter(cds_list,AminoAcidSeq,ActiveSite,InitialAminoAcid):
	for cds in cds_list:
		translation = cds.qualifiers['translation'][0]
		if AminoAcidSeq==translation and translation[ActiveSite-1]==InitialAminoAcid: # and cds.qualifiers['protein_id'][0][0] =='N':
			#print("All good")
			return cds
		else:
			continue
	#print("None of the CDSs matched")
	print cds_list[0]
	return 'NA'

# finds bp in dna seq based off of exon annotations 
def findBase(cds,ActiveSite):
	exons = cds.location
	basePairCount=1
	bpToActiveSite = 3*(ActiveSite-1)
	if exons.strand ==1:
		for exon in exons:
			if basePairCount != bpToActiveSite:
				basePairCount +=1 
			else:
				return exon+2
	elif exons.strand ==-1:
		dnaLoci= list()
		for exon in exons:
			dnaLoci.append(exon)
		dnaLoci.sort()
		dnaLoci = dnaLoci[::-1]
		for loc in dnaLoci:
			if basePairCount != bpToActiveSite:
				basePairCount +=1
			else:
				return loc


#Amino Acid to Codon for most common instance in humans
# may want to optimize to minimize changes from initial codon.
def reverseTranslator(AminoAcid):
	if AminoAcid=='M':
		base='ATG'
	elif AminoAcid == 'F':
		base ='TTC'
	elif AminoAcid == 'L':
		base ='CTG'
	elif AminoAcid == 'I':
		base ='ATC'
	elif AminoAcid == 'V':
		base ='GTG'
	elif AminoAcid == 'S':
		base ='AGC'
	elif AminoAcid == 'P':
		base ='CCC'
	elif AminoAcid == 'T':
		base ='ACC'
	elif AminoAcid == 'A':
		base ='GCC'
	elif AminoAcid == 'Y':
		base ='TAC'
	elif AminoAcid == '*':
		base ='TGA'
	elif AminoAcid == 'H':
		base ='CAC'
	elif AminoAcid == 'Q':
		base ='CAG'
	elif AminoAcid == 'N':
		base ='AAC'
	elif AminoAcid == 'K':
		base ='AAG'
	elif AminoAcid == 'D':
		base ='GAC'
	elif AminoAcid == 'E':
		base ='GAG'
	elif AminoAcid == 'C':
		base ='TGC'
	elif AminoAcid == 'W':
		base ='TGG'
	elif AminoAcid == 'R':
		base ='AGA'
	elif AminoAcid == 'G':
		base ='GGC'
	else:
		base =''
	return base

def find_nth(str1, substr1):
   return str1.find(substr1,str1.find(substr1)+1)


def concatenateNametoNum(name):
	start = name.find('_')
	stop = name.rfind('_')
	return name[start+1:stop]

# test that changed amino acid is the one we are talking about and that sites given map to correct bases
def test(dna,AA,codeline,iAA,fAA,cds,strand):
	firstBreak = codeline.find(',')
	secondBreak = codeline.find(',',firstBreak+1)
	thirdBreak = codeline.find(',',secondBreak+1)
	fourthBreak = codeline.rfind(',')

	codonStartSite = int(codeline[firstBreak+1:secondBreak])
	initialCodon = codeline[secondBreak+1:thirdBreak]
	codonEndSite = int(codeline[thirdBreak+1:fourthBreak])
	finalCodon = codeline[fourthBreak+1:]

	TranslatableInitialCodon = initialCodon
	TranslatableFinalCodon = finalCodon
	if strand ==-1:
		TranslatableInitialCodon = reverse_complement(initialCodon)
		TranslatableFinalCodon = reverse_complement(finalCodon)

	#TEST CASES
	if AA != cds.qualifiers['translation'][0]: #protein seqs match up
		print "AA seqs not equal"
		return False 
	elif dna[codonStartSite-1:codonEndSite] != initialCodon: #codon that is being modified is where its supposed to be
		print dna[codonStartSite-1:codonEndSite] +'!=' + initialCodon + "  :so codeline doesnt match up"
		return False
	elif translate(TranslatableInitialCodon) != iAA: #starting codon is what its supposed to be
		return False
	elif translate(TranslatableFinalCodon) != fAA: #final codon is what its supposed to be
		return False
	else: 
		return True




#######################################Executables#################################################

genbank_file = sys.argv[1] # The Genebank file of proteins is entered first
KinaseFile = sys.argv[2] # The csv file of proteins with order: Name,Protein,Active Site,Initial AA,Final AA,Prot Name,Entrez_GeneID,Anything else...
names, ProteinSeqs, ActiveSites, InitialAminoAcids, FinalAminoAcids = csv_parse(KinaseFile)

gbNames, gbSeqs, gbCDSs =  genbank_parse(genbank_file) #filter through large genbank and return lists of properties for each protein 
gbListLoc = list()
for name in names:
	gbListLoc.append(gbNames.index(concatenateNametoNum(name)))
CDSs = list()
strands = list()
ProtsThatDontWork = list()
GoodProts = list()
titeline = "name"  + "\t" + "ActiveSite " + "\t" + "InitialAminoAcid" + "\t" + "FinalAminoAcid" + "\t" + "ProteinID" + "\t" + "strand"
GoodProts.append(titeline)
n=0
for location in gbListLoc:
	ProteinSeqs[n]  = adjustAminoAcid(ProteinSeqs[n])
	CDSs.append(cdsFilter(gbCDSs[location],ProteinSeqs[n],int(ActiveSites[n]),InitialAminoAcids[n]))

	if CDSs[n] == 'NA':
		strands.append(0)
		ProtsThatDontWork.append([names[n]])
		n+=1
		continue
	else:
		strands.append(CDSs[n].location.strand)
		base = findBase(CDSs[n],int(ActiveSites[n]))
		PID = CDSs[n].qualifiers['protein_id'][0]
		GoodProts.append(names[n]+"\t" + str(base) + "\t" + InitialAminoAcids[n] + "\t" + FinalAminoAcids[n] + "\t" + PID + "\t" + str(strands[n]))# + "\t" + gbSeqs[location])
		output = open(names[n]+".fasta", "w")
		if strands[n] ==1:	
			codeline = (str(names[n])+ ',' + str(base) + ',' + gbSeqs[location][base-1] + gbSeqs[location][base] + 
				gbSeqs[location][base+1] + ',' + str(base+2)+ ',' + reverseTranslator(FinalAminoAcids[n]))
		else:
			codeline = (str(names[n])+ ',' + str(base-2) + ',' + gbSeqs[location][base-3] + gbSeqs[location][base-2] + 
				gbSeqs[location][base-1] + ',' + str(base)+ ',' + reverse_complement(reverseTranslator(FinalAminoAcids[n])))
		fastaInput = SeqRecord(Seq(gbSeqs[location]), id = codeline, description ="")
		print codeline + ' ' + str(strands[n])
		SeqIO.write(fastaInput,output,"fasta")
		output.close()
		if test(gbSeqs[location],ProteinSeqs[n],codeline,InitialAminoAcids[n],FinalAminoAcids[n],CDSs[n],strands[n]):
			print "passed"
		else:
			print "failed"
		n+=1
		continue
output = open('FailedProteins.csv',"wb")
writer = csv.writer(output)
for prot in ProtsThatDontWork:
	writer.writerow(prot)
output.close()

output = open('GoodProteins.csv',"wb")
writer = csv.writer(output)
for prot in GoodProts:
	writer.writerow(prot.split())
output.close()

