#!/usr/bin/env python
import re
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Alphabet import generic_dna, generic_rna
from Bio.Data.CodonTable import unambiguous_dna_by_id
import sys, subprocess, os, tempfile, copy

from Crispr import *

def sg_to_seqfeature( target, sg ):
	location = target.seq.find( sg.protospacer.back_transcribe() )
	if location == -1:
		strand = -1
		location = target.seq.find( sg.protospacer.back_transcribe().reverse_complement() )
	else:
		strand = 1
	if location == -1:
		print location, sg.protospacer.back_transcribe().reverse_complement(), sg.target_seq
		print "Could not find sg in target!"
		return None
	id = sg.protospacer
	quals = { "protospacer": sg.protospacer, "score": sg.score }
	feature = SeqFeature( FeatureLocation( location, location+len(sg.protospacer)), strand=strand, type="sgRNA", id=id, qualifiers=quals )
	return feature

def read_ensembl_tss( fname ):
#Ensembl Gene ID,Ensembl Transcript ID,Associated Gene Name,Chromosome Name,Strand,Transcription Start Site (TSS),APPRIS principal isoform annotation
	handle = open( fname, 'rU')
	tss_dict = {}
	for line in handle:
		if line[0] == "#":
			continue
		geneid, txid, genename, chr, strand, tss, pi = line.rstrip().split(",")
		#do some format conversion
		chr = "chr"+chr
		if strand == '1':
			strand = '+'
		else:
			strand = '-'
		tss = int(tss)
		if pi == '':
			pi = True
		else:
			pi = False

		fake_loc = GenomicLocation( chr, 0, 0, strand )
		fake_exon = Exon( fake_loc, False )
		fake_tx = Transcript( txid, 0, 0, tss, pi, [fake_exon] )
		try:
			tss_dict[txid+"_"+genename] = fake_tx
		except KeyError:
			tss_dict.set_default( txid+"_"+genename, fake_tx )
	handle.close()
	return tss_dict

# start of main script
target_fname = sys.argv[1]
tss_fname = sys.argv[2]
target_handle = open( target_fname, "r" )
genome = SeqIO.to_dict(SeqIO.parse( target_handle, "fasta", alphabet=generic_dna))
#targets = list( rec.upper() for rec in SeqIO.parse( target_handle, "fasta", alphabet=generic_dna))
out_fname = sys.argv[3]

tss_dict = read_ensembl_tss( tss_fname )

weissman_constant = "GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU" # weissman-style constant region for CRISPRi/a
broad_constant = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU"

# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
sense_re=re.compile(r'(?=([ATGCatgc]{20})([ATGCatgc]GG))') #  regex with lookahead to get overlapping sequences. group1 is protospacer, group2 is pam
antisense_re=re.compile(r'(?=(CC[ATGCatgc])([ATGCatgc]{20}))') # group1 is pam, group2 is protospacer

out_fhandle = open( out_fname, mode="w" )
guides = {}
for id,tx in tss_dict.iteritems():
	chromosome = tx.get_chromosome()
	tss = tx.tss
	strand = tx.get_strand()
	try:
		target = genome[chromosome]
	except KeyError:
		print "Chromosome %s not found!" % chromosome
		continue
		
	if strand == "+":
		search_start, search_end = [tss-50,tss+300]
	elif strand == "-": 
		search_start, search_end = [tss-300,tss+50]
	else:
		print "Strand not recognized"
		continue

	# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
	sgs = find_guides( target.seq, constant=weissman_constant, start=search_start, end=search_end )
	print "Found %s potential guides" % len(sgs)
	if len(sgs) == 0:
		out_fhandle.write( "\t".join( [id, 'x'*20, 'x'*3, str(1e6), "\n"] ))
	try:
		guides[id] = sgs
	except KeyError:
		guides.setdefault( id, sgs )

for id, sgs in guides.iteritems():
	# score potential sgRNAs		
	find_offtargets( sgs, genelist="refgene", noncoding=True, mode="searching", bowtie_genome="GRCh38" )
	for sg in sgs:
		sg.calculate_score()
	sgs.sort(key=lambda x: x.score )

	for sg in sgs:		
		score = sg.score
		out_fhandle.write( "\t".join( [id, str(sg.protospacer.back_transcribe()), str(sg.pam), str(score), "\n"] ))

#	outhandle2 = target.id+"_sg.gb"
#	SeqIO.write( target, outhandle2, "gb" )
out_fhandle.close()
