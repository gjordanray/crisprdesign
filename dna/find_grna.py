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

def parse_target_header( sequence_record ):
#reads fasta headers of the format
#>id, start base,premutation sequence,end base,postmutation sequence
#for example
#>gene foo bar,100,ATG,103,TAA
	id, start, seq, end, mutation = sequence_record.id.split(",")
	start = int(start)
	end = int(end)
	assert( end-start+1 ) == len(seq)
	return id, start, end, seq, mutation

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

def hr_to_seqfeature( target, hr ):
	location = target.seq.find( hr.target_seq )
	if location == -1:
		print "Could not find hr template in target!"
		return None
	strand = None
	id = hr.seq
	quals = {"sequence": hr.seq }
	feature = SeqFeature( FeatureLocation( location, location+len(hr.seq)), strand=strand, type="hrtemplate", id=id, qualifiers=quals)
	return feature

def plot_sgrnas( target, target_basen, offset, sites ):
	gd_diagram = GenomeDiagram.Diagram( )
	gd_track_for_features = gd_diagram.new_track(1, name="sgRNAs")
	gd_feature_set = gd_track_for_features.new_set()
	for site in sites:
		record = site[0].split('_')
		target_sequence = record[0]
		target_id = record[1]
		strand = record[2]
		sgrna = ''
		if strand == 'sense':
			sgrna = Seq(target_sequence[:-2], generic_dna)
		elif strand== 'as':
			sgrna = Seq(target_sequence[2:], generic_dna)
			sgrna = sgrna.reverse_complement()
		else:
			print "Strand not recognized!"
		start = target.seq.find( target_sequence )
		if strand == 'sense':
			feature = SeqFeature(FeatureLocation(start,start+23),strand=+1)
			gd_feature_set.add_feature(feature, color="cyan", label=True, name=str(sgrna), label_angle=0, sigil="ARROW")
		else:
			feature = SeqFeature(FeatureLocation(start,start+23),strand=-1)
			gd_feature_set.add_feature(feature, color="blue", label=True, name=str(sgrna), label_angle=0, sigil="ARROW")
		print target_sequence, target_id, strand, sgrna

	feature = SeqFeature(FeatureLocation(target_basestart, target_basestart+1))
	gd_feature_set.add_feature(feature, color="red", label=False)
	
	feature = SeqFeature(FeatureLocation(0, len(target.seq)), strand=+1)
	gd_feature_set.add_feature(feature, color="blue", label=True, name=str(target), label_angle=0, hide=False, label_size=25)
	
	gd_diagram.draw(format="linear", pagesize='LETTER', fragments=1, start=target_basen-offset, end=target_basen+offset)
	gd_diagram.write(target_id+".pdf", "PDF")


# start of main script
target_fname = sys.argv[1]
target_handle = open( target_fname, "rU" )
targets = list( rec.upper() for rec in SeqIO.parse( target_handle, "fasta", alphabet=generic_dna))
out_fname = sys.argv[2]

starting_site_offset = 50 # starting distance around target site for search
max_site_offset = 50 # max distance on each side of target site for search (symmetric)
min_unique_sites = 2 # minimum number of unique sequences to find around the target site

out_fhandle = open( out_fname, mode="w" )
for target in targets:

	# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
	sense_re=re.compile(r'(?=(([ATGCatgc]{20})[ATGCatgc]GG))') # complicated regex with lookahead to get overlapping sequences.
	antisense_re=re.compile(r'(?=(CC[ATGCatgc]([ATGCatgc]{20})))')
	sgs = []
	sense = [SgRna(m.group(2)) for m in sense_re.finditer(str(target.seq))] # find all PAM-containing sequences, but only build sgs from the protospacer itself
	for sg in sense:
		index = target.seq.find( str(sg.protospacer.back_transcribe()) )
		sg.target_seq = target.seq[index-10:index+len(sg.protospacer)+10] # capture 10 bases on each side
		sg.pam = target.seq[index+len(sg.protospacer):index+len(sg.protospacer)+3]
		sgs.append( sg )
	antisense = [ Seq(m.group(2), generic_dna) for m in antisense_re.finditer( str(target.seq)) ] # find all PAM-containing sequences on the other strand
	for seq in antisense:
		index = target.seq.find( seq )
		sg = SgRna( str(seq.reverse_complement()) )
		sg.target_seq = target.seq[index-10:index+len(sg.protospacer)+10] # capture 10 bases on each side
		sg.pam = target.seq[index-3:index].reverse_complement()
		sgs.append( sg )
	print "Found %s potential guides" % len(sgs)

	# score potential sgRNAs		
	find_offtargets( sgs, genelist="refgene", noncoding=True, mode="searching" )
	for sg in sgs:
		sg.calculate_score()
	sgs.sort(key=lambda x: x.score )

	for sg in sgs:		
		score = sg.score
		if score < 100:
			feature = sg_to_seqfeature( target, sg )
			target.features.append( feature )
		out_fhandle.write( "\t".join( [target.id, str(sg.protospacer.back_transcribe()), str(score), "\n"] ))

	outhandle2 = target.id+"_sg.gb"
	SeqIO.write( target, outhandle2, "gb" )
out_fhandle.close()
