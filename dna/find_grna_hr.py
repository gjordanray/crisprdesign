#!/usr/bin/env python
import re
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio.Alphabet import generic_dna, generic_rna
from Bio.Data.CodonTable import unambiguous_dna_by_id
import sys, subprocess, os, tempfile

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

target_fname = sys.argv[1]
target_handle = open( target_fname, "rU" )
targets = list( rec.upper() for rec in SeqIO.parse( target_handle, "fasta", alphabet=generic_dna))
out_fname = sys.argv[2]

starting_site_offset = 50 # starting distance around target site for search
max_site_offset = 50 # max distance around target site for search
min_unique_sites = 2 # minimum number of unique sequences to find around the target site

out_fhandle = open( out_fname, mode="w" )
for target in targets:
	target_id, target_basestart, target_baseend, target_seq, target_mutation = parse_target_header( target )
	print "Targeting bases %i-%i:%s->%s in %s" % (target_basestart, target_baseend, target_seq, target_mutation, target_id)
	if (target_basestart < 0 or target_basestart >= len(target.seq)):
		print "Target site not valid: %d" % target_basestart
		sys.exit(0)

	# Find all overlapping 23mers that either end in GG or start with CC
	sense_re=re.compile(r'(?=(([ATGCatgc]{20})[ATGCatgc]GG))') # complicated regex with lookahead to get overlapping sequences.
	antisense_re=re.compile(r'(?=(CC[ATGCatgc]([ATGCatgc]{20})))')

	site_offset = starting_site_offset
	passing_seqs = []
	sgs = []
#	while len(sgs) < min_unique_sites or site_offset <= max_site_offset:
#		while ( len(sgs) < min_unique_sites ): # keep expanding search until at least min # sgRNAs found (don't care if unique at this stage)
	search_start = target_basestart - site_offset
	search_end = target_basestart + site_offset
	if (search_start < 1 and search_end >= len(target.seq)):
		print "Target site not valid: %d" % target_basestart

	sgs = []
	sense = [SgRna(m.group(2)) for m in sense_re.finditer(str(target.seq), search_start, search_end)] # find all PAM-containing sequences, but only build sgs from the protospacer itself
	for sg in sense:
		index = target.seq.find( str(sg.protospacer.back_transcribe()) )
		sg.target_seq = target.seq[index-10:index+len(sg.protospacer)+10] # capture 10 bases on each side
		sg.pam = target.seq[index+len(sg.protospacer):index+len(sg.protospacer)+3]
		sgs.append( sg )
	#print "--------"
	antisense = [ Seq(m.group(2), generic_dna) for m in antisense_re.finditer( str(target.seq), search_start, search_end ) ] # find all PAM-containing sequences on the other strand
	for seq in antisense:
		index = target.seq.find( seq )
		sg = SgRna( str(seq.reverse_complement()) )
		sg.target_seq = target.seq[index-10:index+len(sg.protospacer)+10] # capture 10 bases on each side
		sg.pam = target.seq[index-3:index].reverse_complement()
		sgs.append( sg )
#		print "          %s" % sg.protospacer
#		print "          %s" % sg.protospacer.reverse_complement()
#		print sg.target_seq
#		print sg.pam
#		print
	find_offtargets( sgs )
	
	scored_sgs = [ (sg, sg.calculate_score()) for sg in sgs ]
		
	hr_sg_list = []
	hr = HrTemplate( target.seq[target_basestart-100:target_basestart-1] + target_mutation + target.seq[target_baseend:target_baseend+100] ) # split sgrna evenly, return 90 bases on either side (200 bases total)
	hr.target_seq = target.seq[target_basestart-100:target_baseend+100]
	hr.find_frame()
	hr_sg_list.append( hr )
	for scored_sg in scored_sgs:
		sg = scored_sg[0]
		if not hr.remove_pam( sg ):
			print "Could not remove pam for %s" % sg.protospacer
			continue
		else:
			print "Found identical translations for guide %s in frame %s" % (sg.protospacer, hr.frame)
			hr_sg_list.append( scored_sg )

	min_score = 1000
#	hr_sg_list.sort( key=lambda x:x[1] )
	hr = hr_sg_list[0]
	feature = hr_to_seqfeature( target, hr )
	target.features.append( feature )

	for scored_sg in hr_sg_list[1:]:
		sg = scored_sg[0]
		score = scored_sg[1]
		out_fhandle.write( "\t".join( [target_id, str(sg.protospacer.back_transcribe()), str(score), str(hr.seq), "\n"] ))
		feature = sg_to_seqfeature( target, sg )
		target.features.append( feature )

	feature = SeqFeature( FeatureLocation( target_basestart, target_baseend), strand=None, type="mutation", id=target_mutation, qualifiers={"mutation": target_mutation} )
	target.features.append( feature )
	outhandle2 = target_id.split(",")[0]+"_sg.gb"
	target.name = target_id.split(",")[0]
	SeqIO.write( target, outhandle2, "gb" )
out_fhandle.close()
