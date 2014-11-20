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

# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
sense_re=re.compile(r'(?=([ATGCatgc]{20})([ATGCatgc]GG))') #  regex with lookahead to get overlapping sequences. group1 is protospacer, group2 is pam
antisense_re=re.compile(r'(?=(CC[ATGCatgc])([ATGCatgc]{20}))') # group1 is pam, group2 is protospacer

for target in targets:

	# parse headers
	target_id, target_basestart, target_baseend, target_seq, target_mutation = parse_target_header( target )
	print "Targeting bases %i-%i:%s->%s in %s" % (target_basestart, target_baseend, target_seq, target_mutation, target_id)
	if (target_basestart < 0 or target_basestart >= len(target.seq)):
		print "Target site not valid: %d" % target_basestart
		sys.exit(0)

	# set search window
	site_offset = starting_site_offset
	passing_seqs = []
	sgs = []
	search_start = target_basestart - site_offset
	search_end = target_basestart + site_offset
	if (search_start < 1 and search_end >= len(target.seq)):
		print "Target site not valid: %d" % target_basestart

	# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
	sgs = find_guides( target.seq, start=search_start, end=search_end )
	print "Found %s potential guides" % len(sgs)

	# score potential sgRNAs		
	find_offtargets( sgs, genelist="refgene", noncoding=True )
	for sg in sgs:
		sg.calculate_score()
	sgs.sort(key=lambda x: x.score )
		
	# make HR templates and try to remove PAMs for each potential sgRNA. only keep a guide if the PAM could be removed
	hr_sg_list = []
	hr = HrTemplate( target.seq[target_basestart-100:target_basestart-1] + target_mutation + target.seq[target_baseend:target_baseend+100] ) # split sgrna evenly, return 90 bases on either side (200 bases total)
	hr.target_seq = target.seq[target_basestart-100:target_baseend+100]
	hr.find_frame()
	for sg in sgs:
		temp_hr = copy.deepcopy(hr) # need to make a deepcopy so that only the pam of the current sgrna is removed
		if not hr.remove_pam( sg ): # used for all-at-once PAM removal
#		if not temp_hr.remove_pam( sg ): # used for one-at-a-time PAM removal
			print "Could not remove pam for %s" % sg.protospacer
			continue
		else:
			# used for all-at-once PAM removal
			print "Found identical translations for guide %s in frame %s" % (sg.protospacer, hr.frame)
			hr_sg_list.append( (hr, sg) )
			# used for one-at-a-time PAM removal
#			print "Found identical translations for guide %s in frame %s" % (sg.protospacer, temp_hr.frame)
#			hr_sg_list.append( (temp_hr, sg) )

	# output
	i = 1
	for hr_sg in hr_sg_list:
		hr = hr_sg[0]
		sg = hr_sg[1]
		feature = hr_to_seqfeature( target, hr )
		target.features.append( feature )
		score = sg.score
		feature = sg_to_seqfeature( target, sg )
		target.features.append( feature )
		out_fhandle.write( "\t".join( [target_id+"_"+str(i), str(sg.protospacer.back_transcribe()), str(sg.build_fullseq()), str(score), str(hr.seq), "\n"] ))
		i+=1

	feature = SeqFeature( FeatureLocation( target_basestart, target_baseend), strand=None, type="mutation", id=target_mutation, qualifiers={"mutation": target_mutation} )
	target.features.append( feature )
	outhandle2 = target_id.split(",")[0]+"_sg.gb"
	target.name = target_id.split(",")[0]
	SeqIO.write( target, outhandle2, "gb" )
out_fhandle.close()
