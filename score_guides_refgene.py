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

guides_fname = sys.argv[1]
out_fname = sys.argv[2]
guides_handle = open( guides_fname, mode="rU" )
guides = []
targets = []
for line in guides_handle:
	l = line.split(" ")
	target = l[0].rstrip()
	guide = l[1].rstrip()
	sg = SgRna( guide, constant_region="GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU" )
				  # weissman constant = "GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU"
				  # broad constant = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU"
	sg.pam = l[2].rstrip()
	guides.append( sg )
	targets.append( target )
find_offtargets( guides, genelist="refgene", noncoding=False)

scores = []
for sg in guides:
	score = sg.score()
	scores.append( score )

out_fhandle = open( out_fname, mode="w" )
for i in range(len(scores)):
	target = targets[i]
	score = scores[i]
	guide = guides[i].build_protospacerpam()
	out_fhandle.write( str(target) + "\t" + str(guide) + "\t" + str(score) + "\n")
#	out_fhandle.write( str(score) + "\n")
	print target, score
#	print score
out_fhandle.close()
