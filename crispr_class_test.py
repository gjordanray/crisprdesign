#!/usr/bin/env python
from crisprdesign import SgRna, HrTemplate, GenomicLocation
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna, generic_rna


sg = SgRna( "GGCATGAGTTGCACAAGAGT" )
sg.target_seq = Seq("CGCTCAGCAAGGCATGAGTTGCACAAGAGTTGGGTATCAGC", generic_dna )
sg.target_site = GenomicLocation( "chr3", 52442063, 52442083, "+" )
sg.find_offtargets()

sg = SgRna("ATGCGTGGTGCAAATGTTGC")
sg.target_site = GenomicLocation( "chrX", 101010, 101030, "+" )
sg.target_seq = Seq("CTATAGCCCTATGCGTGGTGCAAATGTTGCAGGCTAGGAC", generic_dna)
#sg.find_offtargets()
print "protospacer %s" % sg.protospacer
print "full sequence %s" % sg.build_fullseq()
print "model pam %s" % sg.model_pam
print "score %f" % sg.score()


hr = HrTemplate( "GGTGATCGTAGGCTTGATCGGATCGGATAGCATGCCGATCGATGTGCTATAGCCCTATGCGTGGTGCAAATGTTGCAGGCTAGGACAGGGCTCTTCTCGGATATAGCTTAGCTTAGGGATAACGGCTAAA")
hr.frame = 1
print "Frame is %s" % hr.frame
hr.find_frame()
print "Frame is %s" % hr.frame
hr.target_seq = hr.seq
#hr.set_seq( "GGTGATCGTAGGCTTGATCGGATCGGATAGCATGCCGATCGATGTGCTATAGCCCTATGCGTGGTGCAAATGTTGCAGACTAGGACAGGGCTCTTCTCGGATATAGCTTAGCTTAGGGATAACGGCTAAA")
print "template sequence %s" % hr.seq
print "target sequence %s" %  hr.target_seq
if hr.remove_pam(sg):
	print "Yay!"
else:
	print "Boo!"

# now do the reverse complement
sg = SgRna("GCAACATTTGCACCACGCAT")
sg.target_site = GenomicLocation( "chrX", 101010, 101030, "-" )
sg.target_seq = Seq("GTGCTATAGCCCGATGCGTGGTGCAAATGTTGCCTAGGACCAT", generic_dna)
#sg.find_offtargets()
print "protospacer %s" % sg.protospacer
print "full sequence %s" % sg.build_fullseq()
print "model pam %s" % sg.model_pam
print "score %f" % sg.score()


hr = HrTemplate( "GGTGATCGTAGGCTTGATCGGATCGGATAGCATGCCGATCGATGTGCTATAGCCCGATGCGTGGTGCAAATGTTGCCTAGGACAGGGCTCTTCTCGGATATAGCTTAGCTTAGGGATAACGGCTAAA")
hr.frame = -3
print "Frame is %s" % hr.frame
hr.find_frame()
print "Frame is %s" % hr.frame
hr.target_seq = hr.seq
print "template sequence %s" % hr.seq
print "target sequence %s" %  hr.target_seq
if hr.remove_pam(sg):
	print "Yay!"
else:
	print "Boo!"

