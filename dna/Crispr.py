#!/usr/bin/env python
import RNA # Vienna RNA secondary structure prediction
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna, generic_rna
import tempfile
import subprocess
import re
import random
import sys

def hamming_dist( str1, str2 ):
	diffs=0
	for ch1, ch2 in zip( str1, str2 ):
		if ch1 != ch2:
			diffs += 1			
	return diffs

def GCcontent( seq ):
	GC = (seq.count("G") + seq.count("C")) / float(len(seq))
	return GC

def bowtie_search( sgrna_list ): # returns dictionary of {protospacer+pam: GenomicLocations}
#	for sg in sgrna_list:
#		if sg.target_seq == None:
#			print "Must set target_seq to find offtargets!"
#			return []

	# antisense_phredString = 'I4!=======44444++++++++'
	bowtie_constant_options = ['bowtie', '--nomaqround', '--best', '-n 3', '-l 12', '-a', '-e 39', '-p 4', '--suppress', '1,6,7', '--chunkmbs', '256', 'hg19'] #note '--chunkmbs 128' and '--suppress 5,6,7'is one option, but subprocess needs them separated
	temp_bowtiein = tempfile.NamedTemporaryFile(delete=True)
#	temp_bowtiein = open( "bowtie.in", mode="w")

	i=1
	for sg in sgrna_list:
		phred33String = '++++++++44444=======!4I'
		protospacerpam = sg.build_protospacerpam()
		if sg.pam == "": # missing PAM
			phred33String = phred33String[:-3]
		if len(phred33String) > len(protospacerpam):
			while len( phred33String ) > len( protospacerpam ):
				phred33String = phred33String[1:] # try removing from the 5' end
			if (len( phred33String ) < 18) or (len( protospacerpam) != len( phred33String )): # protospacer shorter than a 18mer usually doesn't make sense
				print "Can't figure out phred33String to use! Check manually."
				print "%s %s %s" % (phred33String, protospacerpam, sg.protospacer)
				return None
		elif len(phred33String) < len(protospacerpam): # prepend a 0-penalty character
			while len( phred33String ) < len( protospacerpam ):
				phred33String = "+"+phred33String
			if (len( phred33String ) > 40) or (len( protospacerpam) != len( phred33String )): # protospacer longer than a 40mer usually doesn't make sense
				print "Can't figure out phred33String to use! Check manually."
				print "%s %s %s" % (phred33String, protospacerpam, sg.protospacer)
				return None
		assert( len( protospacerpam) == len( phred33String ) )
		temp_bowtiein.write('@'+str(i)+'\n')
		#print sg.build_protospacerpam()
		temp_bowtiein.write( str(protospacerpam)+'\n')			
		temp_bowtiein.write('+\n')
		temp_bowtiein.write(phred33String+'\n')
	temp_bowtiein.flush() #necessary for bowtie to read files

	temp_bowtieout = tempfile.NamedTemporaryFile(delete=True)
	#temp_bowtieout = "bowtie.out"
	bowtie_cmdline = list(bowtie_constant_options)
	bowtie_cmdline.append(temp_bowtiein.name)
	bowtie_cmdline.append(temp_bowtieout.name)
	subprocess.call( bowtie_cmdline ) #, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb') )
	temp_bowtiein.close() # remove tempfiles

	print "Parsing bowtie output",
	sys.stdout.flush()
	found_locations = {}
	i = 0
	for line in temp_bowtieout:
		if i % 10000 == 0:
			print ".",
			sys.stdout.flush()
		l = line.strip().split('\t')
		s = Seq( l[3], generic_dna )
		loc = GenomicLocation( l[1], long(l[2]), long(l[2])+len(s), l[0] ) # format of bowtie output is strand, chr, start, sequence of read (revcomp if - strand mapped)
		if loc.strand == "-":
			s = s.reverse_complement()
		try:
			found_locations[str(s)].append( loc )
		except KeyError:
			found_locations.setdefault( str(s), [loc] ) # need to convert back to string for proper key referencing
		i+=1
	temp_bowtieout.close()
	#print found_locations
	print "Done!"
	return found_locations

def _read_ccds( fname ):
	fhandle = open( fname, "rU" )
	genes = []
	for line in fhandle:
		split = line.split('\t')
		if split[0] == "#chromosome":
			continue
		# format is    chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from(0-index)	cds_to(0-index)	cds_locations	match_type
		if split[5].split(",")[0] == "Withdrawn":
			continue
		if split[10].rstrip() == "Partial":
			continue
#		print split[2], split[0], split[6], long(split[7])+1, long(split[8])+1
		genes.append( [split[2], split[0], split[6], long(split[7])+1, long(split[8])+1] ) # genes container now holds: geneid, chr, strand, txStart, txEnd (tx sites +1 to conform to UCSC 1-based index)
	fhandle.close()
	return genes

def _read_refgene( fname ):
	fhandle = open( fname, "rU" )
	genes=[]
	for line in fhandle:
		split = line.split('\t')
		# format is    bin name chr strand txStart txEnd cdsStart cdsEnd
		genes.append( [split[1], split[2], split[3], long(split[4]), long(split[5])] ) # genes container now holds: geneid, chr, strand, txStart, txend
	fhandle.close()
	return genes

def find_offtargets( sgrna_list, genelist="refgene" ):
	# run bowtie for protospacer, including PAM
	genomic_sites = bowtie_search( sgrna_list )
	# offtargets are everything that's not the target_site
	print "Parsing reference gene list...",
	if genelist=="refgene":
		genes = _read_refgene( "refGene.hg19.clean.sorted.txt" )
	elif genelist=="ccds":
		genes = _read_ccds( "CCDS.20140807.txt" )
	print "Read %s genes." % len(genes)
	for sg in sgrna_list:
		protospacerpam = str( sg.build_protospacerpam() )
		try:
			offsites = genomic_sites[protospacerpam]
		except KeyError:
			print "Could not find any genomic sites for %s" % protospacerpam
			continue # we didn't find this (only happens if you called a protospacer that doesn't have a genomic match)
		print "Off targets and associated genes for %s:" % protospacerpam
		if not sg.target_site:
			sg.target_site = offsites[0]
			#print "target %s" % sg.target_site
		offsites.remove( sg.target_site ) # first one found should be the target
		gene_dict = {}
		for site in offsites:
			print site,
			gene_dict.setdefault(site,[])
			# read refGene, find chromosome, is this index in a gene?
			for gene in genes:
				rg_gene, rg_chr, rg_start, rg_end = gene[0], gene[1], gene[3], gene[4]
				if rg_chr == site.chr:
					if (rg_start <= site.start <= rg_end) or (rg_start <= site.end <= rg_end):
						try:
							gene_dict[site].append( rg_gene )
							print rg_gene,
						except KeyError:
							continue
							#print site, rg_gene
			print
		sg.offtarget_sites = gene_dict
		print
#		for k,v in gene_dict.items():
#			print k, v
#		print

class GenomicLocation:
	def __init__(self, chr, start, end, strand):
		#instance variables
		self.chr = chr
		self.start = start # start and end are 1-indexed (as per UCSC)
		self.end = end
		self.strand = strand
	def __eq__( self, other):
		return( (self.chr,self.start,self.end,self.strand) == (other.chr,other.start,other.end,other.strand) )
	def __ne__( self, other):
		return not self == other
	def __hash__( self ):
		s = str(self.chr)+str(self.start)+str(self.end)+str(self.strand)
		hash = int(''.join(str(ord(c)) for c in s))
		return hash
	def __str__( self ):
		return "%s %s:%s%s" % (self.chr, self.start, self.end, self.strand)

class HrTemplate:
	def __init__(self, seq):
		# instance variables
		self.seq = Seq(str(seq), generic_dna) # sequence of the hr template itself
		self.target_site = None # GenomicLocation object
		self.target_seq = None # genomic sequence of the region this template targets (not necessarily identical to seq of template)
		self.frame = 1 # can be +1, +2, +3, -1, -2, -3, or 0 (non-coding)
	def __eq__(self, other):
		return (( self.seq, self.target_location, self.target_seq, self.frame ) == (other.seq, other.target_location, other.target_seq, other.frame ))
	def __ne__(self, other):
		return not self == other

	def find_frame(self): #assume that protein seq with fewest number of stops is correct
		revcomp = self.seq.reverse_complement()
		forward_stops = []
		reverse_stops = []
		min_stops = 1000
		for i in range(3):
			#print self.seq[i:].translate()
			#print revcomp[i:].translate()
			#print
			forward_stops = self.seq[i:].translate().count("*")
			reverse_stops = revcomp[i:].translate().count("*")
			if forward_stops < min_stops:
				self.frame = i+1
				min_stops = forward_stops
			if reverse_stops < min_stops:
				self.frame = (-1*i)-1
				min_stops = reverse_stops
				
	def remove_pam( self, sgrna ):
		# find location based on genomic sequence we're targeting
		# then perform mutation based on existing HR template
		# this allows mutations in HR template, serial removal of multiple PAMs, etc.
		possible_bases3 = ["A", "T", "C"] # 3rd base of pam
		possible_bases2 = ["T", "C" ] # 2nd base of pam
		protospacer = sgrna.protospacer.back_transcribe()
		sgrna_template_index = self.target_seq.find( protospacer )
		if sgrna_template_index != -1:
			strand = "+"
		else:
			strand = "-"
			sgrna_template_index = self.target_seq.find(protospacer.reverse_complement())
		if sgrna_template_index == -1:
			print "Could not find protospacer in HR target sequence!"
			# print self.target_seq, protospacer.reverse_complement()
			return False
		# print "Strand %s" % strand
		if strand == "-":
			pos3, pos2 = range(sgrna_template_index-3, sgrna_template_index-1) # only get the last two bases of the PAM
		else:
			pos2,pos3 = range(sgrna_template_index + len(protospacer)+1, sgrna_template_index + len(protospacer)+3) # only get the last 2 bases of the PAM
		# from here on, whether + or - strand, pam_indices order is [pos2,pos3]
		# print pam_indices
		if self.frame == 0: # if we're in a noncoding region, replace with a randomly allowed base
			base = random.choice(possible_bases3)
			seq_copy = self.seq[:pos3] + base + self.seq[pos3+1:]
			self.seq = seq_copy
			return True
		else:	# if we're in a coding region, possible_bases is shuffled and we iterate
			random.shuffle( possible_bases3 )
			random.shuffle( possible_bases3 )
			#print i, self.seq[i]
			possible_bases2.insert(0, self.seq[ pos2 ] ) # add existing 2nd PAM base to the list. that way we start with identity at pos2 and vary pos3, then move on to varying both
			for base2 in possible_bases2:
				for base3 in possible_bases3:
					seq_copy = self.seq[:pos2] + base2 + base3 + self.seq[pos3+1:]
					if self.frame < 1: 
						test_translation = seq_copy.reverse_complement()[abs(self.frame)-1:].translate()
						ref_translation = self.seq.reverse_complement()[abs(self.frame)-1:].translate()
					else:
						test_translation = seq_copy[self.frame-1:].translate()
						ref_translation = self.seq[self.frame-1:].translate()
					#print i
					#print test_translation
					#print ref_translation
					#print
					if str(test_translation) == str(ref_translation):
						print "Found identical translations for frame %s" % self.frame
						#print self.seq
						#print seq_copy
						self.seq = seq_copy
						return True
		#print "Could not remove PAM!"
		return False			

class SgRna:
	# instance variables
	def __init__(self, seq, constant_region="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU" ):
		# turn DNA input into RNA
		seq_copy = seq
		if seq_copy.find("T"): # convert to RNA if not already done
			seq_copy = seq.replace( "T", "U" )
		#print seq, seq_copy
		self.protospacer = Seq(seq_copy, generic_rna) # sequence sans constant portion. can only set protospacer on initialization. always stored as RNA
		self.target_site = None # will eventually become a GenomicLocation. strand (+ means sgrna seq same as + strand, - means sgrna seq same as - strand), 
		self.target_seq = "" # 10 bases on either side of site
		self.offtarget_sites = {} # dict, format = {GenomicLocation: [gene1, gene2, gene3...]}
		self.pam = ""
		self.constant_region = Seq( constant_region, generic_rna )
	def __eq__( self, other ):
		return( (self.protospacer, self.target_site, self.target_seq, self.offtarget_sites, self.constant_region) == (other.protospacer, other.target_site, other.target_seq, other.offtarget_sites, other.constant_region))
	def __ne__( self, other ):
		return not self == other
			
	def build_fullseq( self ):
		fullseq = self.protospacer + self.constant_region
		if fullseq[0] != "G":
			fullseq = "G" + fullseq
		return fullseq
	
	def build_protospacerpam( self ): # returns DNA
		p = self.protospacer.back_transcribe()
		if self.pam != "":
			return Seq( str(p)+str(self.pam), generic_dna )
		if self.target_seq == "":
			print "Can't build PAM without sequence context. Defaulting to protospacer %s" % self.protospacer.back_transcribe()
			return self.protospacer.back_transcribe()
		else:
			temp_target = self.target_seq
			index = temp_target.find( p )
			if index == -1:
				temp_target = temp_target.reverse_complement()
				index = temp_target.find( p )
			if index == -1:
				print "Can't find protospacer in target sequence"
				# print p, temp_target
				return ""
			if self.pam != "":
				pam = self.pam
			else:
				pam = temp_target[ index+len(p):index+len(p)+len(self.pam) ]
				self.pam = pam
			return Seq( str(p)+str(self.pam), generic_dna )
	
	def score(self):
		print "Scoring sgrna %s" % self.protospacer
		score = 0
		# lower score is better
		# +25 if GC content
		# +100 for each U homopolymer (3+) (sliding window)
		# +10 for homoX sequences other than U
		# +25 if 3' base of protospacer is a U
		# +X^2, where X is stability of secondary structure in the protospacer 
		#	(squaring gives reasonably bigger penalties for stability and takes care of negatives)
		# (below are optional, if genomic context and offtarget sites are defined)
		# +10 for each 3' G in a row after the PAM
		# +25 for each offtarget site that passed bowtie filters
		# 	or +100 for each gene (refGene) that offtarget site was in
		# +10 if last base before PAM is C (Doench et al Nat Biotech 2014)
		# -5 if last base before PAM is G (Doench et al Nat Biotech 2014)
		# +20 if no microhomology found TODO
		# ----
		# penalize secondary structure in the protospacer sequence
		secstruct, stability = RNA.fold( str(self.protospacer) )
		print "Protospacer secstruct & stability %s %f kcal/mol" % (secstruct, stability)
		score += stability**2
		# mildly penalize bad secondary structure introduced to the constant region
		full_secstruct, full_stability = RNA.fold( str(self.build_fullseq()))
		constant_secstruct, constant_stability = RNA.fold( str( self.constant_region ))
		if full_secstruct[-len(constant_secstruct):] != constant_secstruct:
			print "constant  secstruct %s" % constant_secstruct
			print "+guide    secstruct %s" % full_secstruct[-len(constant_secstruct):]
			score += hamming_dist( full_secstruct[-len(constant_secstruct):], constant_secstruct )
		# penalize for sequences adjacent PAM
		if self.protospacer[:-1] == "C":
			score += 10
		elif self.protospacer[:-1] == "G":
			score -= 5
		# penalize homoX (3+), especially homoU (pol III terminator)
		matches = re.finditer(r'(?=((\w)\2{2}))', str(self.protospacer)) #
		for match in matches:
			print "Found %s  " % match.group(1)
			if match.group(1) == "UUU" or match.group(1) == "TTT":
				score += 100
			else:
				score += 10
		# remove bad GC sequences
		GC = GCcontent(self.protospacer)
		if (GC < 0.4 or GC > 0.8):
			print "GC content is %f" % GC
			score += 10
		# penalize 3' Us
		for i in [-1, -2, -4]:
			if self.protospacer[:i] == "U" or self.protospacer[:i] == "T": #see Wang et al for sgRNA loading data
				print "3' U found",
				score += 5
		# penalize offtarget sites in genome
		if len(self.offtarget_sites) > 0:
			print "Found %s offtarget sites" % len( self.offtarget_sites.keys() )
			for offsite in self.offtarget_sites.keys():
				num_offtarget_genes = len(self.offtarget_sites[offsite])
				if num_offtarget_genes > 0:
					print "Found offtarget sites in the following genes: ",
					for offtgt in self.offtarget_sites[offsite]:
						print "%s " % offtgt,
					print
					score += num_offtarget_genes * 1000
				else:
					score += 25
		# penalize based on facts  about the site targeted in the genome
		if self.target_seq:
			# PAM profile
			p = self.protospacer.back_transcribe()
			if self.target_seq.find( p ) != -1:
				downstream_seq = self.target_seq
			elif self.target_seq.reverse_complement().find( p ) != -1:
				downstream_seq = self.target_seq.reverse_complement()
			else:
				print "Could not find protospacer exact match in target site!"
				return score
			grna_local_index = downstream_seq.find( p )
			downstream_seq = downstream_seq[grna_local_index + len(self.protospacer):] # everything after the sgrna
#			pam = downstream_seq[:3]
#			print "model PAM %s | PAM %s" % (self.model_pam, self.pam)
#			if pam[1] != self.model_pam[1]:
#				score += 20
#			if pam[2] != self.model_pam[2]:
#				score += 10000 # something is wrong... last base should always match
			# penalize G homopolymers after the PAM
			i=0
			for base in downstream_seq[len(self.pam):]:
				if base == "G":
					i+=1
					score += 10
				else:
					break
			if i>0: print "Found %i Gs following PAM" % i
		print "Total score %f" % score
		print
		return score