#!/usr/bin/env python
try:
	import RNA # Vienna RNA secondary structure prediction
	vienna_loaded = True
except ImportError:
	vienna_loaded = False
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna, generic_rna
from Bio.SeqFeature import SeqFeature, FeatureLocation
import tempfile
import subprocess
import re
import random
import sys
import multiprocessing

# requires: ViennaRNA python, bowtie with hg38 genome and *NO* mitochondrial sequences ("hg38_noM"), cleaned refseq file
# Author: Jacob Corn. jcorn@berkeley.edu

def find_guides( in_seq, sense=True, antisense=True, constant="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU", start=0, end=0 ):
	"""Find potential sgRNAs within a sequence, using a sliding window.
	Args:
		in_seq(Bio.Seq): input sequence
		sense(bool): search sense strand
		antisense(bool): search antisense strand
		constant(string): sequence of constant region to set for guides
		start(int): index in the sequence to start finding guides
		end(int): index in the sequence to end finding guides
	"""

	# default constant region is Broad-style for cutting
	# Find potential sgRNAs (defined as 23-mers ending in NGG or starting in CCN) on both plus and minus strands
	if start == end: # find guides across the entire sequence by default
		end = len( in_seq )
	sgs = []
	if sense:
		sense_re=re.compile(r'(?=([ATGCatgc]{20})([ATGCatgc]GG))') #  regex with lookahead to get overlapping sequences. group1 is protospacer, group2 is pam
		for m in sense_re.finditer(str(in_seq), start, end):
			protospacer = m.group(1)
			foundpam = m.group(2)
			sgstart=m.start(1)
			sgend=m.end(1)
			#print "%s %s forward" % (seq, seqpam)
			sg = SgRna(protospacer, constant_region=constant, target_seq=in_seq[sgstart-10:sgend+10], pam=foundpam)
			sgs.append( sg )
	if antisense:
		antisense_re=re.compile(r'(?=(CC[ATGCatgc])([ATGCatgc]{20}))') # group1 is pam, group2 is protospacer
		for m in antisense_re.finditer( str(in_seq), start, end):
			protospacer = str(Seq(m.group(2), generic_dna).reverse_complement())
			foundpam=Seq(m.group(1), generic_dna).reverse_complement()
			sgstart=m.start(2)
			sgend=m.end(2)
			#print "%s %s reverse" % (seq, seqpam)
			sg = SgRna( protospacer, constant_region=constant, target_seq=in_seq[sgstart-10:sgend+10], pam=foundpam )
			sgs.append( sg )
	return sgs

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

def hamming_dist( str1, str2 ):
	"""Find hamming distance of two sequences
	Args:
		str1(string)
		str2(string)
	"""
	diffs=0
	for ch1, ch2 in zip( str1, str2 ):
		if ch1 != ch2:
			diffs += 1			
	return diffs

def GCcontent( seq ):
	"""Find fractional GC content of sequence
	Args:
		seq(Bio.Seq): input sequence
	"""
	GC = (seq.count("G") + seq.count("C")) / float(len(seq))
	return GC

def bowtie_search( sgrna_list, genome, max_matches ): # returns dictionary of {protospacer+pam: GenomicLocations}
	"""Use bowtie to find where guides map to a genome. Called by find_offtargets()
	Args:
		sgrna_list([SgRna1, SgRna2,...]): guides for which to find offtargets
		max_matches(int): maximum number of times a guide can match the genome before bowtie triages failure. use 0 to find all matches.
		genome(string): bowtie index to use for finding offtargets.
	"""

#	for sg in sgrna_list:
#		if sg.target_seq == None:
#			print "Must set target_seq to find offtargets!"
#			return []

	# use all cpus except for 1
	ncpus = multiprocessing.cpu_count()
	if ncpus > 1:
		ncpus = ncpus - 1
#	phred33String = '++++++++44444=======!4I'
	antisense_phred33String = 'I4!=======44444++++++++' # use antisense pattern so bowtie 5' seed usage corresponds to stringent part of PAM+protospacer
	phred33String = antisense_phred33String
	bowtie_constant_options = ['bowtie', '--nomaqround', '-a', '--best', '-n 3', '-l 12', '-e 39', '-p '+str(ncpus), '--suppress', '1,6,7', '--chunkmbs', '256', genome ] #note '--chunkmbs 128' and '--suppress 5,6,7'is one option, but subprocess needs them separated. -l 12 -n 3 since phred33 scoring also means can have no more than 3 mismatches in 5' 12 bases
	if max_matches > 0:
		bowtie_constant_options.extend( ['-m', str(max_matches)] )
	
	temp_bowtiein = tempfile.NamedTemporaryFile(delete=True)
#	temp_bowtiein = open( "bowtie.in", mode="w")

	i=1
	for sg in sgrna_list:
		protospacerpam = sg.build_protospacerpam()
		if sg.pam == "": # missing PAM
			phred33String = phred33String[:-3]
		if len(phred33String) > len(protospacerpam):
			while len( phred33String ) > len( protospacerpam ):
				# phred33String = phred33String[1:] # try removing from the 5' end
				phred33String = phred33String[:-1] # try removing from the 3' end
			if (len( phred33String ) < 18) or (len( protospacerpam) != len( phred33String )): # protospacer shorter than a 18mer usually doesn't make sense
				print "Can't figure out phred33String to use! Check manually."
				print "%s %s %s" % (phred33String, protospacerpam, sg.protospacer)
				return None
		elif len(phred33String) < len(protospacerpam): # prepend a 0-penalty character
			while len( phred33String ) < len( protospacerpam ):
				#phred33String = "+"+phred33String # add to 5' end of phred33 pattern
				phred33String = phred33String+"+" # add to 3' end of phred33 pattern
			if (len( phred33String ) > 40) or (len( protospacerpam) != len( phred33String )): # protospacer longer than a 40mer usually doesn't make sense
				print "Can't figure out phred33String to use! Check manually."
				print "%s %s %s" % (phred33String, protospacerpam, sg.protospacer)
				return None
		assert( len( protospacerpam) == len( phred33String ) )
		temp_bowtiein.write('@'+str(i)+'\n')
		#print sg.build_protospacerpam()

		# protospacerpam stored as 5'-3'protospacer followed by pam
		# reverse everything so that seed sequence is located in the PAM+3' end of protospacer
		temp_bowtiein.write( str(protospacerpam.reverse_complement())+'\n')			
		temp_bowtiein.write('+\n')
		temp_bowtiein.write(phred33String+'\n')
	temp_bowtiein.flush() #necessary for bowtie to read files

	temp_bowtieout = tempfile.NamedTemporaryFile(delete=True)
	bowtie_cmdline = list(bowtie_constant_options)
	bowtie_cmdline.append(temp_bowtiein.name)
	bowtie_cmdline.append(temp_bowtieout.name)
	print "Running bowtie to find genomic offtarget sites."
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
		# previously revcomped sequence for phred33 pattern. 
		# if found sequence is on + strand, reverse it back now to match protospacer+pam
		# if found sequence is on - strand then already matches
		# all found sequences should now be protospacer+pam, ending in NGG
		if loc.strand == "+":
			s = s.reverse_complement()
		try:
			found_locations[str(s)].append( loc )
		except KeyError:
			found_locations.setdefault( str(s), [loc] ) # need to convert back to string for proper key referencing
		i+=1
	temp_bowtieout.close()
	print "Done!"
	return found_locations

def _read_ccds( fname ):
	"""read ccds-format reference gene list
	Args:
		fname(string): ccds format is    chromosome	nc_accession	gene	gene_id	ccds_id	ccds_status	cds_strand	cds_from(0-index)	cds_to(0-index)	cds_locations	match_type
	"""
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
		exons = split[9]
		exons = exons.lstrip('[').rstrip(']')
		exons = [(long(x.split('-')[0])+1,long(x.split('-')[1])+1) for x in exons.split(',')] # turn into start/end tuples with 1-indexing		
		genes.append( [split[2], split[0], split[6], long(split[7])+1, long(split[8])+1, exons] ) # genes container now holds: geneid, chr, strand, txStart, txEnd (tx sites +1 to conform to UCSC 1-based index)
	fhandle.close()
	return genes

def _read_refgene( fname ):
	"""read refgene-format reference gene list
	Args:
		fname(string): refgene format is    bin name chr strand txStart txEnd cdsStart(multiples) cdsEnd(multiples)
	"""
	fhandle = open( fname, "rU" )
	genes=[]
	for line in fhandle:
		if line[0] == "#":
			continue
		split = line.split('\t')
		# format is    bin name chr strand txStart txEnd cdsStart(multiples) cdsEnd(multiples)
		starts = [long(x) for x in split[6].rstrip(',').split(',')] # need to remove trailing comma
		ends = [long(x) for x in split[7].rstrip(',').split(',')]
		exons = zip(starts,ends)
		genes.append( [split[1], split[2], split[3], long(split[4]), long(split[5]), exons] ) # genes container now holds: geneid, chr, strand, txStart, txend, exons
	fhandle.close()
	return genes

def find_offtargets( sgrna_list, genelist="refgene", noncoding=True, max_matches=10, bowtie_genome="hg38_noM", mode="cut" ):
	"""Finds guideRNA offtargets within a genome
	Args:
		sgrna_list([SgRna1, SgRna2,...]): guides for which to find offtargets
		genelist(string): data file to use for reference list of gene locations. can be "refgene" or "ccds". TODO: change this to use new ensembl data.
		noncoding(bool): use noncoding sequences when considering whether a guide falls within a gene?
		max_matches(int): maximum number of times a guide can match the genome before bowtie triages failure. use 0 to find all matches.
		bowtie_genome(string): bowtie index to use for finding offtargets.
		mode(string): how to determine offtargets. can be "cut" = look within entire transcribed region, "inhibit" = look -50 to +300 of TSS, "activate" = look -400 to -50 of TSS
	"""
	# run bowtie for protospacer, including PAM
	genomic_sites = bowtie_search( sgrna_list, genome=bowtie_genome, max_matches=max_matches )
	# offtargets are everything that's not the target_site
	print "Parsing reference gene list...",
	if genelist=="refgene":
		genes = _read_refgene( "refGene.hg38.clean.sorted.txt" )
	elif genelist=="ccds":
		genes = _read_ccds( "CCDS.20140807.txt" )
	print "Read %s genes." % len(genes)
	print "Finding offtargets including,"
	if noncoding:
		print "both coding and noncoding sequences."
	else:
		print "only coding sequences."
	for sg in sgrna_list:
		protospacerpam = str( sg.build_protospacerpam() )
		print "Off targets and associated genes for %s:" % protospacerpam
		try:
			offsites = genomic_sites[protospacerpam]
		except KeyError:
			print "Could not find any genomic sites for %s" % protospacerpam
			continue # we didn't find this (only happens if you called a protospacer that doesn't have a genomic match)
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
				windows = [] # container used for search windows. does a guide lie within any of these ranges?
				rg_gene, rg_chr, rg_strand, rg_start, rg_end = gene[0], gene[1], gene[2], gene[3], gene[4]
				if rg_chr == site.chr: # guide lies on same chromosome as gene
					if mode == "inhibit": # look -50 to +300 of TSS
						if ( rg_strand ) == "+":
							windows.append( (rg_start-50, rg_start+300) )
						else:
							windows.append( (rg_start+50, rg_start-300) )
					elif mode == "activate": # look -400 to -50 of TSS
						if ( rg_strand ) == "+":
							windows.append( (rg_start-400, rg_start-50) )
						else:
							windows.append( (rg_start+400, rg_start+50) )
					elif mode == "cut": # look through entire transcribed region
						if noncoding:
							windows.append( ( rg_start, rg_end ) )
						else:
							for exon in gene[5]:
								windows.append( (exon[0], exon[1]) )
					for window in windows:
						if (window[0] <= site.start <= window[1]) or (window[0] <= site.end <= window[1]):
							try:
								gene_dict[site].append( rg_gene )
								print rg_gene,
							except KeyError:
								continue
									#print site, rg_gene
									#print site, rg_gene
			print
		sg.offtarget_sites = gene_dict
		print
#		for k,v in gene_dict.items():
#			print k, v
#		print

class GenomicLocation:
	"""Holds information about a location in a genome.
	Args:
		chr(str): chromosome
		start(int): location start (1-indexed)
		end(int): location end (1-indexed)
		strand(str): - or +
	"""
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

class Gene:
	"""Holds information about a gene
	Args:
		id(str): some kind of identification (arbitrary)
		start(int): gene start
		end(int): gene end
		transcripts([Transcripts]): Transcript objects encoded by this gene
	"""
	def __init__(self, id, start, end, transcripts):
		self.id = id
		self.start = start
		self.end = end
		self.transcripts = transcripts
	def __str__( self ):
		return "Gene %s %s %s %s" % (self.id, self.start, self.end, self.transcripts)
	def get_chromosome( self ):
		chrs = set(t.get_chromosome() for t in self.transcripts)
		if len( chrs ) == 1:
			return chrs.pop()
		else:
			print "Found mixed transcript chromosomes!"
			return 0
	def get_strand( self ):
		strands = set( t.get_strand() for t in self.transcripts)
		if len( strands ) == 1:
			return strands.pop()
		else:
			print "Found mixed transcript strands!"
			return 0
	
class Transcript:
	"""Holds information about a transcript.
	Args:
		id(str): some kind of identification (arbitrary)
		start(int): transcript start
		end(int): transcript end
		tss(int): transcription start site (redundant with start?)
		prin_iso(bool): is this the principal isoform?
		exons([Exons]): Exon objects contained within this transcript 
	"""
	def __init__(self, id, start, end, tss, principal_isoform, exons):
		self.id = id
		self.start = start
		self.end = end
		self.tss = tss
		self.prin_iso = principal_isoform
		self.exons = exons
	def __str__( self ):
		return "Transcript %s %s %s" % (self.id, self.start, self.end, self.tss, self.prin_iso, self.exons)
	def get_chromosome( self ):
		chrs = set(exon.location.chr for exon in self.exons)
		if len( chrs ) == 1:
			return chrs.pop()
		else:
			print "Found mixed exon chromosomes!"
			return 0
	def get_strand( self ):
		strands = set( e.location.strand for e in self.exons )
		if len( strands ) == 1:
			return strands.pop()
		else:
			print "Found mixed exon strands!"
			return 0
class Exon:
	"""Holds information about an exon.
	Args:
		location(GenomicLocation): where is this exon located?
		constitutive(bool): is this a constitutive exon?
	"""
	def __init__(self, location, constitutive):
		self.location = location
		self.constitutive = constitutive
	def __str__( self ):
		return "Exon %s %s" % (self.location, self.constitutive)

class HrTemplate:
	"""Holds information about a template for homologous recombination.
	Args:
		seq(Bio.Seq): sequence of the HR template
		target_site(GenomicLocation): genomic location the template is targeting
		target_seq(Bio.Seq): sequence of the targeted site
		frame(int): can be +1, +2, +3, -1, -2, -3, or 0 (non-coding)	
	"""
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
			possible_bases3 = ["A", "T", "G"] # 3rd base of pam
			possible_bases2 = ["A", "G" ] # 2nd base of pam
			pos3, pos2 = range(sgrna_template_index-3, sgrna_template_index-1) # only get the last two bases of the PAM
		else:
			possible_bases3 = ["A", "T", "C"] # 3rd base of pam
			possible_bases2 = ["T", "C" ] # 2nd base of pam
			pos2,pos3 = range(sgrna_template_index + len(protospacer)+1, sgrna_template_index + len(protospacer)+3) # only get the last 2 bases of the PAM
		#print pos2, pos3
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
			random.shuffle( possible_bases2 )
			possible_bases2.insert(0, self.seq[ pos2 ] ) # add existing 2nd PAM base to the list. that way we start with identity at pos2 and vary pos3, then move on to varying both
			for base2 in possible_bases2:
				for base3 in possible_bases3:
					seq_copy = self.seq[:min(pos2,pos3)] + base2 + base3 + self.seq[max(pos2,pos3)+1:] # copy of the sequence with PAM bases replaced
					if self.frame < 1: 
						test_translation = seq_copy.reverse_complement()[abs(self.frame)-1:].translate()
						ref_translation = self.seq.reverse_complement()[abs(self.frame)-1:].translate()
					else:
						test_translation = seq_copy[self.frame-1:].translate()
						ref_translation = self.seq[self.frame-1:].translate()
#					print protospacer
#					print test_translation
#					print ref_translation
#					print
					if str(test_translation) == str(ref_translation):
						#print self.seq
						#print seq_copy
						self.seq = seq_copy
						return True
		#print "Could not remove PAM!"
		return False			

class SgRna:
	"""Holds information about a single guide RNA.
	Args:
		protospacer(Bio.Seq): sequence of the protospacer (sans constant portion). Can only be set on initialization.
		target_site(GenomicLocation): location targeted by the protospacer
		target_seq(str): sequence window +/- 10 bases around protospacer (can be used to find PAM)
		offtarget_sites{GenomicLocation: [geneid1, geneid2,...]}: holds info about potential offtarget sites found in genome of interest, including if those offtargets fall within genes
		pam(str): protospacer adjacent motif for this guide
		constant_region(Bio.Seq): constant region associated with this guide. Used to calculate secondary structure.
		score(float): score of this guide
	"""
	# instance variables
	def __init__(self, seq, constant_region="GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU", target_site=None, target_seq="", pam="" ):
					  # weissman constant = "GUUUAAGAGCUAAGCUGGAAACAGCAUAGCAAGUUUAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUUU"
					     # broad constant = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU"
		# turn DNA input into RNA
		seq_copy = seq
		if seq_copy.find("T"): # convert to RNA if not already done
			seq_copy = seq.replace( "T", "U" )
		#print seq, seq_copy
		self.protospacer = Seq(seq_copy, generic_rna) # sequence sans constant portion. can only set protospacer on initialization. always stored as RNA
		self.target_site = target_site # will eventually become a GenomicLocation. strand (+ means sgrna seq same as + strand, - means sgrna seq same as - strand), 
		self.target_seq = target_seq # 10 bases on either side of site
		self.offtarget_sites = {} # dict, format = {GenomicLocation: [gene1, gene2, gene3...]}
		self.pam = pam
		self.constant_region = Seq( constant_region, generic_rna )
		self.score = 0
	def __eq__( self, other ):
		return( (self.protospacer, self.target_site, self.target_seq, self.offtarget_sites, self.constant_region, self.score) == (other.protospacer, other.target_site, other.target_seq, other.offtarget_sites, other.constant_region, other.score))
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
	
	def calculate_score(self):
		print "Scoring sgrna %s" % self.protospacer
		score = 0
		# lower score is better
		# +1000000 if target site is undefined (not found in genome)
		# +X^2, where X is stability of secondary structure in the protospacer 
		#	(squaring gives reasonably bigger penalties for stability and takes care of negatives)
		# +10 if last base before PAM is C (Doench et al Nat Biotech 2014)
		# -5 if last base before PAM is G (Doench et al Nat Biotech 2014)
		# -5 if PAM is NGG (vs NNG)
		# +1000 for each U homopolymer (4+) (sliding window)
		# +10 for homoX sequences other than U
		# +5 if 3' -1, -2, -4 bases of protospacer is a U (see Wang et al for sgRNA loading data)
		# +25 if GC content 0.4 < x > 0.8
		# +100 for each offtarget site that passed bowtie filters
		# 	or +10000 for each gene (refGene) that offtarget site was in
		# +10 for each 3' G in a row after the PAM
		# +20 if no microhomology found TODO
		# ----
		# major penalty if couldn't find target site in the genome:
		# this can also happen if bowtie search failed due to too MANY matches (-m filter exceeded)
		if not self.target_site:
			print "Target site for %s not found in the genome!" % self.protospacer
			score += 1000000
		if vienna_loaded:
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
		# penalize for sequences adjacent PAM (see Doench et al Nat Biotech 2014)
		if self.protospacer[:-1] == "C":
			score += 10
		elif self.protospacer[:-1] == "G":
			score -= 5
		# penalize homoX (4+), especially homoU (pol III terminator)
		matches = re.finditer(r'(?=((\w)\2{3}))', str(self.protospacer)) #
		for match in matches:
			print "Found %s  " % match.group(1)
			if match.group(1) == "UUUU" or match.group(1) == "TTTT":
				score += 1000
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
					score += num_offtarget_genes * 10000
				else:
					score += 100
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
		self.score = score
		print
		return score