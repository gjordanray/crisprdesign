#!/usr/bin/env python
import RNA # Vienna RNA secondary structure prediction
from Bio.Seq import Seq, MutableSeq
from Bio.Alphabet import generic_dna, generic_rna
import tempfile
import subprocess
import re
import random

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
		

class HrTemplate:
	def __init__(self, seq):
		# instance variables
		self.seq = Seq(seq, generic_dna) # sequence of the hr template itself
		self.target_site = None # GenomicLocation object
		self.target_seq = None # genomic sequence of the region this template targets (not necessarily identical to seq of template)
		self.frame = 1 # can be +1, +2, +3, -1, -2, -3, or 0 (non-coding)
		self.model_pam = "NGG"
	def __eq__(self, other):
		return (( self.seq, self.target_location, self.target_seq, self.frame, self.model_pam ) == (other.seq, other.target_location, other.target_seq, other.frame, other.model_pam))
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
		possible_bases = ["A", "T", "C"]
		protospacer = sgrna.protospacer.back_transcribe()
		sgrna_template_index = self.seq.find( protospacer )
		if sgrna_template_index != -1:
			strand = "+"
		else:
			strand = "-"
			sgrna_template_index = self.seq.find(protospacer.reverse_complement())
		if sgrna_template_index == -1:
			print "Could not find protospacer in HR sequence!"
			return False
		print "Strand %s" % strand
		if strand == "-":
			pam_indices = range(sgrna_template_index-3, sgrna_template_index-1) # only get the last two bases of the PAM
		else:
			pam_indices = range(sgrna_template_index + len(protospacer)+1, sgrna_template_index + len(protospacer)+3) # only get the last 2 bases of the PAM
		if strand == "+":
			pam_indices.reverse() # we want to start with the last base
		# print pam_indices
		if self.frame == 0: # if we're in a noncoding region, replace with a randomly allowed base
			base = random.choice(possible_bases)
			seq_copy = self.seq[:pam_indices[0]] + base + self.seq[pam_indices[0]+1:]
			self.seq = seq_copy
			return True
		else:	# if we're in a coding region, possible_bases is shuffled and we iterate
			random.shuffle( possible_bases )
			random.shuffle( possible_bases )
			for i in pam_indices:
				#print i, self.seq[i]
				for base in possible_bases:
					if base != self.seq[i]:
						seq_copy = self.seq[:i] + base + self.seq[i+1:]
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
							print self.seq
							print seq_copy
							self.seq = seq_copy
							return True
		print "Could not remove PAM!"
		return False
#		print pam_indices
#	 	for index in pam_indices:
#	 		print self.seq[index]
#		pam = pam_region[2:4]
#		codon = Seq.translate( pam )
#		while pam == self.model_pam:
#			if self.frame == "+3" or "-1": # can't change the PAM, since it's the first base of a codon
#				return False				

class SgRna:
	# instance variables
	def __init__(self, seq):
		# turn DNA input into RNA
		seq_copy = seq
		if seq_copy.find("T"): # convert to RNA if not already done
			seq_copy = seq.replace( "T", "U" )
		self.protospacer = Seq(seq_copy, generic_rna) # sequence sans constant portion. can only set protospacer on initialization. always stored as RNA
		self.target_site = None # will eventually become a GenomicLocation. strand (+ means sgrna seq same as + strand, - means sgrna seq same as - strand), 
		self.target_seq = "" # 10 bases on either side of site
		self.offtarget_sites = [] # list, format = [GenomicLocation, [gene1, gene2, gene3...]]
		self.model_pam = "NGG"
		self.constant_region = Seq( "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUUUU", generic_rna )
		self.constant_secstruct, self.constant_stability = RNA.fold( str( self.constant_region ))
	def __eq__( self, other ):
		return( (self.protospacer, self.target_site, self.target_seq, self.offtarget_sites, self.model_pam, self.constant_region, self.constant_secstruct, self.constant_stability) == (other.protospacer, other.target_site, other.target_seq, other.offtarget_sites, other.model_pam, other.constant_region, other.constant_secstruct, other.constant_stability))
	def __ne__( self, other ):
		return not self == other
			
	def build_fullseq( self ):
		fullseq = self.protospacer + self.constant_region
		if fullseq[0] != "G":
			fullseq = "G" + fullseq
		return fullseq
	
	def _n_in_range( self, n, range1, range2 ):
		if (n >= range1) and (n <= range2):
			return True
		else:
			return False
				
	def _bowtie_search( self ):
		if self.target_seq == None:
			print "Must set target_seq to find offtargets!"
			return []
		phred33String = '++++++++44444=======!4I'
		# antisense_phredString = 'I4!=======44444++++++++'
		bowtie_constant_options = ['bowtie', '--nomaqround', '--best', '-n 3', '-l 12', '-k 3', '-e 39', '-p 4', '--suppress', '5,6,7', '--chunkmbs', '128', 'hg19'] #note '--chunkmbs 128' and '--suppress 5,6,7'is one option, but subprocess needs them separated
		temp_bowtein = tempfile.NamedTemporaryFile(delete=True)

		i=1
		temp_bowtein.write('@'+str(i)+'\n')
		start = self.target_seq.find( self.protospacer.back_transcribe() )
		if start != -1: # found protospacer in top strand of target_seq
			end = start+len(self.protospacer)+3
		else:
			start = self.target_seq.find( self.protospacer.back_transcribe() )-3
		if start != -1: # found protospacer in bottom strand of target_seq
			end = start+len(self.protospacer)+3
		else:
			print "Could not find protospacer in target_seq!"
			return []
		print self.target_seq[start:end]
		temp_bowtein.write( str(self.target_seq[start:end])+'\n')			
		temp_bowtein.write('+\n')
		temp_bowtein.write(phred33String+'\n')
		temp_bowtein.flush() #necessary for bowtie to read files

		temp_bowtieout = tempfile.NamedTemporaryFile(delete=True)
#		bowtie_outfile = "bowtie.out"

		bowtie_cmdline = list(bowtie_constant_options)
		bowtie_cmdline.append(temp_bowtein.name)
		bowtie_cmdline.append(temp_bowtieout.name)
		subprocess.call( bowtie_cmdline ) #, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb') )
		temp_bowtein.close() # remove tempfiles

		found_locations = []
		for line in temp_bowtieout:
			l = line.strip().split('\t')
			loc = GenomicLocation( l[2], long(l[3]), long(l[3])+len(self.protospacer), l[1] ) # format of bowtie output is i, strand, chr, start
			found_locations.append(loc)
		temp_bowtieout.close()
		return found_locations

	def _read_refgene( self, fname ):
		fhandle = open( fname, "rU" )
		genes=[]
		for line in fhandle:
			split = line.split('\t')
			# format is    bin name chr strand txStart txEnd cdsStart cdsEnd
		#	genes.setdefault( split[1], (split[2], split[3], long(split[4]), long(split[5])) # genes container now holds:  chr strand txStart txEnd
			genes.append( [split[1], split[2], split[3], long(split[4]), long(split[5])] ) # genes container now holds: geneid, chr, strand, txStart, txend
		fhandle.close()
		return genes

	def find_offtargets(self):
		# run bowtie for protospacer, including PAM
		genomic_sites = self._bowtie_search() # clear any existing data
		# offtargets are everything that's not the target_site
		for site in genomic_sites: # need to skip the target site...
			if site != self.target_site:
				offsite = []
				offsite.append( site )
				# read refGene, find chromosome, is this index in a gene?
				refgene = self._read_refgene( "refGene.hg38.clean.sorted.txt" )
				for gene in refgene:
					rg_gene, rg_chr, rg_start, rg_end = gene[0], gene[1], gene[3], gene[4]
					if rg_chr == site.chr:
						if self._n_in_range( site.start, rg_start, rg_end) or self._n_in_range( site.end, rg_start, rg_end):
							offsite.append( rg_gene )
				self.offtarget_sites.append( offsite )
#			else:
#				print "Found target site!"

	def _GCcontent(self, seq):
		GC = (seq.count("G") + seq.count("C")) / float(len(seq))
		return GC

	def _hamming_dist( self, str1, str2 ):
		diffs=0
		for ch1, ch2 in zip( str1, str2 ):
			if ch1 != ch2:
				diffs += 1			
		return diffs

	def score(self):
		print "Scoring sgrna..."
		score = 0
		# lower score is better
		# +25 if GC content
		# +100 for each U homopolymer (sliding window)
		# +25 if 3' base of protospacer is a U
		# +X^2, where X is stability of secondary structure in the protospacer 
		#	(squaring gives reasonably bigger penalties for stability and takes care of negatives)
		# (below are optional, if genomic context and offtarget sites are defined)
		# +20 if 2nd base in PAM is A instead of G
		# +10 for each 3' G in a row after the PAM
		# +25 if not in DNAse hypersensitive site (TODO cutoff)
		# +25 for each offtarget site that passed bowtie filters
		# 	or +100 for each gene (refGene) that offtarget site was in
		# +20 if no microhomology found
		# ----
		# penalize secondary structure in the protospacer sequence
		secstruct, stability = RNA.fold( str(self.protospacer) )
		print "Protospacer secstruct & stability %s %f kcal/mol" % (secstruct, stability)
		score += stability**2
		# mildly penalize bad secondary structure introduced to the constant region
		full_secstruct, full_stability = RNA.fold( str(self.build_fullseq()))
		if full_secstruct[-len(self.constant_secstruct):] != self.constant_secstruct:
			print "constant  secstruct %s" % self.constant_secstruct
			print "+guide    secstruct %s" % full_secstruct[-len(self.constant_secstruct):]
			score += self._hamming_dist( full_secstruct[-len(self.constant_secstruct):], self.constant_secstruct )
		# penalize homoX (3+), especially homoU (pol III terminator)
		matches = re.finditer(r'(?=((\w)\2{2}))', str(self.protospacer)) #
		for match in matches:
			print "Found %s  " % match.group(1),
			if match.group(1) == "UUU" or match.group(1) == "uuu":
				score += 50
			else:
				score += 10
		print
		# remove bad GC sequences
		GC = self._GCcontent(self.protospacer)
		if (GC < 0.4 or GC > 0.8):
			print "GC content is %f" % GC,
			score += 25
		# penalize 3' Us
		if self.protospacer[:-1] == "U": #TODO, also add -4, -2, see Wang et al for sgRNA loading
			print "3' U found",
			score += 10
		# penalize offtarget sites in genome, also lack of DNAse hypersensitivity
		if self.offtarget_sites:
			print "Found %s offtarget sites" % len( self.offtarget_sites )
			for offsite in self.offtarget_sites:
				num_offtarget_genes = len(offsite) - 1 # all entries > 1 in offsites list are offtarget genes
				if num_offtarget_genes > 0:
					print "Found offtarget sites in the following genes: ,"
					for offtgt in offsite[1:]:
						print offtgt,
					print
					score += num_offtarget_genes * 100
				else:
					score += 25
		# penalize based on facts  about the site targeted in the genome
		if self.target_seq:
			# PAM profile
			if self.target_site.strand == "+":
				downstream_seq = self.target_seq
			else:
				downstream_seq = self.target_seq.reverse_complement()
			grna_local_index = downstream_seq.find( self.protospacer.back_transcribe() )
			if grna_local_index == -1:
				print "Could not find protospacer exact match in target site!"
				return score
			downstream_seq = downstream_seq[grna_local_index + len(self.protospacer):] # everything after the sgrna
			pam = downstream_seq[:3]
			print "model PAM %s | PAM %s" % (self.model_pam, pam)
			if pam[1] != self.model_pam[1]:
				score += 20
			if pam[2] != self.model_pam[2]:
				score += 10000 # something is wrong... last base should always match
			# penalize G homopolymers after the PAM
			i=0
			for base in downstream_seq[len(pam):]:
				if base == "G":
					i+=1
					score += 20
				else:
					break
			if i>0: print "Found %i Gs following PAM" % i
		return score