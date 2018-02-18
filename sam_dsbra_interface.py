################################################
# sam_dsbra_interface.py -r <reference.fa> -q <sample.fastq|.fq> -b <break_index>,<margin>,<last_margin> -o <output_filename.txt> [-v] [-c]
#
################################################

__version__ = '0.0.9'

## Characteristics of Scar: Num_Mismatch, Location, Transitions, Transversions, Num_Insertion_Events, Location, Sequences, Num_Deletion_Events, Location, Length, is_microhomologous		
		
import sys
import re
import subprocess
import pysam
import cPickle
import os
import numpy as np
import csv
from scipy.stats import norm
import time
from datetime import datetime

CIGAR_RE = re.compile('\d+[MIDNSHP]')
BOWTIE2_FOLDER = ''
BWA_FOLDER = ''
SAMTOOLS_FOLDER = ''

RE_FASTA = re.compile(r'\>(.*)[\n|\r]+([ACTGactg\r\n]+)')
RUN_NAME = 'test_run'

SAM_FILENAME = ''

VERBOSE = False

class RepairSample:
	def __init__(self):
		self.total_samples

class RepairPattern:
	def __init__(self,read_name):
		self.read_name = read_name
		self.repair_size = 0
	
		self.ref_begin = -1
		self.ref_end = -1
		self.is_filtered = False
		self.repair_sequence = ''
	
		self.mismatch_locations = [ ]
		self.new_mismatch_bases = [ ]
		self.old_mismatch_bases = [ ]
		self.num_transitions = 0
		self.num_transversions = 0
		
		self.num_insertion_events = 0
		self.insertion_locations = [ ]
		self.inserted_sequences = [ ]
		
		self.num_deletions = 0
		self.deletion_locations = [ ]
		self.deletion_lengths = [ ]
		self.deletion_is_microhomologous = [ ]
		self.deletion_micro = [ ]

	def add_insertion(self, ins_loc, ins_str):
		self.num_insertion_events += 1
		self.insertion_locations.append(ins_loc)
		self.inserted_sequences.append(ins_str)

	def add_deletion(self, del_loc, del_size, is_micro, micro_seq):
		self.num_deletions += 1
		self.deletion_locations.append(del_loc)
		self.deletion_lengths.append(del_size)
		self.deletion_is_microhomologous.append(is_micro)
		self.deletion_micro.append(micro_seq)

	def add_mismatch(self, mismatch_loc, new_base, old_base):
		self.mismatch_locations.append(mismatch_loc)
		self.new_mismatch_bases.append(new_base)
		self.old_mismatch_bases.append(old_base)

		if old_base == 'A':
			if new_base == 'G':
				self.num_transitions += 1
			else:
				self.num_transversions += 1

		elif old_base == 'C':
			if new_base == 'T':
				self.num_transitions += 1
			else:
				self.num_transversions += 1

		elif old_base == 'G':
			if new_base == 'A':
				self.num_transitions += 1
			else:
				self.num_transversions += 1

		else:
			if new_base == 'C':
				self.num_transitions += 1
			else:
				self.num_transversions += 1

	def filter_pattern(self,break_index,index_margin,outside_margin,ref_seq):
		if not self.is_filtered:
			#make list of start/end of events
			event_locations = list(self.insertion_locations)
			event_locations += self.mismatch_locations
			for i in range(len(self.deletion_locations)):
				for j in range(self.deletion_lengths[i]):
					event_locations += [self.deletion_locations[i]+j]
			
			event_locations = (sorted(set(event_locations)))
			bottom_range = -1
			top_range = -1

			for i in range(len(event_locations)):
				#first event within range
				if (event_locations[i] >= break_index-index_margin) and (event_locations[i] <= break_index+index_margin):
					#backtrack
					j = i
					top_range = bottom_range = event_locations[i]
					while (j > 0):
						if event_locations[j]-event_locations[j-1] <= outside_margin:
							bottom_range = event_locations[j-1]
							j -= 1
						else:
							break

					top_range = bottom_range

					j = i
					while (j < len(event_locations)-1):
						if (event_locations[j+1]-event_locations[j] <= outside_margin) or (event_locations[j+1] < break_index + index_margin):
							top_range = event_locations[j+1]
							j += 1
						else:
							break
							
					break


			#TODO: Filter + reconstruct scar
			self.ref_begin = bottom_range
			self.ref_end = top_range

			ins_to_pop = []
			for i in self.insertion_locations:
				if i < bottom_range or i > top_range:
					ins_to_pop.append(i)
					
			for i in ins_to_pop:
				self.num_insertion_events -= 1
				self.inserted_sequences.pop(self.insertion_locations.index(i))
				self.insertion_locations.pop(self.insertion_locations.index(i))

			del_to_pop = []
			for i in self.deletion_locations:
				if i < bottom_range or i > top_range:
					del_to_pop.append(i)

			for i in del_to_pop:
				self.num_deletions -= 1
				self.deletion_is_microhomologous.pop((self.deletion_locations.index(i)))
				self.deletion_lengths.pop(self.deletion_locations.index(i))
				self.deletion_locations.pop(self.deletion_locations.index(i))

			#TODO: figure out how to handle mismatch only
			"""mismatch_to_pop = []
			for i in self.mismatch_locations:
				if i < bottom_range or i > top_
					mismatch_to_pop.append(i)

			for i in mismatch_to_pop:
				old_base = self.old_mismatch_bases[self.mismatch_locations.index(i)]
				new_base = self.new_mismatch_bases[self.mismatch_locations.index(i)]

				if old_base == 'A':
					if new_base == 'G':
						self.num_transitions -= 1
					else:
						self.num_transversions -= 1

				elif old_base == 'C':
					if new_base == 'T':
						self.num_transitions -= 1
					else:
						self.num_transversions -= 1

				elif old_base == 'G':
					if new_base == 'A':
						self.num_transitions -= 1
					else:
						self.num_transversions -= 1

				else:
					if new_base == 'C':
						self.num_transitions -= 1
					else:
						self.num_transversions -= 1
					
				self.new_mismatch_bases.pop(self.mismatch_locations.index(i))
				self.old_mismatch_bases.pop(self.mismatch_locations.index(i))
				self.mismatch_locations.pop(self.mismatch_locations.index(i))"""

			#Reconstructing scar -- should be filtered already
			working_sequence = ref_seq[self.ref_begin:self.ref_end+1]
			
			#this is kind of a patch -- fix this?
			last_event = True
		
			ins_to_pop = []
			for i in range(self.ref_begin,self.ref_end+1)[::-1]:
				if i in self.insertion_locations:
					working_sequence = working_sequence[:i-self.ref_begin] + '[' + self.inserted_sequences[self.insertion_locations.index(i)] + ']' + working_sequence[i-self.ref_begin:]
					#pop last letter here
					if last_event == True:
						working_sequence = working_sequence[:-1]
						last_event = False

				elif i in self.deletion_locations:
					del_length = self.deletion_lengths[self.deletion_locations.index(i)]
					working_sequence = working_sequence[:i-self.ref_begin] + '^(' + working_sequence[i-self.ref_begin:(i-self.ref_begin)+del_length] + ')' + working_sequence[i-self.ref_begin+del_length:]
					
					if last_event == True:
						last_event = False
					
				elif i in self.mismatch_locations:
					working_sequence = working_sequence[:i-self.ref_begin] + '*' + self.new_mismatch_bases[self.mismatch_locations.index(i)] + '*' + working_sequence[i-self.ref_begin+1:]
					
					if last_event == True:
						last_event = False
					
			self.repair_sequence = working_sequence
		
if __name__ == '__main__':
	# INITIALIZATION SECTION -- reads all values from command line and #
	# assigns values to appropriate variables for use downstream       #
	break_index,margin,last_margin = sys.argv[sys.argv.index('-b')+1].split(',')
	break_index = int(break_index)
	margin = int(margin)
	last_margin = int(last_margin)
	fastq_name = sys.argv[sys.argv.index('-q')+1]
	ref_fa = sys.argv[sys.argv.index('-r')+1]
	output_filename = sys.argv[sys.argv.index('-o')+1]
	
	VERBOSE = '-v' in sys.argv
	
	SAM_FILENAME = fastq_name[:fastq_name.rfind('.')]+'.sam'

	#Read in reference fasta, write to new file after checking filetype
	ref_fasta = open(ref_fa,'r')
	ref_fasta_str = ref_fasta.read()

	index_ref = open(RUN_NAME+'_ref.fa','w')
	first_seq = RE_FASTA.match(ref_fasta_str)

	ref_name = ''
	ref_seq = ''

	if first_seq:
		ref_name, ref_seq = first_seq.groups(0)
		ref_seq = re.sub(r'[\n\r]','',ref_seq).upper()
		
		index_ref.write('>'+ref_name+'\n'+ref_seq+'\n')
		
	else:
		raise EOFError('FASTA sequence not detected before end of file ' + str(ref_fa))

	ref_fasta.close()
	index_ref.close()
	
	run_metadata = {'break_index':break_index,'margin':margin,'last_margin':last_margin,'ref_filename':ref_fa,'fastq_file':fastq_name,'ref_seq':ref_seq,'timestamp':str(datetime.now()) }

	# Check if BAM alignment file exists and create new one if not the case.
	# TODO: If BAM exists, verify that reference sequence used is the same
	# sequence used to create alignment
	if True:
		with open(os.devnull, 'wb') as out:
			if VERBOSE:
				print "Creating index from file %s"%ref_fa

			# Build index from reference.
			# TODO: Error handling here?
			if not '-bwa' in sys.argv:
				subprocess.call(BOWTIE2_FOLDER+'bowtie2-build -f '+RUN_NAME+'_ref.fa '+RUN_NAME+'_index',shell=True, stdout = out)
			else:
				subprocess.call(BWA_FOLDER+'bwa index '+RUN_NAME+'_ref.fa',shell=True)
				

			if VERBOSE:
				print "Index creation COMPLETE.\nAligning sample fastq file to index."

			# Align fastq file using BOWTIE2 
			# TODO: Check for fastq errors
			if not '-bwa' in sys.argv:
				run_metadata['alignment_settings'] = "BOWTIE2_FOLDER+'bowtie2 --local -p 8 -x '+RUN_NAME+'_index -q '+fastq_name+' -S '+SAM_FILENAME"
				subprocess.call(BOWTIE2_FOLDER+'bowtie2 --local -p 8 -x '+RUN_NAME+'_index -q '+fastq_name+' -S '+SAM_FILENAME,shell=True)
			else:
				run_metadata['alignment_settings'] = "BWA_FOLDER+'bwa mem '+RUN_NAME+'_ref.fa '+fastq_name+' > '+SAM_FILENAME"
				subprocess.call(BWA_FOLDER+'bwa mem '+RUN_NAME+'_ref.fa '+fastq_name+' > '+SAM_FILENAME,shell=True)
				
			subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+SAM_FILENAME+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 - '\
				+SAM_FILENAME+'.sorted', shell=True, stdout = out)
			subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ SAM_FILENAME + '.sorted.bam ' + SAM_FILENAME + '.sorted.bai',\
				shell=True, stdout = out)
				
	bamfile = pysam.AlignmentFile(SAM_FILENAME+'.sorted.bam','rb')
	
	#VERSION 0.0.4 -- Mutational database now incorporates 
	#First pass on SAM/BAM file creates secondary file with all reads requiring further evaluation (non-WT within margin)
	mut_reads_list = [ ]
	MUT_BAMFILE_NAME = "temp_mut_reads.sam"
	
	if os.path.exists(MUT_BAMFILE_NAME):
		print "WARNING: File \"%s\" exists."%MUT_BAMFILE_NAME
		sys.exit()
		
	mut_bamfile = pysam.AlignmentFile(MUT_BAMFILE_NAME,"wh", template=bamfile, reference_names=bamfile.references, reference_lengths=bamfile.lengths)
	
	print "Beginning first scan of SAM file..."
	start_time = time.time()
	
	#Find average coverage within margin
	coverage_list = []
	
	for pileupcolumn in bamfile.pileup(bamfile.references[0], start=break_index-margin, end=break_index+margin, truncate=True, max_depth=250000):
		coverage_list.append(pileupcolumn.n)
		
		for pileupread in pileupcolumn.pileups:			
			qname = pileupread.alignment.query_name
			
			if qname in mut_reads_list:
				continue
			
			mutated = False
			#check if mismatch occurring
			if pileupread.indel != 0 or pileupread.is_del:
				mutated = True
			#TODO: Figure out what to do with mismatch only
			elif ('-m' in sys.argv) and pileupread.indel == 0 and not pileupread.is_del:
				if ref_seq[pileupcolumn.pos] != pileupread.alignment.query_sequence[pileupread.query_position]:
					mutated = True
					
			if mutated:
				mut_reads_list.append(qname)
				mut_bamfile.write(pileupread.alignment)
					
	run_metadata['coverage'] = int(np.mean(coverage_list))
					
	#Close and sort/index new bamfile
	mut_bamfile.close()
	with open(os.devnull, 'wb') as out:
		subprocess.call(SAMTOOLS_FOLDER+'samtools view -@ 8 -bS '+MUT_BAMFILE_NAME+' | '+SAMTOOLS_FOLDER+'samtools sort -@ 8 - '\
			+MUT_BAMFILE_NAME+'.sorted', shell=True, stdout = out)
		subprocess.call(SAMTOOLS_FOLDER+'samtools index '+ MUT_BAMFILE_NAME + '.sorted.bam ' + MUT_BAMFILE_NAME + '.sorted.bai',\
			shell=True, stdout = out)
	bamfile.close()
	
	print "First scan and sorting ended. Total time: " + str(time.time()-start_time)
	
	#Now we do the dirty work -- probably can parallelize this -- kick identified reads to "analyzer"
	bamfile = pysam.AlignmentFile(MUT_BAMFILE_NAME+".sorted.bam","rb")
	
	#TODO: fix this for all references?
	#TODO: make sure adequate coverage on either side of break
	rp_list = [ ]
	
	for read in bamfile.fetch(bamfile.references[0]):
		rp = RepairPattern(read.query_name)
		lastq = -1
		lastr = -1
		indel = 0
		ins_str = ''

		started = False
		# for (x,y) x = query position, y = ref position
		for pair in read.get_aligned_pairs():
			if not started:
				if pair[0] and pair[1]:
					lastq = pair[0]
					lastr = pair[1]
					started = True
				else:
					continue

			#check for indels
			if not pair[0]: #deletion
				indel -= 1
			elif not pair[1]: #insertion
				indel += 1
				ins_str += read.query_sequence[pair[0]]
			else: #matched reads
				if indel != 0:
					#deletion
					if indel < 0:
						#checking for microhomology
						del_loc = pair[1] + indel
						outside_loc = pair[1]
						
						microhomology = ''
						while (ref_seq[del_loc] == ref_seq[outside_loc]) and (len(microhomology) < abs(indel)):
							microhomology += ref_seq[del_loc]
							del_loc += 1
							outside_loc += 1
							
						del_loc = pair[1]-1
						outside_loc = pair[1] - 1 + indel
						
						rev_microhomology = ''
						while ref_seq[del_loc] == ref_seq[outside_loc] and (len(rev_microhomology) < abs(indel)):
							rev_microhomology += ref_seq[del_loc]
							del_loc -= 1
							outside_loc -= 1
						
						if len(microhomology) < len(rev_microhomology):
							microhomology = rev_microhomology
							
						is_mmej = False
						if len(microhomology) >= 2:
							is_mmej = True
						
						rp.add_deletion(pair[1]+indel,abs(indel),is_mmej,microhomology)
					#insertion
					else:
						rp.add_insertion(pair[1],ins_str)
						ins_str = ''
					
					indel = 0

				#now check mismatch
				if read.query_sequence[pair[0]] != ref_seq[pair[1]]:
					rp.add_mismatch(pair[1],read.query_sequence[pair[0]],ref_seq[pair[1]])

		#TODO: Filter outside of window
		rp.filter_pattern(break_index,margin,last_margin,ref_seq)
		rp_list.append(rp)
	
	#Save all this info to database
	#TODO: have this run concurrently to make life easier/faster
	collection_name = fastq_name[:fastq_name.rfind('.')]+'_'+str(break_index)+'_'+str(margin)+'_'+str(last_margin)
	collection = []
	for rp in rp_list:
		rp_entry = {'read_name':rp.read_name,
					'ref_begin': rp.ref_begin,\
					'repair_size': rp.repair_size,\
					'mismatch_locations': rp.mismatch_locations,\
					'new_mismatch_bases': rp.new_mismatch_bases,\
					'num_mismatch': (rp.num_transitions + rp.num_transversions),\
					'num_transitions': rp.num_transitions,\
					'num_transversions': rp.num_transversions,\
					'num_insertions': rp.num_insertion_events,\
					'insertion_locs': rp.insertion_locations,\
					'insertion_seqs': rp.inserted_sequences,\
					'num_deletions': rp.num_deletions,\
					'deletion_locs': rp.deletion_locations,\
					'deletion_lens': rp.deletion_lengths,\
					'deletion_micro': rp.deletion_is_microhomologous,\
					'repair_sequence': rp.repair_sequence,\
					'deletion_is_micro': rp.deletion_is_microhomologous,\
					'micro_seq': rp.deletion_micro}
					
		collection.append(rp_entry)
		
	with open(output_filename,'wb') as op:
			rp_list = []
			rp_list_count = []
			
			for rp in collection:
				#remove read name to promote compression
				del rp['read_name']
				if not rp in rp_list:
					rp_list.append(rp)
					rp_list_count.append(1)
				else:
					rp_list_count[rp_list.index(rp)] += 1
			
			#Write to file
			op.write('#'+str(int(np.mean(coverage_list)))+'\n')
			
			for key in rp_list[0].keys():
				op.write(key)
				op.write('\t')
			op.write('count\t')
			op.write('mut_freq\t')
			op.write('\n')
			
			num_mut_seqs = np.sum(rp_list_count)
			
			for i in range(len(rp_list)):
				for x in rp_list[i]:
					op.write(str(rp_list[i][x])+'\t')
				op.write(str(rp_list_count[i])+'\t')
				op.write(str(float(rp_list_count[i])/num_mut_seqs))
				op.write('\n')