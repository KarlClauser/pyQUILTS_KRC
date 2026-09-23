# proteome-genes.txt: ENSP00000000233 ENST00000000233 ENSG00000004059 ARF5
# proteome-descriptions.txt: ENSP00000000233-ARF5    ADP ribosylation factor 5 [Source:HGNC Symbol;Acc:HGNC:658]
# proteome.bed: chr7	127588498	127591299	ENSP00000000233	1000	+	127588498	127591299	0	6	67,81,110,72,126,87	0,584,986,1567,2464,2714
# proteome.bed: chr1	924431	939291	ENSP00000411579	1000	+	924431	939291	0	7	517,92,182,51,125,90,17	0,1490,5723,6607,11340,14608,14843

import sys
import os
import gzip
from subprocess import call, check_call, CalledProcessError

sREAD_CHR_BED_EEXECUTABLE = 'read_chr_bed_win64.exe'	#Windows
#sREAD_CHR_BED_EEXECUTABLE = 'read_chr_bed'				#Unix

class Gene:
	def __init__(self, chromosome, gene_id, protein_id, transcript_id, gene_name, gene_desc, strand):
		self.chromosome = 'chr'+chromosome
		self.gene_id = gene_id
		self.protein_id = protein_id
		self.transcript_id = transcript_id
		self.gene_name = gene_name
		self.gene_desc = gene_desc
		self.strand = '+' if strand=='1' else '-'
		self.trans_exons = {}
		self.prot_exons = {}
	
	def add_exon(self, exon_no, prot_start, prot_end, trans_start, trans_end):
		self.trans_exons[int(exon_no)] = [int(trans_start)-1, int(trans_end)]
		if prot_start != '' and prot_end != '':
			self.prot_exons[int(exon_no)] = [int(prot_start)-1, int(prot_end)]

	def print_bed_line(self):
		#exon_nos = self.prot_exons.keys()
		#exon_nos.sort()							#KRC 1/24/2024 python 2
		exon_nos = sorted(self.prot_exons.keys())	#KRC 1/24/2024 python 3
		#k = d.keys(); k.sort(). Use k = sorted(d)
		total_length = 0
		# End phase dealings - there can be an end phase on the final codon, which throws off 
		# calculations, especially on negative strands. Cutting off the end phase if the length
		# of the NT sequence is not divisible by 3
		for i in exon_nos:
			if i in self.prot_exons.keys():
				total_length += self.prot_exons[i][1]-self.prot_exons[i][0]
		endphase_modifier = total_length%3
		if self.strand == '+':
			self.prot_exons[max(exon_nos)][1] -= endphase_modifier
			start = self.prot_exons[min(exon_nos)][0]
			end = self.prot_exons[max(exon_nos)][1]
		else:
			exon_nos.reverse()
			self.prot_exons[max(exon_nos)][0] += endphase_modifier
			end = self.prot_exons[min(exon_nos)][1]
			start = self.prot_exons[max(exon_nos)][0]
		lengths = []
		starts = []
		for i in exon_nos:
			if i in self.prot_exons.keys():
				ex_start, ex_end = self.prot_exons[i][0], self.prot_exons[i][1]
				lengths.append(str(ex_end-ex_start))
				starts.append(str(ex_start-start))
		print_str = "%s\t%d\t%d\t%s\t1000\t%s\t%d\t%d\t0\t%d\t%s\t%s\n" % (self.chromosome, start, end, self.protein_id, self.strand, start, end, len(self.prot_exons), ','.join(lengths), ','.join(starts))
		# Negative strand!?
		return print_str
	
	def print_trans_line(self):
		#exon_nos = self.trans_exons.keys()
		#exon_nos.sort()							#KRC 1/24/2024 python 2
		exon_nos = sorted(self.trans_exons.keys())	#KRC 1/24/2024 python 3
		#k = d.keys(); k.sort(). Use k = sorted(d)

		if self.strand == '+':
			start = self.trans_exons[min(exon_nos)][0]
			end = self.trans_exons[max(exon_nos)][1]
		else:
			exon_nos.reverse()
			end = self.trans_exons[min(exon_nos)][1]
			start = self.trans_exons[max(exon_nos)][0]
		lengths = []
		starts = []
		for i in exon_nos:
			if i in self.trans_exons.keys():
				ex_start, ex_end = self.trans_exons[i][0], self.trans_exons[i][1]
				lengths.append(str(ex_end-ex_start))
				starts.append(str(ex_start-start))
		print_str = "%s\t%d\t%d\t%s\t1000\t%s\t%d\t%d\t0\t%d\t%s\t%s\n" % (self.chromosome, start, end, self.transcript_id, self.strand, start, end, len(self.trans_exons), ','.join(lengths), ','.join(starts))
		# Negative strand!?
		return print_str

bMakeBedDNAonly = False
iNumArgs = len(sys.argv)

if iNumArgs > 2 and bMakeBedDNAonly == False:
	# Clean up the proteome FASTA - remove everything that starts with an X or contains a stop codon
	sCleanedProteomeFilename = 'proteome.fasta'
	#w = open('cleaned_proteome.fasta','w')
	w = open(sCleanedProteomeFilename,'w')		#KRC 1/24/2024 avoid later unix mv command by giving final name now
	header = ''
	seq = ''
	
	#KRC 1/24/2024 allow gzip fasta file
	#f = open(sys.argv[2],'r')
	#line = f.readline()
	# while line:
	#     if line[0] == '>':
	#         if seq != '' and seq[0] != 'X' and '*' not in seq:
	#             w.write(header)
	#             w.write(seq+'\n')
	#         header = line
	#         seq = ''
	#    else:
	#        seq += line.rstrip()
	#    line = f.readline()
	fil = sys.argv[2]
	if not fil.endswith('.gz'):
		f = open(fil, 'r')
		Llines = f.read().splitlines()	#list of lines with \n removed
	else:
		with gzip.open(fil, 'rt') as f:
			Llines = f.read().splitlines()	#list of lines with \n removed
	
	iRemovedEntries = 0
	for line in Llines:
		if line[0] == '>':
			if seq != '' and seq[0] != 'X' and '*' not in seq:
				w.write(header+'\n') #KRC 1/27/2024 added \n
				w.write(seq+'\n')
			else:
				iRemovedEntries += 1
			header = line
			seq = ''
		else:
			seq += line.rstrip()
	
	print('Removed entries that start with an X or contain a stop codon *:', iRemovedEntries)
	print('Wrote cleaned proteome file:', sCleanedProteomeFilename)
	
	if seq[0] != 'X' and '*' not in seq: # last one
	    w.write(header)
	    w.write(seq+'\n')
	f.close()
	w.close()
	#os.system('mv %s orig_proteome.fasta' % sys.argv[2])	#KRC 1/24/2024 keep original file name
	#os.system('mv cleaned_proteome.fasta proteome.fasta')	#KRC 1/24/2024 name already given at open

if iNumArgs > 1 and bMakeBedDNAonly == False:
	# Finding where each column lives in the input file
	f = open(sys.argv[1], 'r')
	header = f.readline().rstrip().split('\t')
	chrm_pos = header.index("Chromosome/scaffold name")
	gene_id_pos = header.index("Gene stable ID")
	transcript_id_pos = header.index("Transcript stable ID")
	protein_id_pos = header.index("Protein stable ID")
	gene_name_pos = header.index("Gene name")
	gene_desc_pos = header.index("Gene description")
	strand_pos = header.index("Strand")
	exon_no_pos = header.index("Exon rank in transcript")
	exon_start_pos = header.index("Exon region start (bp)")
	exon_end_pos = header.index("Exon region end (bp)")
	cds_start_pos = header.index("Genomic coding start")
	cds_end_pos = header.index("Genomic coding end")
	genes = {}
	
	# Parsing the input file
	line = f.readline()
	while line:
		spline = line.rstrip('\n').split('\t')
		gene_name = spline[transcript_id_pos]
		if gene_name not in genes:
			new_gene = Gene(spline[chrm_pos], spline[gene_id_pos], spline[protein_id_pos], spline[transcript_id_pos], spline[gene_name_pos], spline[gene_desc_pos], spline[strand_pos])
			genes[gene_name] = new_gene
		genes[gene_name].add_exon(spline[exon_no_pos], spline[cds_start_pos], spline[cds_end_pos], spline[exon_start_pos], spline[exon_end_pos])
		line = f.readline()	
	f.close()
	
	wgenes = open('proteome-genes.txt','w')
	wdesc = open('proteome-descriptions.txt','w')
	wbed = open('proteome.bed','w')
	wtrans = open('transcriptome.bed','w')
	for g in genes:
		gene = genes[g]
		if gene.trans_exons != {}:
			wgenes.write('%s\t%s\t%s\n' % (gene.protein_id,gene.transcript_id,gene.gene_name))
			wdesc.write('%s\t%s\n' % (gene.protein_id, gene.gene_desc))
			wtrans.write(gene.print_trans_line())
			if gene.prot_exons != {}:
				wbed.write(gene.print_bed_line())
	wgenes.close()
	wdesc.close()
	wtrans.close()
	wbed.close()

if iNumArgs > 3:
	#KRC 1/28/2024 make a central copy proteome.bed.dna file for the proteome/genome combination
	sProteomeBed = 'proteome.bed'
	if os.path.isfile(sProteomeBed):
		print ('Making proteome.bed.dna for central use with this proteome genome pair')
		try:
			script_dir = '..\\..\\..\\pyQUILTS_KRC\\' #expect runnning this script running in QUILTS\references\proteome\<proteome>
			#check_call("%s/read_chr_bed %s/log/proteome.bed %s" % (script_dir, results_folder, args.genome), shell=True)
			#check_call("%s%s %s/log/proteome.bed %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, results_folder, args.genome), shell=True)
			check_call("%s%s %s %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, sProteomeBed, sys.argv[3]), shell=True)
		except CalledProcessError:
			raise SystemExit("ERROR: read_chr_bed didn't work - now we don't have a proteome.bed.dna file. Try recompiling read_chr_bed.c.\nAborting program.")
		print ('Made proteome.bed.dna file.')
	else:
		print('ERROR: File not found, should have been made earlier by this script',  sProteomeBed)
