#Filename:	quiltsWrapper.py
#Purpose:	Replacement for Quilts main function
#Author:	Karl Clauser
#Created:	March 30, 20223

import sys
import os
import pandas as pd
import shutil #for file copying
import gzip
import argparse
from subprocess import call, check_call, CalledProcessError
import multiprocessing
import time
import warnings
from datetime import datetime
import quilts as qu

sREAD_CHR_BED_EEXECUTABLE = 'read_chr_bed_win64.exe'	#Windows
#sREAD_CHR_BED_EEXECUTABLE = 'read_chr_bed'				#Unix

"""
Overview
------------------------------------------------------------------------------
QuiltsWrapper.py replaces the quilts __main __ function to treat quilts.py as a library and call Emily’s processing functions
o	Split out code by type: junctions, variants, fusion for separate processing of each data type into its own file
o	Input data for all samples in one directory
o	Output data for all samples in one directory
o	Create intermediate results dir for running, then delete un-needed intermediate files when done.

A major motivating goal was to restructure the inputs and outputs to be more friendly to processing cohorts
instead of individual samples. Native Quilts expects a directory for each sample with inputs, outputs, and 
intermediates buried within. For a cohort one would prefer to have an input directory containing all samples
, an output directory containing  all samples, and the intermediates should be in a separate directory tree
so they can easily be deleted. Thus, one need not spend time burying the inputs and unburying the outputs.
Also one would like to easily process the data for each data type (somatic variants, germline variants,
splice junctions) separately to deal with the issues of each data type without the need to re-process the 
other data types each time. Quilts outputs are combined by a Spectrum Mill utility after running Quilts,
rather than have Quilts do it.

QUILTS native function			Replacement Here
------------------------------------------------------------------------------
qu.__main __					fMainMulti()
qu.set_up_output_dir()			fSetUpGlobalPaths()
qu.combine_output_fastas()		fConcatenateIntermediateFastas()
qu.pull_*_files()				fListFilesOfType()
qu.convert_star_to_mapsplice()	fConvertStarToMapsplice()
qu.parse_input_arguments()		fFill_argparseNamespace()
"""

###############################################################################
############################# Static Variables ####### ########################
###############################################################################
CODON_MAP = {"TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
	"ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TCT":"S",
	"TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P","ACT":"T","ACC":"T",
	"ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","TAA":"*",
	"TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","AAT":"N","AAC":"N","AAA":"K","AAG":"K",
	"GAT":"D","GAC":"D","GAA":"E","GAG":"E","TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R",
	"CGC":"R","CGA":"R","CGG":"R","AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G",
	"GGA":"G","GGG":"G"}

###############################################################################
############################# fSetUpGlobalPaths ###############################
###############################################################################
#KRC 3/30/2023 Revised qu.set_up_output_dir to act locally here, main in this script will set the qu.Globals
def fSetUpGlobalPaths(output_dir, args, sResultsPrefix):
	''' Sets up the results folder (creates it, sets up a log file and status file).
		Because this will probably confuse people (me) later, "output directory" refers
		to the directory the user specifies in which the "results folder" will reside. The
		results folder contains the actual results.'''
		
	# I hate this global variable stuff but is there another way to do it
	# without having to pass the logfile/statusfile to all the functions so they can be used?
# 	global logfile
# 	global referrorfile
# 	global statusfile
# 	global results_folder
	
	# Will the directory they give us be where we put the results, or where we put a folder
	# containing the results?
	
	# Checks to make sure the output directory exists.
	# Currently aborts if this is the case, but I could also just create the folder 
	# in the working directory? Maybe allow user input for this? Probably not though.
	# Wouldn't be good if they were running it through qsub or whatever.
	if not os.path.isdir(output_dir):
		#raise SystemExit("ERROR: Output directory does not exist.\nAborting program.")
		os.makedirs(output_dir)
	
	# Makes the results folder. Calls it results_(date and time). Looks ugly, but I'm cool with that.
	# Lets us avoid the problem of dealing with multiple results folders in the same output directory.
	today = datetime.today()
	day_time_string = str(today.year)+str(today.month).zfill(2)+str(today.day).zfill(2)+'.'+str(today.hour).zfill(2)+str(today.minute).zfill(2)+str(today.second).zfill(2)
	results_folder = output_dir+'/' + sResultsPrefix + '_IntermediateResults_'+day_time_string
	os.makedirs(results_folder)
	
	# Gives us addresses for the logfile and statusfile,
	# and writes the starting date and time to them.
	logfile = results_folder+'/log.txt'	
	statusfile = results_folder+'/status.txt'
	referrorfile = results_folder+'/reference_mismatches.txt'
	qu.write_to_log("Logfile created: "+str(today), logfile)

	#The globals are not yet created
	#qu.write_to_status("Status file created: "+str(today))
	#qu.write_to_log("Reference mismatches file created: "+str(today), referrorfile)
	
	# Creates some folders within the results folder
	# Currently just copying what's in the Perl version.
	# Likely to change as I start building out other functionality.
	os.makedirs(results_folder+"/log")
	if args.junction:
		os.makedirs(results_folder+"/junction")
	os.makedirs(results_folder+"/fasta")
	os.makedirs(results_folder+"/fasta/parts")
	
	# Moves reference proteome to the working area.
	try:
		#shutil.copy(args.proteome+"/proteome.fasta",results_folder+"/fasta/"+args.proteome.split("/")[-1]+".fasta")
		shutil.copy(args.proteome+"/proteome.fasta",results_folder+"/fasta/reference_proteome.fasta")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .fasta file not found at "+args.proteome+"/proteome.fasta.\nAborting program.")	
	try:
		shutil.copy(args.proteome+"/proteome.bed",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome .bed file not found at "+args.proteome+"/proteome.bed.\nAborting program.")
	# Not sure if these four should be quite so important.
	try:
		shutil.copy(args.proteome+"/proteome-descriptions.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene descriptions file not found at "+args.proteome+"/proteome-descriptions.txt.\nAborting program.")
	try:
		shutil.copy(args.proteome+"/proteome-genes.txt",results_folder+"/log/")
	except IOError:
		raise SystemExit("ERROR: Reference proteome gene names file not found at "+args.proteome+"/proteome-genes.txt.\nAborting program.")
	if args.junction:
		try:
			shutil.copy(args.proteome+"/proteome-descriptions.txt",results_folder+"/log/alternative-descriptions.txt")
		except IOError:
			raise SystemExit("ERROR: Reference proteome gene descriptions file not found at "+args.proteome+"/proteome-descriptions.txt.\nAborting program.")
		try:
			shutil.copy(args.proteome+"/proteome-genes.txt",results_folder+"/log/alternative-genes.txt")
		except IOError:
			raise SystemExit("ERROR: Reference proteome gene names file not found at "+args.proteome+"/proteome-genes.txt.\nAborting program.")
	
	return (results_folder, logfile, referrorfile, statusfile)


###############################################################################
############################## fProcessVariants ###############################
###############################################################################
def fProcessVariants (args, script_dir, results_folder, logfile, ref_prot, LvariantFiles = []):

	# Call sREAD_CHR_BED_EEXECUTABLE, which takes the reference genome and proteome.bed file as input and produces
	# a fasta file of exomes.

	#KRC 1/24/2024 No need to make a new local bed.dna file for every sample, but do so if lacking a central copy
	sPreprocessedProteome = os.path.join(args.proteome, 'proteome.bed.dna')
	if os.path.isfile(sPreprocessedProteome):
		#shutil.copy(sPreprocessedProteome, os.path.join(results_folder, 'log', 'proteome.bed.dna'))
		#print ('Copied proteome.bed.dna file from:', sPreprocessedProteome)
		sProteomeBedDNA = sPreprocessedProteome
	else:
		print ('Making proteome.bed.dna file local for each input file.')
		print ('Skip this step by using pyQuilts_KRC/prepare_proteome script to make it.')
		sProteomeBed = results_folder + '/log/proteome.bed'
		sProteomeBedDNA = sProteomeBed + '.dna'
		try:
			#check_call("%s/read_chr_bed %s/log/proteome.bed %s" % (script_dir, results_folder, args.genome), shell=True)
			#check_call("%s%s %s/log/proteome.bed %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, results_folder, args.genome), shell=True)
			check_call("%s%s %s %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, sProteomeBed, args.genome), shell=True)
		except CalledProcessError:
			raise SystemExit("ERROR: read_chr_bed didn't work - now we don't have a proteome.bed.dna file. Try recompiling read_chr_bed.c.\nAborting program.")
		print ('Made proteome.bed.dna file.')
	# Commented the above out for speed - it's slow, so for current testing purposes I'm just copying it from elsewhere
	#shutil.copy('/ifs/data/proteomics/tcga/scripts/quilts/pyquilts/proteome.bed.dna', results_folder+"/log/")

	# Time to merge and quality-threshold the variant files!
	if args.somatic:
		#som_flag = qu.merge_and_qual_filter(args.somatic, args.variant_quality_threshold)
		som_flag = qu.merge_and_qual_filter_vcf(LvariantFiles, args.somatic, results_folder, args.variant_quality_threshold) #KRC 1/22/2024 separate directories for inputs and intermediates
		if som_flag:
			args.somatic = None
	if args.germline:
		#germ_flag = qu.merge_and_qual_filter(args.germline, args.variant_quality_threshold)
		germ_flag = qu.merge_and_qual_filter_vcf(LvariantFiles, args.germline, results_folder, args.variant_quality_threshold) #KRC 1/22/2024 separate directories for inputs and intermediates
		if germ_flag:
			args.germline = None
	qu.quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
	qu.write_to_status("Merge and qual filter finished")
		
	# Now let's remove everything in the somatic file that is duplicated in the germline file
	# since if it shows up in both, it's a germline variant.
	if args.somatic and args.germline:
		qu.remove_somatic_duplicates(args.germline, args.somatic)
	qu.write_to_status("Somatic duplicates removed")

	# Next, create a proteome.bed file containing only variants...probably.
	# I could probably combine them more prettily, but for now I'll just concatenate the files.
	if args.somatic:
		#qu.get_variants(args.somatic+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S")
		qu.get_variants(results_folder+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "S") #KRC 1/22/2024 separate directories for inputs and intermediates
	if args.germline:
		#qu.get_variants(args.germline+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G")
		qu.get_variants(results_folder+"/merged_pytest/merged.vcf", results_folder+"/log/proteome.bed", "G") #KRC 1/22/2024 separate directories for inputs and intermediates

	qu.write_to_status("Get variants completed, proteome.bed file written")
	
	# Combine them, not very prettily:
	if args.somatic and args.germline:
		dest = open(results_folder+"/log/proteome.bed.var",'w')
		for f in [results_folder+"/log/proteome.bed.S.var",results_folder+"/log/proteome.bed.G.var"]:
			with open(f,'r') as src:
				shutil.copyfileobj(src, dest)
		dest.close()
	elif args.somatic:
		shutil.copy(results_folder+"/log/proteome.bed.S.var", results_folder+"/log/proteome.bed.var")
	elif args.germline:
		shutil.copy(results_folder+"/log/proteome.bed.G.var", results_folder+"/log/proteome.bed.var")
	qu.write_to_status("Somatic and germline combined")

	# Finish out the variants, if there are any	
	if args.somatic or args.germline:
		# Combine (some more) and sort variants.
		# Also variants are sorted by what they do to the sequence (add/remove stop, change AA, etc)
		# I will continue to do this, probably? I'll just ignore more later on
		dfproteomeBed = pd.read_table(results_folder+"/log/proteome.bed", delimiter= '\t', comment='#', header=0 ) #for .csv use delimiter = ','
		print('Sorting variants amongst',  dfproteomeBed.shape[0], 'genes in\n', sProteomeBedDNA)

		#KRC 1/28/2024 allow separate location of reference proteome and intermediate personalized files
		#qu.sort_variants(results_folder+"/log/proteome.bed.dna", results_folder+"/log/proteome.bed.var", "proteome", ref_prot)
		qu.sort_variants(sProteomeBedDNA, results_folder+'/log/proteome.bed.var', "proteome", ref_prot, results_folder+'/log/')

		qu.write_to_status("Variants sorted")

		# Translate the variant sequences into a fasta file.
		#translate(results_folder+"/log/", "proteome.aa.var.bed.dna", logfile, 'aa')
		qu.translate_saavs(results_folder+"/log/", "proteome.aa.var.bed.dna", logfile, ref_prot)
		qu.write_to_status("Translated SAAVs")
		qu.translate(results_folder+"/log/", "proteome.indel.var.bed.dna", logfile, 'indel')
		qu.write_to_status("Translated indels")
		shutil.copy(results_folder+"/log/proteome.aa.var.bed.dna.fasta", results_folder+"/fasta/parts/proteome.aa.fasta")
		shutil.copy(results_folder+"/log/proteome.indel.var.bed.dna.fasta", results_folder+"/fasta/parts/proteome.indel.fasta")

	return(0)

###############################################################################
############################# fProcessJunctions ###############################
###############################################################################
def fProcessJunctions (args, script_dir, results_folder, logfile, ref_prot, LjuncFiles = []):

	qu.write_to_status("Starting on junctions now.")
	junc_flag = None
	# If MapSplice was used instead of Tophat, copy its junctions.txt file to 
	# junctions.bed so the merge_junction_files function will pick up on it
	if args.junction_file_type == 'tophat' or args.junction_file_type == 'star':
		if args.junction_file_type == 'tophat':
			junc_flag = qu.convert_tophat_to_mapsplice(args.junction)
		else:
			#junc_flag = qu.convert_star_to_mapsplice(args.junction)
			junc_flag = fConvertStarToMapsplice(args.junction, results_folder+"/junction", LjuncFiles) #KRC 3/30/2023 /junction created in fSetUpGlobalPaths
		if junc_flag:
			args.junction = None
			qu.write_to_log("Could not convert your junction file - check junction file location and suffix", logfile)
			warnings.warn("Could not convert your junction file - check junction file location and suffix. Skipping..." % args.junction)
			qu.quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file

	# Merge junction files found in junction folder.
	if not junc_flag:
		if args.junction_file_type == 'star':
			junc_flag = qu.merge_junction_files(results_folder+"/junction", results_folder+'/log') #KRC 3/30/2023 input is now within the intermediate folder results_folder+"/junction"
		else:
			junc_flag = qu.merge_junction_files(args.junction, results_folder+'/log')
	if junc_flag:
		args.junction = None
		qu.quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
	else:
		qu.filter_known_transcripts(args.proteome+'/transcriptome.bed', results_folder+'/log', logfile)
		qu.filter_alternative_splices(results_folder+'/log/', args.threshB, args.threshD, args.threshN, logfile, ref_prot)
		qu.write_to_status("Filtered alternative splices into the appropriate types.")
	
		# Make a fasta out of the alternative splices with conserved exon boundaries
		qu.write_to_status("About to do a read_chr_bed")
		try:
			check_call("%s%s %s/log/merged-junctions.alt.filtered.bed %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE,  results_folder, args.genome), shell=True)
			# Don't know why this copies instead of moving. If I never use merged-junctions.filter.A.bed.dna again, just move it or have read_chr_bed output the alternative.bed.dna file instead.
			shutil.copy(results_folder+'/log/merged-junctions.alt.filtered.bed.dna', results_folder+'/log/alternative.bed.dna')
		except CalledProcessError:
			warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.alt.filtered.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries. Try recompiling read_chr_bed.c.")
		qu.write_to_status("Done with read_chr_bed to create alternative.bed.dna")
		
		# Make a fasta out of the alternative splices with conserved donor boundaries
		qu.write_to_status("About to do a read_chr_bed")
		try:
			check_call("%s%s %s/log/merged-junctions.donor.filtered.bed %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, results_folder, args.genome), shell=True)
			# Don't know why this copies instead of moving.
			shutil.copy(results_folder+'/log/merged-junctions.donor.filtered.bed.dna', results_folder+'/log/donor.bed.dna')
		except CalledProcessError:
			warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.donor.filtered.bed.dna file. Will not have a fasta file of alternative splices with conserved exon boundaries. Try recompiling read_chr_bed.c.")
		qu.write_to_status("Done with read_chr_bed to create donor.bed.dna")
		
		# Now to tackle the novels...
		if not args.junction_skipNovels == True:
			qu.write_to_status("About to do a read_chr_bed")
			try:
				check_call("%s%s %s/log/merged-junctions.novel.filtered.bed %s" % (script_dir, sREAD_CHR_BED_EEXECUTABLE, results_folder, args.genome), shell=True)
				shutil.copy(results_folder+'/log/merged-junctions.novel.filtered.bed.dna', results_folder+'/log/novel.bed.dna')
			except CalledProcessError:
				warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a merged-junctions.novel.filtered.bed.dna file. Will not have a fasta file of novel spliceforms. Try recompiling read_chr_bed.c.")
			qu.write_to_status("Done with read_chr_bed to create novel.bed.dna")

		# Translating the junctions. Looks like it requires a slightly different function than the old translation function.
		# Basically, though, can use the indel translation function (keep the exon with the variant and everything that comes after) for all of them except the ones where the exon boundaries are both new, in which case we need to do all six reading frames. And we saved those...where? notA?
		# Single frame translations
		qu.translate(results_folder+"/log/", "alternative.bed.dna", logfile, 'juncA')
		qu.translate(results_folder+"/log/", "donor.bed.dna", logfile, 'juncAN')

		# Six-frame translations
		if not args.junction_skipNovels == True:
			qu.translate_novels(results_folder+"/log/", "novel.bed.dna", logfile)
	
		# Move junctions to results folder
		shutil.copy(results_folder+"/log/alternative.bed.dna.fasta", results_folder+"/fasta/parts/proteome.alternative_splices.fasta")
		with open(results_folder+"/log/donor.bed.dna.fasta", 'r') as src:
			with open(results_folder+"/fasta/parts/proteome.alternative_splices.fasta",'a', newline='\n') as dest: #write unix LF when running on Windows KRC 7/26/2022
				shutil.copyfileobj(src, dest)
		if not args.junction_skipNovels == True:
			shutil.copy(results_folder+"/log/novel.bed.dna.fasta", results_folder+"/fasta/parts/proteome.novel_splices.fasta")
		qu.write_to_status("Translated alternative-splice junctions")

	return(0)
###############################################################################
############################## fProcessFusions ################################
###############################################################################
def fProcessFusions(args, script_dir, results_folder):

	qu.write_to_status("Starting on fusions now.")
	fusion_flag = qu.merge_and_filter_fusion_files(args.fusion, results_folder)
	if fusion_flag:
		args.fusion = None
		qu.quit_if_no_variant_files(args) # Check to make sure we still have at least one variant file
	else:
		#fusion time
		qu.create_fusion_bed(results_folder)
		qu.write_to_status("About to do a read_chr_bed")
		try:
			check_call("%ssREAD_CHR_BED_EEXECUTABLE %s/log/fusions.bed %s" % (script_dir, results_folder, args.genome), shell=True)
		except CalledProcessError:
			warnings.warn("WARNING: read_chr_bed didn't work - now we don't have a fusions.bed.dna file. Will not be able to perform fusions. Try recompiling read_chr_bed.c.")
		qu.write_to_status("Done with read_chr_bed to fusions.bed.dna")
		qu.translate_fusions(results_folder)
		shutil.copy(results_folder+"/log/fusions.fasta", results_folder+"/fasta/parts/proteome.fusions.fasta")

	return(0)


###############################################################################
###################### fConcatenateIntermediateFastas #########################
###############################################################################
### replaces quilts function: combine_output_fastas
### so that the final result is not buried in an individual directory tree with the intermediate files
def fConcatenateIntermediateFastas(sIntermediateDir, sOutputFinalDir, sPrefixFinalOutputFile):
	'''Combines all of the fastas in the output directory into one'''
	files = os.listdir( sIntermediateDir +'/fasta/parts')
	for f in files:
		if f.endswith('.fasta'):
			with open(sIntermediateDir + '/fasta/parts/' + f, 'r') as src:
				with open(os.path.join(sOutputFinalDir, sPrefixFinalOutputFile + '_variant_proteome.fasta'), 'a', newline='\n') as dest: #write unix LF when running on Windows KRC 7/26/2022
					shutil.copyfileobj(src, dest)
	return(0)

###############################################################################
############################# fListFilesOfType ################################
###############################################################################
### replaces quilts functions: pull_*_files
###		so that a file suffix sType is an argument:
###		input .gz files can be used
def fListFilesOfType(junc_dir, sType):
	'''Finds all *.tab files in the directory.'''
	files = os.listdir(junc_dir)
	junc_files = []
	for f in files:
		#if f.endswith('.tab'):
		if f.endswith(sType):
			junc_files.append(f)
	return junc_files

###############################################################################
########################## fConvertStarToMapsplice ############################
###############################################################################
### replaces quilts function: convert_star_to_mapsplice
### so that:
###		output can be located with the intermediate files instead of the input files
###		input .gz files can be used w/o prior unzipping
#def fConvertStarToMapsplice(junc_dir):
def fConvertStarToMapsplice(sInDir, sOutDir, junc_files = []):
	'''Converts STAR junction files into mapsplice format for further processing'''
	if not os.path.isdir(sInDir):
		# No junction files are going to be found. Gotta leave the function.
		warnings.warn("Unable to open junction folder %s. Giving up on finding splice junctions." % sInDir)
		return 1

	w = open(sOutDir+'/junctions.txt','w', newline='\n') #write unix LF when running on Windows KRC 7/23/2022
	count = 0
	for fil in junc_files:
		if not fil.endswith('.gz'):
			f = open(sInDir + '/' + fil, 'r')
			Llines = f.read().splitlines()	#list of lines with \n removed
		else:
			with gzip.open(sInDir + '/' + fil, 'rt') as f:
				Llines = f.read().splitlines()	#list of lines with \n removed

# 		f = open(sInDir+'/'+fil,'r')
# 		for line in f.readlines():
		for line in Llines:
			spline = line.split()
			chr = spline[0]
			#STAR defines the junction start/end as intronic bases, while many other software define them as exonic bases.
			#
			start = int(spline[1])-1	#first base of the intron (1-based), change to last of exon
			end = int(spline[2])+1		#last base of the intron (1-based), change to first of exon
			if spline[3] == '2':		#strand (0: undefined, 1: +, 2: -)
				strand = '-'
			else:
				strand = '+'
			junc_name = "JUNC_%d" % count
			count += 1
			reads = spline[6]
			new_len = end-start+50+1
			w.write("%s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t255,0,0\t2\t50,50\t0,%d\n" % (chr, start, end, junc_name, reads, strand, start, end, new_len))
		f.close()
	w.close()

	return(0)

def draftGzip(sDataPathFrom):
	if not sDataPathFrom.endswith('.gz'):
		print('Hi')
	else:
		gzip.open(sDataPathFrom, 'rb')
###############################################################################
############################### fMainSingle ###################################
###############################################################################

# Quilts original flow after	Split code into functions by type junctions, variants, fusions
def fMainSingle(args): #args is an argparse.Namespace object
	# Parse input, make sure we have at least one input type provided
	args = qu.parse_input_arguments()

	#script_dir = os.path.dirname(os.path.realpath(__file__)) # can this really be the best way to do this!?
	script_dir = '' #default is to run from within quilts dir KRC 7/24/2022

	# Set up log/status files
	output_dir = args.output_dir
	(results_folder, logfile, referrorfile, statusfile) = fSetUpGlobalPaths(output_dir, args)
	qu.results_folder	= results_folder	#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?
	qu.statusfile		= statusfile
	qu.referrorfile		= referrorfile
	qu.logfile			= logfile			#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?

	qu.write_to_status("Started")
	qu.write_to_log("Version Python.0", logfile)
	qu.write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)
	
	# Set up codon map
	qu.codon_map = CODON_MAP
	
	# Prep a map of the reference proteome for quality checks
	ref_prot = qu.save_ref_prot(results_folder+"/fasta/reference_proteome.fasta")

	#process the inputs
	if args.somatic or args.germline:
		fProcessVariants(args, script_dir, results_folder, logfile)

	if args.junction:
		fProcessJunctions (args, script_dir, results_folder, logfile, ref_prot)

	if args.fusion:
		fProcessFusions(args, script_dir, results_folder)


	if args.germline or args.somatic or args.junction:
		qu.tabulate_ref_mismatches()

	# Combine all of the proteome parts into one.
	#qu.combine_output_fastas(results_folder+'/fasta/')
	sPrefixFinalOutputFile = 'test'
	fConcatenateIntermediateFastas(results_folder, output_dir, sPrefixFinalOutputFile)
	qu.write_to_status("DONE")

	return(0)
###############################################################################
############################# fManifestDestiny ################################
###############################################################################
#given an input filename and the dataframe for a sample manifest return a value from another column in the sample manifest
def fManifestDestiny(dfSM, sInputFile, sInputFileColumnHeader, sDestinyColumnHeader):

	if not sInputFile in dfSM[sInputFileColumnHeader].values:
		return(sInputFile)

	#retrieve the alternate sample name from the manifest for the given RNA-seq junction data filename
	sDestinyValue = dfSM.loc[dfSM[sInputFileColumnHeader] == sInputFile][sDestinyColumnHeader]
	if isinstance(sDestinyValue, pd.Series):
		#print('Sample Manifest Alternate name:', sDestinyValue , 'Series of length:', len(sDestinyValue))
		if len(sDestinyValue) == 1:
			sDestinyValue = sDestinyValue.iloc[0] #get the first element of the series by position, the index will be it's row number in the dataframe
		else:
			print('Error Multiple files found in manifest:', sDestinyValue)
			sys.exit(0)

	return(str(sDestinyValue))
###############################################################################
################################ fMainMulti ###################################
###############################################################################
# Quilts revised flow: multiple input files yielding multiple output files 1-1 correspondence
# sWorkingDir should be P:\Projects\QUILTS\pyQUILTS_KRC (drive letter may be different)
#		intendd to be from QUILTS\processing-scripts\QuiltsPreprocessWrapper.py
def fMainMulti(args, sSampleManifestPath = '', sWorkingDir = '', iNumProcessesorsToUse = 1):

	sSaveCWD = os.getcwd()
	if sWorkingDir != '':
		os.chdir(sWorkingDir)

	# Check that we have a somatic and/or germline and/or junction file. Abort if not.
	if not args.somatic and not args.germline and not args.junction and not args.fusion:
		sys.exit("ERROR: Must have at least one variant file (somatic, germline, fusion, and/or junctions).\nAborting program.")

	#script_dir = os.path.dirname(os.path.realpath(__file__)) # can this really be the best way to do this!?
	script_dir = '' #default is to run from within quilts dir KRC 7/24/2022

	# Set up codon map
	qu.codon_map = CODON_MAP

	#get list of input files
	if not args.sample_manifest_junction == '':
		dfSMJ = pd.read_table(args.sample_manifest_junction, delimiter= ',', comment='#', header=0 ) #for .csv use delimiter = ','
	elif args.junction and args.junction_file_type == 'star':
		#allow .tab and .tab.gz
		#LinputFiles = [os.path.join(args.junction, f) for f in os.listdir(args.junction) if '.tab' in f] #full path
		LinputFiles = fListFilesOfType(args.junction, '.tab')		#just the files
		LinputFiles += fListFilesOfType(args.junction, '.tab.gz')	#just the files
		print('\n', LinputFiles)
		#if len(LinputFiles) ==0:
		#	LinputFiles = fListFilesOfType(args.junction, '.tab.gz')
		if len(LinputFiles) > 0:
			print ('STAR junction files to process:', len(LinputFiles))
		else:
			print ('Error could not find STAR junction files of type .tab or .tab.gz')
	elif args.somatic:
		LinputFiles = fListFilesOfType(args.somatic, '.vcf')
		#if len(LinputFiles) > 0:
		print (len(LinputFiles), 'Somatic variant vcf files to process')
	elif args.germline:
		LinputFiles = fListFilesOfType(args.germline, '.vcf')
		#if len(LinputFiles) > 0:
		print (len(LinputFiles), 'Germline variant vcf files to process')

	if iNumProcessesorsToUse > 1:
		fQueue_RunQUILTSmulti (args, LinputFiles, iNumProcessesorsToUse)
	else:
		#run Quilts once for each input file
		#the output Location is separate from all the intermediate file Quilts creates, so the intermediate files can be easily deleted
		iFileCount = 0
		for sInputFile in LinputFiles:
			iFileCount += 1
			print('\n', iFileCount, 'Processing: ', sInputFile, )
			fRunQUILTSsingle(args, sInputFile, iFileCount, script_dir = script_dir)
			"""
			sInputFilePrefix = os.path.splitext(sInputFile)[0]
			sOutputFilePrefix = sInputFilePrefix
			if not args.sample_manifest_junction == '':
				sOutputFilePrefix = fManifestDestiny(dfSMJ, sInputFile, 'junctionFilename', 'QUILTSprefix') #get an alternate sample name for the output
			print('To output prefix: ', sOutputFilePrefix)

			# Create a subdir based on sInputFilePrefix + date/time to become the intermediate results directory, and Set up log/status files within
			output_dir = args.output_dir
			(results_folder, logfile, referrorfile, statusfile) = fSetUpGlobalPaths(output_dir, args, sInputFilePrefix)
			qu.results_folder	= results_folder	#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?
			qu.statusfile		= statusfile
			qu.referrorfile		= referrorfile
			qu.logfile			= logfile			#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?
	
			qu.write_to_status("Started")
			qu.write_to_log("Version Python.0", logfile)
			qu.write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)

			# Prep a map of the reference proteome for quality checks
			bTroubleshoot = False
			if iFileCount == 1:
				bTroubleshoot = True
			ref_prot = qu.save_ref_prot(results_folder+"/fasta/reference_proteome.fasta", bTroubleshoot)
	
			#process the inputs
			if args.somatic or args.germline:
				fProcessVariants(args, script_dir, results_folder, logfile, ref_prot, [sInputFile])
	
			if args.junction:
				fProcessJunctions (args, script_dir, results_folder, logfile, ref_prot, [sInputFile]) #KRC 3/30/2023 create a list with only 1 file in it
	
			if args.fusion:
				fProcessFusions(args, script_dir, results_folder)
	
	
			if args.germline or args.somatic or args.junction:
				qu.tabulate_ref_mismatches()
	
			# Combine all of the proteome parts into one.
			fConcatenateIntermediateFastas(results_folder, output_dir, sOutputFilePrefix)
			qu.write_to_status("DONE")
			"""

	print('Finished all QUILTS processing')
	os.chdir(sSaveCWD)

	return(0)

##############################################################################
############################## fRunQUILTSsingle ##############################
##############################################################################
def fRunQUILTSsingle(args, sInputFile, iFileCount, script_dir = ''):

	sInputFilePrefix = os.path.splitext(sInputFile)[0]
	sOutputFilePrefix = sInputFilePrefix
	if not args.sample_manifest_junction == '':
		sOutputFilePrefix = fManifestDestiny(dfSMJ, sInputFile, 'junctionFilename', 'QUILTSprefix') #get an alternate sample name for the output
	print('To output prefix: ', sOutputFilePrefix)

	# Create a subdir based on sInputFilePrefix + date/time to become the intermediate results directory, and Set up log/status files within
	output_dir = args.output_dir
	(results_folder, logfile, referrorfile, statusfile) = fSetUpGlobalPaths(output_dir, args, sInputFilePrefix)
	qu.results_folder	= results_folder	#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?
	qu.statusfile		= statusfile
	qu.referrorfile		= referrorfile
	qu.logfile			= logfile			#KRC 3/30/2023 In Quilts I think a local variable is mostly/always passed?
	
	qu.write_to_status("Started")
	qu.write_to_log("Version Python.0", logfile)
	qu.write_to_log("Reference DB used: "+args.proteome.split("/")[-1], logfile)

	# Prep a map of the reference proteome for quality checks
	bTroubleshoot = False
	if iFileCount == 1:
		bTroubleshoot = True
	ref_prot = qu.save_ref_prot(results_folder+"/fasta/reference_proteome.fasta", bTroubleshoot)
	
	qu.codon_map = CODON_MAP
	#process the inputs
	if args.somatic or args.germline:
		fProcessVariants(args, script_dir, results_folder, logfile, ref_prot, [sInputFile])
	
	if args.junction:
		fProcessJunctions (args, script_dir, results_folder, logfile, ref_prot, [sInputFile]) #KRC 3/30/2023 create a list with only 1 file in it
	
	if args.fusion:
		fProcessFusions(args, script_dir, results_folder)
	
	
	if args.germline or args.somatic or args.junction:
		qu.tabulate_ref_mismatches()
	
	# Combine all of the proteome parts into one.
	fConcatenateIntermediateFastas(results_folder, output_dir, sOutputFilePrefix)
	qu.write_to_status("DONE")

	return(0);

##############################################################################
########################### fWorker_RunQUILTSsingle ##########################
##############################################################################
def fWorker_RunQUILTSsingle(Qtasks):
	"""
	Worker function that retrieves commands from the queue and executes them.
	"""
	while True:
		(iProcess_id, args, sInputFile) = Qtasks.get() # tuple of args for the subprocess
		if iProcess_id is None:  # Sentinel value to signal termination
			break
		try:
			print(f'Processing: {iProcess_id} {sInputFile}')
			fRunQUILTSsingle(args, sInputFile, iProcess_id)
		except Exception as e:
			print(f"Process {os.getpid()} encountered an error with {sInputFile}: {e}")
		#wait to check for free processors
		time.sleep(1) # Check every 1 seconds

##############################################################################
########################### fQueue_RunQUILTSmulti ############################
##############################################################################
def fQueue_RunQUILTSmulti (args, LinputFiles, iNumProcesses):

	m_processors = multiprocessing.cpu_count()  # Number of available processors
	print(m_processors, ' CPU available. For this queue using a max of: ',  iNumProcesses)
	m_processors = iNumProcesses
	
	# Create a queue to hold the commands
	print(f'Files to add to Queue_RunQUILTSmulti:  {len(LinputFiles)}')
	Qtasks = multiprocessing.Queue()

	#fill the queue for the input files
	iProcess_id = 1
	for sInputFile in LinputFiles:
		print(f'Queueing: {iProcess_id} {sInputFile}')
		# num args to put here must be same as number put below for sentinel values
		Qtasks.put((iProcess_id, args, sInputFile)) # tuple of args for the worker function
		iProcess_id += 1
	print(f'Total num samples queued: {iProcess_id}' )

	# Create and start worker 1st batch of processes using limited num: m_processors
	Lprocesses = []
	for _ in range(m_processors):
		p = multiprocessing.Process(target=fWorker_RunQUILTSsingle, args=(Qtasks,))
		Lprocesses.append(p)
		p.start()

	# Add sentinel values to the queue to signal workers to terminate
	for _ in range(m_processors):
		# num args to put here must be same as number put above for data
		Qtasks.put((None, None,None))

	# Wait for all processes to complete
	for p in Lprocesses:
		p.join()

	print("All commands executed.")

###############################################################################
########################## fFill_ArgparseNamespace ############################
###############################################################################
# Provides a programmatic alternative to CLI of providing the input arguments
# to bypass native QUILTS function parse_input_arguments().
# Handles defaults for options not provided as input.
def fFill_ArgparseNamespace( Doptions):

	# Check that we have a somatic and/or germline and/or junction file. Abort if not.
	if 'somatic' not in Doptions and 'junction' not in Doptions and 'germline' not in Doptions and 'fusion' not in Doptions:
		sys.exit("ERROR: Options must specify at least one variant directory (somatic, germline, junctions, and/or fusion). fFill_ArgparseNamespace() \nAborting.")
		# parser.add_argument('--germline',type=str, help="full path to folder containing germline variant VCF file(s)")
		# parser.add_argument('--somatic',type=str, help="full path to folder containing somatic variant VCF file(s)")
		# parser.add_argument('--junction', type=str, help="full path to folder containing junction file(s)")
	if 'output_dir' not in Doptions:
		sys.exit("ERROR: Options must specify output_dir full path to output folder (defaults to .). fFill_ArgparseNamespace() \nAborting.")
		# parser.add_argument('--output_dir', type=str, default=".", help="full path to output folder (defaults to .)")
	#if 'junction' in Doptions and 'junction_file_type' not in Doptions:
	#	sys.exit("ERROR: Options must specify a supported junction_file_type (star, mapsplice, or tophat). fFill_ArgparseNamespace() \nAborting.")
	#set up defaults
	DoptionsDefault = { 'junction' : '', 'somatic' : '', 'germline' : '', 'fusion' : '',
										# relative path to folder containing reference proteome
		'proteome'					:	r'..\references\proteome\Gencode42\\', #need trailing slash
										# relative path to folder containing reference genome
		'genome'					:	r'..\references\genome\ensembl_hg38_KRC\\', #need trailing slash
		'variant_quality_threshold'	: 0.0,	# Quality threshold for SAAVs somatic or germline (default=0.0)
		'junction_skipNovels'		: True,	# type=bool, if true generate only junction types bothConserved and donorConserved 
		'junction_file_type'		: 'star', #star, mapsplice, or tophat supported 
		'threshB'					:	2,	# Minimum # reads for junctions with both exon boundaries annotated.
		'threshD'					:	3,	# Minimum # reads for junctions with only donor exon boundary annotated.
		'threshN'					:	3,	# Minimum # reads for junctions, novel, without donor exon boundary annotated.
		'sample_manifest_junction'	:	""	# sample manifest for mapping input file to output file
	}
	# parser.add_argument('--sample_manifest_junction', type=str, default="", help="sample manifest for mapping input file to output file")
	# parser.add_argument('--no_missed_cleavage', action='store_true', default=False, help="Tryptic peptide fasta by default allows for a single missed cleavage; adding this argument will tell the virtual trypsinizer to assume perfect cleavage")
	# # Pull out the arguments
	# args = parser.parse_args()

	#Overwrite defaults with input options for input keys present (allows defaults to not be provided).
	Doptions = DoptionsDefault | Doptions

	#Create a new Namespace object for a dict
	args = argparse.Namespace(**Doptions)

	return(args)

	"""Unused namespace creation and modification operations left here for future reference
	#argparse.Namespace object
	# Manually create the Namespace object with desired attribute values
	Anamespace = argparse.Namespace(test_option='my_value', positional_arg=123)

	# Add new items to namespace manually
	args.new_item = 42
	args.another_item = "hello"

	# Update the existing namespace's internal dictionary
	# Get the dictionary and add a new key-value pair
	Dargs = vars(args)
	Dargs['new_item'] = 42
	print(f"After adding new items: {args}")
	# You can still access it using dot notation:
	print(f"Accessing via dot notation: {args.new_item}")

	# Update the existing namespace's internal dictionary with a dictionary of new items
	DnewItems = {
		'new_item_1': 'value1',
		'new_item_2': 123
	}
	# Update the namespace's internal dictionary
	args.__dict__.update(DnewItems)
	"""

###############################################################################
################################# __main__ ####################################
###############################################################################
if __name__ == "__main__":
	# Parse input, make sure we have at least one input type provided
	args = qu.parse_input_arguments()

	#fMainSingle(args)
	fMainMulti(args)
	sys.exit(0)