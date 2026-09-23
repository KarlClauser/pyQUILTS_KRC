#Filename:	gunzip_dir.py
#Purpose:	Replacement for gunzip_dir.pl/compress_chr_dir.pl that does not require the external executables gunzip, compress_chr
#Author:	Karl Clauser
#Created:	June 7, 2023

import os
import glob			#for reading list of files in directory
import argparse
import gzip
import shutil
import sys

###############################################################################
################################### fGunzip ###################################
###############################################################################
def fGunzip(args):
	#get all the .gz files in the directory and ungzip them
	sWildcardPath = args.dir + "/" + '*.gz'
	Lfiles = [f for f in glob.glob(sWildcardPath, recursive=False)]
	for sDataPathFrom in Lfiles:
		#sDataPathFrom	= os.path.join(args.dir, sFile)
		with gzip.open(sDataPathFrom, 'rb') as f_in:
			with open(sDataPathFrom.replace('.gz', ''), 'wb') as f_out:	#trim .gz off the filename of the unzip result
				shutil.copyfileobj(f_in, f_out)
		print ("Unzipped and copied: %s" % sDataPathFrom)
	return(0)
###############################################################################
################################## fCompress ##################################
###############################################################################
def fCompress(args):
	#get all the .fa files in the directory
	sWildcardPath = args.dir + "/" + '*.fa'
	Lfiles = [f for f in glob.glob(sWildcardPath, recursive=False)]

	#compress by removing header lines and converting all sequence lines into a single line
	#output is 1 file per fasta entry
	for sFile in Lfiles:
		with open(sFile, 'r') as file:
			FileOut = ''
			sPreviousFileOut = ''
			# Iterate over the lines of the file
			for sLine in file:		# Iterate over the lines of the file
				if sLine.startswith('>'): #header line
					#>MT dna:chromosome chromosome:GRCh37:MT:1:16569:1 REF
					sChr = sLine.split(' ')[0].lstrip('>') #get the first field before a space
					sFileOut = args.dir + '/' + 'chr' + sChr + '.fa.cmp1'
					if len(sPreviousFileOut) > 0:
						FileOut.close()					#close the previous one
					FileOut = open(sFileOut,'w')
					print ("Writing: %s" % sFileOut)
					sPreviousFileOut = sFileOut
				else:
					FileOut.write(sLine.strip().upper()) #no line ending tacked on

		print ("Compressed: %s" % sFile)

	return()
###############################################################################
################################## fGunzip_Compress ##################################
###############################################################################
def fGunzip_Compress(args):
	#get all the .gz files in the directory and ungzip them
	sWildcardPath = args.dir + "/" + '*.gz'
	Lfiles = [f for f in glob.glob(sWildcardPath, recursive=False)]

	#compress by removing header lines and converting all sequence lines into a single line
	#output is 1 file per fasta entry

	for sFile in Lfiles:
		print ("Reading: %s" % sFile)
		with gzip.open(sFile, 'rt') as FileIn:
			Llines = FileIn.read().splitlines()	#list of lines with \n removed
			FileOut = ''
			sPreviousFileOut = ''
			# Iterate over the lines of the file
			for sLine in Llines:		# Iterate over the lines of the file
				if sLine.startswith('>'): #header line
					#>MT dna:chromosome chromosome:GRCh37:MT:1:16569:1 REF
					sChr = sLine.split(' ')[0].lstrip('>') #get the first field before a space
					sFileOut = args.dir + '/' + 'chr' + sChr + '.fa.cmp1'
					if len(sPreviousFileOut) > 0:
						FileOut.close()					#close the previous one
					FileOut = open(sFileOut,'w')
					print ("Writing: %s" % sFileOut)
					sPreviousFileOut = sFileOut
				else:
					FileOut.write(sLine.upper()) #no line ending tacked on
			FileOut.close() #close the last one

		print ("Compressed: %s" % sFile)

	return()

###############################################################################
################################# __main__ ####################################
###############################################################################
if __name__ == "__main__":

	# Pull out the CLI arguments
	parser = argparse.ArgumentParser(description="gunzip_dir.py")
	parser.add_argument('--dir', type=str, default="", help="Dir where .gz files are")
	args = parser.parse_args() # Pull out the arguments
	#fGunzip(args)
	#fCompress(args)

	fGunzip_Compress(args) #skip the intermediate .fa

	sys.exit(0)