#Filename:	gunzip_dir.py
#Purpose:	Replacement for gunzip_dir.pl that does not require the external executable gunzip
#Author:	Karl Clauser
#Created:	June 7, 2023

import os
import glob			#for reading list of files in directory
import argparse
import gzip
import shutil

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

	#compress by removing header line and converting all sequence lines into a single line

	return()

###############################################################################
################################# __main__ ####################################
###############################################################################
if __name__ == "__main__":

	# Pull out the CLI arguments
	parser = argparse.ArgumentParser(description="gunzip_dir.py")
	parser.add_argument('--dir', type=str, default="", help="Dir where .gz files are")
	args = parser.parse_args() # Pull out the arguments
	fGunzip(args)
	fCompress(args)

	sys.exit(0)