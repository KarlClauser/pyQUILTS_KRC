# pyQUILTS_KRC
This fork of Emily Kawaler's pyQUILTS had 2 major motivations
1. Run on Windows with Python 3 instead of v2.7
2. Restructure the inputs and outputs to be more friendly to processing cohorts
  instead of individual samples. Native Quilts expects a directory for each sample with inputs, outputs, and 
  intermediates buried within. For a cohort one would prefer to have an input directory containing all samples,
  an output directory containing  all samples, and the intermediates should be in a separate directory tree
  so they can easily be deleted. Thus, one need not spend time burying the inputs and unburying the outputs.
  Also one would like to easily process the data for each data type (somatic variants, germline variants,
  splice junctions) separately to deal with the issues of each data type without the need to re-process the 
  other data types each time. Quilts outputs are expected to be combined by a Spectrum Mill utility after running Quilts,
  rather than have Quilts do it.

Major Changes
  •	Revised python code to run with Python 3 instead of v2.7
  •	Revised to run on Windows
  •	Recompiled C program for Windows, read_chr_bed_win64.exe, write LF not CRLF
  •	Quilts.py, write LF not CRLF
  •	Rewrote prepare genome perl scripts into a python script
  •	Created QuiltsWrapper.py which replaces the quilts __main __ function to treat quilts.py as a library and call Emily’s processing functions
  •	Split out code by type: junctions, variants, fusion for separate processing of each data type into its own file
  •	Input data for all samples in one directory
  •	Output data for all samples in one directory
  •	Create intermediate results dir for running, then delete un-needed intermediate files when done.
  •	Added option: to Skip novel junctions (5’ junction not conserved) to eliminate later need to remove them.
  •	Added option: to support sample manifest for mapping input file to output file.
  •	Accelerated SAAV processing from 1hr to 1 min/patient by adding check for matching protein acc number before entering function process_gene.
  •	Diminished storage needed for inputs. Revised to read STAR .SJ.out.tab.gz files w/o ungzipping

