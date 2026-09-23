// read_chr_bed_win64.cpp

//Karl Clauser
//July 24, 2022
//Compile in Microsoft Visual Studio 2019

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

#include <stdio.h>
#include <string.h>

//#ifdef _MSC_VER
// 
//KRC 7/24/2022 Added preprocessor definition _CRT_SECURE_NO_WARNINGS to eliminate errors about strcpy, strcat, fopen 
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#include <stdlib.h>		//KRC added for atoi
#include <ctype.h>		//KRC added for isdigit
#include <stddef.h>		// KRC added for size_t
#include <sys/types.h>	// KRC added for ssize_t
#include <iostream>		// KRC added for getline
#include <sys/stat.h>
#include <fcntl.h>

//https://stackoverflow.com/questions/735126/are-there-alternate-implementations-of-gnu-getline-interface/735472#735472
#include <errno.h>
#include <stdint.h>

	  // if typedef doesn't exist (msvc, blah)
typedef intptr_t ssize_t;

ssize_t getline(char **lineptr, size_t *n, FILE *stream); //Windows version of getline
//#endif


int main ( int argc, char *argv[] )
{
	char filename[512];
	char filename_out[512];
	char chr[256];
	char chr_dir[2000];
	char name[2000];
	char qual[256];
	char strand[256];
	char lengths[200000];
	char offsets[200000];
	char buf[500000];

	//std::ifstream fi;
	//fi.open((LPCTSTR)szFilename);
	//if (!fi.is_open())
	//{
	//	return;
	//}

	//char szline[10240];
	//while (fi.getline(szline, 10240))

	if (argc>=2)
	{
		strcpy(filename,argv[1]);
		strcpy(filename_out,argv[1]);
		if (argc>=3) { strcpy(chr_dir,argv[2]); } else { strcpy(chr_dir,"/ifs/data/proteomics/tcga/databases/genome_human"); } 
		strcat(filename_out,".dna");
		FILE *fp = fopen (filename, "r");
		FILE *fop = fopen (filename_out, "wb"); //KRC 7/26/2022 b added so Unix LF are generated, otherwise Windows CRLF
		char * line = NULL;
		size_t len = 0;
		ssize_t read;
		while ((read = getline(&line, &len, fp)) != -1) 
		//while ((read = std::basic_istream::getline(&line, &len, fp)) != -1) 
		//while ((read = std::istream::getline(&line, &len, fp)) != -1) 
		//while ((read = std::cin.getline(&line, &len, fp)) != -1) 
		//while ((read = std::stdio::getline(&line, &len, fp)) != -1) 
		//while ( fgets(line, len, fp) != NULL ) 
		{
			int i=0;
			int j=0;
			for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { chr[j]=line[i]; }
			chr[j]='\0';
			if (line[i]!='\0')
			{
				i++;
				int start=atoi(line+i);
				for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
				if (line[i]!='\0')
				{
					i++;
					for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
					if (line[i]!='\0')
					{
						i++;
						for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { name[j]=line[i]; }
						name[j]='\0';
						if (line[i]!='\0')
						{
							i++;
							for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { qual[j]=line[i]; }
							qual[j]='\0';
							if (line[i]!='\0')
							{
								i++;
								for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { strand[j]=line[i]; }
								strand[j]='\0';
								if (line[i]!='\0')
								{
									i++;
									for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
									if (line[i]!='\0')
									{
										i++;
										for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
										if (line[i]!='\0')
										{
											i++;
											for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
											if (line[i]!='\0')
											{
												i++;
												for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { ; }
												if (line[i]!='\0')
												{
													i++;
													for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { lengths[j]=line[i]; }
													lengths[j]='\0';
													if (line[i]!='\0')
													{
														i++;
														for(j=0;line[i]!='\n' && line[i]!='\r' && line[i]!='\t' && line[i]!='\0';i++,j++) { offsets[j]=line[i]; }
														offsets[j]='\0';
														if (line[i]!='\0')
														{
															i++;
															fprintf(fop,"%s",line);
															fprintf(fop,">%s (MAP:%s:%d%s %s %s)\n",name,chr,start,strand,lengths,offsets);
															int k=0;
															int l=0;
															int m=0;
															char filename_chr[256];
															int padding=600;
															strcpy(filename_chr,chr_dir);
															//strcat(filename_chr,"/"); KRC 9/27/2022 now included in inputpath
															strcat(filename_chr,chr);
															strcat(filename_chr,".fa.cmp1");
															FILE *fp_chr = fopen(filename_chr, "r");
															if (fp_chr == NULL) {
																sprintf(filename_chr,"%schr%s.fa.cmp1",chr_dir,chr);
																//strcpy(filename_chr,chr_dir);
																//strcat(filename_chr,"chr"); //KRC 9/27/2022 handles when chr is is just the integer without chr prefix
																//strcat(filename_chr,chr);
																//strcat(filename_chr,".fa.cmp1");
																fp_chr = fopen(filename_chr, "r");
															}
															if (fp_chr!=NULL)
															{
																buf[0]='\0';
																fseek (fp_chr, start-padding, SEEK_SET);
																fgets(buf, padding+1, fp_chr);
																if (buf!=NULL)
																{
																	fprintf(fop,"%s\t-1\t%s\t%d\t%d\t%s\t0\t%s\n",name,chr,start-padding,padding,strand,buf);
																}
																int start_=0;
																int length=0;
																for(k=0,l=0,m=0;offsets[l]!='\0';k++)
																{
																	if (isdigit(offsets[l]))
																	{
																		start_=start+atoi(offsets+l);
																		if (isdigit(lengths[m]))
																		{
																			length=atoi(lengths+m);
																			buf[0]='\0';
																			fseek (fp_chr, start_, SEEK_SET);
																			if (length<500000)
																			{
																				fgets(buf, length+1, fp_chr);
																				if (buf!=NULL)
																				{
																					fprintf(fop,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\n",name,k,chr,start_,length,strand,qual,buf);
																				}
																			} else { fprintf(fop,"%s\t%d\t%s\t%d\t%d\t%s\t%s\tError Reading Sequence\n",name,k,chr,start_,length,strand,qual); }
																		}
																	}
																	for(;offsets[l]!=',' && offsets[l]!='\0';l++) { ; }
																	if (offsets[l]==',') { l++; }
																	for(;lengths[m]!=',' && lengths[m]!='\0';m++) { ; }
																	if (lengths[m]==',') { m++; }
																}
																buf[0]='\0';
																fseek (fp_chr, start_+length, SEEK_SET);
																fgets(buf, padding+1, fp_chr);
																if (buf!=NULL)
																{
																	fprintf(fop,"%s\t+1\t%s\t%d\t%d\t%s\t0\t%s\n",name,chr,start_+length,padding,strand,buf);
																}
																fclose(fp_chr);
															} else { printf("Error opening: %s\n",filename_chr); }
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return 0;
}

//#ifdef _MSC_VER

//https://stackoverflow.com/questions/735126/are-there-alternate-implementations-of-gnu-getline-interface/735472#735472

/* The original code is public domain -- Will Hartung 4/9/09 */
/* Modifications, public domain as well, by Antti Haapala, 11/10/17
- Switched to getc on 5/23/19 */

//KRC 7/26/2022 moved to top
//#include <stdio.h>
//#include <stdlib.h>
//#include <errno.h>
//#include <stdint.h>
//
//// if typedef doesn't exist (msvc, blah)
//typedef intptr_t ssize_t;

ssize_t getline(char **lineptr, size_t *n, FILE *stream) {
	size_t pos;
	int c;

	if (lineptr == NULL || stream == NULL || n == NULL) {
		errno = EINVAL;
		return -1;
	}

	c = getc(stream);
	if (c == EOF) {
		return -1;
	}

	if (*lineptr == NULL) {
		//*lineptr = malloc(128);	//KRC 7/26/2022 compiler error cannot convert void* to char*
		*lineptr = (char*)malloc(128); 
		if (*lineptr == NULL) {
			return -1;
		}
		*n = 128;
	}

	pos = 0;
	while(c != EOF) {
		if (pos + 1 >= *n) {
			size_t new_size = *n + (*n >> 2);
			if (new_size < 128) {
				new_size = 128;
			}
			//char *new_ptr = realloc(*lineptr, new_size); 	//KRC 7/26/2022 compiler error cannot convert void* to char*
			char *new_ptr = (char*)realloc(*lineptr, new_size);
			if (new_ptr == NULL) {
				return -1;
			}
			*n = new_size;
			*lineptr = new_ptr;
		}

		((unsigned char *)(*lineptr))[pos ++] = c;
		if (c == '\n') {
			break;
		}
		c = getc(stream);
	}

	(*lineptr)[pos] = '\0';
	return pos;
}
//#endif