#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

guideseq_v2n.py
===============

Rewrite of GuideSeq2 (from github.com/tsailabSJ/guideseq), altered to accept 
single-end sequencing such as that from Nanopore. Also attempts to compensate
for a lack of UMIs with aligned read size comparison.

"""

# import regular python classes
import os
import sys
import yaml
import argparse
import traceback
import subprocess
from tabulate import tabulate
import pandas as pd

# Set up logger
import log
logger = log.createCustomLogger('root',"guideseq_V2")

# import custom python scripts
#	If 'Nano' is in the name, changes have been made from tsailabsSJ version of GuideSeq2
from filterBackgroundSites import filterBackgroundSites,filterBlackList
from visualizationNano import visualizeOfftargets
import identifyOfftargetSitesNano as identifyOfftargetSites
import validation

# Color setup
GREEN = "\033[32m"
BLUE = "\033[34m"
YELLOW = "\033[33m"
RED = "\033[31m"
RESET = "\033[0m"

DEFAULT_YAML = os.path.dirname(os.path.realpath(__file__)) + "/default.yaml"

# Setup and parse command line arguments
def parse_args():
	parser = argparse.ArgumentParser()

	subparsers = parser.add_subparsers(description='Main run mode: main',
									   help='Use this to run individual steps of the pipeline',
									   dest='command')

	main_parser = subparsers.add_parser('main', help='Run a single step or a series of steps for one or all samples')
	main_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True)
	main_parser.add_argument('--sample', '-s', help='Specify sample to process (default is all)', default='all')
	main_parser.add_argument('--step', help='Specify steps: filter, align, identify, visualize, tableize, seperated by +', default='filter+align+identify+visualize+tableize')
	main_parser.add_argument('--overwrite', help='overwrite specifications in the yaml file', default=None,type=yaml.load)
	
	init_parser = subparsers.add_parser('init', help='initialize a default yaml in the current folder')

	return parser.parse_args()

# Retrieve settings from manifest. If a value isn't given, check default.yaml
def get_parameters(manifest_data):

	# init
	default_refseqName = os.path.dirname(os.path.realpath(__file__)) + "/refseq_gene_name.tsv"
	with open(DEFAULT_YAML, 'r') as f:
		default = yaml.safe_load(f)
	with open(manifest_data, 'r') as f:
		return_dict = yaml.safe_load(f) # this is user input YAML
	default['analysis_folder'] = os.getcwd()
	default['refseq_names'] = default_refseqName
	default['Manhattan.R'] = os.path.dirname(os.path.realpath(__file__)) + "/Manhattan.R"
	validation.validateManifest(return_dict)
	
	# assign default if user YAML missing some parameters
	return_dict = assign_default(return_dict,default)
	return_dict['control_coord_chr'],tmp = return_dict['control_primer_coord'].split(":")
	return_dict['control_coord_start'],return_dict['control_coord_end'] = [int(x) for x in tmp.split("-")]

	return return_dict

# Sub-routine for get_parameters()
def assign_default(return_dict,default):
	
	section_list = ['demultiplex']
	for k in section_list:
		if k in return_dict:
			for p in default[k]:
				if not p in return_dict[k]:
					return_dict[k][p] = default[k][p]
		else:
			return_dict[k] = default[k]

	# default global variables
	for p in default:
		if not p in section_list+['samples']:
			if not p in return_dict:
				return_dict[p] = default[p]
	
	# set abs path for input and output
	for k in section_list:
		return_dict[k]['out_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['out_dir'])
		return_dict[k]['input_dir'] = os.path.join(return_dict['analysis_folder'], return_dict[k]['input_dir'])

	return return_dict


class GuideSeq:

	def __init__(self):
		pass

	# Process the run manifest and display settings at command line
	def parseManifest(self, manifest_path,sample='all',user_dict=None):
		logger.info(f'{GREEN}Loading manifest...{RESET}')

		try:
			manifest_data = get_parameters(manifest_path)
			self.samples = {}
			self.parameters = manifest_data
			if user_dict == "None":
				user_dict = None
			if user_dict:
				for p in user_dict:
					self.parameters[p] = user_dict[p]
			if sample != "all":
				self.samples[sample] = manifest_data['samples'][sample]
			else:
				self.samples = manifest_data['samples']
			del self.parameters['samples']
			logger.info("\n"+tabulate([(k,v) for k,v in self.parameters.items()])) 
			logger.info("\n"+tabulate([(k,v) for k,v in self.samples[list(self.samples.keys())[0]].items()])) 

		except Exception as e:
			logger.error(
				f'{RED}Manifest initiation failed. Possibly due to incorrect or malformed manifest file. Please ensure your manifest contains all required fields.{RESET}')
			logger.error(traceback.format_exc())
			sys.exit()


	# Filter reads to prepare them for alignment. 
	#	- Removes RAPID sequences to allow for better alignment.
	#	- Aligns dsODN sequences on 5' end of read, as expected by downstream processes
	#	- Generates read count and length statistics for sample comparisons
	# Outputs a '[sampleName]_filtered.fastq' file in the 'filtured' directory for use by alignment step
	def filterReads(self):
		logger.info(f'{GREEN}Filtering reads to only those containing the dsODN sequence and trimming RAPID primer sequences...{RESET}')

		try:
			filtered_dir = os.path.join(self.parameters['analysis_folder'], 'filtered')
			os.makedirs(filtered_dir, exist_ok=True)

			seqkit = self.parameters['seqkit']
			RAPIDall = self.parameters['RAPIDall']
			RAPIDprimer = self.parameters['RAPIDprimer']
			dsODN1 = self.parameters['dsODN_primer']
			dsODN2 = self.parameters['dsODN_primer_revcomp']

			for sample in self.samples:
				try:
					logger.info(f"Filtering sample {sample}...")
					# Confirm the fastq file listed in the manifest 
					if 'read' not in self.samples[sample]:
						raise KeyError(f"'read' not found for sample '{sample}' in manifest. Please specify read path for each sample.")	

					# Confirm the fastq file exists before starting seqkit filtering.
					if not os.path.isfile(self.samples[sample]['read']):
						raise OSError(f"FASTQ file for sample '{sample}' not found. Please make sure manifest 'read:' parameters are correct.")

					# Naming of revised fastq files
					fastq_filter_1 = os.path.join(filtered_dir, sample + '_filterRAPIDwhole.fastq')
					fastq_filter_2 = os.path.join(filtered_dir, sample + '_filterRAPIDall.fastq')
					fastq_filter_3 = os.path.join(filtered_dir, sample + '_filterRAPIDall_dsODN1.fastq')
					fastq_filter_4 = os.path.join(filtered_dir, sample + '_filterRAPIDall_dsODN2.fastq')
					fastq_filter_5 = os.path.join(filtered_dir, sample + '_filtered.fastq')

					# Begin collecting statistics on each filtering step
					seqkit_output = []
					seqkit_output = subprocess.run(
						f"{seqkit} stats {self.samples[sample]['read']}", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n')

					# Remove full RAPID sequence from reads. 
					#   If the RAPID sequence is in multiple places in the read, this will result in multiple fragments being output.
					seqkit_output.append(subprocess.run(
						f"cat {self.samples[sample]['read']} | {seqkit} amplicon -j 15 -u -m 2 -F {RAPIDall} -r {str(len(RAPIDall)+1)}:99999 --quiet > {fastq_filter_1} && {seqkit} stats {fastq_filter_1} | tail -n 1", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n'))

					# Remove RAPID primer sequence from reads.
					#   This is to accommodate places where the RAPID primer bound nonspecifically, which are still usable by the process.
					seqkit_output.append(subprocess.run(
						f"cat {fastq_filter_1} | {seqkit} amplicon -j 15 -u -m 2 -F {RAPIDprimer} -r {str(len(RAPIDprimer)+1)}:99999 --quiet > {fastq_filter_2} && {seqkit} stats {fastq_filter_2} | tail -n 1", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n'))

					# Search for dsODN1 sequence. When found, trim and/or reverse-complement to put it on the 5' end of the read.
					#    Afterwards, remove any reads smaller than 70nt (including 34nt of dsODN1), as reads smaller than this won't align anyway.
					seqkit_output.append(subprocess.run(
						f"cat {fastq_filter_2} | {seqkit} amplicon -j 15 -m 0 -F {dsODN1} -r 1:99999 --quiet | seqkit seq -m 70 -g --quiet > {fastq_filter_3} && seqkit stats {fastq_filter_3} | tail -n 1", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n'))

					# Search for dsODN2 sequence. When found, trim and/or reverse-complement to put it on the 5' end of the read.
					#    Afterwards, remove any reads smaller than 70nt (including 34nt of dsODN2), as reads smaller than this won't align anyway.
					seqkit_output.append(subprocess.run(
						f"cat {fastq_filter_2} | {seqkit} amplicon -j 15 -m 0 -F {dsODN2} -r 1:99999 --quiet | seqkit seq -m 70 -g --quiet > {fastq_filter_4} && seqkit stats {fastq_filter_4} | tail -n 1", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n'))

					# Combine dsODN1 and dsODN2 fastq files. File to be passed to analyze step.
					seqkit_output.append(subprocess.run(
						f"cat {fastq_filter_3} {fastq_filter_4} > {fastq_filter_5} && {seqkit} stats {fastq_filter_5} | tail -n 1", shell=True, capture_output=True, text=True
					).stdout.strip().split('\n'))
					
					# Save path to filtered fastq file
					self.samples[sample]['filtered'] = os.path.join(filtered_dir, sample + '_filtered.fastq')

					# Convert strings to a clean table and output:
					seqkit_output.insert(1, ['---- ------ ---- -------- ------- ------- ------- -------'])
					seqkit_matrix = ["".join(row).split() for row in seqkit_output]
					seqkit_col_widths = [max(len(cell) for cell in col) for col in zip(*seqkit_matrix)]
					for row in seqkit_matrix:	    
					    format_str = "  ".join(f"{{:<{width}}}" for width in seqkit_col_widths)
					    logger.info(f'   {format_str.format(*row)}')
					print("")
				
				except:
					logger.error(f'{RED}Filtering failed for sample {sample}{RESET}')
					logger.error(traceback.format_exc())
				
		except Exception as e:
			logger.error(f'{RED}Error filtering reads{RESET}')
			logger.error(traceback.format_exc())
			quit()
		logger.info('Finished filtering reads.')


	# Align (filtered) reads against a reference genome
	#	- Uses minimap2 to align reads to genome
	#	- Remove read duplicates
	# Final output is in aligned directory with .dedup.sam extension. This file gets passed to the identify step.
	def alignReads(self):
		logger.info(f'{GREEN}Aligning reads and collapsing reads based on alignment and size...{RESET}')

		try:
			filtered_dir = os.path.join(self.parameters['analysis_folder'], 'filtered')
			aligned_dir = os.path.join(self.parameters['analysis_folder'], 'aligned')
			os.makedirs(aligned_dir, exist_ok=True)

			for sample in self.samples:
				try:
					# Check for filtered read file. (./filtered/[sample]_filtered.fastq)
					if 'filtered' not in self.samples[sample]:
						self.samples[sample]['filtered'] = os.path.join(filtered_dir, sample + '_filtered.fastq')
					# Confirm the fastq file exists before starting minimap2ing.
					if not os.path.isfile(self.samples[sample]['filtered']):
						raise OSError(f"Filtered FASTQ file for sample '{sample}' not found. Please filter step before analyzing.")

					# Align treatment sample with minimap2
					sample_alignment_path = os.path.join(aligned_dir, sample + '.dedup.sam')
					self._minimap2_align(
						self.parameters['reference_genome'],
						self.samples[sample]['filtered'],
						sample_alignment_path,
					)
					self.samples[sample]['aligned'] = sample_alignment_path

				except:
					logger.error(f'Failed for sample {sample}')
					logger.error(traceback.format_exc())
				logger.info('Finished aligning reads to genome.')

		except Exception as e:
			logger.error('Error aligning')
			logger.error(traceback.format_exc())
			quit()

	def _minimap2_align(self, reference, read, outfile):
		"""Align a single FASTQ file against the reference with minimap2 and
		remove PCR duplicates by position, producing a SAM output file.

		Parameters
		----------
		reference : str
			Path to the reference genome FASTA (or pre-built minimap2 index).
		read1 : str
			Path to the single-end FASTQ file.
		outfile : str
			Destination SAM path (must end in .dedup.sam).

		PCR duplicate definition
		------------------------
		Reads are considered duplicates if they share the same strand, the same
		starting genomic position (± 1 nt) and the same CIGAR-derived aligned
		length (± 1 nt).  Among duplicates, the read with the highest mapping
		quality (MAPQ) is kept; ties are broken by keeping the first encountered.
		"""
		minimap2 = self.parameters.get('minimap2', 'minimap2')
		samtools = self.parameters['samtools']
		njobs    = self.parameters['njobs']

		sort_bam = outfile.replace('.dedup.sam', '.sort.bam')
		raw_sam  = outfile.replace('.dedup.sam', '.raw.sam')

		# 1. Map with minimap2 and convert/sort to BAM
		map_cmd = (
			f"{minimap2} -a -x map-ont -t {njobs} {reference} {read} "
			f"| {samtools} view -bS -@ {njobs} - "
			f"| {samtools} sort -@ {njobs} -o {sort_bam}"
		)
		logger.info(f"minimap2 align: {map_cmd}")
		subprocess.check_call(map_cmd, shell=True)

		# 2. Index the sorted BAM
		subprocess.check_call(f"{samtools} index {sort_bam}", shell=True)

		# 3. Convert sorted BAM to SAM for duplicate removal
		subprocess.check_call(
			f"{samtools} view -h -o {raw_sam} {sort_bam}", shell=True
		)

		# 4. Position-based PCR duplicate removal
		logger.info(f"Removing PCR duplicates by position: {raw_sam}")
		self._remove_positional_duplicates(raw_sam, outfile)

		# 5. Sort dedup SAM to BAM and index
		dedup_bam = outfile.replace('.dedup.sam', '.dedup.bam')
		subprocess.check_call(
			f"{samtools} sort -@ {njobs} -o {dedup_bam} {outfile}", shell=True
		)
		subprocess.check_call(f"{samtools} index {dedup_bam}", shell=True)
		logger.info(f"Alignment complete: {outfile}")

	@staticmethod
	def _remove_positional_duplicates(in_sam, out_sam):
		"""Remove PCR duplicates from a SAM file based on genomic position.

		For each group of reads sharing the same strand, start position (± 4 nt)
		and CIGAR-derived aligned length (± 4 nt), retain only the read with the
		highest MAPQ score.

		Parameters
		----------
		in_sam : str
			Path to input SAM file (must be coordinate-sorted).
		out_sam : str
			Path to write the deduplicated SAM file.
		"""
		import re

		def cigar_ref_length(cigar):
			"""Sum reference-consuming CIGAR operations (M, D, N, =, X)."""
			if cigar == '*':
				return 0
			return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)
					   if op in 'MDN=X')

		def duplicate_key(pos, ref_len, strand):
			"""Round position and length to nearest 4 nt to allow ± 2 nt tolerance."""
			return (round(pos / 4) * 4, round(ref_len / 4) * 4, strand)

		headers = []
		reads   = []  # list of (flag_int, pos_int, mapq_int, ref_len, raw_line)

		with open(in_sam, 'r') as f:
			for line in f:
				if line.startswith('@'):
					headers.append(line)
					continue
				fields  = line.split('\t')
				flag    = int(fields[1])
				pos     = int(fields[3])
				mapq    = int(fields[4])
				cigar   = fields[5]
				ref_len = cigar_ref_length(cigar)
				strand  = '-' if flag & 16 else '+'
				reads.append((flag, pos, mapq, ref_len, strand, line))

		# Build duplicate groups: key → best (mapq, line)
		best = {}
		for flag, pos, mapq, ref_len, strand, line in reads:
			# Discard unmapped, secondary, and supplementary alignments
			if flag & 4 or flag & 256 or flag & 2048:
				continue
			key = duplicate_key(pos, ref_len, strand)
			if key not in best or mapq > best[key][0]:
				best[key] = (mapq, line)

		with open(out_sam, 'w') as f:
			for h in headers:
				f.write(h)
			for mapq, line in best.values():
				f.write(line)

		n_in  = len(reads)
		n_out = len(best)
		logger.info(f"PCR dedup: {n_in} reads in → {n_out} reads out "
					f"({n_in - n_out} duplicates/unmapped removed)")

	def identifyOfftargetSites(self):
		logger.info(f'{GREEN}Identifying offtarget sites...{RESET}')

		try:
			self.identified = {}

			# Identify offtarget sites for each sample
			for sample in self.samples:
				try:
					# default input
					if not "aligned" in self.samples[sample]:
						self.samples[sample]['aligned'] = os.path.join(self.parameters['analysis_folder'], 'aligned', sample + '.dedup.sam')				
					# Prepare sample annotations
					sample_data = self.samples[sample]
					annotations = {}
					annotations['Description'] = sample_data['description']
					annotations['Targetsite'] = sample

					annotations['Sequence'] = sample_data['target']

					self.identified[sample] = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.txt')

					identifyOfftargetSites.analyze(self.samples[sample]['aligned'], self.parameters['reference_genome'], self.identified[sample], annotations,
												   self.parameters['window_size'], self.parameters['max_score'], self.parameters['control_primer'],self.parameters)

				except:
					logger.error(f'{RED}Failed for sample {sample}{RESET}')   
					logger.error(traceback.format_exc())

			logger.info('Finished identifying offtarget sites.')

		except Exception as e:
			logger.error(f'{RED}Error identifying offtarget sites.{RESET}')
			logger.error(traceback.format_exc())
			quit()

	def filterBackgroundSites(self):
		logger.info(f'{GREEN}Filtering background sites{RESET}')

		# self.filtered = {}
		for sample in self.samples:
			try:
				out = os.path.join(self.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
				filterBlackList(self.parameters['bedtools'], self.identified[sample], self.parameters['blacklist'], out)
				self.identified[sample] = out
				
				logger.info('Finished background filtering for {0} sample'.format(sample))

			except Exception as e:
				logger.error('Error filtering background sites: %s'%(sample))
				logger.error(traceback.format_exc())

	def visualize(self):
		logger.info(f'{GREEN}Visualizing off-target sites{RESET}')

		for sample in self.samples:
			try:
				infile = self.identified[sample]
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '_offtargets')
				try:
					self.parameters['PAM']
					visualizeOfftargets(infile, outfile, title=sample,PAM=self.parameters['PAM'],genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				except:
					visualizeOfftargets(infile, outfile, title=sample,PAM="NGG",genome=self.parameters['genome'],refseq_names=self.parameters['refseq_names'])
				
				# Manhattan plot
				outfile = os.path.join(self.parameters['analysis_folder'], 'visualization', sample + '.Manhattan.pdf')
				# command = f"{self.parameters['Rscript']} {self.parameters['Manhattan.R']} {self.parameters['on_target_reference_sequence']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				command = f"{self.parameters['Rscript']} --vanilla {self.parameters['Manhattan.R']} {self.samples[sample]['target']} {infile.replace('.txt','.rm_no_match.annot.tsv')} {outfile}"
				logger.info(command)
				subprocess.call(command,shell=True)
				
			except Exception as e:
				logger.error('Error visualizing off-target sites: %s'%(sample))
				logger.error(traceback.format_exc())
		logger.info('Finished visualizing off-target sites')

	# Generates a tab-seperated table of 'visualize' output of all samples
	# in the manifest (in manifest order) in order to better compare them.
	#	- Loads each _offtargets.tsv file as a pandas dataframe.
	#	- Performs an outer-join on these dataframes
	#	- After all joins, already sorted alphabetically by Coordinates, but then sorts by mismatches (low to high),
	#	  and then re-sorts by total reads across samples (high to low), which should put the best hits at the top of the list
	#	- Outputs this joined table as a new TSV with the same name as the manifest
	def tableize(self, manifest_path):
		logger.info(f'{GREEN}Combining all sample off-targets into one table...{RESET}')

		try:
			visualize_dir = os.path.join(self.parameters['analysis_folder'], 'visualization')
			table_path = manifest_path.replace('.yaml', '.tsv')
			print(table_path)

			# Generate output pandas dataframe
			output_df = pd.DataFrame()
			common_cols = ['Mismatches', 'Annotation', 'SubstitutionsOnly.Sequence', 'GapsAllowed.Sequence']

			for sample in self.samples:
				try:
					# Establish location of visualization output tsv file. (./visualization/[sample]_offtargets.tsv)
					self.samples[sample]['visualize'] = os.path.join(visualize_dir, sample + '_offtargets.tsv')

					# Load TSV data into temporary dataframe
					# Int64 used instead of int to avoid int column being converted to float when 'NaN' appears in a column after a merge
					sample_df = pd.read_csv(self.samples[sample]['visualize'], sep='\t', header=0, dtype={'Mismatches': 'Int64'})
					sample_df[f'{sample}'] = sample_df.pop('Reads')
	
					# Merge sample dataframe with output dataframe
					#	- If this is the first round of merging, just transfer sample into output dataframe
					#	- Spot-checks for mismatches in annotations / guide sequence, just in case something went wrong.
					#	- Merges 'common' columns between datasets
					#	- After merge, moves sample-specific reads columns to end for readability
					if output_df.empty:
						output_df = sample_df
					else:	
						output_df = pd.merge(output_df, sample_df, on='Coordinates', how='outer')

						for col in common_cols:
							mismatch = output_df[f'{col}_x'].notna() & output_df[f'{col}_y'].notna() & (output_df[f'{col}_x'] != output_df[f'{col}_y'])
							if mismatch.any():
								print(f"Warning: {mismatch.sum()} disagreeing rows in {col} between final output table and sample {sample}")

							output_df[col] = output_df[f'{col}_x'].combine_first(output_df[f'{col}_y'])
							output_df.drop(columns=[f'{col}_x', f'{col}_y'], inplace=True)

					unique_cols = [c for c in output_df.columns if c not in common_cols + ['Coordinates']]
					output_df = output_df[['Coordinates'] + common_cols + unique_cols]
					output_df[unique_cols] = output_df[unique_cols].fillna(0).astype('Int64')
				
				except:
					logger.error(f'Failed for sample {sample}')
					logger.error(traceback.format_exc())

			# Reorder final table for ease of analysis
			output_df['Sample_sum'] = output_df[unique_cols].sum(axis=1, skipna=True)
			output_df = output_df.sort_values(by=['Sample_sum','Mismatches'], ascending=[False, True], na_position='last')

			logger.info('Finished generating final table.')
			clean_output = output_df.to_csv(sep='\t', index=False)
			logger.info(clean_output)
			table_file = open(table_path, 'w')
			table_file.write(clean_output)

		except Exception as e:
			logger.error('Error generating final table.')
			logger.error(traceback.format_exc())
			quit()

def main():
	# Pull out arguments from command line. 
	args = parse_args()

	# Use 'init' command-line argument to get a copy of default.yaml for editing.
	if args.command == "init":
		subprocess.call(f"cp {DEFAULT_YAML} Guideseq.yaml",shell=True)
		print("Guideseq.yaml has been initialize. \nPlease add sample info and then run GUIDE-seq analysis using the 'main' subcommand.")
	
	# Use 'main' command-line argument to run analysis'
	# Individual steps can be run with the --step command plus a combination of:
	#		filter 		-	Remove RAPID primer sequences and put dsODN sequence on 5' end of read.
	#		align 		-	Align filtered reads to genome of choice and deduplicate reads based on alignment.
	#		identify 	-	Utilize HOMER to find sites that are similar to guide sequence. Predetermined background sites are also filtered at this step.
	#		visualize 	-	Generate graphical and text tables of results for further analysis.	
	#		tableize	-	Puts all samples in a tsv table for simple comparison.
	elif args.command == "main":
		logger.info("User command: %s"%(" ".join(sys.argv)))
		g = GuideSeq()
		g.parseManifest(args.manifest,args.sample,args.overwrite)
		steps  = args.step.split("+")
		if "filter" in steps:
			g.filterReads()
		if "align" in steps:
			g.alignReads()
		if "identify" in steps:
			g.identifyOfftargetSites()
			g.filterBackgroundSites()
		if "visualize" in steps:
			try:
				g.identified
			except:
				g.identified = {}
				for sample in g.samples:
					if 'identified' in g.samples[sample]:
						g.identified[sample] = g.samples[sample]['identified']
						g.identified["Control_"+sample] = g.samples[sample]['controlidentified']
					else:
						g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
						sample = "Control_"+sample
						g.identified[sample] = os.path.join(g.parameters['analysis_folder'], 'identified', sample + '_identifiedOfftargets.rmblck.txt')
			g.visualize()
		if "tableize" in steps:
			g.tableize(args.manifest)
 
if __name__ == '__main__':
	main()
