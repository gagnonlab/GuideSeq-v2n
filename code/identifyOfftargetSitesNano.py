# IdentifyOffTargetSiteSequences.py
#
# 2015-10-05 Replaced swalign with regex matching.
# 2017-05-31 Replaced nwalign with a explicit search of realignments which uses regex.search.
# 2017-06-03 Output the best offtarget sequences with and/or without bulges, if any.
# 2026-07-01 Adapted for single-end, UMI-free ONT reads (no paired-end logic, no barcode dedup).

from __future__ import print_function

import argparse
import collections
import numpy
import os
import string
import operator
import pyfaidx
import re
import regex
import logging
from Levenshtein import distance
import pandas as pd
import dill

logger = logging.getLogger('root')


class chromosomePosition():
	#Accumulates per-position read counts from single-end, UMI-free alignments.

	def __init__(self, reference_genome):
		self.chromosome_dict = {}
		self.chromosome_position_dict = {}   # replaces chromosome_barcode_dict
		self.position_summary = []
		self.index_stack = {}
		self.genome = pyfaidx.Fasta(reference_genome)

	def addPosition(self, chromosome, position, strand, primer, count):
		"""Record a single read at (chromosome, position, strand) with primer label."""
		if chromosome not in self.chromosome_position_dict:
			self.chromosome_position_dict[chromosome] = {}
		if position not in self.chromosome_position_dict[chromosome]:
			d = self.chromosome_position_dict[chromosome][position] = {}
			for s in ('+', '-'):
				d[s + 'primer1_total']           = 0
				d[s + 'primer2_total']           = 0
				d[s + 'primer1_mispriming_total'] = 0
				d[s + 'primer2_mispriming_total'] = 0
				d[s + 'nomatch_total']            = 0

		d = self.chromosome_position_dict[chromosome][position]
		d[strand + primer + '_total'] += count

	def getSequence(self, genome, chromosome, start, end, strand="+"):
		if strand == "+":
			seq = self.genome[chromosome][int(start):int(end)]
		elif strand == "-":
			seq = self.genome[chromosome][int(start):int(end)].reverse.complement
		return seq

	def SummarizePositions(self):
		"""Build a flat list of per-position count rows (no UMI dimension)."""
		self.position_summary = [
			[chromosome, position,
			 d['+primer1_total'] + d['+primer2_total'] + d['+primer1_mispriming_total'] + d['+primer2_mispriming_total'] + d['+nomatch_total'],  # +_total
			 d['-primer1_total'] + d['-primer2_total'] + d['-primer1_mispriming_total'] + d['-primer2_mispriming_total'] + d['-nomatch_total'],  # -_total
			 d['+primer1_total'], d['+primer2_total'],
			 d['-primer1_total'], d['-primer2_total'],
			 d['+primer1_mispriming_total'], d['+primer2_mispriming_total'],
			 d['-primer1_mispriming_total'], d['-primer2_mispriming_total'],
			]
			for chromosome in sorted(self.chromosome_position_dict)
			for position, d in sorted(self.chromosome_position_dict[chromosome].items())
		]
		self.chr_dataframe_dict = {
			c: pd.DataFrame(list(self.chromosome_position_dict[c].keys()))
			for c in self.chromosome_position_dict
		}
		return self.position_summary

	def SummarizeIndex(self, windowsize):
		"""Collapse positions within windowsize bp into windows; compute per-window totals."""
		last_chromosome, last_position, window_index = 0, 0, 0
		index_summary = []

		for row in self.position_summary:
			chromosome, position = row[0], row[1]
			if chromosome != last_chromosome or abs(position - last_position) > 25:
				window_index += 1
			last_chromosome, last_position = chromosome, position
			if window_index not in self.index_stack:
				self.index_stack[window_index] = []
			self.index_stack[window_index].append(row)

		for index in self.index_stack:
			rows = self.index_stack[index]
			# Sort by total read count (plus + minus)
			rows_sorted = sorted(rows, key=lambda r: r[2] + r[3])

			chromosomes    = [r[0] for r in rows_sorted]
			positions      = [r[1] for r in rows_sorted]
			plus_totals    = [r[2] for r in rows_sorted]
			minus_totals   = [r[3] for r in rows_sorted]

			total_plus   = sum(plus_totals)
			total_minus  = sum(minus_totals)
			total_sum    = total_plus + total_minus
			total_geom   = min(total_plus, total_minus)

			plus_primer1   = sum(r[4] for r in rows_sorted)
			plus_primer2   = sum(r[5] for r in rows_sorted)
			minus_primer1  = sum(r[6] for r in rows_sorted)
			minus_primer2  = sum(r[7] for r in rows_sorted)
			plus_p1_mis    = sum(r[8] for r in rows_sorted)
			plus_p2_mis    = sum(r[9] for r in rows_sorted)
			minus_p1_mis   = sum(r[10] for r in rows_sorted)
			minus_p2_mis   = sum(r[11] for r in rows_sorted)

			primer1             = plus_primer1 + minus_primer1
			primer2             = plus_primer2 + minus_primer2
			primer1_mispriming  = plus_p1_mis  + minus_p1_mis
			primer2_mispriming  = plus_p2_mis  + minus_p2_mis
			primer_geom         = (primer1 * primer2) ** 0.5
			position_std        = numpy.std(positions)
			min_position        = min(positions)
			max_position        = max(positions)

			most_frequent_chromosome = rows_sorted[-1][0]
			most_frequent_position   = rows_sorted[-1][1]
			BED_chromosome           = most_frequent_chromosome
			BED_name                 = f"{BED_chromosome}:{min_position}-{max_position}"

			offtarget_sequence = self.getSequence(
				self.genome, most_frequent_chromosome,
				most_frequent_position - windowsize,
				most_frequent_position + windowsize
			)

			# +.mi, -.mi, bi.sum.mi, bi.geometric_mean.mi have no UMI equivalent and remain 0.
			# primer1/2 and mispriming columns now carry real raw read counts.
			summary_list = [str(x) for x in [
				index, most_frequent_chromosome, most_frequent_position, offtarget_sequence,
				BED_chromosome, min_position, max_position, BED_name,
				0, 0, 0, 0,                  # +.mi, -.mi, bi.sum.mi, bi.geometric_mean.mi
				total_plus, total_minus, total_sum, total_geom,
				primer1, primer2, primer1_mispriming, primer2_mispriming,
				primer_geom, position_std
			]]

			# Keep windows with reads from at least one strand
			if total_sum > 0:
				index_summary.append(summary_list)

		return index_summary


# ---------------------------------------------------------------------------
# Sequence alignment helpers (unchanged from original)
# ---------------------------------------------------------------------------

def regexFromSequence(seq, lookahead=True, indels=1, errors=7):
	seq = seq.upper()
	from Bio.Data import IUPACData
	IUPAC_notation_regex = IUPACData.ambiguous_dna_values
	pattern = ''
	for c in seq:
		pattern += "[" + IUPAC_notation_regex[c] + "]"
	if lookahead:
		pattern = '(?b:' + pattern + ')'
	pattern_standard = pattern + '{{s<={0}}}'.format(errors)
	pattern_gap      = pattern + '{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors)
	return pattern_standard, pattern_gap


def extendedPattern(seq, indels=1, errors=7):
	IUPAC_notation_regex_extended = {
		'N': '[ATCGN]', '-': '[ATCGN]', 'B': '[TCGB]', 'D': '[ATGD]',
		'H': '[ATCH]',  'V': '[ACGV]',  'Y': '[CTY]',  'R': '[AGR]',
		'W': '[ATW]',   'S': '[CGS]',   'K': '[TGK]',  'M': '[ACM]',
		'A': 'A',		'T': 'T',		'C': 'C',      'G': 'G'
	}
	realign_pattern = ''
	for c in seq:
		realign_pattern += IUPAC_notation_regex_extended[c]
	return ('(?b:' + realign_pattern + ')' +
			'{{i<={0},d<={0},s<={1},3i+3d+1s<={1}}}'.format(indels, errors))


def realignedSequences(targetsite_sequence, chosen_alignment, errors=7, PAM="NGG"):
	match_sequence = chosen_alignment.group()
	substitutions, insertions, deletions = chosen_alignment.fuzzy_counts
	realigned_fuzzy = (substitutions, max(0, insertions - 1), max(0, deletions - 1))

	if insertions:
		if targetsite_sequence.index(PAM) > len(targetsite_sequence) / 2:
			targetsite_realignments = [
				targetsite_sequence[:i + 1] + '-' + targetsite_sequence[i + 1:]
				for i in range(targetsite_sequence.index(PAM) + 1)
			]
		else:
			targetsite_realignments = [
				targetsite_sequence[:i] + '-' + targetsite_sequence[i:]
				for i in range(targetsite_sequence.index(PAM), len(targetsite_sequence))
			]
	else:
		targetsite_realignments = [targetsite_sequence]

	realigned_target_sequence, realigned_offtarget_sequence = None, ''

	for seq in targetsite_realignments:
		if deletions:
			match_realignments = [match_sequence[:i + 1] + '-' + match_sequence[i + 1:]
								  for i in range(len(match_sequence) - 1)]
			match_pattern = [match_sequence[:i + 1] + seq[i + 1] + match_sequence[i + 1:]
							 for i in range(len(match_sequence) - 1)]
		else:
			match_realignments = match_pattern = [match_sequence]

		x = extendedPattern(seq, errors)
		for y_pattern, y_alignment in zip(match_pattern, match_realignments):
			m = regex.search(x, y_pattern, regex.BESTMATCH)
			if m and m.fuzzy_counts == realigned_fuzzy:
				realigned_target_sequence, realigned_offtarget_sequence = seq, y_alignment

	return realigned_target_sequence, realigned_offtarget_sequence


def alignSequences(targetsite_sequence, window_sequence, max_score=7, PAM="NGG"):

	print(targetsite_sequence)
	print(window_sequence)
	print(max_score)
	print(PAM)

	window_sequence = window_sequence.upper()
	query_regex_standard, query_regex_gap = regexFromSequence(targetsite_sequence, errors=max_score)

	alignments_mm    = []
	alignments_bulge = []
	alignments_mm.append(('+', regex.search(query_regex_standard, window_sequence, regex.BESTMATCH)))
	alignments_mm.append(('-', regex.search(query_regex_standard, reverseComplement(window_sequence), regex.BESTMATCH)))
	alignments_bulge.append(('+', regex.search(query_regex_gap, window_sequence, regex.BESTMATCH)))
	alignments_bulge.append(('-', regex.search(query_regex_gap, reverseComplement(window_sequence), regex.BESTMATCH)))

	lowest_distance_score, lowest_mismatch = 100, max_score + 1
	chosen_alignment_b, chosen_alignment_m = None, None
	chosen_alignment_strand_b, chosen_alignment_strand_m = '', ''

	for strand_m, match_m in alignments_mm:
		if match_m is not None:
			mismatches = match_m.fuzzy_counts[0]
			if mismatches < lowest_mismatch:
				chosen_alignment_m      = match_m
				chosen_alignment_strand_m = strand_m
				lowest_mismatch         = mismatches

	for strand_b, match_b in alignments_bulge:
		if match_b is not None:
			substitutions, insertions, deletions = match_b.fuzzy_counts
			if insertions or deletions:
				distance_score = substitutions + (insertions + deletions) * 3
				edistance      = substitutions + insertions + deletions
				if distance_score < lowest_distance_score and edistance < lowest_mismatch:
					chosen_alignment_b      = match_b
					chosen_alignment_strand_b = strand_b
					lowest_distance_score   = distance_score

	if chosen_alignment_m:
		offtarget_sequence_no_bulge  = chosen_alignment_m.group()
		mismatches                   = chosen_alignment_m.fuzzy_counts[0]
		start_no_bulge               = chosen_alignment_m.start()
		end_no_bulge                 = chosen_alignment_m.end()
	else:
		offtarget_sequence_no_bulge, mismatches, start_no_bulge, end_no_bulge, chosen_alignment_strand_m = '', '', '', '', ''

	bulged_offtarget_sequence = score = length = substitutions = insertions = deletions = ''
	bulged_start = bulged_end = realigned_target = ''
	realigned_target = 'none'

	if chosen_alignment_b:
		realigned_target, bulged_offtarget_sequence = realignedSequences(
			targetsite_sequence, chosen_alignment_b, max_score, PAM
		)
		if bulged_offtarget_sequence:
			length         = len(chosen_alignment_b.group())
			substitutions, insertions, deletions = chosen_alignment_b.fuzzy_counts
			score          = substitutions + (insertions + deletions) * 3
			bulged_start   = chosen_alignment_b.start()
			bulged_end     = chosen_alignment_b.end()
		else:
			chosen_alignment_strand_b = ''

	return [offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m,
			start_no_bulge, end_no_bulge,
			bulged_offtarget_sequence, length, score,
			substitutions, insertions, deletions,
			chosen_alignment_strand_b, bulged_start, bulged_end, realigned_target]


# ---------------------------------------------------------------------------
# Primer / control-primer matching (unchanged logic, kept for dsODN detection)
# ---------------------------------------------------------------------------

def assignPrimerstoReads(read_sequence, sam_flag, dsODN_dict=None):
	if int(sam_flag) & 16:
		read_sequence = reverseComplement(read_sequence)
	flag, seq1, myDistance1, distance2 = match_dsODN(read_sequence, dsODN_dict['i7-'], dsODN_dict['i7-_match_distance'])
	if flag:
		extend_sequence  = read_sequence[len(dsODN_dict['i7-']):len(dsODN_dict['dsODN_primer'])][:dsODN_dict['i7-_mispriming_length']]
		correct_sequence = dsODN_dict['dsODN_primer'][len(dsODN_dict['i7-']):len(dsODN_dict['dsODN_primer'])][:dsODN_dict['i7-_mispriming_length']]
		if distance(extend_sequence, correct_sequence) >= dsODN_dict['i7-_mispriming_distance']:
			return "primer1", False, seq1, myDistance1, distance2
		else:
			return "primer1", True, seq1, myDistance1, distance2
	flag, seq2, myDistance2, distance2 = match_dsODN(read_sequence, dsODN_dict['i7+'], dsODN_dict['i7+_match_distance'])
	if flag:
		extend_sequence  = read_sequence[len(dsODN_dict['i7+']):len(dsODN_dict['dsODN_primer_revcomp'])][:dsODN_dict['i7+_mispriming_length']]
		correct_sequence = dsODN_dict['dsODN_primer_revcomp'][len(dsODN_dict['i7+']):len(dsODN_dict['dsODN_primer_revcomp'])][:dsODN_dict['i7+_mispriming_length']]
		if distance(extend_sequence, correct_sequence) >= dsODN_dict['i7+_mispriming_distance']:
			return "primer2", False, seq2, myDistance2, distance2
		else:
			return "primer2", True, seq2, myDistance2, distance2
	if myDistance1 < myDistance2:
		return "nomatch", False, seq1, myDistance1, distance2
	return "nomatch", False, seq2, myDistance2, distance2


def match_dsODN(read_sequence, primer, cutoff):
	read_start_sequence = read_sequence[:len(primer)]
	myDistance          = distance(read_start_sequence, primer)
	distance2           = myDistance - read_start_sequence.count("N")
	if myDistance <= cutoff:
		return True,  read_start_sequence, myDistance, distance2
	else:
		return False, read_start_sequence, myDistance, distance2


def contain_control_primer(read_sequence, sam_flag, myDict=None):
	control_primer = myDict['control_primer']
	if control_primer == "":
		return "nomatch"
	if int(sam_flag) & 16:
		read_sequence = reverseComplement(read_sequence)
	myDistance = distance(read_sequence[:len(control_primer)], control_primer)
	if myDistance <= myDict['control_primer_match_distance']:
		return "control_rev" if int(sam_flag) & 16 else "control_fwd"
	return "nomatch"


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def is_control(chr, pos, ref_chr, ref_start, ref_end):
	if chr != ref_chr:
		return False
	return ref_start <= pos <= ref_end


def get_num_pos_given_pos(df, pos):
	out = []
	for i in [500, 1000, 2000]:
		half = int(i / 2)
		out.append(df[0].between(pos - half, pos + half).sum())
	return out


def py2min(myList):
	return min(i for i in myList if i != "")


def py2max(myList):
	return max(i for i in myList if i != "")


def reverseComplement(sequence):
	import sys
	if sys.version_info[0] < 3:
		tab = string.maketrans("ACGTacgt", "TGCATGCA")
	else:
		tab = str.maketrans("ACGTacgt", "TGCATGCA")
	return sequence.translate(tab)[::-1]


def save_object(obj, out):
	try:
		with open(out, "wb") as f:
			dill.dump(obj, f)
	except Exception as ex:
		print("Error during pickling object (Possibly unsupported):", ex)


def parseReadName(read_name):
	"""Single-end reads have no UMI; treat the full read name as the identifier."""
	return read_name, 1


def cigar_read_length(cigar):
	"""Return the number of reference-consuming bases described by a CIGAR string.
	Counts M, D, N, =, X operations (reference-consuming); excludes S, I, H, P.
	"""
	return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar) if op in 'MDN=X')


# ---------------------------------------------------------------------------
# Main analysis function
# ---------------------------------------------------------------------------

def analyze(sam_filename, reference_genome, outfile, annotations, windowsize,
			max_score, control_primer, myDict=None):
	"""Identify CRISPR off-target sites from a single-end, UMI-free SAM file.

	Key differences from the paired-end UMI version:
	- Reads are not filtered by pair-mate flag (flag & 128).  Every mapped,
	  non-supplementary read above the MAPQ threshold is counted.
	- There is no barcode/UMI deduplication; raw read counts are used directly.
	- Integration position is taken from the read's own mapping coordinate
	  (forward reads: leftmost position; reverse reads: rightmost position),
	  since there is no mate to anchor the cut site.
	- The geometric-mean filter on UMI counts is replaced by a simple read-
	  count threshold (total_sum > 0); downstream filtering by
	  filterBackgroundSites should apply any desired count cutoffs.
	"""

	output_folder = os.path.dirname(outfile)
	if output_folder and not os.path.exists(output_folder):
		os.makedirs(output_folder)

	primer_log  = open(outfile + ".primer.tsv",    'w')
	tl_filter   = open(outfile + ".tl_filter.tsv", 'w')

	logger.info("Processing SAM file %s", sam_filename)

	chromosome_position = chromosomePosition(reference_genome)
	control_primer_count     = 0
	control_counts           = 0
	total_dsODN              = 0
	control_primer_count_dict = {}
	total_dsODN_count_dict    = {}

	with open(sam_filename, 'r') as f:
		for line in f:
			if line.startswith('@'):
				continue                       # skip SAM header lines
			fields = line.split('\t')
			if len(fields) < 11:
				continue

			(full_read_name, sam_flag, chromosome, position, mapq,
			 cigar, name_of_mate, position_of_mate,
			 template_length, read_sequence, read_quality) = fields[:11]

			sam_flag_int      = int(sam_flag)
			mapq_int          = int(mapq)
			position_int      = int(position)
			template_length_int = int(template_length)

			# Drop supplementary / secondary alignments
			if sam_flag_int & 2048 or sam_flag_int & 256:
				continue

			# Drop reads below MAPQ threshold
			if mapq_int < myDict['mapq_threshold']:
				continue

			# Drop reads with implausibly large template length
			if abs(template_length_int) > 50000:
				print(full_read_name, read_sequence, sam_flag,
					  chromosome, template_length, sep="\t", file=tl_filter)
				continue

			# ----- control primer check -----
			control_read = contain_control_primer(read_sequence, sam_flag, myDict)
			control_primer_count_dict[control_read] = \
				control_primer_count_dict.get(control_read, 0) + 1
			if control_read != "nomatch":
				control_primer_count += 1

			# ----- dsODN primer assignment -----
			primer, flag, seq, myDistance, distance2 = assignPrimerstoReads(
				read_sequence, sam_flag, dsODN_dict=myDict
			)
			if primer != "nomatch":
				if flag:
					total_dsODN += 1
				else:
					primer = primer + "_mispriming"
			total_dsODN_count_dict[primer] = \
				total_dsODN_count_dict.get(primer, 0) + 1

			# ----- integration position (single-end) -----
			# Input FASTQs always have the dsODN on the 5' end of the read.
			# For reads aligning to the forward strand the dsODN is at the
			# leftmost mapped base (position_int).
			# For reads aligning to the reverse strand the read is reverse-
			# complemented on the genome, so the 5' dsODN end maps to the
			# rightmost base, estimated as position_int + read_length - 1.
			if sam_flag_int & 16:   # reverse strand — dsODN at rightmost mapped base
				read_position = position_int + cigar_read_length(cigar) - 1
				strand        = "-"
			else:                   # forward strand — dsODN at leftmost base
				read_position = position_int
				strand        = "+"

			chromosome_position.addPosition(chromosome, read_position, strand, primer, 1)

			print(full_read_name, read_sequence, sam_flag, chromosome,
				  read_position, seq, myDistance, distance2, primer,
				  sep="\t", file=primer_log)

	primer_log.close()
	tl_filter.close()

	# ----- Summarise positions into windows -----
	chromosome_position.SummarizePositions()

	if control_primer_count == 0:
		control_primer_count = -1

	with open(outfile, 'w') as f:
		print('#BED_Chromosome', 'BED_Min.Position', 'BED_Max.Position', 'BED_Name', 'Filename',
			  'WindowIndex', 'WindowChromosome', 'Position', 'WindowSequence',
			  '+.mi', '-.mi', 'bi.sum.mi', 'bi.geometric_mean.mi',
			  '+.total', '-.total', 'total.sum', 'total.geometric_mean',
			  'primer1', 'primer2', 'primer1_mispriming', 'primer2_mispriming',
			  'primer.geometric_mean', 'position.stdev',
			  'BED_Site_Name', 'BED_Score', 'BED_Site_Chromosome',
			  'Site_SubstitutionsOnly.Sequence', 'Site_SubstitutionsOnly.NumSubstitutions',
			  'Site_SubstitutionsOnly.Strand', 'Site_SubstitutionsOnly.Start', 'Site_SubstitutionsOnly.End',
			  'Site_GapsAllowed.Sequence', 'Site_GapsAllowed.Length', 'Site_GapsAllowed.Score',
			  'Site_GapsAllowed.Substitutions', 'Site_GapsAllowed.Insertions', 'Site_GapsAllowed.Deletions',
			  'Site_GapsAllowed.Strand', 'Site_GapsAllowed.Start', 'Site_GapsAllowed.End',
			  'Cell', 'Targetsite', 'TargetSequence', 'RealignedTargetSequence',
			  'control_primer_reads', 'control_counts', 'normalization_ratio',
			  '#pos_500bp', '#pos_1kb', '#pos_2kb',
			  sep='\t', file=f)

		summary        = chromosome_position.SummarizeIndex(windowsize)
		target_sequence = annotations["Sequence"]
		annotation      = [annotations['Description'], annotations['Targetsite'], annotations['Sequence']]
		output_dict     = {}

		for row in summary:
			window_sequence, window_chromosome, window_start, window_end, BED_name = row[3:8]

			non_bulged_target_start_absolute = bulged_target_start_absolute = ''

			if target_sequence:
				(offtarget_sequence_no_bulge, mismatches, chosen_alignment_strand_m,
				 start_no_bulge, end_no_bulge,
				 bulged_offtarget_sequence, length, dist, substitutions,
				 insertions, deletions, chosen_alignment_strand_b,
				 bulged_start, bulged_end, realigned_target_sequence) = \
					alignSequences(target_sequence, window_sequence, max_score, myDict['PAM'])

				BED_score      = 1
				BED_chromosome = window_chromosome

				if chosen_alignment_strand_m == "+":
					non_bulged_target_start_absolute = start_no_bulge + int(row[2]) - windowsize
					non_bulged_target_end_absolute   = end_no_bulge   + int(row[2]) - windowsize
				elif chosen_alignment_strand_m == "-":
					non_bulged_target_start_absolute = int(row[2]) + windowsize - end_no_bulge
					non_bulged_target_end_absolute   = int(row[2]) + windowsize - start_no_bulge
				else:
					non_bulged_target_start_absolute, non_bulged_target_end_absolute = '', ''

				if chosen_alignment_strand_b == "+":
					bulged_target_start_absolute = bulged_start + int(row[2]) - windowsize
					bulged_target_end_absolute   = bulged_end   + int(row[2]) - windowsize
				elif chosen_alignment_strand_b == "-":
					bulged_target_start_absolute = int(row[2]) + windowsize - bulged_end
					bulged_target_end_absolute   = int(row[2]) + windowsize - bulged_start
				else:
					bulged_target_start_absolute, bulged_target_end_absolute = '', ''

				if not (chosen_alignment_strand_m or chosen_alignment_strand_b):
					BED_chromosome, BED_score, BED_name = '', '', ''

				output_row = (row[4:8] + [os.path.basename(sam_filename)] +
							  row[0:4] + row[8:] +
							  [str(x) for x in [
								  BED_name, BED_score, BED_chromosome,
								  offtarget_sequence_no_bulge, mismatches,
								  chosen_alignment_strand_m,
								  non_bulged_target_start_absolute, non_bulged_target_end_absolute,
								  bulged_offtarget_sequence, length, dist,
								  substitutions, insertions, deletions,
								  chosen_alignment_strand_b,
								  bulged_target_start_absolute, bulged_target_end_absolute]] +
							  [str(x) for x in annotation] + [realigned_target_sequence])
			else:
				output_row = [str(x) for x in
							  row[4:8] + [os.path.basename(sam_filename)] +
							  row[0:4] + row[8:] + [""] * 17 + annotation + ['none']]

			if non_bulged_target_start_absolute != '' or bulged_target_start_absolute != '':
				output_row_key = '{0}_{1}_{2}'.format(
					window_chromosome,
					py2min([non_bulged_target_start_absolute, bulged_target_start_absolute]),
					py2max([non_bulged_target_end_absolute,   bulged_target_end_absolute])
				)
			else:
				output_row_key = '{0}_{1}_{2}'.format(window_chromosome, window_start, window_end)

			if output_row_key in output_dict:
				# output_row[15] = total.sum (column index accounts for the 4 inserted UMI zero-columns)
				output_dict[output_row_key][15] = str(int(output_row[15]) + int(output_dict[output_row_key][15]))
			else:
				output_dict[output_row_key] = output_row

		# Determine control-region read count for normalisation
		for key in sorted(output_dict.keys()):
			current_pos = int(output_dict[key][7])
			chr_        = output_dict[key][0]
			if is_control(chr_, current_pos,
						  myDict['control_coord_chr'],
						  myDict['control_coord_start'],
						  myDict['control_coord_end']):
				control_counts = int(float(output_dict[key][15]))
				break
		if control_counts == 0:
			control_counts = -1

		for key in sorted(output_dict.keys()):
			current_pos = int(output_dict[key][7])
			chr_        = output_dict[key][0]
			norm_ratio  = float(output_dict[key][15]) / float(control_counts)
			print(*output_dict[key] +
				  [str(control_primer_count),
				   str(control_counts),
				   str(norm_ratio)] +
				  get_num_pos_given_pos(chromosome_position.chr_dataframe_dict[chr_], current_pos),
				  sep='\t', file=f)

	if myDict.get('save_pickle'):
		save_object(chromosome_position, outfile + ".pkl")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
	parser = argparse.ArgumentParser(
		description='Identify off-target candidates from single-end, UMI-free sequencing data.')
	parser.add_argument('--ref',           help='Reference Genome Fasta', required=True)
	parser.add_argument('--samfile',       help='SAM file', nargs='*')
	parser.add_argument('--outfile',       help='File to output identified sites to.', required=True)
	parser.add_argument('--window',        help='Window around breakpoint to search for off-target',
						type=int, default=25)
	parser.add_argument('--control_primer', help='Control primer forward sequence',
						type=str, default=None)
	parser.add_argument('--max_score',     help='Score threshold', type=int, default=7)
	parser.add_argument('--target',        default='')

	args = parser.parse_args()
	annotations = {'Description': 'test description', 'Targetsite': 'dummy targetsite',
				   'Sequence': args.target}
	analyze(args.samfile[0], args.ref, args.outfile, annotations,
			args.window, args.max_score, args.control_primer)


if __name__ == "__main__":
	main()
