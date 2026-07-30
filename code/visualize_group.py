from __future__ import print_function
import svgwrite
import sys
import os
import logging
import argparse
import pandas as pd
import subprocess
from Bio import SeqUtils

logger = logging.getLogger('root')
logger.propagate = False

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3', 'R': '#B3B3B3', '-': '#B3B3B3'}
for c in ['Y','S','W','K','M','B','D','H','V','.']:
    colors[c] = "#B3B3B3"

def parseSitesFile(infile):
    offtargets = []
    total_seq = 0
    with open(infile, 'r') as f:
        header_line = f.readline().rstrip('\n')
        header_items = header_line.split('\t')
        # Everything from column index 8 onward is a sample-derived count column;
        # the header gives us the name of each one, and there can be any number of them.
        sample_names = header_items[8:]

        for line in f:
            line = line.rstrip('\n')
            line_items = line.split('\t')
            
            # File construction:
            # 0     Coordinates
            # 1     Mismatches
            # 2     Name
            # 3     Annotation
            # 4     SubstitutionsOnly.Sequence
            # 5     GapsAllowed_Sequence
            # 6     Target Sequence
            # 7     Realigned Target Sequence
            # 8+    Sample-derived counts (variable number, named in the header row)

            no_bulge_offtarget_sequence = line_items[4]
            bulge_offtarget_sequence = line_items[5]
            target_seq = line_items[6]
            realigned_target_seq = line_items[7]
            # print (offtarget_reads,control_primer)

            # Grab the sample-derived values for this row, padding with '' if the
            # row happens to be shorter than the header (ragged TSV).
            sample_values = line_items[8:]
            if len(sample_values) < len(sample_names):
                sample_values = sample_values + [''] * (len(sample_names) - len(sample_values))

            if no_bulge_offtarget_sequence != '' or bulge_offtarget_sequence != '':
                if bulge_offtarget_sequence:
                    total_seq += 1
                if no_bulge_offtarget_sequence:
                    total_seq += 1
                offtargets.append({'seq': no_bulge_offtarget_sequence.strip(),
                                   'bulged_seq': bulge_offtarget_sequence.strip(),
                                   'sample_values': [v.strip() for v in sample_values],
                                   'coord': line_items[0],
                                   'annot': line_items[3],
                                   'num_mismatch': str(line_items[1]),
                                   'off_name': line_items[2],
                                   'target_seq': target_seq.strip(),
                                   'realigned_target_seq': realigned_target_seq.strip()
                                   })
    return offtargets, target_seq, total_seq, sample_names

def find_PAM(seq,PAM):
    try:
        PAM_index = seq.index(PAM)                                                          
    except:
        # PAM on the left
        left_search = SeqUtils.nt_search(seq[:len(PAM)], PAM)
        if len(left_search)>1:
            PAM_index = left_search[1]
        else:
            right_search = SeqUtils.nt_search(seq[-len(PAM):], PAM)
            if len(right_search)>1:
                PAM_index = len(seq)-len(PAM)
            else:
                print ("PAM: %s not found in %s. Set PAM index to 20"%(PAM,seq))
                PAM_index=20
    return PAM_index

def visualizeOfftargets(infile, outfile, title, PAM, genome=None,refseq_names=None):
    output_folder = os.path.dirname(outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, target_seq, total_seq, sample_names = parseSitesFile(infile)

    # How far apart each dynamic sample column is placed.
    sample_col_spacing = 80
                                                                                       
    # Initiate canvas (width is 100% so it naturally accommodates however many
    # sample columns there are).
    dwg = svgwrite.Drawing(outfile, profile='full', size=(u'100%', 100 + total_seq*(box_size + 1)))

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # x-position where the first sample column starts (mirrors the old fixed
    # 'Reads' column position).
    samples_x_start = x_offset + box_size * len(target_seq) + 900

    ## Find PAM
    tick_locations = []
    tick_legend = []
    PAM_index = find_PAM(target_seq,PAM)
    count = 0
    for i in range(PAM_index,0,-1):
        count = count+1
        if count % 10 == 0:
            tick_legend.append(count)
            # print (count,i)
            tick_locations.append(i)
    if len(PAM)>=3:
        tick_legend+=['P', 'A', 'M']+['-']*(len(PAM)-3)
    else:
        tick_legend+=["PAM"]+['-']*(len(PAM)-3)
    tick_locations+=range(PAM_index+1,len(target_seq)+1)
    if PAM_index == 0:
        tick_legend = []
        tick_locations = []
        tick_legend+=['P', 'A', 'M']+['-']*(len(PAM)-3)
        tick_locations+=range(1,len(PAM)+1)
        count = 0
        for i in range(len(PAM)+1,len(target_seq)+1):
            count = count+1
            if count % 10 == 0 or count == 1:
                tick_legend.append(count)
                # print (count,i)
                tick_locations.append(i)
    # print (zip(tick_locations, tick_legend))
    for x,y in zip(tick_locations, tick_legend):
        dwg.add(dwg.text(y, insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(target_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))
    
    dwg.add(dwg.text('Mismatches', insert=(box_size * (len(target_seq) + 1) + 20, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Name', insert=(box_size * (len(target_seq) + 1) + 140, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Coordinates', insert=(box_size * (len(target_seq) + 1) + 300, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    dwg.add(dwg.text('Annotation', insert=(box_size * (len(target_seq) + 1) + 550, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    for s_i, sample_name in enumerate(sample_names):
        dwg.add(dwg.text(sample_name, insert=(samples_x_start + s_i * sample_col_spacing, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    
    # Draw aligned sequence rows
    y_offset += 1  # leave some extra space after the reference row
    line_number = 0  # keep track of plotted sequences
                                        
    for j, seq in enumerate(offtargets):
        realigned_target_seq = offtargets[j]['realigned_target_seq']
        no_bulge_offtarget_sequence = offtargets[j]['seq']
        bulge_offtarget_sequence = offtargets[j]['bulged_seq']

        if no_bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(no_bulge_offtarget_sequence, target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
        if bulge_offtarget_sequence != '':
            k = 0
            line_number += 1
            y = y_offset + line_number * box_size
            for i, (c, r) in enumerate(zip(bulge_offtarget_sequence, realigned_target_seq)):
                x = x_offset + k * box_size
                if r == '-':
                    if 0 < k < len(realigned_target_seq):
                        x = x_offset + (k - 0.25) * box_size
                        dwg.add(dwg.rect((x, box_size * 1.4 + y), (box_size*0.6, box_size*0.6), fill=colors[c]))
                        dwg.add(dwg.text(c, insert=(x+1, 2 * box_size + y - 2), fill='black', style="font-size:10px; font-family:Courier"))
                elif c == r:
                    dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
                    k += 1
                elif r == 'N':
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1
                else:
                    dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                    dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
                    k += 1

        if no_bulge_offtarget_sequence == '' or bulge_offtarget_sequence == '':
            # Mismatch column
            mismatch_text = dwg.text(seq['num_mismatch'], insert=(box_size * (len(target_seq) + 1) + 22, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(mismatch_text)

            # Offtarget Name
            name_text = dwg.text(seq['off_name'], insert=(box_size * (len(target_seq) + 1) + 140, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(name_text)

            # Coordinates
            coord_text = dwg.text(seq['coord'], insert=(box_size * (len(target_seq) + 1) + 300, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(coord_text)
            
            # Annotations
            annot_text = dwg.text(seq['annot'], insert=(box_size * (len(target_seq) + 1) + 550, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(annot_text)
            
            # Sample-derived count columns
            for s_i, sample_value in enumerate(seq['sample_values']):
                sample_text = dwg.text(str(sample_value), insert=(box_size * (len(target_seq) + 1) + 900 + s_i * sample_col_spacing, y_offset + box_size * (line_number + 2) - 2),
                                  fill='black', style="font-size:15px; font-family:Courier")
                dwg.add(sample_text)


        else:
            # Bulge }
            reads_text02 = dwg.text(u"\u007D", insert=(box_size * (len(target_seq) + 1) + 2, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:23px; font-family:Courier")
            dwg.add(reads_text02)

            # Mismatch column
            mismatch_text = dwg.text(seq['num_mismatch'], insert=(box_size * (len(target_seq) + 1) + 22, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(mismatch_text)
            
            # Offtarget Name
            name_text = dwg.text(seq['off_name'], insert=(box_size * (len(target_seq) + 1) + 140, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(name_text)
            
            # Coordinates
            coord_text = dwg.text(seq['coord'], insert=(box_size * (len(target_seq) + 1) + 300, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(coord_text)
            
            # Annotations
            annot_text = dwg.text(seq['annot'], insert=(box_size * (len(target_seq) + 1) + 550, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
            dwg.add(annot_text)
            
            # Sample-derived count columns
            for s_i, sample_value in enumerate(seq['sample_values']):
                sample_text = dwg.text(str(sample_value), insert=(box_size * (len(target_seq) + 1) + 900 + s_i * sample_col_spacing, y_offset + box_size * (line_number + 1) + 5),
                                  fill='black', style="font-size:15px; font-family:Courier")
                dwg.add(sample_text)
            
    dwg.save()

def main():
    parser = argparse.ArgumentParser(description='Plot visualization plots for re-aligned reads.')
    parser.add_argument("-f","--infile", help="TSV file derived from previous steps", required=True)
    parser.add_argument("-o","--outfile", help="Name of output svgwrite", required=True)
    parser.add_argument("-t","--title", help="Plot title", default="Untitled")
    parser.add_argument("--PAM", help="PAM sequence", default="NGG")    
    args = parser.parse_args()

    print(args)

    visualizeOfftargets(args.infile, args.outfile, args.title, args.PAM)

if __name__ == "__main__":

    main()