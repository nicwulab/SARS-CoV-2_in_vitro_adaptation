#!/usr/bin/python
import sys
import pysam
import logomaker
import matplotlib.pyplot as plt
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
    "---":"-"}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def extract_MBCS(bamfile):
	read_count = 0
	MBCS_list  = []
	start = 23605
	end = 23619
	length = end-start+1
	for aln in pysam.Samfile(bamfile).fetch('BavtPat1.fa',start,end):
		read_count += 1
		MBCS = ''
		for pos in aln.get_aligned_pairs(matches_only=False, with_seq=True):
			query_pos, ref_pos, ref_base = pos
			if ref_pos == None:
				continue
			if start <= ref_pos <= end:
				if query_pos == None:
					MBCS += '-'
				else:
					query_base = aln.query_sequence[query_pos]
					MBCS += query_base
		if len(MBCS) == length:
			MBCS_list.append(MBCS)
	return MBCS_list

def filter_MBCS_by_freq(MBCS_list, min_freq):
	MBCS_dict = Counter(MBCS_list)
	MBCS_high_freq = [MBCS for MBCS in MBCS_dict.keys() if MBCS_dict[MBCS]/len(MBCS_list) > min_freq]
	return (MBCS_high_freq)

def make_sequence_logo(sequence_list):
	height_per_row = .8
	width_per_col = 1.5
	num_cols = 4
	num_rows = 1
	seqlogo_matrix = logomaker.alignment_to_matrix(sequence_list)
	seqlogo = logomaker.Logo(seqlogo_matrix, font_name="Arial", color_scheme="weblogo_protein", width=1)
	seqlogo.style_spines(visible=False)
	seqlogo.ax.set_xticks([])
	seqlogo.ax.set_yticks([])
	plt.savefig("graph/MBCS_seqlogo.png")
	plt.close()

def write_freq(dict_freq, outfilename):
	outfile = open(outfilename,'w')
	outfile.write("\t".join(['MBCS', 'freq'])+"\n")
	for MBCS in sorted(dict_freq.keys(),key=lambda x:dict_freq[x]):
		freq = "{0:.3g}".format(dict_freq[MBCS]/sum(dict_freq.values()))
		outfile.write("\t".join(map(str,[MBCS, freq]))+"\n")
	outfile.close()

def main():
	bamfiles = ['sorted.bam'] 
	outfilename = "results/MBCS_freq.tsv"
	min_freq = 0.01
	for bamfile in bamfiles:
            MBCS_list = extract_MBCS(bamfile)
            MBCS_high_freq = filter_MBCS_by_freq(MBCS_list, min_freq)
            MBCS_list_high_freq = [translation(MBCS) for MBCS in MBCS_list if MBCS in MBCS_high_freq]
            MBCS_dict_high_freq = Counter(MBCS_list_high_freq)
            make_sequence_logo(MBCS_list_high_freq)
            write_freq(MBCS_dict_high_freq, outfilename)

if __name__ == "__main__":
	main()
