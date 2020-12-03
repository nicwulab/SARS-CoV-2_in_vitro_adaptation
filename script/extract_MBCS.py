#!/usr/bin/python
import pysam

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

read_count = 0
for aln in pysam.Samfile('bam/sorted.bam').fetch('BavtPat1.fa',23605,23616):
	read_count += 1
	MBCS = ''
	for pos in aln.get_aligned_pairs(matches_only=False, with_seq=True):
		query_pos, ref_pos, ref_base = pos
		if ref_pos == None:
			continue
		if 23605 <= ref_pos <= 23616:
			if query_pos == None:
				MBCS += '-'
			else:
				query_base = aln.query_sequence[query_pos]
				MBCS += query_base
	if len(MBCS) == 12:
		if MBCS.count('-')%3 == 0:
			print ("\t".join(map(str,[read_count, MBCS, translation(MBCS.upper())])))
