#!/usr/bin/python3
import glob
from collections import defaultdict

def MBCS_freq_to_dict(freq_files):
  MBCS_dict = {}
  for freq_file in freq_files:
    sampleID = freq_file.rsplit('/')[1].rsplit('_')[0]
    infile = open(freq_file,'r') 
    MBCS_dict[sampleID] = defaultdict(float)
    for line in infile.readlines():
      if "freq" in line: continue
      line = line.rstrip().rsplit("\t")
      MBCS = line[0]
      freq = float(line[1])
      MBCS_dict[sampleID][MBCS] = freq
    infile.close()
  return MBCS_dict

def write_MBCS_freq(MBCS_list, MBCS_dict, out_filename):
  print ("writing: %s" % out_filename)
  outfile = open(out_filename,'w')
  outfile.write("\t".join(["sampleID", "MBCS", "freq"])+"\n")
  for sampleID in sorted(MBCS_dict.keys(), key=lambda x:int(x[1::])):
    for MBCS in MBCS_list:
      freq = MBCS_dict[sampleID][MBCS]
      if MBCS == "-----": MBCS = "deletion"
      outfile.write("\t".join(map(str,[sampleID, MBCS, freq]))+"\n")
  outfile.close()
   
def main():
  freq_files = glob.glob("Cov2_fastq/*/MBCS_freq.tsv")
  out_filename = "results/all_MBCS_freq.tsv"
  MBCS_dict = MBCS_freq_to_dict(freq_files)
  MBCS_list = []
  [MBCS_list.extend(MBCS_dict[sampleID].keys()) for sampleID in MBCS_dict.keys()]
  MBCS_list = list(set(MBCS_list))
  write_MBCS_freq(MBCS_list, MBCS_dict, out_filename)

if __name__ == "__main__":
    main()
