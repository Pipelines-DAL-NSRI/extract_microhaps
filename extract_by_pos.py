import argparse
import pandas as pd
import os
import subprocess
import glob

from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()

parser.add_argument("--vcf", default=None, help="the input file in VCF format.")
parser.add_argument("--markers", default=None, help="file (csv, xlsx, or ods) containing marker information from microhaplotypes to be extracted. \
    File should contain microhaplotype name, GRCh37/38 position, chromosome number, and rsID or marker name.")
parser.add_argument("--sep", default=';', help="the separator of snps and pos in the marker file.")
parser.add_argument("--col_markers", nargs="+", default=None, help="The name of columns to be used (e.g. MH,Chrom,Pos,ID). \
    Required columns are indicated above.")
# Requires Python > v3.9
parser.add_argument("--splitmh", action = argparse.BooleanOptionalAction, default = False, help="group output by microhaplotype.")
parser.add_argument("--outfile", default="merged_data", help="the prefix of the filename of the merged dataset (will be in VCF format).")
parser.add_argument("--dir", default=None, help="the directory where all the files will be downloaded")
parser.add_argument("--clean", action = argparse.BooleanOptionalAction, default = True, help="remove intermediate files.")

args = parser.parse_args()

if args.dir:
    directory = args.dir
else:
    directory = "."

BASH_EXTRACT = os.path.abspath('./src/extract.sh')
BASH_ADD_ID = os.path.abspath('./src/add_id.sh')
BASH_MERGE = os.path.abspath('./src/merge.sh')
BASH_UPDATE = os.path.abspath('./src/update_id.sh')
BASH_EXTRACT_MH = os.path.abspath('./src/extract_by_snp.sh')
SNPS_LIST = defaultdict(list)
RENAME_RSID = 'rename_files.txt'
MERGE_FILES = 'merge_files.txt'
merged_file = f'{directory}/merged_dataset.vcf'
final_file = f'{directory}/{args.outfile}'

if os.path.exists(RENAME_RSID):
    os.remove(RENAME_RSID)

if os.path.exists(MERGE_FILES):
    os.remove(MERGE_FILES)

if os.path.exists(merged_file):
    os.remove(merged_file)

out_extracted_list = []
base_extracted = []
chrom_pos = []
rsid_list = []

def pairs(seq1, seq2):
   return [
   (seq1[idx], seq2[idx])
   for idx in range(len(seq1))
   ]

cols = [str(x) for x in args.col_markers[0].split(',')]
if len(cols) < 4:
    raise ValueError("Marker file should contain at least four columns.")

# Read marker file
print("READING MARKER FILE")
marker_file = Path(args.markers)
marker_ext = marker_file.suffix

if marker_ext == ".ods":
    marker_meta = pd.read_excel(marker_file, engine = 'odf', usecols = cols)
if marker_ext == ".xlsx":
    marker_meta = pd.read_excel(marker_file, usecols = cols)
if marker_ext == ".csv":
    marker_meta = pd.read_csv(marker_file, names = cols)

for idx, conts in marker_meta.iterrows():
   chrom = str(conts[1])[3:]
   pos_list = {}
   pos = str(conts[2]).split(args.sep)
   ids = str(conts[3]).split(args.sep)
   if len(pos) != len(ids):
      raise ValueError("Length of POS and rsID for the microhaplotype is not equal.")
   else:
      pos_list[conts[0]] = pairs(pos, ids)
      SNPS_LIST[chrom].append(pos_list)

print("EXTRACTING SNPS")
for key, val in SNPS_LIST.items(): # key is name is chrom
   for pair in val: # gets the microhaplotypes and pairs of pos, rsID
      for microhaps, pos_pairs in pair.items(): 
         for pos, rsid in pos_pairs:
            out_extracted = f'{directory}/chrom{key}_{rsid}'
            out_extracted_list.append(out_extracted)
            chrom_pos.append(f'{key}:{pos}')
            rsid_list.append(rsid)
            base_extracted.append(f'chrom{key}_{rsid}')            
            # Extract the SNPs
            subprocess.run([BASH_EXTRACT, args.vcf, key, pos, pos, out_extracted])

#            print("SETTING MISSING IDS TO CHR:POS")
#            # Set missing var ids to change ID
#            subprocess.run([BASH_ADD_ID, out_extracted, out_renamed])
#            out_renamed_v2 = out_renamed + ".vcf"
#            update = (f'{key}:{pos}',) + (f'{rsid}',)
#            with open(RENAME_RSID, 'a') as file:
#                file.write("\t".join(update) + "\n")
#            with open(MERGE_FILES, 'a') as file:
#                file.write(out_renamed_v2 + "\n")

extracted_files = [str(file) for file in Path(directory).glob("*.psam") if 'temporary' not in str(file)]

print("SETTING MISSING IDS TO CHR:POS")
for idx in range(len(out_extracted_list)):
    e_path = out_extracted_list[idx]
    updated = chrom_pos[idx]
    rsid = rsid_list[idx]
    base_file = base_extracted[idx]

    for file in extracted_files:
        f = file.strip()[:-5]
        if base_file == f:
            print(base_file, f)
            snp_update = (f'{updated}',) + (f'{rsid}',)
            out_renamed = f'{e_path}_updated'
            subprocess.run([BASH_ADD_ID, e_path, out_renamed])
            out_renamed_v2 = out_renamed + ".vcf"

            with open(RENAME_RSID, 'a') as file:
                file.write("\t".join(snp_update) + "\n")
            with open(MERGE_FILES, 'a') as file:
                file.write(out_renamed_v2 + "\n")

print("MERGING FILES")
subprocess.run([BASH_MERGE, MERGE_FILES, merged_file])

print("ADDING RSID TO MARKERS")
subprocess.run([BASH_UPDATE, merged_file, RENAME_RSID, final_file])

if args.splitmh:
    file = final_file + ".vcf"
    snps_list = {}
    for idx, conts in marker_meta.iterrows():
        chrom = str(conts[1])[3:]
        pos_list = {}
        ids = str(conts[3]).split(args.sep)
        # pair pos and ids first
        pos_list[conts[0]] = ids
        if chrom not in snps_list.keys():
            snps_list[chrom] = []
        else:
            snps_list[chrom].append(pos_list)

    for key, val in snps_list.items(): # key is name is chrom
       for pair in val: # gets the microhaplotypes and pairs of pos, rsID
          for microhaps, pos_pairs in pair.items(): 
             snps = ",".join(pos_pairs)
             name = f'{directory}/{microhaps}'
             subprocess.run([BASH_EXTRACT_MH, file, snps, name])

if args.clean:
    if os.path.exists(merged_file):
        os.remove(merged_file)
    for file in Path(directory).glob('*.pgen'):
        file.unlink()
    for file in Path(directory).glob('*_updated.vcf'):
        file.unlink()
    for file in Path(directory).glob('*.psam'):
        file.unlink()
    for file in Path(directory).glob('*.pvar'):
        file.unlink()
    for file in Path(directory).glob('*.log'):
        file.unlink()


for file in Path(directory).glob('*-temporary*'):
    file.unlink()