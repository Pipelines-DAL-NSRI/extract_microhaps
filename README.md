# POS-based extraction of Microhaplotypes

# Requirements
- Python > v3.9
- Pip
- Install function requirements by doing pip install requirements.txt

# Parameters
```--vcf```: The VCF file  
```--markers```: File containing markers to be extracted. Should be in an acceptable format (.csv, .xlsx, or .ods). Requires at least four columns: microhaplotype name, GRCh37/38 position, chromosome number, and rsID. The positions and rsIDs within a microhaplotype are expected to be in a single cell.  
```--sep```: Separator of pos and ids in the marker file. Default is ';'.  
```--col_markers```: The name of columns to be used, should follow the order of microhaplotype name, chromosome number, position, and rsID/ID.  
```--splitmh```: Separate output by microhaplotypes. Use flag to set to True. Default is False. Files will be in VCF format.  
```--outfile```: The prefix of the merged output file. All microhaplotypes will be merged into a single VCF file.  
```--dir```: Directory to save all files. Default is working directory where python file is located.  
```--no-clean```: Remove intermediate files. Use flag to set to False.  

# Usage
```
python extract_by_pos.py --vcf input_file.vcf --markers list_of_microhaplotypes.csv --sep ; --col_names Name,Chrom,Pos,ID --splitmh --outfile merged_files
```

