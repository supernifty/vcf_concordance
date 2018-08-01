# vcf_concordance

Calculate vcf concordance at different levels of AF and DP and generate a heatmap. This is useful if you have replicates and want to determine a good AF and/or DP threshold, or assign a confidence level to reported variants.

## Installation
```
python3 -m venv venv
. ./venv/bin/activate
pip install -r requirements.txt
```

## Usage
Specify your VCF files, and the corresponding sample names:
```
python concordance.py --vcfs vcf_file_1.vcf vcf_file_2.vcf --samples sample_1 sample_2 
```

The software looks for the AF and DP values in the provided VCF. These are required.

There are numerous additional filtering options:
* SNVs only
* filter on bed file

To learn about all command options:
```
python concordance.py --help
```
