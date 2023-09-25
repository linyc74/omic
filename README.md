## OMIC CLI tools

```commandline
git clone https://github.com/linyc74/omic.git
git clone https://github.com/linyc74/omic.git

python omic [COMMAND] [ARGS]
```

Help messages:

```commandline
python omic variant-filtering -h
python omic variant-picking -h
python omic vcf2csv -h
```

### Variant Filtering

The `variant-filering` command flags variants and then remove variants based on flags:

```commandline
python omic variant-filtering \
    --input-vcf input.vcf \
    --output-vcf output.vcf \
    --variant-flagging-criteria "LOW_DP: DP<20, HIGH_MQ: MQ>=30" \
    --variant-removal-flags panel_of_normal,LOW_DP
```

### Variant Picking

The `variant-picking` command picks variants from multiple vcfs:

```commandline
python omic variant-picking \
    --ref-fa hg38.fa \     
    --mutect2 mutect2.vcf \
    --muse muse.vcf \      
    --lofreq lofreq.vcf \  
    --output-vcf output.vcf \
    --min-snv-callers 2 \
    --min-indel-callers 1
```

Available caller intputs for `variant-picking` include:
- `--mutect2`
- `--haplotype-caller`
- `--muse`
- `--varscan`
- `--vardict`
- `--lofreq`
- `--somatic-sniper`

### VCF to CSV

The `vcf2csv` command parses VCF file into CSV format:

```commandline
python omic vcf2csv \
    --input-vcf input.vcf \
    --output-csv output.csv
```
