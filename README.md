## CLI tools for variant calling pipeline

```commandline
git clone https://github.com/linyc74/variant.git

python variant [COMMAND] [ARGS]
```

Help messages:

```commandline
python variant filtering -h
python variant picking -h
python variant vcf2csv -h
```

### Variant Filtering

The `filtering` mode flags variants and then remove variants based on flags:

```commandline
python variant filtering \
    --input-vcf input.vcf \
    --output-vcf output.vcf \
    --variant-flagging-criteria "LOW_DP: DP<20, HIGH_MQ: MQ>=30" \
    --variant-removal-flags panel_of_normal,LOW_DP
```

### Variant Picking

The `picking` mode picks variants from multiple vcfs:

```commandline
python variant picking \
    --ref-fa hg38.fa \     
    --mutect2 mutect2.vcf \
    --muse muse.vcf \      
    --lofreq lofreq.vcf \  
    --output-vcf output.vcf \
    --min-snv-callers 2 \
    --min-indel-callers 1
```

Available caller intputs for the `picking` mode include:
- `--mutect2`
- `--haplotype-caller`
- `--muse`
- `--varscan`
- `--vardict`
- `--lofreq`
- `--somatic-sniper`

### VCF to CSV

The `vcf2csv` mode parses VCF file into CSV format:

```commandline
python variant vcf2csv \
    --input-vcf input.vcf \
    --output-csv output.csv \
```
