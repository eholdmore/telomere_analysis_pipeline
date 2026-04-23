# Telomere Maintenance Gene Variant Pipeline

Snakemake + Slurm pipeline for calling, annotating, and classifying variants
in **TERT, TERC, POT1, and DKC1** from tumor panel sequencing data (hg19).

---

## Overview

```
transferred_samples.csv  (manifest)
        │
        ▼
1. collect_paths          ──► file_paths_mapping.csv
        │
        ├──► 2. run_haplotypecaller (per sample, Slurm array)
        │           │
        │    3. filter_hc (per sample)
        │           │
        ├──► 4. parse_hs_metrics
        │           │
        └──────────────────────────────────────┐
                                               ▼
                                   5. hc_pre_vep
                                      (HC TSVs + M2 cache)
                                               │
                                   6. run_vep (novel variants only)
                                               │
                                   7. hc_post_vep
                                      (merge VEP + M2 + HS metrics)
                                               │
                                   8. annotate_cn_and_classify
                                      (RobustCNV + anchor VAF)
                                               │
                                   9. qc_and_finalize
                                      ┌────────┴────────┐
                                  table          figures + summaries
```

---

## Directory layout

```
telomere_pipeline/
├── Snakefile
├── config/
│   └── config.yaml          ← edit paths here
├── scripts/
│   ├── collect_paths.py
│   ├── hs_metrics_summary.py
│   ├── hc_pre_vep.py
│   ├── hc_post_vep.py
│   ├── get_cn_info.py
│   └── qc_and_finalize.py
├── workflow/
│   └── slurm_profile.yaml   ← copy to ~/.config/snakemake/slurm/config.yaml
└── setup_venv.sh
```

---

## Quick-start

### 1. Install the Python virtual environment
```bash
bash setup_venv.sh            # creates ~/telomere_numpy_venv
```

### 2. Configure the pipeline
Edit `config/config.yaml` — at minimum:

| Key | Description |
|-----|-------------|
| `work_dir` | Where all intermediates and outputs are written |
| `data_root` | Root of the `/data` filesystem containing sample subdirectories |
| `manifest` | Path to `transferred_samples.csv` |
| `mutect2_cache` | Path to `variant_classification_table_v3.csv` |
| `gatk_sif` | Path to the GATK 4.5 Singularity image |
| `vep_sif` | Path to the VEP Singularity image |
| `vep_cache_dir` / `vep_plugin_dir` / `vep_plugin_data` | VEP data paths |

### 3. Install the Slurm profile
```bash
mkdir -p ~/.config/snakemake/slurm
cp workflow/slurm_profile.yaml ~/.config/snakemake/slurm/config.yaml
```

### 4. Dry-run
```bash
source ~/telomere_numpy_venv/bin/activate
snakemake --profile slurm --configfile config/config.yaml -n
```

### 5. Submit
```bash
snakemake --profile slurm --configfile config/config.yaml
```

---

## Outputs

All outputs land in `{work_dir}/output/`:

| File | Description |
|------|-------------|
| `final_variant_table.csv` / `.parquet` | Filtered, annotated variant table |
| `qc_summary.csv` | Per-sample QC flags and coverage metrics |
| `variant_summary.csv` | Counts by gene × classification × consequence |
| `figures/vaf_by_classification.png` | VAF histogram coloured by Germline/Somatic/Ambiguous |
| `figures/cn_distribution.png` | Gene copy number (Ct) histogram |
| `figures/anchor_vaf.png` | Per-sample anchor VAF distribution |
| `figures/coverage_vs_20x.png` | Mean coverage vs. % bases at 20× (QC scatter) |
| `figures/qd_vs_qual.png` | QD by HC variant quality status |

### Final table columns
```
SAMPLE_ID, SAMPLE_QC_STATUS, HC_VARIANT_QUALITY, CHROM, POS, REF, ALT,
SYMBOL.PRIMARY, CALLER, QUAL, QD, MQRS, RPRS, AD, DP, DISEASE_GROUP,
VARIANT_ID, VAF, Consequence.PRIMARY, HGVSp.PRIMARY, HGVSc.PRIMARY,
Feature.PRIMARY, CDS_position.PRIMARY, Protein_position.PRIMARY, Codons.PRIMARY,
IMPACT.PRIMARY, STRAND.PRIMARY, HGNC_ID.PRIMARY, EXON.PRIMARY,
gnomad_non_cancer_AF_popmax, cosmic_CNT, dbSNP, maxMAF,
Ct, sample_anchor, Multiplicity,
Classification_0.2, Classification_0.25, Classification_0.3,
MEAN_TARGET_COVERAGE, PCT_TARGET_BASES_20X, PCT_TARGET_BASES_30X,
PCT_EXC_DUPE, MEDIAN_TARGET_COVERAGE,
FLAG_LOW_COVG, FLAG_LOW_20X, FLAG_HIGH_DUPE, NUM_SAMPLE_FLAGS,
FLAG_HC_LOW_QD, FLAG_HC_LOW_QUAL, FLAG_HC_MQRS, FLAG_HC_RPRS, NUM_FLAGS,
GERMLINE_PON, GERMLINE_RESCUE,
[REVEL / AlphaMissense / SpliceAI columns from VEP plugins]
```

---

## QC thresholds

### Sample-level
| Flag | Threshold |
|------|-----------|
| `FLAG_LOW_COVG` | Mean target coverage < 100× |
| `FLAG_LOW_20X` | % bases at 20× < 95 % |
| `FLAG_HIGH_DUPE` | Duplicate rate > 80 % |

`SAMPLE_QC_STATUS`: PASS (0 flags) / WARN (1 flag) / FAIL (≥ 2 flags)

### Variant-level
| Flag | Threshold |
|------|-----------|
| `FLAG_HC_LOW_QD` | QD < 3 |
| `FLAG_HC_LOW_QUAL` | QUAL < 30 |
| `FLAG_HC_MQRS` | MQRankSum < −4 |
| `FLAG_HC_RPRS` | ReadPosRankSum < −3 |

`HC_VARIANT_QUALITY`: HighConf (0 flags) / LowConf (≥ 1 flag)

### Hard filter applied to final table
- `SAMPLE_QC_STATUS == PASS`
- `HC_VARIANT_QUALITY == HighConf`
- `GERMLINE_PON == FALSE`
- Alternate allele read count `> 9`

---

## Classification model

For each variant, the ploidy-corrected expected somatic VAF is:

```
p_c_vaf = sample_anchor / (Ct / 2)
```

where `sample_anchor` is the median VAF of common germline heterozygous SNPs
(gnomAD non-cancer popmax > 0.01 %, VAF 0.30–0.70) and `Ct` is the
RobustCNV gene-level copy number.

Variants are classified at three thresholds (0.20 / 0.25 / 0.30):
- **Germline** — VAF > 0.85, or |VAF − p_c_vaf| ≤ 0.12
- **Somatic**  — |VAF − p_c_vaf| > threshold AND VAF < sample_anchor
- **Ambiguous** — everything else

---

## Resource requirements (approximate, 30 k samples)

| Step | CPUs | RAM | Wall time |
|------|------|-----|-----------|
| collect_paths | 1 | 2 GB | 5 min |
| run_haplotypecaller | 4 | 16 GB | 12 h |
| filter_hc | 1 | 4 GB | 1 h |
| parse_hs_metrics | 1 | 4 GB | 30 min |
| hc_pre_vep | 4 | 32 GB | 2 h |
| run_vep | 4 | 8 GB | 2 h |
| hc_post_vep | 4 | 32 GB | 2 h |
| annotate_cn_and_classify | 4 | 32 GB | 3 h |
| qc_and_finalize | 4 | 16 GB | 1 h |

HaplotypeCaller and filter_hc run in parallel across all samples (up to 500
concurrent Slurm jobs by default — tune `jobs:` in the profile as needed).
