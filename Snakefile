"""
Telomere Maintenance Gene Variant Pipeline
==========================================
Processes tumor panel sequencing data to produce an annotated table
of variants in TERT, TERC, POT1, and DKC1.

Usage:
    snakemake --profile slurm --configfile config/config.yaml

Steps:
    1.  collect_paths      - Build file paths mapping from sample manifest
    2.  run_haplotypecaller - Run GATK HaplotypeCaller per sample (Slurm array)
    3.  filter_hc          - Filter HC VCFs to target gene regions
    4.  hs_metrics_summary - Parse Picard HS metrics per sample
    5.  hc_pre_vep         - Merge HC + Mutect2 cache, prep novel VCF for VEP
    6.  run_vep            - Run VEP with REVEL/SpliceAI/AlphaMissense plugins
    7.  hc_post_vep        - Merge VEP annotations back, build master call table
    8.  annotate_cn        - Add copy number + anchor VAF, classify variants
    9.  qc_and_finalize    - Apply QC flags, filters, generate outputs
"""

import os
import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ── Derived paths ──────────────────────────────────────────────────────────────
WORK_DIR    = config["work_dir"]
DATA_ROOT   = config["data_root"]
SCRATCH_DIR = config["scratch_dir"]
REF_DIR     = config["ref_dir"]
SCRIPTS_DIR = config["scripts_dir"]

MANIFEST     = config["manifest"]
M2_CACHE     = config["mutect2_cache"]
VENV_ACTIVATE = config["venv_activate"]

# ── Target genes / regions (hg19) ─────────────────────────────────────────────
TARGET_GENES   = ["TERT", "TERC", "POT1", "DKC1"]
TARGET_REGIONS = config["target_regions"]  # e.g. "5:1251147-1297184,..."

# ── Load sample manifest (lazy - only resolved at dag build time) ──────────────
def get_samples():
    df = pd.read_csv(MANIFEST)
    return df["cbio_sample_id"].tolist()

SAMPLES = get_samples()

# ── Final outputs ──────────────────────────────────────────────────────────────
rule all:
    input:
        os.path.join(WORK_DIR, "intermediates", "all_variants_qc_flagged.parquet"),
        os.path.join(WORK_DIR, "output", "final_variant_table.csv"),
        os.path.join(WORK_DIR, "output", "final_variant_table.parquet"),
        os.path.join(WORK_DIR, "output", "qc_summary.csv"),
        os.path.join(WORK_DIR, "output", "figures", "vaf_by_classification.png"),
        os.path.join(WORK_DIR, "output", "figures", "cn_distribution.png"),
        os.path.join(WORK_DIR, "output", "figures", "anchor_vaf.png"),
        os.path.join(WORK_DIR, "output", "figures", "coverage_vs_20x.png"),
        os.path.join(WORK_DIR, "output", "figures", "qd_vs_qual.png"),
        os.path.join(WORK_DIR, "output", "variant_summary.csv"),


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Collect file paths for all samples
# ══════════════════════════════════════════════════════════════════════════════
rule collect_paths:
    input:
        manifest = MANIFEST
    output:
        paths_csv = os.path.join(WORK_DIR, "intermediates", "file_paths_mapping.csv")
    params:
        data_root = DATA_ROOT,
        scripts   = SCRIPTS_DIR
    log:
        os.path.join(WORK_DIR, "logs", "collect_paths.log")
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/collect_paths.py \
            --manifest {input.manifest} \
            --data-root {params.data_root} \
            --output {output.paths_csv} \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Run GATK HaplotypeCaller per sample
# ══════════════════════════════════════════════════════════════════════════════
rule run_haplotypecaller:
    input:
        paths_csv = os.path.join(WORK_DIR, "intermediates", "file_paths_mapping.csv"),
        ref       = os.path.join(REF_DIR, "Homo_sapiens_assembly19.fasta"),
        intervals = os.path.join(REF_DIR, config["intervals_file"]),
        dbsnp     = os.path.join(REF_DIR, "All_20180423.vcf.gz"),
    output:
        vcf = os.path.join(WORK_DIR, "haplotypecaller", "{sample}", "{sample}.vcf.gz"),
        tbi = os.path.join(WORK_DIR, "haplotypecaller", "{sample}", "{sample}.vcf.gz.tbi"),
    params:
        gatk_sif = config["gatk_sif"],
        out_dir  = os.path.join(WORK_DIR, "haplotypecaller", "{sample}"),
        scratch  = SCRATCH_DIR,
        phshome  = os.path.join(os.environ.get("HOME", "/PHShome/eh011")),
    log:
        os.path.join(WORK_DIR, "logs", "haplotypecaller", "{sample}.log")
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 720,   # 12 h
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        set -euo pipefail
        mkdir -p {params.out_dir}

        # Resolve BAM path for this sample from the mapping CSV
        BAM=$(python -c "
import csv, sys
sample = '{wildcards.sample}'
with open('{input.paths_csv}') as f:
    for row in csv.DictReader(f):
        if row['cbio_sample_id'] == sample:
            print(row['bam_path'])
            sys.exit(0)
sys.exit(1)
")

        module load singularity/latest
        singularity exec --bind /data,{params.scratch},{params.phshome} \
            {params.gatk_sif} \
            gatk --java-options "-Xmx12G" HaplotypeCaller \
            --reference {input.ref} \
            --input "$BAM" \
            --output {output.vcf} \
            --intervals {input.intervals} \
            --dbsnp {input.dbsnp} \
            --standard-min-confidence-threshold-for-calling 30 \
            --read-validation-stringency LENIENT \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Filter HC VCFs to target regions and extract per-variant metrics
# ══════════════════════════════════════════════════════════════════════════════
rule filter_hc:
    input:
        vcf       = os.path.join(WORK_DIR, "haplotypecaller", "{sample}", "{sample}.vcf.gz"),
        paths_csv = os.path.join(WORK_DIR, "intermediates", "file_paths_mapping.csv"),
    output:
        tsv = os.path.join(WORK_DIR, "filtered_hc", "{sample}_filtered.tsv"),
    params:
        regions  = TARGET_REGIONS,
        out_dir  = os.path.join(WORK_DIR, "filtered_hc"),
    log:
        os.path.join(WORK_DIR, "logs", "filter_hc", "{sample}.log")
    threads: 1
    resources:
        mem_mb  = 4000,
        runtime = 60,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        set -euo pipefail
        mkdir -p {params.out_dir}

        BRIEFCASE=$(python -c "
import csv, sys
with open('{input.paths_csv}') as f:
    for row in csv.DictReader(f):
        if row['cbio_sample_id'] == '{wildcards.sample}':
            print(row['briefcase'])
            sys.exit(0)
sys.exit(1)
")

        module load bcftools/1.11
        bcftools view -r "{params.regions}" "{input.vcf}" | \
        bcftools norm -m-both | \
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%QUAL\t%INFO/QD\t%INFO/MQRankSum\t%INFO/ReadPosRankSum\t[%AD]\t[%DP]\n' | \
        awk -v sid="{wildcards.sample}" -v bc="$BRIEFCASE" '
        BEGIN{{OFS="\t"}}
        {{
            split($10, ad, ",");
            print $1, $2, $3, $4, $5, $6, $7, $8, $9, ad[2], $11, sid, bc
        }}' > "{output.tsv}" 2> "{log}"
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — Parse HS metrics for each sample
# ══════════════════════════════════════════════════════════════════════════════
rule parse_hs_metrics:
    input:
        paths_csv = os.path.join(WORK_DIR, "intermediates", "file_paths_mapping.csv"),
    output:
        metrics_csv = os.path.join(WORK_DIR, "intermediates", "hs_metrics_summary.csv"),
    params:
        scripts = SCRIPTS_DIR,
    log:
        os.path.join(WORK_DIR, "logs", "hs_metrics_summary.log")
    resources:
        mem_mb  = 4000,
        runtime = 30,
        slurm_partition = config.get("login_partition", "login"),
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/hs_metrics_summary.py \
            --paths-csv {input.paths_csv} \
            --output {output.metrics_csv} \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — Merge HC TSVs + Mutect2 cache, prep novel variants for VEP
# ══════════════════════════════════════════════════════════════════════════════
rule hc_pre_vep:
    input:
        hc_tsvs  = expand(
            os.path.join(WORK_DIR, "filtered_hc", "{sample}_filtered.tsv"),
            sample=SAMPLES
        ),
        m2_cache = M2_CACHE,
    output:
        known_parquet  = os.path.join(WORK_DIR, "intermediates", "hc_obs_known.parquet"),
        novel_parquet  = os.path.join(WORK_DIR, "intermediates", "hc_obs_novel_waiting.parquet"),
        novel_vcf      = os.path.join(WORK_DIR, "intermediates", "hc_novel_variants_to_annotate.vcf"),
    params:
        scripts  = SCRIPTS_DIR,
        tsv_glob = os.path.join(WORK_DIR, "filtered_hc", "*.tsv"),
    log:
        os.path.join(WORK_DIR, "logs", "hc_pre_vep.log")
    resources:
        mem_mb  = 32000,
        runtime = 120,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/hc_pre_vep.py \
            --tsv-glob "{params.tsv_glob}" \
            --m2-cache {input.m2_cache} \
            --output-known {output.known_parquet} \
            --output-novel {output.novel_parquet} \
            --output-vcf {output.novel_vcf} \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6 — Run VEP on novel variants
# ══════════════════════════════════════════════════════════════════════════════
rule run_vep:
    input:
        vcf     = os.path.join(WORK_DIR, "intermediates", "hc_novel_variants_to_annotate.vcf"),
        ref     = os.path.join(REF_DIR, "Homo_sapiens_assembly19.fasta"),
    output:
        vep_out = os.path.join(WORK_DIR, "intermediates", "hc_novel_variants_annotated.txt"),
    params:
        vep_sif        = config["vep_sif"],
        vep_cache_dir  = config["vep_cache_dir"],
        vep_plugin_dir = config["vep_plugin_dir"],
        vep_plugin_data= config["vep_plugin_data"],
        work_dir_abs   = WORK_DIR,
        ref_dir_abs    = REF_DIR,
    log:
        os.path.join(WORK_DIR, "logs", "run_vep.log")
    threads: 4
    resources:
        mem_mb  = 8000,
        runtime = 120,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        set -euo pipefail
        module load singularity/latest
        singularity exec \
            --bind {params.ref_dir_abs}:/root/reference \
            --bind {params.work_dir_abs}/intermediates:/root/work \
            {params.vep_sif} \
            vep \
            --input_file  /root/work/hc_novel_variants_to_annotate.vcf \
            --output_file /root/work/hc_novel_variants_annotated.txt \
            --format vcf \
            --tab \
            --force_overwrite \
            --no_stats \
            --offline \
            --cache \
            --dir_cache    /root/reference/vep_cache \
            --dir_plugins  {params.vep_plugin_dir} \
            --species      homo_sapiens \
            --assembly     GRCh37 \
            --fasta        /root/reference/Homo_sapiens_assembly19.fasta \
            --fork         {threads} \
            --hgvs \
            --symbol \
            --numbers \
            --pick \
            --pick_order   biotype,rank,tsl,canonical,length \
            --no_escape \
            --plugin REVEL,file={params.vep_plugin_data}/revel_with_transcript_ids \
            --plugin AlphaMissense,file={params.vep_plugin_data}/AlphaMissense_hg19.tsv.gz \
            --plugin SpliceAI,file={params.vep_plugin_data}/spliceai_scores.raw.snv.hg19.vcf.gz \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 — Merge VEP annotations back; build master call table
# ══════════════════════════════════════════════════════════════════════════════
rule hc_post_vep:
    input:
        known_parquet  = os.path.join(WORK_DIR, "intermediates", "hc_obs_known.parquet"),
        novel_parquet  = os.path.join(WORK_DIR, "intermediates", "hc_obs_novel_waiting.parquet"),
        vep_out        = os.path.join(WORK_DIR, "intermediates", "hc_novel_variants_annotated.txt"),
        m2_cache       = M2_CACHE,
        metrics_csv    = os.path.join(WORK_DIR, "intermediates", "hs_metrics_summary.csv"),
    output:
        master_parquet = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_final.parquet"),
        master_csv     = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_final.csv"),
    params:
        scripts = SCRIPTS_DIR,
    log:
        os.path.join(WORK_DIR, "logs", "hc_post_vep.log")
    resources:
        mem_mb  = 32000,
        runtime = 120,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/hc_post_vep.py \
            --known-parquet  {input.known_parquet} \
            --novel-parquet  {input.novel_parquet} \
            --vep-output     {input.vep_out} \
            --m2-cache       {input.m2_cache} \
            --metrics-csv    {input.metrics_csv} \
            --output-parquet {output.master_parquet} \
            --output-csv     {output.master_csv} \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 8 — Add copy number, anchor VAF, and classify variants
# ══════════════════════════════════════════════════════════════════════════════
rule annotate_cn_and_classify:
    input:
        master_parquet = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_final.parquet"),
        paths_csv      = os.path.join(WORK_DIR, "intermediates", "file_paths_mapping.csv"),
    output:
        classified_parquet = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_classified.parquet"),
        classified_csv     = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_classified.csv"),
    params:
        scripts = SCRIPTS_DIR,
    log:
        os.path.join(WORK_DIR, "logs", "annotate_cn.log")
    resources:
        mem_mb  = 32000,
        runtime = 180,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/get_cn_info.py \
            --input-parquet  {input.master_parquet} \
            --paths-csv      {input.paths_csv} \
            --output-parquet {output.classified_parquet} \
            --output-csv     {output.classified_csv} \
            > {log} 2>&1
        """


# ══════════════════════════════════════════════════════════════════════════════
# STEP 9 — QC flags, hard filters, final outputs and visualizations
# ══════════════════════════════════════════════════════════════════════════════
rule qc_and_finalize:
    input:
        classified_parquet = os.path.join(WORK_DIR, "intermediates", "master_telomere_calls_classified.parquet"),
    output:
        prefilter_parquet = os.path.join(WORK_DIR, "intermediates", "all_variants_qc_flagged.parquet"),
        final_csv         = os.path.join(WORK_DIR, "output", "final_variant_table.csv"),
        final_parquet     = os.path.join(WORK_DIR, "output", "final_variant_table.parquet"),
        qc_summary        = os.path.join(WORK_DIR, "output", "qc_summary.csv"),
        var_summary       = os.path.join(WORK_DIR, "output", "variant_summary.csv"),
        fig_vaf           = os.path.join(WORK_DIR, "output", "figures", "vaf_by_classification.png"),
        fig_cn            = os.path.join(WORK_DIR, "output", "figures", "cn_distribution.png"),
        fig_anchor        = os.path.join(WORK_DIR, "output", "figures", "anchor_vaf.png"),
        fig_covg          = os.path.join(WORK_DIR, "output", "figures", "coverage_vs_20x.png"),
        fig_qd            = os.path.join(WORK_DIR, "output", "figures", "qd_vs_qual.png"),
    params:
        scripts    = SCRIPTS_DIR,
        output_dir = os.path.join(WORK_DIR, "output"),
    log:
        os.path.join(WORK_DIR, "logs", "qc_and_finalize.log")
    resources:
        mem_mb  = 16000,
        runtime = 60,
        slurm_partition = config.get("partition", "compute"),
    shell:
        """
        source {VENV_ACTIVATE}
        python {params.scripts}/qc_and_finalize.py \
            --input-parquet     {input.classified_parquet} \
            --output-dir        {params.output_dir} \
            --prefilter-parquet {output.prefilter_parquet} \
            > {log} 2>&1
        """
