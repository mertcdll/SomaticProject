# Somatic Variant Calling Pipeline

This repository provides a configuration-driven somatic variant calling workflow built around GATK Mutect2. It supports optional paired tumor-normal calling, scattering intervals for parallel execution, duplicate marking/removal, contamination estimation, read-orientation modeling, and alignment-artifact filtering. The pipeline also includes additional scripts for Panel of Normals (PoN) generation and CNV analysis.

## Features

* Reference indexing with BWA
* Alignment with BWA-MEM or DRAGEN-OS
* SAM to BAM conversion, duplicate marking/removal (sambamba or GATK MarkDuplicatesSpark)
* Interval scattering using GATK SplitIntervals
* Parallel Mutect2 execution across scattered intervals
* Optional contamination estimation with GetPileupSummaries and CalculateContamination
* Optional LearnReadOrientationModel (F1R2)
* Merge scattered VCFs and Mutect2 stats
* FilterMutectCalls and FilterAlignmentArtifacts
* Final PASS-only VCF export and tabix indexing

## Requirements

Command-line tools available in PATH:

* bwa
* samtools
* sambamba (unless using GATK MarkDuplicatesSpark)
* gatk (GATK4)
* bcftools
* tabix

Optional:

* dragen-os and a DRAGEN hash table (if use_dragen is enabled)

Python:

* Python 3.8+
* No external Python packages required beyond the standard library

## Input

The pipeline is driven by a single JSON configuration file. It expects tumor samples (and optionally normal samples) with FASTQ pairs, reference resources, and runtime settings.

Count of tumor_samples and normal_samples should match by index when running paired mode.

## Usage

Run:

python somaticpip.py --config config.json

## Configuration

Minimal configuration outline:

* ref_fasta: reference FASTA
* ref_dict: reference sequence dictionary
* algotype: bwa index algorithm (e.g. bwtsw, is)
* output_folder: base output directory
* run_name: run identifier used for output structuring
* aligner_threads: thread count for alignment and processing
* index_fasta: boolean; if true, run bwa index
* interval_list: target intervals (BED/interval list as used by GATK)
* scatter_count: number of scattered interval shards
* tumor_samples: list of tumor samples with sample_id, r1, r2
* normal_samples: optional list of normal samples with sample_id, r1, r2
* germline_resource: gnomAD or similar germline resource VCF for Mutect2
* index_image: BWA-MEM index image required by FilterAlignmentArtifacts
* panel_of_normals: optional PoN VCF for Mutect2
* variants_for_contamination: optional common variants VCF used by GetPileupSummaries
* downsampling_stride: optional Mutect2 parameter
* max_reads_per_alignment_start: optional Mutect2 parameter
* max_suspicious_reads_per_alignment_start: optional Mutect2 parameter
* max_population_af: optional Mutect2 parameter
* lrom: boolean; enable LearnReadOrientationModel
* interval_padding: integer; interval padding passed to Mutect2
* use_dragen: boolean; align with dragen-os instead of bwa
* remove_all_duplicates: boolean; remove duplicates during deduplication
* remove_sequencing_duplicates: boolean; applicable to MarkDuplicatesSpark
* use_gatk_mark_duplicates: boolean; use GATK MarkDuplicatesSpark instead of sambamba

## Outputs

For each tumor sample (and paired normal if provided), outputs are written under:

<output_folder>/<run_name>/

Alignment and preprocessing:

* tumors/<tumor_id>/sorted_filtered_<tumor_id>.bam and .bai
* normals/<normal_id>/sorted_filtered_<normal_id>.bam and .bai (if provided)
* <sample_dir>/<sample_id>_dup_metrics.txt (duplicate metrics)

Variant calling (per tumor sample):

* scattered_intervals/ (interval shards from SplitIntervals)
* mutect2_outputs/shard_*/output.vcf
* mutect2_outputs/shard_*/f1r2.tar.gz (if enabled by Mutect2 options)
* mutect2_outputs/shard_*/normal-pileups.table (if normal + variants_for_contamination)
* mutect2_outputs/shard_*/tumor-pileups.table (if variants_for_contamination)

Merged outputs:

* merged_outputs/<tumor_name>_clean-unfiltered.vcf.gz
* merged_outputs/merged.stats
* merged_outputs/artifact-priors.tar.gz (if lrom enabled)
* merged_outputs/<tumor_name>_clean.tsv (tumor pileups, if variants_for_contamination)
* merged_outputs/<normal_name>_clean.tsv (normal pileups, if matched normal)
* merged_outputs/contamination.table and segments.table (if variants_for_contamination)
* merged_outputs/<tumor_name>_clean-flagged.vcf.gz
* merged_outputs/<tumor_name>_clean-filtered.vcf.gz
* merged_outputs/<tumor_name>_clean-further-filtered.vcf.gz and .tbi (final PASS-only VCF)

## Workflow Summary

1. Optional: BWA index reference
2. Align tumor FASTQs (BWA-MEM or DRAGEN-OS)
3. Convert SAM to BAM, mark/remove duplicates, sort and index
4. Optional: align and preprocess matched normal
5. Split target intervals into shards
6. Run Mutect2 per shard in parallel
7. Sort shard VCFs
8. Optional: GetPileupSummaries for tumor (and normal if present)
9. Optional: LearnReadOrientationModel from per-shard F1R2 archives
10. Merge shard VCFs and stats
11. Optional: gather pileup summaries and calculate contamination
12. FilterMutectCalls to flag calls
13. FilterAlignmentArtifacts and export final PASS-only VCF

## Panel of Normals and CNV Scripts

This repository also includes additional scripts for:

* Panel of Normals (PoN) generation
* CNV analysis

Refer to the corresponding script files for their specific configuration and usage.

## Notes

* FilterAlignmentArtifacts requires a BWA-MEM index image. Ensure index_image points to a valid image created for the same reference.
* When using paired mode, normal_samples must be aligned by index to tumor_samples.
* For DRAGEN alignment, a valid DRAGEN hash table must be provided and use_dragen enabled.
