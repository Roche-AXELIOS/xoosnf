# Profile Feature Comparison Tables

Each profile explicitly sets all parameters listed below unless marked with `—`, which indicates the parameter is not set by the profile and inherits the global default from `nextflow.config`.

## Germline Duplex, Somatic & cfDNA

| Parameter                                      | germline_wgs_duplex | germline_wgs_duplex_benchmark | somatic_wgs_tn_duplex | somatic_wgs_tn_duplex_benchmark | cfdna_wgs_duplex   | cfdna_wgs_duplex_benchmark |
| ---------------------------------------------- | ------------------- | ----------------------------- | --------------------- | ------------------------------- | ------------------ | -------------------------- |
| **Mode & Alignment**                           |                     |                               |                       |                                 |                    |                            |
| `variant_calling_mode`                         | `germline`          | `germline`                    | `somatic_tn`          | `somatic_tn`                    | `none`             | `none`                     |
| `aligner`                                      | `pb_giraffe`        | `pb_giraffe`                  | `pb_fq2bam`           | `pb_fq2bam`                     | `pb_fq2bam`        | `pb_fq2bam`                |
| `enable_pangenome_personalization`             | ✅                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| **Demultiplex**                                |                     |                               |                       |                                 |                    |                            |
| `exclude_partial_reads`                        | ✅                  | ✅                            | ✅                    | ✅                              | ✅                 | ✅                         |
| `enable_demux_strand_detect`                   | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| **Deduplication**                              |                     |                               |                       |                                 |                    |                            |
| `dedup_strategy`                               | `markdup`           | `markdup`                     | `markdup`             | `markdup`                       | `consensus`        | `consensus`                |
| `dedup_aligner`                                | —                   | —                             | —                     | —                               | —                  | —                          |
| `read_collapser_preset`                        | `wgs-duplex`        | `wgs-duplex`                  | `wgs-duplex`          | `wgs-duplex`                    | `wgs-duplex-cfdna` | `wgs-duplex-cfdna`         |
| `read_collapser_cluster_by_umi`                | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| **Subsampling**                                |                     |                               |                       |                                 |                    |                            |
| `enable_subsample_bam`                         | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `subsample_bam_coverage_target`                | `null`              | `30`                          | `null`                | `null`                          | `null`             | `null`                     |
| **Metrics**                                    |                     |                               |                       |                                 |                    |                            |
| `enable_alignment_metrics_autosomes`           | ✅                  | ✅                            | ✅                    | ✅                              | ✅                 | ✅                         |
| `enable_alignment_metrics_targets`             | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_alignment_metrics_high_confidence`     | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ✅                         |
| `enable_hp_metrics`                            | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `enable_gc_bias_metrics`                       | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `enable_hs_metrics`                            | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_pre_dedup_metrics`                     | ❌                  | ❌                            | ❌                    | ❌                              | ✅                 | ✅                         |
| `enable_pre_subsample_metrics`                 | ❌                  | ❌                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `alignment_metrics_min_base_quality_pre_dedup` | `30`                | `30`                          | —                     | —                               | —                  | —                          |
| `alignment_metrics_min_base_quality`           | `30`                | `30`                          | `30`                  | `30`                            | `30`               | `30`                       |
| **Contamination**                              |                     |                               |                       |                                 |                    |                            |
| `enable_contamination_estimation`              | ✅                  | ✅                            | ✅                    | ✅                              | ✅                 | ✅                         |
| `enable_tfe`                                   | ❌                  | ❌                            | ❌                    | ❌                              | ✅                 | ✅                         |
| **Tumor-normal mis-match Verification**        |                     |                               |                       |                                 |                    |                            |
| `enable_somatic_tn_pair_verification`          | ❌                  | ❌                            | —                     | —                               | ❌                 | ❌                         |
| **Variant Calling**                            |                     |                               |                       |                                 |                    |                            |
| `enable_small_variant_caller`                  | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `enable_small_variant_caller_retrain`          | ❌                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_copy_number_caller`                    | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `enable_mitochondrial_variant_caller`          | —                   | —                             | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_str_caller`                            | ✅                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_msi_caller`                            | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_merge_vcf`                             | ✅                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_vcfeval`                               | —                   | ✅                            | —                     | —                               | —                  | —                          |
| **Variant Calling Partitioning**               |                     |                               |                       |                                 |                    |                            |
| `enable_variant_calling_partitioning`          | ✅                  | ✅                            | ✅                    | ✅                              | —                  | —                          |
| `cpu_variant_calling_partition_count`          | `20`                | `20`                          | `40`                  | `40`                            | —                  | —                          |
| `gpu_variant_calling_partition_count`          | `1`                 | `1`                           | `4`                   | `4`                             | —                  | —                          |

MSI calling is currently opt-in and is not enabled by any production profile. When enabled for `somatic_tn`, the optional `msi_sites_thresholds_tsv` resource can be supplied; when enabled for `somatic_to`, the optional `msi_background_threshold_json` resource can be supplied. The infrequently changed MSI thresholds are handled internally by the module rather than exposed as profile parameters.

## Somatic Tumor-Only Duplex

| Parameter                                  | somatic_wgs_to_duplex | somatic_wgs_to_duplex_benchmark |
| ------------------------------------------ | --------------------- | ------------------------------- |
| **Mode & Alignment**                       |                       |                                 |
| `variant_calling_mode`                     | `somatic_to`          | `somatic_to`                    |
| `aligner`                                  | `pb_fq2bam`           | `pb_fq2bam`                     |
| `enable_pangenome_personalization`         | ❌                    | ❌                              |
| **Demultiplex**                            |                       |                                 |
| `exclude_partial_reads`                    | ✅                    | ✅                              |
| `enable_demux_strand_detect`               | ❌                    | ❌                              |
| **Deduplication**                          |                       |                                 |
| `dedup_strategy`                           | `markdup`             | `markdup`                       |
| `read_collapser_preset`                    | `wgs-duplex`          | `wgs-duplex`                    |
| `read_collapser_cluster_by_umi`            | ❌                    | ❌                              |
| **Subsampling**                            |                       |                                 |
| `enable_subsample_bam`                     | ❌                    | ✅                              |
| `subsample_bam_coverage_target`            | `null`                | `100`                           |
| **Metrics**                                |                       |                                 |
| `enable_alignment_metrics_autosomes`       | ✅                    | ✅                              |
| `enable_alignment_metrics_targets`         | ❌                    | ❌                              |
| `enable_alignment_metrics_high_confidence` | ✅                    | ✅                              |
| `enable_hp_metrics`                        | ❌                    | ✅                              |
| `enable_gc_bias_metrics`                   | ❌                    | ✅                              |
| `enable_hs_metrics`                        | ❌                    | ❌                              |
| `enable_pre_dedup_metrics`                 | ❌                    | ❌                              |
| `enable_pre_subsample_metrics`             | ❌                    | ✅                              |
| `alignment_metrics_min_base_quality`       | `30`                  | `30`                            |
| **Contamination**                          |                       |                                 |
| `enable_contamination_estimation`          | ✅                    | ✅                              |
| **Variant Calling**                        |                       |                                 |
| `enable_small_variant_caller`              | ✅                    | ✅                              |
| `enable_small_variant_caller_retrain`      | ❌                    | ❌                              |
| `enable_copy_number_caller`                | ✅                    | ✅                              |
| `enable_germline_filtering`                | ✅                    | ✅                              |
| `enable_str_caller`                        | ❌                    | ❌                              |
| `enable_merge_vcf`                         | ✅                    | ✅                              |
| **Variant Calling Partitioning**           |                       |                                 |
| `enable_variant_calling_partitioning`      | ✅                    | ✅                              |
| `cpu_variant_calling_partition_count`      | `40`                  | `40`                            |
| `gpu_variant_calling_partition_count`      | `4`                   | `4`                             |

## WGS Simplex

| Parameter                                      | wgs_simplex   | wgs_simplex_benchmark | germline_wgs_simplex | germline_wgs_simplex_benchmark |
| ---------------------------------------------- | ------------- | --------------------- | -------------------- | ------------------------------ |
| **Mode & Alignment**                           |               |                       |                      |                                |
| `variant_calling_mode`                         | `none`        | `none`                | `germline`           | `germline`                     |
| `aligner`                                      | `pb_fq2bam`   | `pb_fq2bam`           | `pb_giraffe`         | `pb_giraffe`                   |
| `enable_pangenome_personalization`             | ❌            | ❌                    | ✅                   | ✅                             |
| **Demultiplex**                                |               |                       |                      |                                |
| `exclude_partial_reads`                        | ✅            | ✅                    | ❌                   | ❌                             |
| `enable_demux_strand_detect`                   | ❌            | ❌                    | ❌                   | ❌                             |
| **Deduplication**                              |               |                       |                      |                                |
| `dedup_strategy`                               | `markdup`     | `markdup`             | `markdup`            | `markdup`                      |
| `dedup_aligner`                                | —             | —                     | —                    | —                              |
| `read_collapser_preset`                        | `wgs-simplex` | `wgs-simplex`         | `wgs-simplex`        | `wgs-simplex`                  |
| `read_collapser_cluster_by_umi`                | ❌            | ❌                    | ❌                   | ❌                             |
| **Subsampling**                                |               |                       |                      |                                |
| `enable_subsample_bam`                         | ❌            | —                     | ❌                   | ✅                             |
| `subsample_bam_coverage_target`                | `null`        | `null`                | `null`               | `40`                           |
| **Metrics**                                    |               |                       |                      |                                |
| `enable_alignment_metrics_autosomes`           | ✅            | ✅                    | ✅                   | ✅                             |
| `enable_alignment_metrics_targets`             | ❌            | ❌                    | ❌                   | ❌                             |
| `enable_alignment_metrics_high_confidence`     | ✅            | ✅                    | ✅                   | ✅                             |
| `enable_hp_metrics`                            | ❌            | ✅                    | ❌                   | ✅                             |
| `enable_gc_bias_metrics`                       | ❌            | ✅                    | ❌                   | ✅                             |
| `enable_hs_metrics`                            | ❌            | ❌                    | ❌                   | ❌                             |
| `enable_pre_dedup_metrics`                     | ❌            | ❌                    | ❌                   | ❌                             |
| `enable_pre_subsample_metrics`                 | ❌            | ❌                    | ❌                   | ❌                             |
| `alignment_metrics_min_base_quality_pre_dedup` | `0`           | `0`                   | `0`                  | `0`                            |
| `alignment_metrics_min_base_quality`           | `0`           | `0`                   | `0`                  | `0`                            |
| **Contamination**                              |               |                       |                      |                                |
| `enable_contamination_estimation`              | ❌            | ❌                    | ✅                   | ✅                             |
| **Variant Calling**                            |               |                       |                      |                                |
| `enable_small_variant_caller`                  | ❌            | ❌                    | ✅                   | ✅                             |
| `enable_small_variant_caller_retrain`          | ❌            | ❌                    | ❌                   | ✅                             |
| `enable_copy_number_caller`                    | ❌            | ❌                    | ✅                   | ✅                             |
| `enable_str_caller`                            | ❌            | ❌                    | ❌                   | ❌                             |
| `enable_merge_vcf`                             | ❌            | ❌                    | ✅                   | ✅                             |
| `enable_vcfeval`                               | —             | —                     | —                    | ✅                             |
| `small_variant_caller_sequencing_protocol`     | —             | —                     | `simplex`            | `simplex`                      |

SVC Germline WGS `filter_variants` uses the inherited
`small_variant_caller_germline_model_type` value (`large` by default).
SVC Germline WGS model retraining always uses the upstream `small` model family.
The germline simplex profiles retain partial reads and use XOOS read_collapser
MARKDUP.
They use `pb_giraffe` with pangenome personalization and therefore require
`unfiltered_giraffe_index_dir` and `unfiltered_giraffe_reference_linearized`.
When copy number calling is enabled, Parabricks Giraffe emits up to three
multimapping alignments for secondary-alignment rescue.

## TE Duplex

| Parameter                                      | te_duplex   | te_duplex_benchmark | te_duplex_umi | te_duplex_umi_benchmark |
| ---------------------------------------------- | ----------- | ------------------- | ------------- | ----------------------- |
| **Mode & Alignment**                           |             |                     |               |                         |
| `variant_calling_mode`                         | `none`      | `none`              | `none`        | `none`                  |
| `aligner`                                      | `pb_fq2bam` | `pb_fq2bam`         | `pb_fq2bam`   | `pb_fq2bam`             |
| `enable_pangenome_personalization`             | ❌          | ❌                  | ❌            | ❌                      |
| **Demultiplex**                                |             |                     |               |                         |
| `exclude_partial_reads`                        | ❌          | ❌                  | ❌            | ❌                      |
| `enable_demux_strand_detect`                   | ❌          | ❌                  | ❌            | ❌                      |
| **Deduplication**                              |             |                     |               |                         |
| `dedup_strategy`                               | `consensus` | `consensus`         | `consensus`   | `consensus`             |
| `dedup_aligner`                                | —           | —                   | —             | —                       |
| `read_collapser_preset`                        | `te-duplex` | `te-duplex`         | `te-duplex`   | `te-duplex`             |
| `read_collapser_wiggle_room`                   | `0`         | `0`                 | —             | —                       |
| `read_collapser_cluster_by_umi`                | ❌          | ❌                  | ✅            | ✅                      |
| **Subsampling**                                |             |                     |               |                         |
| `enable_subsample_bam`                         | ❌          | ✅                  | ❌            | ✅                      |
| `subsample_bam_coverage_target`                | `null`      | `null`              | `null`        | `null`                  |
| **Metrics**                                    |             |                     |               |                         |
| `enable_alignment_metrics_autosomes`           | ❌          | ❌                  | ❌            | ❌                      |
| `enable_alignment_metrics_targets`             | ✅          | ✅                  | ✅            | ✅                      |
| `enable_alignment_metrics_high_confidence`     | ❌          | ❌                  | ❌            | ❌                      |
| `enable_target_enrichment_metrics`             | ❌          | ✅                  | ❌            | ✅                      |
| `enable_hp_metrics`                            | ❌          | ✅                  | ❌            | ✅                      |
| `enable_gc_bias_metrics`                       | ❌          | ✅                  | ❌            | ✅                      |
| `enable_hs_metrics`                            | ❌          | ✅                  | ❌            | ✅                      |
| `enable_pre_dedup_metrics`                     | ✅          | ✅                  | ✅            | ✅                      |
| `enable_pre_subsample_metrics`                 | ❌          | ❌                  | ❌            | ❌                      |
| `alignment_metrics_min_base_quality_pre_dedup` | `30`        | `30`                | `30`          | `30`                    |
| `alignment_metrics_min_base_quality`           | `30`        | `30`                | `30`          | `30`                    |
| **Contamination**                              |             |                     |               |                         |
| `enable_contamination_estimation`              | ❌          | ❌                  | ❌            | ❌                      |
| **Variant Calling**                            |             |                     |               |                         |
| `enable_small_variant_caller`                  | ❌          | ❌                  | ❌            | ❌                      |
| `enable_small_variant_caller_retrain`          | ❌          | ❌                  | ❌            | ❌                      |
| `enable_copy_number_caller`                    | ❌          | ❌                  | ❌            | ❌                      |
| `enable_str_caller`                            | ❌          | ❌                  | ❌            | ❌                      |
| `enable_merge_vcf`                             | ❌          | ❌                  | ❌            | ❌                      |

## TE Simplex

| Parameter                                      | te_simplex   | te_simplex_benchmark | te_simplex_umi | te_simplex_umi_benchmark |
| ---------------------------------------------- | ------------ | -------------------- | -------------- | ------------------------ |
| **Mode & Alignment**                           |              |                      |                |                          |
| `variant_calling_mode`                         | `none`       | `none`               | `none`         | `none`                   |
| `aligner`                                      | `pb_fq2bam`  | `pb_fq2bam`          | `pb_fq2bam`    | `pb_fq2bam`              |
| `enable_pangenome_personalization`             | ❌           | ❌                   | ❌             | ❌                       |
| **Demultiplex**                                |              |                      |                |                          |
| `exclude_partial_reads`                        | ❌           | ❌                   | ❌             | ❌                       |
| `enable_demux_strand_detect`                   | ❌           | ❌                   | ❌             | ❌                       |
| **Deduplication**                              |              |                      |                |                          |
| `dedup_strategy`                               | `consensus`  | `consensus`          | `consensus`    | `consensus`              |
| `dedup_aligner`                                | —            | —                    | —              | —                        |
| `read_collapser_preset`                        | `te-simplex` | `te-simplex`         | `te-simplex`   | `te-simplex`             |
| `read_collapser_cluster_by_umi`                | ❌           | ❌                   | ✅             | ✅                       |
| **Subsampling**                                |              |                      |                |                          |
| `enable_subsample_bam`                         | ❌           | ✅                   | ❌             | ✅                       |
| `subsample_bam_coverage_target`                | `null`       | `null`               | `null`         | `null`                   |
| **Metrics**                                    |              |                      |                |                          |
| `enable_alignment_metrics_autosomes`           | ❌           | ❌                   | ❌             | ❌                       |
| `enable_alignment_metrics_targets`             | ✅           | ✅                   | ✅             | ✅                       |
| `enable_alignment_metrics_high_confidence`     | ❌           | ❌                   | ❌             | ❌                       |
| `enable_target_enrichment_metrics`             | ❌           | ✅                   | ❌             | ✅                       |
| `enable_hp_metrics`                            | ❌           | ✅                   | ❌             | ✅                       |
| `enable_gc_bias_metrics`                       | ❌           | ✅                   | ❌             | ✅                       |
| `enable_hs_metrics`                            | ❌           | ✅                   | ❌             | ✅                       |
| `enable_pre_dedup_metrics`                     | ✅           | ✅                   | ✅             | ✅                       |
| `enable_pre_subsample_metrics`                 | ❌           | ❌                   | ❌             | ❌                       |
| `alignment_metrics_min_base_quality_pre_dedup` | `0`          | `0`                  | `0`            | `0`                      |
| `alignment_metrics_min_base_quality`           | `30`         | `30`                 | `30`           | `30`                     |
| **Contamination**                              |              |                      |                |                          |
| `enable_contamination_estimation`              | ❌           | ❌                   | ❌             | ❌                       |
| **Variant Calling**                            |              |                      |                |                          |
| `enable_small_variant_caller`                  | ❌           | ❌                   | ❌             | ❌                       |
| `enable_small_variant_caller_retrain`          | ❌           | ❌                   | ❌             | ❌                       |
| `enable_copy_number_caller`                    | ❌           | ❌                   | ❌             | ❌                       |
| `enable_str_caller`                            | ❌           | ❌                   | ❌             | ❌                       |
| `enable_merge_vcf`                             | ❌           | ❌                   | ❌             | ❌                       |
