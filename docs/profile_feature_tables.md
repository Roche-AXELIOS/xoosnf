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
| `read_collapser_preset`                        | `wgs-duplex`        | `wgs-duplex`                  | `wgs-duplex`          | `wgs-duplex`                    | `wgs-duplex-cfdna` | `wgs-duplex-cfdna`         |
| `enable_read_collapser_dedup`                  | —                   | —                             | ✅                    | ✅                              | —                  | —                          |
| `read_collapser_cluster_by_umi`                | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| **Subsampling**                                |                     |                               |                       |                                 |                    |                            |
| `enable_subsample_bam`                         | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `subsample_bam_coverage_target`                | `null`              | `30`                          | `null`                | `null`                          | `null`             | `null`                     |
| **Metrics**                                    |                     |                               |                       |                                 |                    |                            |
| `enable_alignment_metrics_autosomes`           | ✅                  | ✅                            | ✅                    | ✅                              | ✅                 | ✅                         |
| `enable_alignment_metrics_targets`             | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_alignment_metrics_high_confidence`     | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `enable_hp_metrics`                            | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `enable_gc_bias_metrics`                       | ❌                  | ✅                            | ❌                    | ✅                              | ❌                 | ✅                         |
| `enable_hs_metrics`                            | ❌                  | ❌                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_pre_dedup_metrics`                     | ❌                  | ❌                            | ❌                    | ❌                              | ✅                 | ✅                         |
| `enable_pre_subsample_metrics`                 | ❌                  | ❌                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `alignment_metrics_min_base_quality_pre_dedup` | `30`                | `30`                          | —                     | —                               | —                  | —                          |
| `alignment_metrics_min_base_quality`           | `30`                | `30`                          | `30`                  | `30`                            | `30`               | `30`                       |
| **Contamination**                              |                     |                               |                       |                                 |                    |                            |
| `enable_contamination_estimation`              | ✅                  | ✅                            | ✅                    | ✅                              | ✅                 | ✅                         |
| **Variant Calling**                            |                     |                               |                       |                                 |                    |                            |
| `enable_small_variant_caller`                  | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `enable_small_variant_caller_retrain`          | ❌                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_copy_number_caller`                    | ✅                  | ✅                            | ✅                    | ✅                              | ❌                 | ❌                         |
| `enable_str_caller`                            | ✅                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_merge_vcf`                             | ✅                  | ✅                            | ❌                    | ❌                              | ❌                 | ❌                         |
| `enable_vcfeval`                               | —                   | ✅                            | —                     | —                               | —                  | —                          |
| **Variant Calling Partitioning**               |                     |                               |                       |                                 |                    |                            |
| `enable_variant_calling_partitioning`          | ✅                  | ✅                            | ✅                    | ✅                              | —                  | —                          |
| `cpu_variant_calling_partition_count`          | `20`                | `20`                          | `40`                  | `40`                            | —                  | —                          |
| `gpu_variant_calling_partition_count`          | `1`                 | `1`                           | `4`                   | `4`                             | —                  | —                          |

## WGS Simplex

| Parameter                                      | wgs_simplex   | wgs_simplex_benchmark |
| ---------------------------------------------- | ------------- | --------------------- |
| **Mode & Alignment**                           |               |                       |
| `variant_calling_mode`                         | `none`        | `none`                |
| `aligner`                                      | `pb_fq2bam`   | `pb_fq2bam`           |
| `enable_pangenome_personalization`             | ❌            | ❌                    |
| **Demultiplex**                                |               |                       |
| `exclude_partial_reads`                        | ✅            | ✅                    |
| `enable_demux_strand_detect`                   | ❌            | ❌                    |
| **Deduplication**                              |               |                       |
| `dedup_strategy`                               | `markdup`     | `markdup`             |
| `read_collapser_preset`                        | `wgs-simplex` | `wgs-simplex`         |
| `enable_read_collapser_dedup`                  | —             | —                     |
| `read_collapser_cluster_by_umi`                | ❌            | ❌                    |
| **Subsampling**                                |               |                       |
| `enable_subsample_bam`                         | ❌            | ✅                    |
| `subsample_bam_coverage_target`                | `null`        | `null`                |
| **Metrics**                                    |               |                       |
| `enable_alignment_metrics_autosomes`           | ✅            | ✅                    |
| `enable_alignment_metrics_targets`             | ❌            | ❌                    |
| `enable_alignment_metrics_high_confidence`     | ✅            | ✅                    |
| `enable_hp_metrics`                            | ❌            | ✅                    |
| `enable_gc_bias_metrics`                       | ❌            | ✅                    |
| `enable_hs_metrics`                            | ❌            | ❌                    |
| `enable_pre_dedup_metrics`                     | ❌            | ❌                    |
| `enable_pre_subsample_metrics`                 | ❌            | ❌                    |
| `alignment_metrics_min_base_quality_pre_dedup` | `22`          | `22`                  |
| `alignment_metrics_min_base_quality`           | `22`          | `22`                  |
| **Contamination**                              |               |                       |
| `enable_contamination_estimation`              | ❌            | ❌                    |
| **Variant Calling**                            |               |                       |
| `enable_small_variant_caller`                  | ❌            | ❌                    |
| `enable_small_variant_caller_retrain`          | ❌            | ❌                    |
| `enable_copy_number_caller`                    | ❌            | ❌                    |
| `enable_str_caller`                            | ❌            | ❌                    |
| `enable_merge_vcf`                             | ❌            | ❌                    |

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
| `read_collapser_preset`                        | `te-simplex` | `te-simplex`         | `te-simplex`   | `te-simplex`             |
| `enable_read_collapser_dedup`                  | —            | —                    | —              | —                        |
| `read_collapser_cluster_by_umi`                | ❌           | ❌                   | ✅             | ✅                       |
| **Subsampling**                                |              |                      |                |                          |
| `enable_subsample_bam`                         | ❌           | ✅                   | ❌             | ✅                       |
| `subsample_bam_coverage_target`                | `null`       | `null`               | `null`         | `null`                   |
| **Metrics**                                    |              |                      |                |                          |
| `enable_alignment_metrics_autosomes`           | ❌           | ❌                   | ❌             | ❌                       |
| `enable_alignment_metrics_targets`             | ✅           | ✅                   | ✅             | ✅                       |
| `enable_alignment_metrics_high_confidence`     | ❌           | ❌                   | ❌             | ❌                       |
| `enable_hp_metrics`                            | ❌           | ✅                   | ❌             | ✅                       |
| `enable_gc_bias_metrics`                       | ❌           | ✅                   | ❌             | ✅                       |
| `enable_hs_metrics`                            | ❌           | ❌                   | ❌             | ❌                       |
| `enable_pre_dedup_metrics`                     | ✅           | ✅                   | ✅             | ✅                       |
| `enable_pre_subsample_metrics`                 | ❌           | ❌                   | ❌             | ❌                       |
| `alignment_metrics_min_base_quality_pre_dedup` | `29`         | `29`                 | `29`           | `29`                     |
| `alignment_metrics_min_base_quality`           | `29`         | `29`                 | `29`           | `29`                     |
| **Contamination**                              |              |                      |                |                          |
| `enable_contamination_estimation`              | ❌           | ❌                   | ❌             | ❌                       |
| **Variant Calling**                            |              |                      |                |                          |
| `enable_small_variant_caller`                  | ❌           | ❌                   | ❌             | ❌                       |
| `enable_small_variant_caller_retrain`          | ❌           | ❌                   | ❌             | ❌                       |
| `enable_copy_number_caller`                    | ❌           | ❌                   | ❌             | ❌                       |
| `enable_str_caller`                            | ❌           | ❌                   | ❌             | ❌                       |
| `enable_merge_vcf`                             | ❌           | ❌                   | ❌             | ❌                       |
