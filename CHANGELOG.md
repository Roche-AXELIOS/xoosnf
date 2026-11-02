# Roche-AXELIOS/xoosnf: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.2.0]

### `Added`

- Include `filter_variants` small-variant caller QC metrics from germline (pretrained/retrained), tumor-normal, tumor-only, and mitochondrial workflows in the merged metrics report and optional metrics summary.
- Added `enable_target_enrichment_metrics` to pass `--enable-target-enrichment-metrics` to target-region alignment metrics and include its scalar coverage-uniformity and on-target read metrics in the report. It is enabled in the `te_simplex_benchmark`, `te_simplex_umi_benchmark`, `te_duplex_benchmark`, and `te_duplex_umi_benchmark` profiles.
- Variant phasing for somatic tumor-normal and tumor-only workflows. The filtered small-variant VCF is phased with the corresponding Mutect2 bamout BAM before downstream evaluation, copy-number calling, and germline filtering; partitioned Mutect2 BAMs are merged before phasing.
- Added `te_duplex`, `te_duplex_benchmark`, `te_duplex_umi`, and `te_duplex_umi_benchmark` profiles for Target Enrichment duplex workflows. The UMI profiles enable UMI clustering; the benchmark profiles enable HP, GC-bias, and Hybrid Selection metrics.
- Microsatellite instability (MSI) calling for tumor/normal pairs (beta / experimental), enabled via `enable_msi_caller` (default `false`, `variant_calling_mode = somatic_tn` only). Runs `msi_caller tumor-normal` on each verified tumor/normal pair and publishes the results directory under `<run>/<sample>/msi_caller/`. An optional `msi_sites_thresholds_tsv` override supplies per-site thresholds; when omitted, the bundled threshold-bearing default is used. The caller uses fixed defaults for its infrequently changed thresholds.
- Tumor-only microsatellite instability (MSI) calling (beta / experimental) is now available with `enable_msi_caller` in `variant_calling_mode = somatic_to`. It runs `msi_caller tumor-only` for each tumor BAM and uses the bundled panel-of-normals background thresholds by default; `msi_background_threshold_json` can override that resource.
- New options `alignment_metrics_additional_regions` and `alignment_metrics_additional_region_beds` compute alignment metrics over additional stratification regions. Region labels are paired with BED paths, merged into a single labeled BED, evaluated in one label-stratified `alignment_metrics all` run, and reported per region.
- Somatic tumor-only WGS variant calling (`variant_calling_mode = somatic_to`): Mutect2 tumor-only calling with small-variant filtering, tumor-only copy number calling, and optional germline filtering (`enable_germline_filtering`). Small-variant filtering (`filter_variants tumor-only-wgs`) uses the tumor-only-wgs model shipped in the small variant caller container by default; supply `svc_model_somatic_tumor_only` to override. An optional panel of normals VCF (`svc_somatic_to_panel_of_normals_vcf` / `svc_somatic_to_panel_of_normals_tbi`) can be provided for PON-based features. Germline filtering consumes the small-variant VCF and the copy-number caller's purity/ploidy (SEG) output. Small-variant and copy-number VCFs are merged into a single VCF by default (`enable_merge_vcf = true`); when germline filtering is enabled, its filtered VCF is merged in place of the raw small-variant VCF. Adds the `somatic_wgs_to_duplex` and `somatic_wgs_to_duplex_benchmark` profiles, the `CALL_SOMATIC_VARIANTS_TO` subworkflow, and the `svc_model_somatic_tumor_only`, `svc_somatic_to_{truth_vcf,eval_bed}_dir`, `cnc_somatic_to_blocklist_bed`, `cnc_somatic_to_snp_blocklist_bed`, and `svc_somatic_to_panel_of_normals_{vcf,tbi}` resource parameters.
- Optional indexed panel-of-normals VCF support for somatic small-variant filtering through the reusable `svc_panel_of_normals_vcf` and `svc_panel_of_normals_tbi` parameters.
- Mitochondrial BAM post-processing now consumes the Giraffe graph components emitted by alignment, allowing germline chrM calling to run with either static or personalized pangenome graphs.
- New `dedup_aligner` parameter to run a different aligner for the post-read-collapser re-alignment in the consensus deduplication path. When unset (default), the re-alignment reuses `aligner`, preserving existing behaviour. Giraffe and non-Giraffe aligners cannot be mixed across the read-collapser boundary because they target different reference spaces.
- New `germline_wgs_simplex_benchmark` profile for benchmarking germline WGS simplex variant calling (adds HP/GC-bias metrics, SVC retraining, and vcfeval).
- Somalier-based Tumor/Normal pair verification (`enable_somatic_tn_pair_verification`) to confirm paired samples originate from the same individual before somatic variant calling, requires `variant_calling_mode` set to `somatic_tn`.
- `somatic_tn_pair_mismatch_action` controls the response to a detected mismatch: `warn` (log and continue), `abort` (stop the pipeline), or `skip` (exclude the mismatched pair from downstream variant calling).
- Additional somalier options: `somalier_sites_vcf` (known polymorphic sites for genotyping) and `somatic_tn_pair_concordance_cutoff` (mismatch threshold, default `0.9`).
- Support for pre-aligned CRAM input via the `cram` `file_type`. Requires `reference_fasta` and `reference_fai` for the CRAM-to-BAM conversion.
- Germline mitochondrial variant calling, enabled via the new `enable_mitochondrial_variant_caller` parameter (default `false`, germline mode only). Calls variants on chrM with GATK Mutect2 in mitochondria-mode followed by SVC germline-mitochondria threshold filtering, detecting heteroplasmy down to AF >= 1%.
- New option `enable_downsampled_bwa_metrics` to run downsampled BWA metrics subworkflow: raw demuxed FASTQs are subsampled, re-aligned with the configured aligner (`bwa` or `pb_fq2bam`), and full alignment metrics are computed on the re-aligned BAM. Controlled by `downsampled_bwa_metrics_aligner`, `downsampled_bwa_metrics_read_count`, and `downsampled_bwa_metrics_seed`.
- The germline copy number caller now accepts a blocklist BED via `cnv_blocklist_bed`, passed to `copy_number_caller germline-wgs --blocklist-bed`. When unset it defaults to the container-packaged `hg38_decoy_intervals.bed`.

### `Changed`

- Enabled tumor-fraction estimation by default for the `cfdna_wgs_duplex` and `cfdna_wgs_duplex_benchmark` profiles. Samples that provide a `reporter_vcf` now produce tumor-fraction output without requiring `--enable_tfe true`.
- Hybrid Selection metrics now passes `alignment_metrics_min_mapping_quality` to Picard `CollectHsMetrics` as `--MINIMUM_MAPPING_QUALITY` and publishes its summary, per-target, and per-base coverage outputs.
- Demux 1.1 is no longer supported: xoosnf uses the newer Demux CLI options `--reads`, `--output-dir`, `--strand-index`, and `--min-concordant-duplex-bases-for-early-stop`, and requires the canonical `metric_name` and `value` metrics headers.
- `read_collapser_cluster_by_umi` now passes `--cluster-by-umi=true|false` when explicitly configured and otherwise leaves the tool default unchanged. Added the `read_collapser_wiggle_room` parameter for consensus clustering.
- Germline copy-number calling now maps pipeline `run_type` values to the copy-number caller's `--library-type` option (`sbx-d`, `sbx-fast`, or `sbx-s`) and fails clearly for unsupported run types.
- Copy number caller invocations now preserve the tool's `metrics/`, `tracks/`, and `supplemental/` hierarchy and let report/MultiQC consumers discover the nested metrics file.
- Cross-sample contamination estimation is enabled by default for `germline_wgs_simplex` and `germline_wgs_simplex_benchmark`.
- The `germline_wgs_simplex` and `germline_wgs_simplex_benchmark` profiles now use `pb_giraffe` with pangenome personalization by default. These profiles require the unfiltered Giraffe index and linearized-reference resources.
- The `germline_wgs_simplex_benchmark` profile now subsamples BAMs to 40× before benchmark processing.
- The germline simplex profiles now disable STR calling because `str_caller` does not support simplex data.
- SVC Germline WGS `filter_variants` now uses `small_variant_caller_germline_model_type` (`large` by default; `small` opt-in) for bundled models. Model retraining invokes `train_model germline-wgs --model-type small`.
- `XOOS_MSI_CALLER_TUMOR_NORMAL`/`XOOS_MSI_CALLER_TUMOR_ONLY` no longer force `--min-depth`/`--min-sites`; `msi_caller` applies its own defaults, so `--min-sites` effectively becomes `50` instead of the previous `10`. The caller uses embedded catalogs and thresholds by default, while `msi_sites_thresholds_tsv` and `msi_background_threshold_json` remain available as explicit resource overrides.
- Somatic tumor-normal small variant calling now defaults to automatic sample-type detection when selecting the bundled BWA model. Samplesheets may still explicitly select `ffpe` or `cell-line`.
- Alignment metric **keys** now include the aligner name and use a hierarchical `<stage>/<aligner>/<region>` layout (e.g. `alignment_metrics/dedup/bwa/autosomes/...`, `alignment_metrics/downsampled_metrics/pb_fq2bam/autosomes/...`). Only the in-report metric keys change; output JSON filenames and published directory paths are unchanged.
- Split the subsample FASTQ module into `seqkit_head` (deterministic) and `seqkit_random` (random) variants.
- Enabled high-confidence alignment metrics (`enable_alignment_metrics_high_confidence = true`) in the `cfdna_wgs_duplex_benchmark` profile, aligning it with the other WGS benchmark profiles.
- Renamed the experimental simplex germline profiles to `germline_wgs_simplex` and `germline_wgs_simplex_benchmark`.
- Set `small_variant_caller_sequencing_protocol = simplex` for the germline simplex profiles. These profiles retain partial reads (`exclude_partial_reads = false`) and use XOOS read_collapser MARKDUP.
- Updated germline small-variant-caller invocations to the unified `germline-wgs` workflow. The pipeline passes the explicit metadata-derived library type (`sbx-d`, `sbx-fast`, or `sbx-s`) to feature extraction, filtering, and model training.
- Enabled Hybrid Selection (HS) metrics (`enable_hs_metrics = true`) in the `te_simplex_benchmark` and `te_simplex_umi_benchmark` profiles. When `capture_bed` is not provided, HS metrics use `target_bed` as a fallback and warn that on-bait metrics may be inaccurate.
- Removed the explicit `enable_subsample_bam` setting from the `wgs_simplex_benchmark` profile; it now inherits the global default (`false`).
- The copy number caller now uses `--output-prefix` to name outputs directly instead of renaming files after the run.
- Updated the `read_collapser` invocations for its cross-module CLI/output alignment: `--bam-input`→`--bam` and `--bed-input`→`--bed` (markdup/consensus), `--prefix`→`--output-prefix` with a `_` separator and the default rescue BAM renamed to `<prefix>_rescued_secondary.bam` (rescue-secondary), and all metric TSVs now read from the tool's `metrics/` subdirectory. The markdup and consensus `metrics_dir` outputs now preserve that `metrics/` nesting (previously flattened), matching the rescue-secondary layout and the tool's output hierarchy; metric parsing and MultiQC consumers were updated accordingly. No change to published output paths.
- `max_parabricks_alignment_read_length` sets `--max-read-length` and `--max-read-length-on-gpu` for the Parabricks fq2bam aligner, and its default is `1024`; Parabricks Giraffe uses its own default read length.
- Extended the input samplesheet `file_type` enum to accept `cram`.
- New `VERIFY_CRAM` module that checks each pre-aligned CRAM carries `M5` (reference-checksum) tags on every `@SQ` header record before conversion.
- The Parabricks `minimap2` module now feeds input FASTQs through a named pipe (`--in-fq input.fastq.gz`) instead of expanding every file onto the command line, unifying it with the `fq2bam` and `giraffe` aligners and removing the `ARG_MAX` command-line-length ceiling for samples with many FASTQ shards.
- Alignment metrics now uses giraffe linearized index when giraffe aligners (`pb_giraffe` or `vg_giraffe`) are used.
- Update the STR caller invocation to pass the renamed `--bam` flag (was `--bam-input`).
- Duplicate marking now uses XOOS read_collapser when `dedup_strategy = markdup`. Parabricks fq2bam, Giraffe, and minimap2 always disable their integrated duplicate marking.
- The pipeline now requires that the pipeline's `dedup_strategy` to be compatible with the `bam_type` for any externally deduped BAM/CRAM input.
- The pipeline now correctly skips externally deduplicated BAMs when calculating pre-dedup metrics.
- Copy Number Caller WGS invocations now omit `--mappability-bigwig` when no override is configured, allowing the caller to select its bundled hg38 track from the detected library type.

### `Dependencies`

- Update Parabricks container to `4.7.1-1`.
- Add `somalier` `0.2.19` for Tumor/Normal pair verification.
- Update seqkit to `2.13.0`.

### `Fixed`

- Normalized somatic tumor-only vcfeval metric names to `vcf_eval/small_variant_pretrained/*`, matching germline output and the configured metrics summary patterns.
- Somatic CNC and MSI VCF replacement now preserves each caller's output directory (`copy_number_caller/` or `msi_caller/`) instead of publishing files under a generic `results/` directory.
- The `te_simplex`, `te_simplex_benchmark`, `te_simplex_umi`, and `te_simplex_umi_benchmark` profiles no longer report invalid Q93 pre-dedup accuracy. Simplex pre-dedup bases are capped below Q29, so the previous `alignment_metrics_min_base_quality_pre_dedup = 29` filtered out every base and alignment metrics silently reported zero errors / Phred 93. The alignment metrics thresholds are now `alignment_metrics_min_base_quality_pre_dedup = 0` and `alignment_metrics_min_base_quality = 30`, so pre-dedup metrics are computed over the available simplex bases.
- Accuracy parsing (`base_level_accuracy_summary.tsv`, `hp_errors.tsv`) now reports `NA` instead of the meaningless Phred 93 that `alignment_metrics` emits when no bases are evaluated (undefined error rate). A genuinely error-free sample (defined 0% error rate) still keeps its real Phred value.
- Parabricks Giraffe now emits up to three multimapping alignments when copy number calling is enabled, allowing the secondary-alignment rescue step to supply eligible alignments to the copy number caller.
- Set `alignment_metrics_min_base_quality_pre_dedup` and `alignment_metrics_min_base_quality` to `0` for all WGS simplex profiles so alignment metrics do not filter preserved simplex base-caller quality scores.
- Samplesheet schema-validation errors now show the clean nf-schema message as the headline error instead of a raw Java/Groovy stack trace. Validation is now performed eagerly on the main thread in `ROUTE_INPUTS` before any channel operator closures run.
- Samplesheet schema-validation errors are now reported once, not duplicated across multiple actor threads.
- The Parabricks `fq2bam`, `giraffe`, and `minimap2` modules no longer hang when `pbrun` exits early (e.g. on an argument-parsing error). The background FIFO writer feeding `input.fastq.gz` is now tracked and killed via an `EXIT` trap, so no orphaned process is left blocked on the named pipe keeping the container alive until timeout.
- The Parabricks `minimap2` module now emits a single, sample-level read group instead of one read group per input FASTQ. Feeding `pbrun minimap2` a single concatenated stream (via the same named pipe used by `fq2bam`/`giraffe`) makes it see one input, and `--read-group-sm` is now set so the `SM` tag is populated.
- All aligner modules (`bwa/mem`, `vg/giraffe`, and Parabricks `fq2bam`/`giraffe`/`minimap2`) now fail fast with a clear error when no input FASTQ files are found, instead of hanging (the piped aligners previously ran `zcat` on empty stdin) or running the aligner on an empty stream.
- Read groups across all aligners now set the `PM` (platform model) tag to `AXELIOS_1` and use a structured library name `LB=unknown.<run_type>.<sample_id>`. Because Parabricks has no `--read-group-pm` flag, the `fq2bam` and `minimap2` modules now pass the full read group as an inline `@RG` string on the input argument instead of discrete `--read-group-*` flags.
- HaplotypeCaller and Mutect2 index files (`.tbi`/`.bai`) are now tracked via their producing process output channels instead of being reconstructed from string-built paths. The old paths bypassed Nextflow's DAG, so FilePorter staged them as untracked "foreign files" and downstream tasks intermittently failed with `file does not exist` even when the file was present. Adds an explicit `SAMTOOLS_INDEX` step for the GATK realigned BAMs to ensure the index is DAG-tracked (the GATK-generated `.bai` is not declared as a module output and is therefore invisible to Nextflow's dependency graph).
- Germline SVC retraining no longer aborts the pipeline for samples without matching benchmark resources; such samples are skipped for retraining with a warning.
- Added regression coverage for `XOOS_PARSE_STR_METRICS` locating `metrics.tsv` when staged nested under a `metrics/` directory.
- The pipeline completion handler no longer aborts. `workflow.onComplete` read `workflow.success`, but inside that nested subworkflow closure the implicit `workflow` binding resolves to null, throwing "Cannot get property 'success' on null object" and skipping the completion summary (and failing the driver job). The run's success flag is now captured in the enclosing scope and read from there.
- Germline copy-number results now use the sample name from the BAM `SM` tag consistently with HaplotypeCaller and the small-variant caller, while reheadering CNC VCFs back to the samplesheet sample name.

## [1.1.0] - 2026-07-12

Initial release of Roche-AXELIOS/xoosnf, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Changed`

- Update the STR caller invocation to pass `--aligner giraffe` (for `pb_giraffe`/`vg_giraffe`) instead of the removed boolean `--giraffe-input-bam` flag. For other aligners the flag is omitted so the caller infers the aligner from the BAM `@PG` header (the caller's `--aligner` default). The `--aligner` option accepts only `unknown`, `bwa`, or `giraffe`.

### `Fixed`

### `Dependencies`

### `Deprecated`

[1.1.0]: https://github.com/Roche-AXELIOS/xoosnf/releases/tag/1.1.0
[1.2.0]: https://github.com/Roche-AXELIOS/xoosnf/compare/1.1.0...1.2.0
