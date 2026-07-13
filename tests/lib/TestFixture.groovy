import groovy.json.JsonOutput

/*
 * TestFixture
 *
 * Builds the on-disk fixture tree required by stub (`-stub`) tests, and emits a
 * `resources.json` params file wiring the reference params to the files it created.
 *
 * Stub tests never read file *content* (processes only emit stub outputs), but
 * the pipeline logic still runs: params with `checkIfExists: true` must point at
 * existing files, directory params are globbed at runtime so their internal
 * filenames must match, and samplesheets/JSON are parsed for real.
 *
 * Isolation: build under a *per-test* root — `${launchDir}/test-data` — not a shared
 * location. nf-test hash-isolates `launchDir` per test and recreates it each run, so
 * fixtures never leak between tests. `launchDir` is the parent of `outputDir`, so the
 * tree sits outside `params.outdir` and is never swept into output snapshots.
 *
 * Resources config: each builder records the reference params it implies; `writeResources()`
 * serializes them to `${root}/resources.json`. The test loads them — the nf-test analogue
 * of the real pipeline's `-c resources.config`:
 *
 *   setup {
 *       TestFixture.under("${launchDir}/test-data")
 *           .bwaRef().variantCallingResources().fastqs().samplesheets()
 *           .writeResources()
 *   }
 *   when {
 *       params {
 *           load("${launchDir}/test-data/resources.json")   // references, from the fixture
 *           outdir  = "$outputDir"
 *           input   = "${launchDir}/test-data/inputs/input_samplesheet.csv"
 *           aligner = 'bwa'                                  // test-specific flags only
 *       }
 *   }
 *
 * Each test declares exactly the fixtures it needs via which builders it calls; the
 * generated `resources.json` then contains only those references. Every builder is
 * idempotent (placeholder files are only created if missing).
 * `under(root)` starts a builder; `root()` returns the root path; `writeResources()`
 * writes the params file and returns the root; `all()` materializes everything.
 */
class TestFixture {

    /** Samplesheet rows shared across pipeline tests. */
    static final String SAMPLESHEET = """\
sample_name,sample_sid
HG001,CTATTGAAGGCT
HG002,CCTGCGGTTGAA
"""

    static final String SAMPLESHEET_SOMATIC = """\
sample_name,sample_sid,TN,group
HG001,CTATTGAAGGCT,N,pair1
HG002,CCTGCGGTTGAA,T,pair1
"""

    /** Sample identifiers used across fixtures. */
    private static final List<String> SAMPLES = ['HG001', 'HG002']

    private final String root

    /** Reference params implied by the builders called, written by writeResources(). */
    private final Map<String, String> resourceParams = [:]

    private TestFixture(String root) {
        this.root = root
        new File(root).mkdirs()
    }

    /** Entry point: start a builder rooted at `root`. */
    static TestFixture under(String root) {
        return new TestFixture(root)
    }

    /** The fixture root path (terminates a fluent chain). */
    String root() {
        return root
    }

    /**
     * Write the recorded reference params to `${root}/resources.json` and return the root.
     * Tests load this via `params { load("${root}/resources.json") }` — the nf-test
     * analogue of the real pipeline's `-c resources.config`.
     */
    String writeResources() {
        def asStrings = resourceParams.collectEntries { k, v -> [(k): v.toString()] }
        new File("${root}/resources.json").text = JsonOutput.prettyPrint(JsonOutput.toJson(asStrings)) + '\n'
        return root
    }

    /**
     * Materialize the complete fixture tree, write resources.json, and return the root.
     *
     * Convenience for a test that exercises every code path; most tests should call
     * only the fine-grained builders they need so resources.json stays minimal.
     */
    String all() {
        return bwaRef()
            .giraffeRef()
            .referenceSdf()
            .variantCallingResources()
            .vcfevalResources()
            .fastqs()
            .bams()
            .samplesheets()
            .writeResources()
    }

    // ---- Fine-grained fixtures -------------------------------------------

    /** BWA index directory — globbed by the pipeline as a directory. */
    TestFixture bwaRef() {
        def dir = "${root}/resources/bwa_ref"
        [
            'test.fna', 'test.fna.fai', 'test.fna.amb', 'test.fna.ann',
            'test.fna.bwt', 'test.fna.pac', 'test.fna.sa',
            'test.dict', 'test.autosomes.bed', 'test.bloom', 'test.bw',
            'test.chrX_PAR.bed', 'test.chrY_PAR.bed', 'cnv_seed.bed',
        ].each { emptyFile("${dir}/${it}") }
        resourceParams.putAll([
            reference_fasta          : "${dir}/test.fna",
            reference_fai            : "${dir}/test.fna.fai",
            reference_dict           : "${dir}/test.dict",
            autosomes_bed            : "${dir}/test.autosomes.bed",
            chr_x_par_bed            : "${dir}/test.chrX_PAR.bed",
            chr_y_par_bed            : "${dir}/test.chrY_PAR.bed",
            mappability_bw           : "${dir}/test.bw",
            demux_strand_detect_index: "${dir}/test.bloom",
            bwa_index_dir            : dir,
        ])
        return this
    }

    /** Giraffe index directory — globbed by the pipeline as a directory. */
    TestFixture giraffeRef() {
        def dir = "${root}/resources/giraffe_ref"
        [
            'test.gbz', 'test.min', 'test.dist', 'test.gbwt', 'test.hapl',
            'test.pg', 'test.ri', 'test.vg', 'test.fa', 'test.fa.fai', 'test.dict',
            'test.shortread.withzip.min', 'test.shortread.zipcodes',
        ].each { emptyFile("${dir}/${it}") }
        // Pipeline derives {prefix}.zipcodes from the .gbz basename.
        emptyFile("${dir}/test.zipcodes")
        resourceParams.putAll([
            giraffe_index_dir           : dir,
            giraffe_reference_linearized: "${dir}/test.fa",
        ])
        return this
    }

    /** Pre-built RTG SDF directory for the reference_sdf parameter. */
    TestFixture referenceSdf() {
        def dir = "${root}/resources/reference_sdf"
        [
            'done', 'format.log', 'mainIndex', 'nameIndex0', 'namedata0',
            'namepointer0', 'progress', 'seqdata0', 'seqpointer0',
            'sequenceIndex0', 'suffixIndex0', 'suffixdata0', 'suffixpointer0',
            'summary.txt',
        ].each { emptyFile("${dir}/${it}") }
        resourceParams.putAll([reference_sdf: dir])
        return this
    }

    /** SVC models, population-AF VCF, STR catalog. */
    TestFixture variantCallingResources() {
        def resources = "${root}/resources"
        [
            'svc_model_snv.bin', 'svc_model_indel.bin', 'svc_model_somatic_tn.bin',
        ].each { emptyFile("${resources}/variantCaller/${it}") }

        emptyFile("${resources}/gatk-best-practices/hg38/af-only-gnomad.hg38.vcf.gz")
        emptyFile("${resources}/gatk-best-practices/hg38/af-only-gnomad.hg38.vcf.gz.tbi")

        // Parsed as JSON by the pipeline — must be valid, empty list is fine.
        textFile("${resources}/STR/hg38/full_repeat_catalog.json", '[]\n')
        resourceParams.putAll([
            svc_model_germline_snv        : "${resources}/variantCaller/svc_model_snv.bin",
            svc_model_germline_indel      : "${resources}/variantCaller/svc_model_indel.bin",
            svc_model_somatic_tumor_normal: "${resources}/variantCaller/svc_model_somatic_tn.bin",
            pop_af_vcf                    : "${resources}/gatk-best-practices/hg38/af-only-gnomad.hg38.vcf.gz",
            pop_af_tbi                    : "${resources}/gatk-best-practices/hg38/af-only-gnomad.hg38.vcf.gz.tbi",
            repeat_catalog                : "${resources}/STR/hg38/full_repeat_catalog.json",
        ])
        return this
    }

    /** vcfeval truth VCFs, eval BEDs, and SVC train/test BEDs (sample-named). */
    TestFixture vcfevalResources() {
        def resources = "${root}/resources"
        SAMPLES.each { sample ->
            emptyFile("${resources}/truth_vcf/${sample}_truth.vcf.gz")
            emptyFile("${resources}/truth_vcf/${sample}_truth.vcf.gz.tbi")
            emptyFile("${resources}/eval_bed/${sample}_eval.bed")
            emptyFile("${resources}/train_bed/${sample}_train.bed")
            emptyFile("${resources}/test_bed/${sample}_test.bed")
        }
        resourceParams.putAll([
            germline_small_variant_truth_vcf_dir: "${resources}/truth_vcf",
            germline_small_variant_eval_bed_dir : "${resources}/eval_bed",
        ])
        return this
    }

    /**
     * FASTQs: duplex run_dir, per-sample directories for the align subworkflow,
     * and simplex basecall placeholders.
     */
    TestFixture fastqs() {
        def inputs = "${root}/inputs"
        SAMPLES.each { sample ->
            emptyFile("${inputs}/raw_fastq/duplex/raw_fastq_${sample}.fastq.gz")
            emptyFile("${inputs}/raw_fastq/simplex/raw_fastq_${sample}.fastq.gz")
            // Per-sample directories consumed by align subworkflow tests.
            emptyFile("${inputs}/fastq/${sample}/raw_fastq_${sample}.fastq.gz")
        }
        return this
    }

    /** Empty BAM placeholders + indexes. */
    TestFixture bams() {
        def bamDir = "${root}/inputs/bam"
        SAMPLES.each { sample ->
            emptyFile("${bamDir}/${sample}.bam")
            emptyFile("${bamDir}/${sample}.bam.bai")
        }
        return this
    }

    /** Input samplesheets + per-run sample sheets (parsed by the pipeline). */
    TestFixture samplesheets() {
        def inputs = "${root}/inputs"
        def fastqDir = "${inputs}/raw_fastq/duplex"
        def simplexDir = "${inputs}/raw_fastq/simplex"
        def bamDir = "${inputs}/bam"
        def samplesheet = "${inputs}/samplesheet_stub.csv"
        def samplesheetSomatic = "${inputs}/samplesheet_somatic_tn.csv"

        textFile(samplesheet, SAMPLESHEET)
        textFile(samplesheetSomatic, SAMPLESHEET_SOMATIC)

        textFile("${inputs}/input_samplesheet.csv", """\
run_name,run_type,file_type,run_dir,samplesheet
test_run,SBX-D,demultiplexed_fastq,${fastqDir},${samplesheet}
""")

        textFile("${inputs}/input_samplesheet_bam.csv", """\
run_name,run_type,file_type,run_dir,samplesheet
test_run,SBX-D,bam,${bamDir},${samplesheet}
""")

        textFile("${inputs}/input_samplesheet_somatic_tn.csv", """\
run_name,run_type,file_type,run_dir,samplesheet
test_run,SBX-D,bam,${bamDir},${samplesheetSomatic}
""")

        textFile("${inputs}/input_samplesheet_simplex.csv", """\
run_name,run_type,file_type,run_dir,samplesheet
test_run,YSU,basecall_fastq,${simplexDir},${samplesheet}
""")
        return this
    }

    // ---- Helpers ---------------------------------------------------------

    /** Empty-file helper: creates parent dirs and an empty file if missing. */
    private static void emptyFile(String path) {
        def f = new File(path)
        f.parentFile?.mkdirs()
        if (!f.exists()) {
            f.text = ''
        }
    }

    /** Text-file helper: creates parent dirs and writes content if missing. */
    private static void textFile(String path, String content) {
        def f = new File(path)
        f.parentFile?.mkdirs()
        if (!f.exists()) {
            f.text = content
        }
    }
}
