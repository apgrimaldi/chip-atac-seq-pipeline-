nextflow.enable.dsl=2

include { FASTQC }                   from '../modules/local/fastqc.nf'
include { TRIMGALORE }               from '../modules/local/trimgalore.nf'
include { BOWTIE2_BUILD }            from '../modules/local/bowtie2_build.nf'
include { BOWTIE2 }                  from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }            from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }           from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES }    from '../modules/local/picard_markduplicates.nf'
include { FILTERING }                from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }        from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }         from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }        from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }         from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }      from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }                from '../modules/local/calc_frip.nf'
include { DEEPTOOLS }                from '../modules/local/deeptools.nf'
include { DIFFBIND }                 from '../modules/local/diffbind.nf'
include { MULTIQC }                  from '../modules/local/multiqc.nf'
include { SAMTOOLS_INDEX }           from '../modules/local/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/local/samtools_index.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input 

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    def reference_file   = null
    def gtf_file         = null
    def bowtie2_index    = null
    def blacklist_path   = null
    def m_genome         = params.macs_gsize 

    if (params.genomes && params.genomes.containsKey(params.genome)) {
        def gdata      = params.genomes[params.genome]
        reference_file = params.reference_file ?: gdata.fasta
        gtf_file       = params.gtf_file       ?: gdata.gtf
        bowtie2_index  = params.bowtie2_index  ?: gdata.bowtie2
        blacklist_path = params.blacklist      ?: (gdata.containsKey('blacklist') ? gdata.blacklist : null)
        if (!m_genome) m_genome = gdata.macs_gsize
    } else {
        reference_file = params.reference_file
        gtf_file       = params.gtf_file
        bowtie2_index  = params.bowtie2_index
        blacklist_path = params.blacklist
    }

    if (!m_genome || m_genome == 'custom') {
        m_genome = 'hs'
    }

    ch_index_internal = Channel.empty() 
    if (bowtie2_index) {
        ch_index_internal = Channel.fromPath("${bowtie2_index}/*.bt2*").collect()
    } else if (reference_file) {
        BOWTIE2_BUILD ( file(reference_file) )
        ch_index_internal = BOWTIE2_BUILD.out.index.collect()
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    FASTQC ( ch_input )
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions, TRIMGALORE.out.versions)

    BOWTIE2 ( TRIMGALORE.out.reads, ch_index_internal )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    def ch_ref_path = reference_file ? file(reference_file) : []
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, ch_ref_path, [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    
    if (blacklist_path) {
        FILTERING ( SAMTOOLS_INDEX.out.bam_bai, file(blacklist_path) )
        SAMTOOLS_INDEX_FINAL ( FILTERING.out.bam )
        ch_final_bams = SAMTOOLS_INDEX_FINAL.out.bam_bai
    } else {
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
    }

    SAMTOOLS_STATS ( ch_final_bams.map { meta, bam, bai -> [ meta, bam ] } )
    DEEPTOOLS ( ch_final_bams )

    ch_bams_branched = ch_final_bams
        .branch { meta, bam, bai ->
            control: meta.is_control == true
            ip:      true
        }

    if (params.protocol == 'atac') {
        ch_macs_input = ch_bams_branched.ip.map { meta, bam, bai -> [ meta, bam, [] ] }
    } else {
        ch_macs_input = ch_bams_branched.ip
            .map { meta, bam, bai -> [ meta.control, meta, bam ] }
            .combine(ch_bams_branched.control.map { meta, bam, bai -> [ meta.id, bam ] }, by: 0)
            .map { control_id, meta, ip_bam, control_bam -> [ meta, ip_bam, control_bam ] }
    }

    ch_narrow_peaks = Channel.empty()
    ch_broad_peaks  = Channel.empty()
    ch_macs_logs_mqc = Channel.empty()
    ch_narrow_counts_mqc = Channel.empty()
    ch_broad_counts_mqc  = Channel.empty()

    if (params.protocol == 'atac') {
        MACS3_ATAC_NARROW ( ch_macs_input.map{ meta, ip, ctrl -> [meta, ip] }, m_genome )
        MACS3_ATAC_BROAD ( ch_macs_input.map{ meta, ip, ctrl -> [meta, ip] }, m_genome )
        ch_narrow_peaks = MACS3_ATAC_NARROW.out.peaks
        ch_broad_peaks  = MACS3_ATAC_BROAD.out.peaks
        ch_narrow_counts_mqc = MACS3_ATAC_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_ATAC_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_ATAC_NARROW.out.versions.map{ it[1] }.mix(MACS3_ATAC_BROAD.out.versions.map{ it[1] })
    } else {
        MACS3_CHIP_NARROW ( ch_macs_input, m_genome )
        MACS3_CHIP_BROAD ( ch_macs_input, m_genome )
        ch_narrow_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_broad_peaks  = MACS3_CHIP_BROAD.out.peaks
        ch_narrow_counts_mqc = MACS3_CHIP_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_CHIP_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_CHIP_NARROW.out.xls.map{ it[1] }.mix(MACS3_CHIP_BROAD.out.xls.map{ it[1] })
    }

    ch_frip_input = ch_bams_branched.ip.map { meta, bam, bai -> [ meta, bam ] }.join(ch_narrow_peaks)
    CALC_FRIP ( ch_frip_input )

    ch_homer_mqc = Channel.empty()
    if (!params.skip_homer && reference_file && gtf_file) {
        ch_homer_input = ch_narrow_peaks.mix(ch_broad_peaks)
            .filter { meta, peak -> peak != null && peak.exists() && peak.size() > 0 }
        
        HOMER_ANNOTATEPEAKS ( ch_homer_input, file(reference_file), file(gtf_file) )
        ch_homer_mqc = HOMER_ANNOTATEPEAKS.out.stats_mqc.map{ it[1] }.collect().ifEmpty([])
    }

    ch_diffbind_mqc = Channel.empty()
    if (!params.skip_diffbind) {
        ch_diffbind_prep = ch_bams_branched.ip
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .join( ch_narrow_peaks.map { meta, peak -> [ meta.id, peak ] } )
            .map { id, meta, bam, bai, peak -> [ meta, bam, bai, peak ] }

        ch_db_samplesheet = ch_diffbind_prep
            .map { meta, bam, bai, peak -> "${meta.id},${meta.condition},${meta.replicate},${bam.name},${peak.name},narrow" }
            .collectFile(name: 'samplesheet_diffbind.csv', newLine: true, seed: 'SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller')

        DIFFBIND (
            ch_db_samplesheet,
            ch_final_bams.map{ it[1] }.collect(), 
            ch_final_bams.map{ it[2] }.collect(), 
            ch_narrow_peaks.map{ it[1] }.collect()
        )
        ch_diffbind_mqc = DIFFBIND.out.mqc_html.collect().ifEmpty([])
        ch_versions = ch_versions.mix(DIFFBIND.out.versions)
    }

  ch_summary_mqc = Channel.value("Protocol: ${params.protocol}\nGenome: ${params.genome}")
        .collectFile(name: 'summary.txt')
        .collect()

    ch_versions_mqc = ch_versions
        .unique()
        .collectFile(name: 'collated_versions.yml')
        .collect()
        .ifEmpty([])
    
    ch_fastqc_mqc   = FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([])
    ch_trim_mqc     = TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([])
    ch_bowtie_mqc   = BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([])
    ch_picard_mqc   = PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([])
    ch_stats_mqc    = SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([])
    ch_frip_mqc     = CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([])
    
    ch_macs_mqc     = ch_macs_logs_mqc.collect().ifEmpty([])
    
    ch_counts_mqc   = ch_narrow_counts_mqc.mix(ch_broad_counts_mqc)
        .map{ it[1] }
        .collect()
        .ifEmpty([])
    
    ch_deeptools_mqc = DEEPTOOLS.out.fingerprint_txt.map{ it[1] }
        .mix(DEEPTOOLS.out.fingerprint_metrics.map{ it[1] })
        .collect()
        .ifEmpty([])

    ch_homer_final    = ch_homer_mqc.collect().ifEmpty([])
    ch_diffbind_final = ch_diffbind_mqc.collect().ifEmpty([])

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        ch_summary_mqc,
        ch_fastqc
