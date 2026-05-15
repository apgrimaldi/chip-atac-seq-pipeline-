nextflow.enable.dsl=2

include { FASTQC } from '../modules/local/fastqc.nf'
include { TRIMGALORE } from '../modules/local/trimgalore.nf'
include { BOWTIE2_BUILD } from '../modules/local/bowtie2_build.nf'
include { BOWTIE2 } from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT } from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS } from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES } from '../modules/local/picard_markduplicates.nf'
include { FILTERING } from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW } from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD } from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW } from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD } from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS } from '../modules/local/homer_annotate.nf'
include { CALC_FRIP } from '../modules/local/calc_frip.nf'
include { DEEPTOOLS } from '../modules/local/deeptools.nf'
include { DIFFBIND } from '../modules/local/diffbind.nf'
include { LANCEOTRON } from '../modules/local/lanceotron.nf'
include { MULTIQC } from '../modules/local/multiqc.nf'
include { SAMTOOLS_INDEX } from '../modules/local/samtools_index.nf'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/local/samtools_index.nf'
include { PROFILEPLYR as PROFILEPLYR_LANCE } from '../modules/local/profileplyr.nf'
include { PROFILEPLYR as PROFILEPLYR_MACS } from '../modules/local/profileplyr.nf'

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
    if (!m_genome || m_genome == 'custom') { m_genome = 'hs' }

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

    ch_bams_branched = ch_final_bams.branch { meta, bam, bai ->
        control: meta.is_control == true
        ip:      true
    }

    if (params.protocol == 'atac') {
        ch_lanceotron_input = ch_bams_branched.ip
            .join(DEEPTOOLS.out.bw_lanceotron)
            .map { meta, bam, bai, bw -> [ meta, bam, bw, [], [] ] }
    } else {
        ch_ip_l6n = ch_bams_branched.ip
            .join(DEEPTOOLS.out.bw_lanceotron)
            .map { meta, bam, bai, bw -> [ meta.control, meta, bam, bw ] }

        ch_ctrl_l6n = ch_bams_branched.control
            .join(DEEPTOOLS.out.bw_lanceotron)
            .map { meta, bam, bai, bw -> [ meta.id, bam, bw ] }

        ch_lanceotron_input = ch_ip_l6n
            .combine(ch_ctrl_l6n, by: 0)
            .map { ctrl_id, meta, ip_bam, ip_bw, ctrl_bam, ctrl_bw -> 
                [ meta, ip_bam, ip_bw, ctrl_bam, ctrl_bw ] 
            }
    }

    LANCEOTRON ( ch_lanceotron_input )
    ch_versions = ch_versions.mix(LANCEOTRON.out.versions)

    if (params.protocol == 'atac') {
        ch_macs_input = ch_bams_branched.ip.map { meta, bam, bai -> [ meta, bam, [] ] }
        MACS3_ATAC_NARROW ( ch_macs_input.map{ meta, ip, ctrl -> [meta, ip] }, m_genome )
        MACS3_ATAC_BROAD ( ch_macs_input.map{ meta, ip, ctrl -> [meta, ip] }, m_genome )
        ch_narrow_peaks = MACS3_ATAC_NARROW.out.peaks
        ch_broad_peaks  = MACS3_ATAC_BROAD.out.peaks
        ch_narrow_counts_mqc = MACS3_ATAC_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_ATAC_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_ATAC_NARROW.out.versions.map{ it[1] }.mix(MACS3_ATAC_BROAD.out.versions.map{ it[1] })
    } else {
        ch_ip_m3 = ch_bams_branched.ip.map { meta, bam, bai -> [ meta.control, meta, bam ] }
        ch_ct_m3 = ch_bams_branched.control.map { meta, bam, bai -> [ meta.id, bam ] }
        ch_macs_chip_input = ch_ip_m3.combine(ch_ct_m3, by: 0).map { cid, meta, ip, ct -> [ meta, ip, ct ] }

        MACS3_CHIP_NARROW ( ch_macs_chip_input, m_genome )
        MACS3_CHIP_BROAD ( ch_macs_chip_input, m_genome )
        ch_narrow_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_broad_peaks  = MACS3_CHIP_BROAD.out.peaks
        ch_narrow_counts_mqc = MACS3_CHIP_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_CHIP_BROAD.out.count_broad
        ch_macs_logs_mqc = MACS3_CHIP_NARROW.out.xls.map{ it[1] }.mix(MACS3_CHIP_BROAD.out.xls.map{ it[1] })
    }

    ch_frip_input = ch_bams_branched.ip.map { meta, bam, bai -> [ meta.id, meta, bam ] }
        .join(ch_narrow_peaks.map { meta, peak -> [ meta.id, peak ] })
        .map { id, meta, bam, peak -> [ meta, bam, peak ] }
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
            .map { meta, bam, bai, peak -> 
                def cond = meta.condition ?: meta.antibody
                def repl = meta.replicate ?: "1"
                "${meta.id},${cond},${repl},${bam.name},${peak.name},narrow" 
            }
            .collectFile(name: 'samplesheet_diffbind.csv', newLine: true, seed: 'SampleID,Condition,Replicate,bamReads,Peaks,PeakCaller')

        DIFFBIND ( ch_db_samplesheet, ch_final_bams.map{ it[1] }.collect(), ch_final_bams.map{ it[2] }.collect(), ch_narrow_peaks.map{ it[1] }.collect() )
        ch_diffbind_mqc = DIFFBIND.out.mqc_html.mix(DIFFBIND.out.mqc_txt).collect().ifEmpty([])
        ch_versions = ch_versions.mix(DIFFBIND.out.versions)
    }

    PROFILEPLYR_LANCE ( 
        LANCEOTRON.out.peaks.map{ it[1] }.collect(), 
        DEEPTOOLS.out.bw_display.map{ it[1] }.collect(),
        "lanceotron"
    )

    PROFILEPLYR_MACS ( 
        ch_narrow_peaks.map{ it[1] }.collect(), 
        DEEPTOOLS.out.bw_display.map{ it[1] }.collect(),
        "macs"
    )

    ch_profileplyr_mqc = PROFILEPLYR_LANCE.out.mqc_html
        .mix(PROFILEPLYR_MACS.out.mqc_html)
        .collect()
        .ifEmpty([])

    ch_versions = ch_versions.mix(
        PROFILEPLYR_LANCE.out.versions,
        PROFILEPLYR_MACS.out.versions
    )

    ch_summary_mqc = Channel.value("Protocol: ${params.protocol}\nGenome: ${params.genome}").collectFile(name: 'summary.txt').collect()
    ch_versions_mqc = ch_versions.unique().collectFile(name: 'collated_versions.yml').collect().ifEmpty([])
    
    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        ch_summary_mqc,
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        DEEPTOOLS.out.fingerprint_txt.map{ it[1] }.mix(DEEPTOOLS.out.fingerprint_metrics.map{ it[1] }).collect().ifEmpty([]),
        ch_macs_logs_mqc.collect().ifEmpty([]),
        ch_narrow_counts_mqc.mix(ch_broad_counts_mqc).map{ it[1] }.collect().ifEmpty([]),
        CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([]),
        ch_homer_mqc,
        ch_diffbind_mqc,
        ch_profileplyr_mqc,
        LANCEOTRON.out.counts_mqc.collect().ifEmpty([]),
        ch_versions_mqc
    )
}
