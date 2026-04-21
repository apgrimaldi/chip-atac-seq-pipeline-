// --- INCLUDE DEI MODULI ---
include { FASTQC } from '../modules/local/fastqc.nf'
include { TRIMGALORE } from '../modules/local/trimgalore.nf'
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
include { MULTIQC } from '../modules/local/multiqc.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input
    ch_index

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    // 1. Qualità iniziale
    FASTQC ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // 2. Trimming
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    // 3. Allineamento
    BOWTIE2 ( TRIMGALORE.out.reads, ch_index )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    // 4. Ordinamento
    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    // 5. Duplicati (Emette meta, bam, bai)
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, [], [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    // 6. Filtraggio Blacklist
    def blacklist_val = params.genomes[ params.genome ]?.blacklist ?: null
    if (blacklist_val) {
        // file() assicura il download se è un URL
        ch_blacklist = file(blacklist_val) 
        FILTERING ( PICARD_MARKDUPLICATES.out.bam, ch_blacklist )
        
        // Bedtools distrugge l'indice, quindi creiamo una terna con un placeholder vuoto []
        ch_final_bams = FILTERING.out.bam.map { it -> [ it[0], it[1], [] ] }
        ch_versions = ch_versions.mix(FILTERING.out.versions)
    } else {
        // Se non filtriamo, usiamo l'output di Picard (meta, bam, bai)
        ch_final_bams = PICARD_MARKDUPLICATES.out.bam
    }

    // 7. Statistiche Allineamento
    SAMTOOLS_STATS ( ch_final_bams.map { it -> [ it[0], it[1] ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // 8. DeepTools (Generazione BigWig)
    DEEPTOOLS ( 
        ch_final_bams.map { it -> [ it[0], it[1] ] },
        ch_final_bams.map { it -> it[2] } 
    )

    // 9. Peak Calling
    ch_peaks = Channel.empty()
    ch_frip_peaks = Channel.empty() 

    if (params.protocol == 'atac') {
        ch_macs_input = ch_final_bams.map { it -> [ it[0], it[1] ] }
        
        MACS3_ATAC_NARROW ( ch_macs_input )
        MACS3_ATAC_BROAD  ( ch_macs_input )
        
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks 
        ch_versions = ch_versions.mix(MACS3_ATAC_NARROW.out.versions, MACS3_ATAC_BROAD.out.versions)
    } 
    else if (params.protocol == 'chip') {
        // Separiamo Controlli (Input/IgG) e Campioni (IP)
        ch_control_bams = ch_final_bams
            .filter { it -> 
                def m = it[0]
                m.antibody == 'none' || !m.antibody || m.antibody == 'IgG' || m.antibody == '' 
            }
            .map { it -> [ it[0].id, it[1] ] }

        ch_macs3_chip_input = ch_final_bams
            .filter { it -> 
                def m = it[0]
                m.antibody && m.antibody != 'none' && m.antibody != 'IgG' && m.antibody != '' 
            }
            .map { it -> [ it[0].control, it[0], it[1] ] } 
            .join(ch_control_bams)
            .map { it -> [ it[1], it[2], it[3] ] } // [meta, bam_ip, bam_ctrl]

        MACS3_CHIP_NARROW ( ch_macs3_chip_input )
        MACS3_CHIP_BROAD  ( ch_macs3_chip_input )
        
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_versions = ch_versions.mix(MACS3_CHIP_NARROW.out.versions, MACS3_CHIP_BROAD.out.versions)
    }

    // 10. FRiP (Fraction of Reads in Peaks)
    ch_frip_input = ch_final_bams
        .map { it -> [ it[0], it[1] ] }
        .join(ch_frip_peaks)
    
    CALC_FRIP ( ch_frip_input )
    ch_versions = ch_versions.mix(CALC_FRIP.out.versions)

    // 11. Annotazione con HOMER
    def fasta = params.genomes[ params.genome ]?.fasta ?: null
    def gtf   = params.genomes[ params.genome ]?.gtf   ?: null
    if (fasta && gtf) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta), file(gtf) )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    }

    // 12. MULTIQC
    def summary_info = "Protocol: ${params.protocol.toUpperCase()}\nGenome: ${params.genome}"
    def ch_workflow_summary = Channel.value(summary_info).collectFile(name: 'workflow_summary_mqc.txt')

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        ch_workflow_summary,
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOW_LOG = BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        ch_peaks.map{ it[1] }.collect().ifEmpty([]),
        CALC_FRIP.out.summary.map{ it[1] }.collect().ifEmpty([]),
        ch_versions.unique().collect().ifEmpty([])
    )

    emit:
    bam      = ch_final_bams
    peaks    = ch_peaks
    bigwig   = DEEPTOOLS.out.bw
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
