questo è l'analysis vecchio prima delle ultime modifiche : 
include { FASTQC }                 from '../modules/local/fastqc.nf'
include { TRIMGALORE }             from '../modules/local/trimgalore.nf'
include { BOWTIE2 }                from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }          from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }         from '../modules/local/samtools_stats.nf' 
include { PICARD_MARKDUPLICATES }  from '../modules/local/picard_markduplicates.nf'
include { FILTERING }              from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }      from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }       from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }      from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }       from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }    from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }              from '../modules/local/calc_frip.nf'      
include { DEEPTOOLS }              from '../modules/local/deeptools.nf'      
include { MULTIQC }                from '../modules/local/multiqc.nf'        

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input    // [meta, [reads]]
    ch_index    // path_to_index_folder

    main:
    ch_versions = Channel.empty()

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

    // 5. Duplicati
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, [], [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    // 6. Filtraggio
    def blacklist_path = params.genomes[ params.genome ]?.blacklist ?: null
    if (blacklist_path) {
        ch_blacklist = file(blacklist_path) 
        FILTERING ( PICARD_MARKDUPLICATES.out.bam, ch_blacklist )
        ch_final_bams = FILTERING.out.bam
        ch_versions = ch_versions.mix(FILTERING.out.versions)
    } else {
        ch_final_bams = PICARD_MARKDUPLICATES.out.bam
    }

    // 7. Statistiche Allineamento
    SAMTOOLS_STATS ( ch_final_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // 8. DeepTools (BigWig)
    DEEPTOOLS ( ch_final_bams )

    // 9. Peak Calling
    ch_peaks = Channel.empty()
    ch_frip_peaks = Channel.empty() 

    if (params.protocol == 'atac') {
        MACS3_ATAC_NARROW ( ch_final_bams )
        MACS3_ATAC_BROAD  ( ch_final_bams )
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks 
        ch_versions = ch_versions.mix(MACS3_ATAC_NARROW.out.versions, MACS3_ATAC_BROAD.out.versions)
    } 
    else if (params.protocol == 'chip') {
        ch_control_bams = ch_final_bams
            .filter { meta, bam -> meta.antibody == 'none' || !meta.antibody || meta.antibody == '' }
            .map { meta, bam -> [ meta.id, bam ] }

        ch_macs3_chip_input = ch_final_bams
            .filter { meta, bam -> meta.antibody && meta.antibody != 'none' && meta.antibody != '' }
            .map { meta, bam -> [ meta.control, meta, bam ] } 
            .join(ch_control_bams)
            .map { id_ctrl, meta, bam_ip, bam_ctrl -> [ meta, bam_ip, bam_ctrl ] }

        MACS3_CHIP_NARROW ( ch_macs3_chip_input )
        MACS3_CHIP_BROAD  ( ch_macs3_chip_input )
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_versions = ch_versions.mix(MACS3_CHIP_NARROW.out.versions, MACS3_CHIP_BROAD.out.versions)
    }

    // 10. FRiP
    ch_frip_input = ch_final_bams.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    // 11. Annotazione
    def fasta = params.genomes[ params.genome ]?.fasta ?: null
    def gtf   = params.genomes[ params.genome ]?.gtf   ?: null
    if (fasta && gtf) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta), file(gtf) )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    }

    // --- 12. MULTIQC ---
    
    // Creazione del sommario (Genoma, Protocollo, etc.)
    def summary_info = """
    Pipeline Configuration:
    - Protocol: ${params.protocol.toUpperCase()}
    - Genome:   ${params.genome}
    - Input:    ${params.input}
    - Outdir:   ${params.outdir}
    """.stripIndent()
    
    // Trasforma la stringa in un file che MultiQC può leggere
    def ch_workflow_summary = Channel.value(summary_info).collectFile(name: 'workflow_summary_mqc.txt')

    // Chiamata al modulo MULTIQC
    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        ch_workflow_summary,
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        ch_peaks.map{ it[1] }.collect().ifEmpty([]),            // Logs di MACS3
        CALC_FRIP.out.summary.map{ it[1] }.collect().ifEmpty([]),
        ch_versions.unique().collect().ifEmpty([])             // Versioni software
    )

    emit:
    bam      = ch_final_bams
    peaks    = ch_peaks
    bigwig   = DEEPTOOLS_BAMCOVERAGE.out.bw
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
