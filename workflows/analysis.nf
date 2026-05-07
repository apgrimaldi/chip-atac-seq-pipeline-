nextflow.enable.dsl=2

// --- INCLUDE DEI MODULI ---
include { FASTQC }                 from '../modules/local/fastqc.nf'
include { TRIMGALORE }             from '../modules/local/trimgalore.nf'
include { BOWTIE2_BUILD }          from '../modules/local/bowtie2_build.nf'
include { BOWTIE2 }                from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }          from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }         from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES }  from '../modules/local/picard_markduplicates.nf'
include { FILTERING }              from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }      from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }       from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }      from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }       from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }    from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }              from '../modules/local/calc_frip.nf'
include { DEEPTOOLS }              from '../modules/local/deeptools.nf'
include { MULTIQC }                from '../modules/local/multiqc.nf'
include { SAMTOOLS_INDEX }         from '../modules/local/samtools_index.nf'
// Aggiungiamo un secondo alias per l'indice post-filtraggio se necessario
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_FINAL } from '../modules/local/samtools_index.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input 

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    // --- LOGICA DI ASSEGNAZIONE PARAMETRI GENOMA ---
    def fasta_file     = params.fasta_file
    def gtf_file       = params.gtf_file
    def bowtie2_index  = params.bowtie2_index
    def blacklist_path = params.blacklist ?: params.blacklist_file // Gestisce entrambi i nomi
    def m_genome       = params.macs_gsize ?: params.genome

    // --- GESTIONE INDICE BOWTIE2 ---
    ch_index_internal = Channel.empty()

    if (bowtie2_index) {
        ch_index_internal = Channel.fromPath("${bowtie2_index}/*.bt2*").collect()
    } else if (fasta_file) {
        BOWTIE2_BUILD ( file(fasta_file) )
        ch_index_internal = BOWTIE2_BUILD.out.index.collect()
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    } else {
        error "ERRORE: Fornire --fasta_file o --bowtie2_index"
    }

    // --- START PIPELINE ---

    // 1. QUALITY CONTROL & TRIMMING
    FASTQC ( ch_input )
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions, TRIMGALORE.out.versions)

    // 2. ALIGNMENT
    BOWTIE2 ( TRIMGALORE.out.reads, ch_index_internal )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    // 3. SORTING
    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    // 4. DEDUPLICATION (Picard)
    // Passiamo il fasta come file semplice per evitare problemi di tuple
    def ch_fasta_ref = fasta_file ? file(fasta_file) : []
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, ch_fasta_ref, [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    // 5. INDEXING POST-DEDUP
    SAMTOOLS_INDEX ( PICARD_MARKDUPLICATES.out.bam )
    
    // 6. BLACKLIST FILTERING (Se presente)
    if (blacklist_path) {
        // FILTERING si aspetta [meta, bam, bai]
        FILTERING ( SAMTOOLS_INDEX.out.bam_bai, file(blacklist_path) )
        
        // Dopo il filtraggio serve un nuovo indice perché il BAM è cambiato
        SAMTOOLS_INDEX_FINAL ( FILTERING.out.bam )
        ch_final_bams = SAMTOOLS_INDEX_FINAL.out.bam_bai
        ch_versions = ch_versions.mix(FILTERING.out.versions)
    } else {
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
    }

    // 7. QUALITY METRICS (Samtools Stats & DeepTools)
    SAMTOOLS_STATS ( ch_final_bams.map { meta, bam, bai -> [ meta, bam ] } )
    DEEPTOOLS ( ch_final_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions, DEEPTOOLS.out.versions)

    // 8. PEAK CALLING
    ch_peaks = Channel.empty()
    ch_macs_logs_mqc = Channel.empty()
    ch_frip_peaks = Channel.empty()
    
    // Prepariamo l'input pulito per MACS3 [meta, bam]
    ch_macs_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam ] }

    if (params.protocol == 'atac') {
        MACS3_ATAC_NARROW ( ch_macs_input )
        MACS3_ATAC_BROAD  ( ch_macs_input )
        
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks
        ch_macs_logs_mqc = MACS3_ATAC_NARROW.out.versions.map{ it[1] }.mix(MACS3_ATAC_BROAD.out.versions.map{ it[1] })
    } 
    else if (params.protocol == 'chip') {
        // Per il ChIP-seq passiamo [meta, bam, control_bam] -> qui control è vuoto se non gestito da samplesheet
        ch_macs_chip_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam, [] ] }
        MACS3_CHIP_NARROW ( ch_macs_chip_input )
        MACS3_CHIP_BROAD  ( ch_macs_chip_input )
        
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_macs_logs_mqc = MACS3_CHIP_NARROW.out.xls.map{ it[1] }.mix(MACS3_CHIP_BROAD.out.xls.map{ it[1] })
    }

    // 9. ANNOTATION & FRIP
    // Uniamo i BAM con i picchi per calcolare il FRiP
    ch_frip_input = ch_final_bams.map { meta, bam, bai -> [ meta, bam ] }.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    ch_homer_mqc = Channel.empty()
    if (fasta_file && gtf_file) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta_file), file(gtf_file) )
        ch_homer_mqc = HOMER_ANNOTATEPEAKS.out.stats.map{ it[1] }.collect().ifEmpty([])
    }

    // 10. MULTIQC
    ch_versions_multiqc = ch_versions.unique().collectFile(name: 'collated_versions.yml')

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        Channel.value("Protocol: ${params.protocol}\nGenome Size: ${m_genome}").collectFile(name: 'summary.txt'),
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        DEEPTOOLS.out.fingerprint_txt.map{ it[1] }.collect().ifEmpty([]),
        ch_macs_logs_mqc.collect().ifEmpty([]),
        CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([]),
        ch_homer_mqc,
        ch_versions_multiqc.collect()
    )
}
