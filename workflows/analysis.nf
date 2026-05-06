nextflow.enable.dsl=2

// --- INCLUDE DEI MODULI ---
include { FASTQC }                  from '../modules/local/fastqc.nf'
include { TRIMGALORE }              from '../modules/local/trimgalore.nf'
include { BOWTIE2 }                 from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }           from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }          from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES }   from '../modules/local/picard_markduplicates.nf'
include { FILTERING }               from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }       from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }        from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }       from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }        from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }     from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }               from '../modules/local/calc_frip.nf'
include { DEEPTOOLS }               from '../modules/local/deeptools.nf'
include { MULTIQC }                 from '../modules/local/multiqc.nf'
include { SAMTOOLS_INDEX }          from '../modules/local/samtools_index.nf'

workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input
    ch_index

    main:
    ch_versions = Channel.empty()
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    // --- GESTIONE DINAMICA RISORSE GENOMA ---
    def blacklist_path = params.blacklist_file ?: params.genomes[ params.genome ]?.blacklist ?: null
    def fasta_file     = params.fasta_file ?: params.genomes[ params.genome ]?.fasta ?: null
    def gtf_file       = params.gtf_file ?: params.genomes[ params.genome ]?.gtf ?: null

    // --- INIZIO STEPS ---

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
    
    ch_picard_bams = PICARD_MARKDUPLICATES.out.bam.map { it ->
        it.size() == 3 ? it : [ it[0], it[1], it[2] ?: [] ]
    }

    // 6. Filtraggio Blacklist
    if (blacklist_path) {
        ch_blacklist = file(blacklist_path)
        FILTERING ( ch_picard_bams, ch_blacklist )
        SAMTOOLS_INDEX ( FILTERING.out.bam )
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
        ch_versions = ch_versions.mix(FILTERING.out.versions, SAMTOOLS_INDEX.out.versions)
    } else {
        SAMTOOLS_INDEX ( ch_picard_bams.map { it -> [it[0], it[1]] } )
        ch_final_bams = SAMTOOLS_INDEX.out.bam_bai
    }

    // 7. Statistiche Allineamento
    SAMTOOLS_STATS ( ch_final_bams.map { it -> [ it[0], it[1] ] } )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // 8. DeepTools
    DEEPTOOLS ( ch_final_bams )
    ch_versions = ch_versions.mix(DEEPTOOLS.out.versions)

    // 9. Peak Calling
    ch_peaks = Channel.empty()
    ch_frip_peaks = Channel.empty()
    ch_narrow_counts_mqc = Channel.empty()
    ch_broad_counts_mqc  = Channel.empty()
    ch_macs_logs_mqc = Channel.empty()

    if (params.protocol == 'atac') {
        ch_macs_input = ch_final_bams.map { it -> [ it[0], it[1] ] }
        MACS3_ATAC_NARROW ( ch_macs_input )
        MACS3_ATAC_BROAD  ( ch_macs_input )
        
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks
        ch_narrow_counts_mqc = MACS3_ATAC_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_ATAC_BROAD.out.count_broad
        
        // Per ATAC raccogliamo i log/versioni
        ch_macs_logs_mqc = MACS3_ATAC_NARROW.out.versions.map{ it[1] }
                            .mix(MACS3_ATAC_BROAD.out.versions.map{ it[1] })
    } 
    else if (params.protocol == 'chip') {
        ch_macs3_chip_input = ch_final_bams.map { it -> [ it[0], it[1], [] ] } 
        MACS3_CHIP_NARROW ( ch_macs3_chip_input )
        MACS3_CHIP_BROAD  ( ch_macs3_chip_input )
        
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_narrow_counts_mqc = MACS3_CHIP_NARROW.out.count_narrow
        ch_broad_counts_mqc  = MACS3_CHIP_BROAD.out.count_broad

        // Per ChIP passiamo i file .xls che MultiQC adora per le statistiche
        ch_macs_logs_mqc = MACS3_CHIP_NARROW.out.xls.map{ it[1] }
                            .mix(MACS3_CHIP_BROAD.out.xls.map{ it[1] })
    }

    // 10. FRiP
    ch_frip_input = ch_final_bams.map { it -> [ it[0], it[1] ] }.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    // 11. Annotazione
    ch_homer_mqc = Channel.empty()
    if (fasta_file && gtf_file) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta_file), file(gtf_file) )
        ch_homer_mqc = HOMER_ANNOTATEPEAKS.out.stats.map{ it[1] }.collect().ifEmpty([])
    }

    // 12. MULTIQC
    ch_versions_multiqc = ch_versions.unique().collectFile(name: 'collated_versions.yml')
    ch_all_counts_mqc = ch_narrow_counts_mqc.mix(ch_broad_counts_mqc).collect().ifEmpty([])

    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        Channel.value("Protocol: ${params.protocol}\nGenome: ${params.genome}").collectFile(name: 'summary.txt'),
        FASTQC.out.zip.map{ it[1] }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ it[1] }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ it[1] }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ it[1] }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ it[1] }.collect().ifEmpty([]),
        DEEPTOOLS.out.bw.map{ it instanceof List ? it[1] : it }.collect().ifEmpty([]),
        ch_macs_logs_mqc.collect().ifEmpty([]), 
        ch_all_counts_mqc,                      
        CALC_FRIP.out.frip.map{ it[1] }.collect().ifEmpty([]), 
        ch_homer_mqc,
        ch_versions_multiqc.collect()                    
    )
}
