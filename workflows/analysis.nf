include { FASTQC }                from '../modules/local/fastqc.nf'
include { TRIMGALORE }            from '../modules/local/trimgalore.nf'
include { BOWTIE2 }               from '../modules/local/bowtie2.nf'
include { SAMTOOLS_SORT }         from '../modules/local/samtools_sort.nf'
include { SAMTOOLS_STATS }        from '../modules/local/samtools_stats.nf'
include { PICARD_MARKDUPLICATES } from '../modules/local/picard_markduplicates.nf'
include { FILTERING }             from '../modules/local/filtering.nf'
include { MACS3_ATAC_NARROW }     from '../modules/local/macs3_atac_narrow.nf'
include { MACS3_ATAC_BROAD }      from '../modules/local/macs3_atac_broad.nf'
include { MACS3_CHIP_NARROW }     from '../modules/local/macs3_chip_narrow.nf'
include { MACS3_CHIP_BROAD }      from '../modules/local/macs3_chip_broad.nf'
include { HOMER_ANNOTATEPEAKS }   from '../modules/local/homer_annotate.nf'
include { CALC_FRIP }             from '../modules/local/calc_frip.nf'
include { DEEPTOOLS }             from '../modules/local/deeptools.nf'
include { MULTIQC }               from '../modules/local/multiqc.nf'


workflow ATAC_CHIP_PIPELINE {
    take:
    ch_input    // [meta, [reads]]
    ch_index    // path_to_index_folder

    main:
    ch_versions = Channel.empty()

    // --- Definizione Canali Statici (Assets) ---
    // Correzione: Recuperiamo il percorso corretto per la configurazione MultiQC
    ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
  

    // --- 1. Qualità Iniziale ---
    FASTQC ( ch_input )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // --- 2. Trimming (Rimozione Adattatori) ---
    TRIMGALORE ( ch_input )
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    // --- 3. Allineamento al Genoma ---
    BOWTIE2 ( TRIMGALORE.out.reads, ch_index )
    ch_versions = ch_versions.mix(BOWTIE2.out.versions)

    // --- 4. Ordinamento e Indicizzazione BAM ---
    SAMTOOLS_SORT ( BOWTIE2.out.bam )

    // --- 5. Identificazione e Rimozione Duplicati PCR ---
    PICARD_MARKDUPLICATES ( SAMTOOLS_SORT.out.bam, [], [] )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)

    // --- 6. Filtraggio (Blacklist e Qualità) ---
    def blacklist_path = params.genomes[ params.genome ]?.blacklist ?: null
    if (blacklist_path) {
        ch_blacklist = file(blacklist_path) 
        FILTERING ( PICARD_MARKDUPLICATES.out.bam, ch_blacklist )
        // Correzione: Assegniamo esplicitamente l'output BAM + BAI a ch_final_bams
        ch_final_bams = FILTERING.out.bam
        ch_versions = ch_versions.mix(FILTERING.out.versions)
    } else {
        // Correzione: Se non c'è blacklist, usiamo l'output di Picard (che include BAI)
        ch_final_bams = PICARD_MARKDUPLICATES.out.bam
    }

    // --- 7. Statistiche Allineamento Finali ---
    SAMTOOLS_STATS ( ch_final_bams )
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    // --- 8. DeepTools (Generazione BigWig) ---
    // Correzione Sistematica: Scompattiamo esplicitamente la lista [meta, bam, bai]
    // Questo previene l'errore "Invalid method invocation call".
    DEEPTOOLS (
        ch_final_bams.map { [meta, bam, bai] -> bam }, // 1° Argomento: BAM
        ch_final_bams.map { [meta, bam, bai] -> bai }  // 2° Argomento: BAI
    )

    // --- 9. Peak Calling (Specifico per Protocollo) ---
    ch_peaks = Channel.empty()
  
    if (params.protocol == 'atac') {
        // --- 9.1 ATAC-seq Peak Calling ---
        MACS3_ATAC_NARROW ( ch_final_bams )
        MACS3_ATAC_BROAD  ( ch_final_bams )
        // Uniamo i picchi stretti e larghi
        ch_peaks = MACS3_ATAC_NARROW.out.peaks.mix(MACS3_ATAC_BROAD.out.peaks)
        // Definiamo quale set di picchi usare per il calcolo del FRiP
        ch_frip_peaks = MACS3_ATAC_NARROW.out.peaks 
        ch_versions = ch_versions.mix(MACS3_ATAC_NARROW.out.versions, MACS3_ATAC_BROAD.out.versions)
    } 
    else if (params.protocol == 'chip') {
        // --- 9.2 ChIP-seq Peak Calling (Richiede IP e Controllo) ---
        
        // Correzione: Creiamo il canale dei BAM di controllo (nessun anticorpo specificato)
        // Correzione: Scompattiamo la lista [meta, bam, bai]
        ch_control_bams = ch_final_bams
            .filter { [meta, bam, bai] -> meta.antibody == 'none' || !meta.antibody || meta.antibody == '' }
            .map { [meta, bam, bai] -> [ meta.id, bam ] } // Mappiamo ID e BAM

        // Correzione: Creiamo l'input per MACS3 (unendo IP e relativo Controllo)
        // Correzione: Scompattiamo la lista [meta, bam, bai]
        ch_macs3_chip_input = ch_final_bams
            .filter { [meta, bam, bai] -> meta.antibody && meta.antibody != 'none' && meta.antibody != '' }
            .map { [meta, bam, bai] -> [ meta.control, meta, bam ] } // Mappiamo l'ID del controllo, metadati e BAM IP
            .join(ch_control_bams) // Join basato sull'ID del controllo (meta.control == control.id)
            .map { [id_ctrl, meta, bam_ip, bam_ctrl] -> [ meta, bam_ip, bam_ctrl ] } // Riorganizziamo la tupla

        // Chiamiamo MACS3 per ChIP-seq (Narrow e Broad)
        MACS3_CHIP_NARROW ( ch_macs3_chip_input )
        MACS3_CHIP_BROAD  ( ch_macs3_chip_input )
        // Uniamo i picchi
        ch_peaks = MACS3_CHIP_NARROW.out.peaks.mix(MACS3_CHIP_BROAD.out.peaks)
        // Definiamo i picchi per il FRiP score
        ch_frip_peaks = MACS3_CHIP_NARROW.out.peaks
        ch_versions = ch_versions.mix(MACS3_CHIP_NARROW.out.versions, MACS3_CHIP_BROAD.out.versions)
    }

    // --- 10. Calcolo del FRiP Score ---
    // Join basato sui metadati (primo elemento della tupla)
    // Correzione: Scompattiamo la lista [meta, bam, bai] per assicurarci di prendere meta e bam
    ch_frip_input = ch_final_bams.map{ [meta, bam, bai] -> [meta, bam] }.join(ch_frip_peaks)
    CALC_FRIP ( ch_frip_input )

    // --- 11. Annotazione Funzionale dei Picchi ---
    def fasta = params.genomes[ params.genome ]?.fasta ?: null
    def gtf   = params.genomes[ params.genome ]?.gtf   ?: null
    // Eseguiamo l'annotazione solo se fasta e gtf sono definiti
    if (fasta && gtf) {
        HOMER_ANNOTATEPEAKS ( ch_peaks, file(fasta), file(gtf) )
        ch_versions = ch_versions.mix(HOMER_ANNOTATEPEAKS.out.versions)
    }

    // --- 12. Generazione del Report MultiQC ---
    
    // Generiamo il file di riepilogo della configurazione
    def summary_info = """
    Pipeline Configuration:
    - Protocol: ${params.protocol.toUpperCase()}
    - Genome:   ${params.genome}
    - Input:    ${params.input}
    - Outdir:   ${params.outdir}
    """.stripIndent()
    
    def ch_workflow_summary = Channel.value(summary_info).collectFile(name: 'workflow_summary_mqc.txt')

    // Raccogliamo tutti i log e i risultati dai moduli precedenti
    MULTIQC (
        ch_multiqc_config.collect().ifEmpty([]),
        ch_workflow_summary,
        // Correzione Sistematica: Scompattiamo la lista [meta, path] per prendere solo il file
        FASTQC.out.zip.map{ [meta, path] -> path }.collect().ifEmpty([]),
        TRIMGALORE.out.log.map{ [meta, path] -> path }.collect().ifEmpty([]),
        BOWTIE2.out.log.map{ [meta, path] -> path }.collect().ifEmpty([]),
        PICARD_MARKDUPLICATES.out.metrics.map{ [meta, path] -> path }.collect().ifEmpty([]),
        SAMTOOLS_STATS.out.stats.map{ [meta, path] -> path }.collect().ifEmpty([]),
        ch_peaks.map{ [meta, path] -> path }.collect().ifEmpty([]),            
        CALC_FRIP.out.summary.map{ [meta, path] -> path }.collect().ifEmpty([]),
        ch_versions.unique().collect().ifEmpty([])
    )

    // --- Definizione degli Output del Workflow ---
    emit:
    bam      = ch_final_bams
    peaks    = ch_peaks
    bigwig   = DEEPTOOLS.out.bw
    multiqc  = MULTIQC.out.report
    versions = ch_versions
}
