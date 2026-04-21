nextflow.enable.dsl=2


include { ATAC_CHIP_PIPELINE } from './workflows/analysis.nf'

/**
 * Funzione per parsare le righe del CSV e creare la struttura dati [meta, [reads]]
 */
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id         = row.sample.trim()
    meta.antibody   = (row.antibody && row.antibody.trim() != "") ? row.antibody.trim() : 'none'
    meta.control    = (row.control && row.control.trim() != "") ? row.control.trim() : 'none'
    
    // Verifica Single-End o Paired-End
    meta.single_end = (row.fastq_2 == null || row.fastq_2.trim() == "") ? true : false

    // Definizione file path (Nextflow gestisce URL S3 o percorsi locali)
    def fastq_1 = file(row.fastq_1, checkIfExists: true)
    def fastqs = [ fastq_1 ]
    
    if (!meta.single_end) {
        def fastq_2 = file(row.fastq_2, checkIfExists: true)
        fastqs << fastq_2
    }

    return [ meta, fastqs ]
}

workflow {
    
    // 1. Controllo validità parametri iniziali
    if (!params.input) { error "Errore: Specifica un file samplesheet con --input" }
    log.info "Protocollo selezionato: ${params.protocol.toUpperCase()}"
    log.info "Genoma di riferimento: ${params.genome}"

    // 2. Lettura del Samplesheet (CSV)
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header:true, sep:',')
        .map { row -> create_fastq_channel(row) }

    // 3. Gestione Indice Bowtie2
    // Se non fornito manualmente, lo cerca nel config igenomes basandosi sulla stringa del genoma
    def index_path = params.bowtie2_index ?: params.genomes[ params.genome ]?.bowtie2 ?: null

    if (!index_path) {
        error "Errore: Indice Bowtie2 non trovato per '${params.genome}'. Controlla --genome o --bowtie2_index"
    }

    // Se l'indice è una cartella S3 o locale, usiamo collect() per passarlo a tutti i task
    ch_index = Channel.fromPath(index_path, checkIfExists: true).collect()

    // 4. Lancio della Pipeline
    // Passiamo i campioni e l'indice al workflow principale
    ATAC_CHIP_PIPELINE ( ch_input, ch_index )
    
    // Messaggio finale al completamento
    workflow.onComplete {
        log.info "Pipeline completata con successo!"
        log.info "I risultati sono disponibili in: ${params.outdir}"
    }
}
