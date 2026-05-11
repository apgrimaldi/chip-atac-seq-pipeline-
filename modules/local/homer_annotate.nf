process HOMER_ANNOTATEPEAKS {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/homer:4.11--pl526hc9558a2_3'

    input:
    tuple val(meta), path(peak)
    path fasta
    path gtf

    output:
    tuple val(meta), path("*.annotatePeaks.txt"), emit: txt
    tuple val(meta), path("*.homer_stats.txt")   , emit: stats 
    path "versions.yml"                           , emit: versions

    script:
    // Determiniamo se stiamo lavorando su narrow o broad per il nome del file
    def type = peak.name.contains('narrow') ? 'narrow' : 'broad'
    def prefix = "${meta.id}.${type}"
    """
    # 1. Gestione GTF (HOMER non accetta .gz)
    GTF_FILE=\$(basename ${gtf})
    if [[ \${GTF_FILE} == *.gz ]]; then
        gunzip -c ${gtf} > reference.gtf
        ANNOT_FILE="reference.gtf"
    else
        ANNOT_FILE="${gtf}"
    fi

    # 2. Esecuzione Annotazione
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf \${ANNOT_FILE} \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    # 3. Estrazione statistiche (Colonna 8 di HOMER contiene l'annotazione)
    echo -e "Sample\\tIntergenic\\tTTS\\texon\\tintron\\tpromoter-TSS" > ${prefix}.homer_stats.txt
    
    # Calcolo conteggi (grep -c restituisce 0 grazie a || true se non trova nulla)
    INTERGENIC=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "Intergenic" || true)
    TTS=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "TTS" || true)
    EXON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "exon" || true)
    INTRON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "intron" || true)
    PROMOTER=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -i -c "promoter-TSS" || true)
    
    echo -e "${prefix}\\t\$INTERGENIC\\t\$TTS\\t\$EXON\\t\$INTRON\\t\$PROMOTER" >> ${prefix}.homer_stats.txt

    # 4. Versioni
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: 4.11
    END_VERSIONS
    """
}
