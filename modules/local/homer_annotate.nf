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
    tuple val(meta), path("*.homer_stats.txt") , emit: stats 
    path "versions.yml"                         , emit: versions

    script:
    def type = peak.name.contains('narrow') ? 'narrow' : 'broad'
    def prefix = "${meta.id}.${type}"
    """
    annotatePeaks.pl \\
        $peak \\
        $fasta \\
        -gtf $gtf \\
        -cpu $task.cpus \\
        > ${prefix}.annotatePeaks.txt

    # Estraiamo le statistiche pulite
    # Usiamo 'cut' per isolare la colonna dell'annotazione (solitamente la 8a)
    echo -e "Sample\\tIntergenic\\tTTS\\texon\\tintron\\tpromoter-TSS" > ${prefix}.homer_stats.txt
    
    # Contiamo solo nella colonna specifica per evitare match errati nel nome del picco
    TOTAL=\$(wc -l < ${prefix}.annotatePeaks.txt)
    INTERGENIC=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "Intergenic" || true)
    TTS=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "TTS" || true)
    EXON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "exon" || true)
    INTRON=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -c "intron" || true)
    PROMOTER=\$(cut -f8 ${prefix}.annotatePeaks.txt | grep -i -c "promoter-TSS" || true)
    
    echo -e "${prefix}\\t\$INTERGENIC\\t\$TTS\\t\$EXON\\t\$INTRON\\t\$PROMOTER" >> ${prefix}.homer_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homer: 4.11
    END_VERSIONS
    """
}
