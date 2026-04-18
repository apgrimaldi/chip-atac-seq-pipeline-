process BOWTIE2 {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input:
    tuple val(meta), path(reads)
    path index // La cartella dell'indice (passata come canale collect)

    output:
    tuple val(meta), path("*.raw.bam"), emit: bam
    tuple val(meta), path("*.log")     , emit: log
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rg_args = "--rg-id ${prefix} --rg SM:${prefix} --rg PL:ILLUMINA --rg LB:lib1"
    
    // Logica per gestire Single-End vs Paired-End
    def input_reads = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    
    // Aggiungiamo flag specifici per ATAC se necessario (opzionale, basato sui tuoi params)
    def extra_args = params.protocol == 'atac' ? "--no-mixed --no-discordant" : ""

    """
    # Trova la base dell'indice in modo robusto (cerca file .1.bt2 o .1.bt2l per indici grandi)
    INDEX_BASE=\$(find -L . -name "*.1.bt2*" | sed 's/\\.1\\.bt2.*//' | head -n 1)

    bowtie2 \\
        -x \$INDEX_BASE \\
        $input_reads \\
        -p $task.cpus \\
        $rg_args \\
        --very-sensitive \\
        $extra_args \\
        -X 2000 \\
        2> ${prefix}.bowtie2.log \\
        | samtools view -@ $task.cpus -b > ${prefix}.raw.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
