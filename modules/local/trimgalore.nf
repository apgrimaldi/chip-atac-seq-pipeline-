process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high' 

    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

    publishDir "${params.outdir}/02_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    // Cattura i file con i suffissi standard di TrimGalore
    tuple val(meta), path("*.fq.gz")       , emit: reads
    tuple val(meta), path("*_report.txt")  , emit: log
    path "versions.yml"                    , emit: versions

    script:
    // TrimGalore suggerisce di allocare meno core di quelli totali del task 
    // perché lancia processi paralleli per cutadapt e pigz
    def cores = task.cpus ? Math.max(Math.floor(task.cpus / 2) as int, 1) : 2

    if (meta.single_end) {
        // --- LOGICA SINGLE-END ---
        """
        trim_galore \\
            --cores $cores \\
            --gzip \\
            $reads

        cat <<EOF > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
            cutadapt: \$(cutadapt --version | head -n 1)
        EOF
        """
    } else {
        // --- LOGICA PAIRED-END ---
        """
        trim_galore \\
            --cores $cores \\
            --paired \\
            --gzip \\
            ${reads[0]} \\
            ${reads[1]}

        cat <<EOF > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/ .*\$//')
            cutadapt: \$(cutadapt --version | head -n 1)
        EOF
        """
    }
}
