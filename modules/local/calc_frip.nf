process CALC_FRIP {
    tag "$meta.id"
    label 'process_medium' 
    container "quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61231030a9fcb-0"

    input:
    tuple val(meta), path(bam), path(peak)

    output:
    tuple val(meta), path("*.FRiP.txt"), emit: frip
    path "versions.yml"                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # 1. Conta le reads totali mappate
    TOTAL_MAPPED=\$(samtools view -c -F 4 $bam)

    # 2. Calcola le reads nei picchi
    READS_IN_PEAKS=\$(samtools view -u $bam | bedtools intersect -a stdin -b $peak -u -bed | wc -l)

    # 3. Calcolo finale senza header per il bargraph MultiQC
    awk -v rip="\$READS_IN_PEAKS" -v total="\$TOTAL_MAPPED" -v samp="${prefix}" \\
        'BEGIN {
            OFS="\\t"; 
            if (total > 0) print samp, rip/total; else print samp, "0"
        }' > ${prefix}.FRiP.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
