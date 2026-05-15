process PROFILEPLYR {
    tag "${label}"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    publishDir "${params.outdir}/06_visualization/profileplyr_${label}", mode: 'copy'

    input:
    path peaks
    path bigwigs
    val label

    output:
    path "*.pdf"                       , emit: pdf, optional: true
    path "*.png"                       , emit: png, optional: true
    path "profileplyr_mqc.html"        , emit: mqc_html, optional: true
    path "versions.yml"                , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(profileplyr)
    library(base64enc)
    library(rtracklayer)

    # Identificazione file
    peak_files <- list.files(pattern = "\\\\.(bed|narrowPeak|broadPeak)\$")
    bw_files <- list.files(pattern = "\\\\.(bw|bigWig)\$")

    # Importazione picchi
    peaks_gr <- rtracklayer::import(peak_files[1])

    # Creazione oggetto Profileplyr (Nota la 'P' maiuscola)
    # Usiamo i nomi dei file puliti per la colonna sampleData
    pro_obj <- Profileplyr(
        peaks_gr,
        binsize = 50,
        distance = 2000,
        sampleData = data.frame(bamReads = bw_files, row.names = basename(bw_files))
    )

    # Generazione Heatmap
    png("profile_heatmap.png", width=1000, height=1200, res=150)
    generateEnrichedHeatmap(pro_obj)
    dev.off()

    pdf("profile_heatmap.pdf", width=7, height=9)
    generateEnrichedHeatmap(pro_obj)
    dev.off()

    # Preparazione HTML per MultiQC
    img_64 <- base64encode("profile_heatmap.png")
    cat(paste0(
        "\\n",
        "<div style='text-align: center; padding: 20px;'>\\n",
        "  <h3>Profile Analysis: ", "${label}", "</h3>\\n",
        "  <img src='data:image/png;base64,", img_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
        "</div>"
    ), file="profileplyr_mqc.html")

    # Versioni
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    profileplyr: ", packageVersion("profileplyr")),
        paste0("    rtracklayer: ", packageVersion("rtracklayer"))
    ), "versions.yml")
    """
}
