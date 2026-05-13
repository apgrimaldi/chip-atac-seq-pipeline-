process PROFILEPLYR {
    tag "profileplyr_analysis"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-profileplyr:1.22.0--r44hdfd78af_0'

    input:
    path peaks
    path bigwigs

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

    # 1. Caricamento dei picchi
    # Legge file BED o NarrowPeak (es. quelli in uscita da LanceOtron)
    peaks_gr <- profileplyr::readPeakFile("${peaks}")

    # 2. Generazione dell'oggetto profileplyr
    # Utilizza i BigWig per mappare il segnale intorno ai picchi (+/- 2kb)
    pro_obj <- profileplyr(
        peaks_gr,
        binsize = 50,
        distance = 2000,
        sampleData = data.frame(bamReads = "${bigwigs}")
    )

    # 3. Plotting per MultiQC (Immagine rimpicciolita a 500px)
    png("profile_heatmap.png", width=1000, height=1200, res=150)
    generateEnrichedHeatmap(pro_obj)
    dev.off()

    pdf("profile_heatmap.pdf", width=7, height=9)
    generateEnrichedHeatmap(pro_obj)
    dev.off()

    # 4. Generazione codice HTML per il report
    img_64 <- base64encode("profile_heatmap.png")
    cat(paste0(
        "\\n",
        "<div style='text-align: center; padding: 20px;'>\\n",
        "  <img src='data:image/png;base64,", img_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
        "</div>"
    ), file="profileplyr_mqc.html")

    # Tracking versioni
    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    profileplyr: ", packageVersion("profileplyr"))
    ), "versions.yml")
    """
}
