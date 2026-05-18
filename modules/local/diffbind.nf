process DIFFBIND {
    tag "diffbind_analysis"
    label 'process_high'
    container 'quay.io/biocontainers/bioconductor-diffbind:3.20.0--r45ha27e39d_0'

    input:
    path samplesheet
    path bams
    path bais
    path peaks

    output:
    path "*.pdf"                       , emit: pdf, optional: true
    path "*.csv"                       , emit: csv, optional: true
    path "*_mqc.html"                  , emit: mqc_html, optional: true
    path "diffbind_correlation_mqc.txt", emit: mqc_txt, optional: true
    path "*.png"                       , emit: png, optional: true
    path "versions.yml"                , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(DiffBind)
    library(base64enc)

    samples <- read.csv("${samplesheet}")
    samples\$bamReads <- basename(as.character(samples\$bamReads))
    samples\$Peaks    <- basename(as.character(samples\$Peaks))
    if ("bamControl" %in% colnames(samples)) {
        samples\$bamControl <- basename(as.character(samples\$bamControl))
    }

    db_obj <- dba(sampleSheet=samples)
    
    sample_info <- dba.show(db_obj)
    keep_mask <- as.numeric(sample_info\$Intervals) > 0
    if(sum(keep_mask) < length(keep_mask)) {
        db_obj <- dba(db_obj, mask=keep_mask)
    }

    png("diffbind_correlation.png", width=1000, height=1000, res=150)
    plot(db_obj, margin=20)
    dev.off()

    img_corr_64 <- base64encode("diffbind_correlation.png")
    
    cat(paste0(
        "\\n",
        "<div style='text-align: center; padding: 20px;'>\\n",
        "  <img src='data:image/png;base64,", img_corr_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
        "</div>"
    ), file="diffbind_corr_mqc.html")

    db_obj <- dba.count(db_obj, bParallel=FALSE, bUseSummarizeOverlaps=TRUE)

    try({
        cor_matrix <- dba.overlap(db_obj, mode=DBA_OL_COR)
        write.table(cor_matrix, file="diffbind_correlation_mqc.txt", sep="\t", quote=FALSE, col.names=NA)
    }, silent=TRUE)

    analysis_status <- try({
        contrast_category <- if ("Condition" %in% colnames(samples) && length(unique(samples\$Condition)) > 1) DBA_CONDITION else DBA_ANTIBODY
        db_obj <- dba.contrast(db_obj, categories=contrast_category, minMembers=2)
        db_obj <- dba.analyze(db_obj)
    }, silent=TRUE)

    if (!inherits(analysis_status, "try-error") && !is.null(db_obj\$contrasts)) {
        res_db <- dba.report(db_obj)
        write.csv(as.data.frame(res_db), "diff_bind_results.csv")

        png("diffbind_pca.png", width=1000, height=800, res=150)
        dba.plotPCA(db_obj, attributes=contrast_category, label=DBA_ID)
        dev.off()

        img_pca_64 <- base64encode("diffbind_pca.png")
        
        cat(paste0(
            "\\n",
            "<div style='text-align: center; padding: 20px;'>\\n",
            "  <img src='data:image/png;base64,", img_pca_64, "' style='width: 500px; max-width: 100%; height: auto;'>\\n",
            "</div>"
        ), file="diffbind_pca_mqc.html")
    }

    writeLines(c(
        "\\"${task.process}\\":",
        paste0("    diffbind: ", packageVersion("DiffBind"))
    ), "versions.yml")
    """
}
