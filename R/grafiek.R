#' Import data from Galaxy DESeq2
#'
#' Reads in DESeq2 output from Galaxy into a tibble.
#'
#' @param file Path to the DESeq2 results file.
#' @return A tibble containing the DESeq2 results with columns for gene expression statistics.
#' @export
import_galaxy_deseq <- function(file){
  tidy_degs <- read_delim(file, col_names = c("gene", "base_mean", "avg_log_fc", "std_err", "wald_stats", "p_val", "p_val_adj"))
  return(tidy_degs)
}



#' Generates a volcano plot from differential expression data, highlighting specified gene sets.
#'
#' @param tidy_degs A tibble containing differential expression results.
#' @param label_features A vector of gene names to highlight in the plot.
#' @param gene_column The column name for genes in the deg_file. Default is "gene".
#' @param logfc_column The column name for log fold change values. Default is "avg_log_fc".
#' @param pvalue_column The column name for p-values. Default is "p_val".
#' @param padj_column The column name for adjusted p-values. Default is "p_val_adj".
#' @param label_padj The cutoff for adjusted p-values to label genes. Default is 0.05.
#' @param label_log2fc The log2 fold change cutoff for labeling genes. Default is 0.
#' @param point_alpha Transparency of points in the plot. Default is 1.
#' @param gs_alpha Transparency of gene set points. Default is 1.
#' @param pt_size Size of the points. Default is 1.
#' @param gs_size Size of the points representing gene sets. Default is 2.
#' @param base_color Color for non-significant points. Default is "gray".
#' @param sig_color Color for significant points. Default is "black".
#' @param gs_color Color for points in the gene set. Default is "purple".
#' @param logfc_collapse A threshold for collapsing log fold change values.
#' @param pval_collapse A threshold for collapsing p-values.
#' @param feature_label Whether to label the genes in the gene set. Default is TRUE.
#' @param max_overlaps Maximum number of label overlaps. Default is 10.
#' @return A ggplot object showing the volcano plot.
#' @export
volcano_plot_gs <- function(tidy_degs, label_features, gene_column = "gene", logfc_column = "avg_log_fc",
                            pvalue_column = "p_val", padj_column = "p_val_adj",
                            label_padj = .05, label_log2fc = 0, point_alpha = 1,
                            gs_alpha = 1, pt_size = 1, gs_size = 2,
                            base_color = "gray", sig_color = "black", gs_color = "purple",
                            logfc_collapse = NULL, pval_collapse = NULL,
                            feature_label = TRUE, max_overlaps = 10){
  tryCatch({
    if(!is.null(logfc_collapse)){
      tidy_degs[[logfc_column]][tidy_degs[[logfc_column]] < -logfc_collapse] <- -logfc_collapse
      tidy_degs[[logfc_column]][tidy_degs[[logfc_column]] > logfc_collapse] <- logfc_collapse
    }
    if(!is.null(pval_collapse)){
      tidy_degs[[pvalue_column]][tidy_degs[[pvalue_column]] < pval_collapse] <- pval_collapse
    }
    label_features <- label_features %>% tolower() %>% Hmisc::capitalize()
    missing_features <- label_features[!(label_features %in% tidy_degs[[gene_column]])]
    label_data <- filter(tidy_degs, tidy_degs[[gene_column]] %in% label_features)
    sig_data <- filter(tidy_degs, abs(tidy_degs[[logfc_column]]) > label_log2fc & tidy_degs[[padj_column]] < label_padj)
    p1 <- ggplot(data = tidy_degs, aes(x = tidy_degs[[logfc_column]], y = -log10(tidy_degs[[pvalue_column]]))) +
      geom_point(alpha = point_alpha, color = base_color, size = pt_size) +
      geom_point(alpha = point_alpha, size = pt_size, data = sig_data, color = sig_color, aes(x = sig_data[[logfc_column]], y = -log10(sig_data[[pvalue_column]]))) +
      theme_classic() +
      geom_point(size = gs_size, color = gs_color, alpha = gs_alpha, data = label_data,
                 aes(x = label_data[[logfc_column]], y = -log10(label_data[[pvalue_column]]))) +
      xlab("LogFC") +
      ylab("-Log10(p value)")
    if(feature_label){
      p1 <- p1 + ggrepel::geom_text_repel(segment.alpha = .5, color = gs_color, data = label_data,
                                          aes(x = label_data[[logfc_column]], y = -log10(label_data[[pvalue_column]]),
                                              label = label_data[[gene_column]]), max.overlaps = max_overlaps, fontface = "italic")
    }
    if(length(missing_features) > 0){
      warning(paste0("missing feature: ",  missing_features, " "))
    }
    return(p1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

volcano_plot <- function(tidy_degs, gene_column = "geneid", logfc_column = "log2fc",
                         pvalue_column = "pvalue", padj_column = "padj", color_by = "cutoffs",
                         cutoff_log2fc = 0, cutoff_padj = .05, point_alpha = 1){
  tryCatch({
    if(color_by == "cutoffs"){
      label.data <- filter(tidy_degs, abs(tidy_degs[[logfc_column]]) > cutoff_log2fc & tidy_degs[[padj_column]] < cutoff_padj)
      p1 <- ggplot(data = tidy_degs, aes(x = tidy_degs[[logfc_column]], y = -log10(tidy_degs[[pvalue_column]]))) +
        geom_point(alpha = point_alpha, color = "gray") +
        geom_point(alpha = point_alpha, data = label.data, color = "black", aes(x = label.data[[logfc_column]], y = -log10(label.data[[pvalue_column]]))) +
        geom_text_repel(data = label.data, color = "black", label = label.data[[gene_column]], aes(x = label.data[[logfc_column]], y = -log10(label.data[[pvalue_column]]))) +
        theme_classic() +
        xlab("LogFC") +
        ylab("-Log10(p value)")
    }else{
      p1 <- ggplot(data = tidy_degs, aes(x = tidy_degs[[logfc_column]], y = -log10(tidy_degs[[pvalue_column]]), color = tidy_degs[[color_by]])) +
        geom_point(alpha = point_alpha) +
        scale_color_gradient(name = color_by, trans = "log") +
        theme_classic() +
        xlab("LogFC") +
        ylab("-Log10(p value)")
    }
    return(p1)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#' Convert DESeq2 output from Galaxy to GSEA rank file format
#'
#' Converts the DESeq2 output into a format suitable for GSEA rank files.
#'
#' @param tidy_degs A tibble containing DESeq2 results.
#' @return A named numeric vector where the names are gene IDs and the values are the log fold changes.
#' @export
galaxy_deseq_to_gsea_rnk <- function(tidy_degs){
  temp <- tidy_degs %>% select(gene, avg_log_fc)
  colnames(temp) <- c("id", "t")
  temp <- temp[complete.cases(temp),]
  temp$id <- temp$id %>% toupper()
  temp <- setNames(temp$t, temp$id)
  return(temp)
}

#' Create a bar plot of NES values from FGSEA results
#'
#' Plots normalized enrichment scores (NES) from FGSEA output, highlighting significant pathways.
#'
#' @param fgsea_output A tibble containing the FGSEA results.
#' @param padj_cutoff Adjusted p-value cutoff for highlighting significant pathways. Default is 0.05.
#' @return A ggplot object showing the NES values of the pathways.
#' @export
fgsea_nes_tree <- function(fgsea_output, padj_cutoff = .05){
  sig_data <- filter(fgsea_output, padj <= padj_cutoff)
  fgsea_output %>%
    arrange(nes) %>%
    mutate(pathway = factor(pathway, levels = pathway)) %>%
    ggplot(aes(x = pathway, y = nes)) +
    geom_col(fill = "gray") +
    geom_col(data = sig_data, fill = "black") +
    coord_flip() +
    theme_classic()
}


#' Saves a ggplot object to a file with options to customize font, size, and layout.
#'
#' @param filename The name of the file to save the plot to.
#' @param plot The ggplot object to save.
#' @param font The font family to use. Default is "Arial".
#' @param font_size The font size. Default is 10.
#' @param line_width The width of the axis lines. Default is 0.5.
#' @param ncol The number of columns for the plot layout. Default is 1.
#' @param nrow The number of rows for the plot layout. Default is 1.
#' @param base_height The base height for the plot. Default is 3.71.
#' @param base_asp The aspect ratio of the plot. Default is 1.618.
#' @param base_width Optional base width for the plot.
#' @param strip_text Whether to remove axis text and titles. Default is FALSE.
#' @return No return value. Saves the plot to the specified file.
#' @export
save_ploty <- function(filename, plot, font = "Arial", font_size = 10, line_width = .5,
                       ncol = 1, nrow = 1, base_height = 3.71, base_asp = 1.618, base_width = NULL,
                       strip_text = FALSE){
  plot <- plot + theme(text = element_text(family = font, size = font_size, color = 'black'),
                       axis.text = element_text(family = font, size = font_size, color = 'black'),
                       legend.text = element_text(family = font, size = font_size, color = 'black'),
                       axis.ticks = element_line(size = line_width, color = "black"),
                       axis.line = element_line(size = line_width, color = "black"))
  if(strip_text){
    plot <- plot + theme(axis.text.x = element_blank(), plot.title = element_blank(),
                         axis.text.y = element_blank(), axis.title.x = element_blank(),
                         axis.title.y = element_blank(), legend.position = "none")
  }
  save_plot(filename = filename, plot = plot, ncol = ncol, nrow = nrow, base_height = base_height,
            base_asp = base_asp, base_width = base_width)
}

#' Capitalize the first letter of each word in a string
#'
#' @param x A string to capitalize.
#' @return A string with the first letter of each word capitalized.
#' @export
capitalize <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep="", collapse=" ")
}


