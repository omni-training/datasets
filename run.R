#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(anndataR)
  library(SingleCellExperiment)
})

parser <- ArgumentParser(description = "Benchmarking entrypoint")

# dummy args (old API)
parser$add_argument("--output_dir", dest = "output_dir")
parser$add_argument("--name", dest = "name")

# parameters
parser$add_argument("--dataset_name", dest = "dataset_name", type = "character",
                    help = "name of the dataset", 
                    choices = c("sc-mix", "be1", "cb", "1.3m"), required = TRUE)
parser$add_argument("--output", action = "append", default = list(),
                    help = "Repeatable: key=path")

args <- parser$parse_args()

output <- setNames(sub("^[^=]*=", "", args$output), 
                   sub("=.*$", "", args$output))

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))

cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
for (i in 1:length(output)) {
  cat(sprintf("[output] %s: %s\n", names(output)[i], output[[i]]))
}
cat(sprintf("----------------------------------\n"))

# set outputs

make_qc_df <- function(
  nFeature_min = NA_real_, nFeature_max = NA_real_,
  nCount_min = NA_real_, nCount_max = NA_real_,
  percent_mt_min = NA_real_, percent_mt_max = NA_real_
) {
  data.frame(
    metric = c("nFeature", "nCount", "percent.mt"),
    min = c(nFeature_min, nCount_min, percent_mt_min),
    max = c(nFeature_max, nCount_max, percent_mt_max),
    stringsAsFactors = FALSE
  )
}

if (args$dataset_name == "sc-mix") {
  # download processed scMixology dataset
  options(timeout=180)
  #url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  url <- "https://omnibenchmark.mls.uzh.ch/datasets/sc-mix/sincell_with_class_5cl.RData"
  #url <- "file:///Users/mark/projects/omb/omni-training/sincell_with_class_5cl.RData" # for local
  bn <- basename(url)
  #raw_path <- file.path(args$output_dir, bn)
  raw_path <- file.path(tempdir(), bn)
  if (!file.exists(raw_path)) {
    download.file(url, destfile = raw_path)
  }

  load(raw_path)
  sce <- sce_sc_10x_5cl_qc
  file.remove(raw_path)

  # use provided cell line annotations as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$cell_line

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 6200,
    nCount_max = 60000,
    percent_mt_max = 10
  )
} else if (args$dataset_name == "be1") {
  library(DropletUtils)
  library(GEOquery)
  library(stringr)
  # download GEO files for be1
  gse_id <- "GSE243665"
  getGEOSuppFiles(
    GEO = gse_id,
    #baseDir = args$output_dir
    baseDir = tempdir()
  )

  # collect compressed matrix/feature/barcode files across samples
  #raw_dir <- file.path(args$output_dir, gse_id)
  raw_dir <- file.path(temp_dir(), gse_id)
  files <- list.files(
    raw_dir,
    pattern = "\\.(mtx|tsv)\\.gz$", full.names = TRUE
  )
  samples <- unique(str_match(basename(files), "_(.*?)_")[, 2])

  for (s in samples) {
    # create 10x-style folders and move each sample's files into place
    sample_dir <- file.path(raw_dir, s)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    s_files <- files[grepl(paste0("_", s, "_"), basename(files))]
    file.copy(s_files, file.path(raw_dir, s), overwrite = TRUE)

    for (f in s_files) {
      bn <- basename(f)
      # strip the GEO/sample prefix so filenames match 10x conventions
      new_name <- sub("^[^_]+_[^_]+_", "", bn)
      dest <- file.path(sample_dir, new_name)
      file.rename(f, dest)
    }
  }

  # read each sample and concatenate them into a single SCE object
  sce_list <- lapply(samples, function(sample) {
    read10xCounts(
      samples = file.path(raw_dir, sample),
      sample.names = sample,
      col.names = TRUE
    )
  })
  sce <- do.call(cbind, sce_list)

  unlink(raw_dir, recursive = TRUE)
  metadata(sce) <- list()
  rownames(sce) <- rowData(sce)$Symbol
  rownames(sce) <- make.unique(rownames(sce))

  # use sample names as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$Sample

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 5000,
    nCount_max = 25000,
    percent_mt_max = 5
  )
} else if (args$dataset_name == "cb") {
  library(SingleCellMultiModal)
  # load Cord blood CITEseq data
  sce <- CITEseq(
    DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
    DataClass = "SingleCellExperiment"
  )

  # keep human genes and drop the "HUMAN_" prefix for consistency
  gene_m <- grep("^HUMAN_", rownames(sce), value = TRUE)
  sce <- sce[gene_m, ]
  rownames(sce) <- sub("^HUMAN_", "", rownames(sce))

  # use provided celltype annotations as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$celltype

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 2500,
    nCount_max = 4000,
    percent_mt_max = 5
  )
} else if (args$dataset_name == "1.3m") {
  # fetch 10x 1.3M data
  library(TENxBrainData)
  sce <- TENxBrainData()
  rownames(sce) <- rowData(sce)$Symbol

  # cluster labels are set to the library index for the full pipeline to run
  # the ARIs computed from this dataset do not have any biological significance
  colData(sce)$clusters.truth <- colData(sce)$Library

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 2500,
    nCount_max = 4000,
    percent_mt_max = 5
  )
}

# filter NA annotated cells
sce <- sce[, !is.na(colData(sce)$clusters.truth)]

# write outputs
# h5ad ----------------------------------------------------------------------
write_h5ad(sce, path = output[["data_h5ad"]])
# ground truth clusters -----------------------------------------------------
write.table( data.frame(cell_id = colnames(sce), truths = sce$clusters.truth),
             output[["data_clusters_truth"]],
             sep = "\t", quote = FALSE, row.names = FALSE)
# number of clusters --------------------------------------------------------
writeLines(as.character(length(unique(sce$clusters.truth))), 
           con = output[["data_clusters_truth_num"]])
