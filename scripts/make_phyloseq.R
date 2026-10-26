#!/usr/bin/env Rscript

# Build a phyloseq object from QIIME 2 artifacts that have been exported with
# `qiime tools export`:
#   - feature table  -> feature-table.biom
#   - taxonomy       -> taxonomy.tsv
#   - rooted tree    -> tree.nwk
#
# Usage:
#   Rscript make_phyloseq.R --table feature-table.biom --taxonomy taxonomy.tsv \
#       --tree tree.nwk --metadata metadata.tsv --output physeq.rds

suppressPackageStartupMessages({
    library(phyloseq)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) {
        if (required) {
            stop("missing argument: ", flag, call. = FALSE)
        }
        return(NULL)
    }
    args[[i + 1L]]
}

table_file <- get_arg("--table")
taxonomy_file <- get_arg("--taxonomy")
tree_file <- get_arg("--tree", required = FALSE)
metadata_file <- get_arg("--metadata")
output_file <- get_arg("--output")

# `phyloseq::import_biom` fails for tables with a single sample or feature.
# Fall back to biomformat in that case.
read_otu <- function(path) {
    otu <- tryCatch(
        suppressWarnings(phyloseq::import_biom(path)),
        error = function(e) NULL
    )
    if (!is.null(otu)) {
        return(phyloseq::otu_table(otu))
    }

    biom <- suppressWarnings(biomformat::read_biom(path))
    m <- as.matrix(biomformat::biom_data(biom))
    if (is.null(rownames(m)) || is.null(colnames(m))) {
        stop("could not read feature or sample names from ", path, call. = FALSE)
    }
    phyloseq::otu_table(m, taxa_are_rows = TRUE)
}

read_taxonomy <- function(path) {
    tax <- read.delim(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )
    if (ncol(tax) < 2) {
        stop("taxonomy file has fewer than two columns: ", path, call. = FALSE)
    }

    id_column <- if ("Feature ID" %in% names(tax)) "Feature ID" else names(tax)[[1]]
    tax_column <- if ("Taxon" %in% names(tax)) "Taxon" else names(tax)[[2]]

    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    parts <- strsplit(as.character(tax[[tax_column]]), ";", fixed = TRUE)
    matrix_rows <- vapply(parts, function(x) {
        x <- sub("^[a-z]__", "", x, perl = TRUE)
        x <- trimws(gsub("_", " ", x))
        x <- x[x != "" & !is.na(x)]
        length(x) <- length(ranks)
        x
    }, character(length(ranks)))

    mat <- t(matrix_rows)
    rownames(mat) <- as.character(tax[[id_column]])
    colnames(mat) <- ranks
    mat
}

read_metadata <- function(path) {
    meta <- read.delim(
        path,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE,
        na.strings = c("", "NA")
    )
    if (ncol(meta) < 1) {
        stop("metadata file has no columns: ", path, call. = FALSE)
    }

    ids <- as.character(meta[[1]])
    keep <- !is.na(ids) & ids != "" & ids != "#q2:types"
    meta <- meta[keep, , drop = FALSE]
    ids <- ids[keep]

    if (anyDuplicated(ids)) {
        stop(
            "duplicate sample IDs in metadata: ",
            paste(unique(ids[duplicated(ids)]), collapse = ", "),
            call. = FALSE
        )
    }

    rownames(meta) <- ids
    meta[[1]] <- NULL
    meta
}

otu <- read_otu(table_file)
otu_matrix <- as(phyloseq::otu_table(otu), "matrix")
if (!phyloseq::taxa_are_rows(otu)) {
    otu_matrix <- t(otu_matrix)
}

taxonomy <- read_taxonomy(taxonomy_file)
metadata <- read_metadata(metadata_file)

tree <- NULL
if (!is.null(tree_file) && file.exists(tree_file)) {
    tree <- tryCatch(
        phyloseq::read_tree(tree_file),
        error = function(e) NULL
    )
    if (!is.null(tree) && length(tree$tip.label) < 2) {
        tree <- NULL
    }
}

# Keep the features and samples that are present in every component. QIIME 2
# exports from the same analysis should already agree; this guards against
# partially exported or manually edited artifacts.
features <- rownames(otu_matrix)
if (!is.null(taxonomy)) {
    features <- intersect(features, rownames(taxonomy))
}
if (!is.null(tree)) {
    missing_from_tree <- setdiff(features, tree$tip.label)
    if (length(missing_from_tree) > 0) {
        warning(
            length(missing_from_tree),
            " features are missing from the tree and will be dropped"
        )
    }
    features <- intersect(features, tree$tip.label)
}
if (length(features) == 0) {
    stop("no features are shared by the table, taxonomy and tree", call. = FALSE)
}

samples <- colnames(otu_matrix)
missing_metadata <- setdiff(samples, rownames(metadata))
if (length(missing_metadata) > 0) {
    warning(
        length(missing_metadata),
        " samples are missing from the metadata and will be dropped"
    )
}
samples <- intersect(samples, rownames(metadata))
if (length(samples) == 0) {
    stop("no samples are shared by the table and metadata", call. = FALSE)
}

otu_matrix <- otu_matrix[features, samples, drop = FALSE]
taxonomy <- taxonomy[features, , drop = FALSE]
metadata <- metadata[samples, , drop = FALSE]

if (!is.null(tree)) {
    extra_tips <- setdiff(tree$tip.label, features)
    if (length(extra_tips) > 0) {
        tree <- ape::drop.tip(tree, extra_tips)
    }
    if (length(tree$tip.label) < 2) {
        tree <- NULL
    }
}

components <- list(
    phyloseq::otu_table(otu_matrix, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy),
    phyloseq::sample_data(metadata)
)
if (!is.null(tree)) {
    components <- c(components, list(phyloseq::phy_tree(tree)))
}

physeq <- do.call(phyloseq::phyloseq, components)
saveRDS(physeq, output_file)

cat(sprintf(
    "Wrote %s: %d taxa x %d samples\n",
    output_file, phyloseq::ntaxa(physeq), phyloseq::nsamples(physeq)
))
