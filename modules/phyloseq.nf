/*
===============================================================================
Export QIIME 2 artifacts into a phyloseq object for downstream analysis in R.
===============================================================================
*/

process make_phyloseq {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path tax_vsearch
    path rooted_tree
    path metadata
    path phyloseq_script

    output:
    path "physeq_vsearch.rds", emit: phyloseq_rds

    script:
    """
    mkdir -p phyloseq_inputs
    qiime tools export --input-path $asv_freq --output-path phyloseq_inputs/table
    qiime tools export --input-path $tax_vsearch --output-path phyloseq_inputs/taxonomy
    qiime tools export --input-path $rooted_tree --output-path phyloseq_inputs/tree

    Rscript --vanilla $phyloseq_script \
        --table phyloseq_inputs/table/feature-table.biom \
        --taxonomy phyloseq_inputs/taxonomy/taxonomy.tsv \
        --tree phyloseq_inputs/tree/tree.nwk \
        --metadata $metadata \
        --output physeq_vsearch.rds
    """
}

process make_phyloseq_combined {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/combined", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path tax_vsearch
    path rooted_tree
    path metadata
    path phyloseq_script

    output:
    path "physeq_combined.rds", emit: phyloseq_rds

    script:
    """
    mkdir -p phyloseq_inputs
    qiime tools export --input-path $asv_freq --output-path phyloseq_inputs/table
    qiime tools export --input-path $tax_vsearch --output-path phyloseq_inputs/taxonomy
    qiime tools export --input-path $rooted_tree --output-path phyloseq_inputs/tree

    Rscript --vanilla $phyloseq_script \
        --table phyloseq_inputs/table/feature-table.biom \
        --taxonomy phyloseq_inputs/taxonomy/taxonomy.tsv \
        --tree phyloseq_inputs/tree/tree.nwk \
        --metadata $metadata \
        --output physeq_combined.rds
    """
}
