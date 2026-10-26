/*
===============================================================================
Processes used to combine libraries that belong to the same biological sample
(e.g. the same DNA extract prepared and sequenced in two separate runs).

Libraries are still denoised separately (one `pool` per sequencing run) and
merged by `mergeASV`. When `--group_samples_by` is set, the merged feature
table is then summed per biological sample and downstream diversity/taxonomy
outputs are produced for the combined samples in `results/combined`.
===============================================================================
*/

process combine_runs {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/combined", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path metadata
    path metadata_script
    val group_column

    output:
    path "dada2-ccs_table_combined.qza", emit: grouped_table
    path "metadata_combined.tsv", emit: grouped_metadata
    path "combining_summary.txt"

    script:
    """
    qiime feature-table group \
        --i-table $asv_freq \
        --p-axis sample \
        --p-mode sum \
        --m-metadata-file $metadata \
        --m-metadata-column '$group_column' \
        --o-grouped-table dada2-ccs_table_combined.qza

    python3 $metadata_script \
        --input $metadata \
        --output metadata_combined.tsv \
        --column '$group_column' 2> combining_summary.txt
    """
}

process export_combined {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/combined", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path grouped_table
    path asv_seq
    path tax_vsearch

    output:
    path "vsearch_merged_freq_tax_combined.tsv", emit: freq_tax_tsv
    path "feature-table-tax_vsearch_combined.biom", emit: biom_vsearch
    path "merged_freq_tax_combined.qzv"

    script:
    """
    qiime feature-table transpose --i-table $grouped_table \
        --o-transposed-feature-table transposed-asv.qza

    qiime metadata tabulate --m-input-file $asv_seq \
        --m-input-file $tax_vsearch \
        --m-input-file transposed-asv.qza \
        --o-visualization merged_freq_tax_combined.qzv

    qiime tools export --input-path merged_freq_tax_combined.qzv \
        --output-path merged_freq_tax_tsv
    mv merged_freq_tax_tsv/metadata.tsv vsearch_merged_freq_tax_combined.tsv

    qiime tools export --input-path $tax_vsearch --output-path tax_export
    qiime tools export --input-path $grouped_table --output-path asv_freq/
    sed 's/Feature ID/#OTUID/' tax_export/taxonomy.tsv | sed 's/Taxon/taxonomy/' | \
        sed 's/Consensus/confidence/' > biom-taxonomy.tsv

    biom add-metadata -i asv_freq/feature-table.biom \
        -o feature-table-tax_vsearch_combined.biom \
        --observation-metadata-fp biom-taxonomy.tsv \
        --sc-separated taxonomy
    """
}

process barplot_combined {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/combined", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_tab
    path tax
    path metadata
    val(oname)

    output:
    path "${oname}"

    script:
    """
    qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
        --m-metadata-file $metadata \
        --o-visualization $oname
    """
}

process qiime2_diversity_combined {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results/combined/phylogeny_diversity", mode: params.publish_dir_mode
    label 'cpu8'

    input:
    path grouped_table
    path grouped_metadata
    path rooted_tree
    val(rarefaction_depth)
    val(alpha_depth)

    output:
    path "core-metrics-diversity"
    path "alpha-rarefaction-curves.qzv"
    path "core-metrics-diversity/bray_curtis_distance_matrix.tsv", emit: bray_mat
    path "core-metrics-diversity/weighted_unifrac_distance_matrix.tsv", emit: wunifrac_mat
    path "core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv", emit: unifrac_mat

    script:
    """
    n_group=\$(awk 'NR>1 && \$1 != "#q2:types"' $grouped_metadata | wc -l)
    if [ \${n_group} -gt 1 ];
    then
        qiime diversity core-metrics-phylogenetic \
            --p-n-jobs-or-threads $task.cpus \
            --i-phylogeny $rooted_tree \
            --i-table $grouped_table \
            --m-metadata-file $grouped_metadata \
            --p-sampling-depth $rarefaction_depth \
            --output-dir ./core-metrics-diversity

        qiime tools export --input-path ./core-metrics-diversity/bray_curtis_distance_matrix.qza \
            --output-path ./core-metrics-diversity
        mv ./core-metrics-diversity/distance-matrix.tsv \
            ./core-metrics-diversity/bray_curtis_distance_matrix.tsv
        qiime tools export --input-path ./core-metrics-diversity/weighted_unifrac_distance_matrix.qza \
            --output-path ./core-metrics-diversity
        mv ./core-metrics-diversity/distance-matrix.tsv \
            ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv
        qiime tools export --input-path ./core-metrics-diversity/unweighted_unifrac_distance_matrix.qza \
            --output-path ./core-metrics-diversity
        mv ./core-metrics-diversity/distance-matrix.tsv \
            ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
    else
        # A single combined sample has no meaningful distance matrix
        mkdir ./core-metrics-diversity
        touch ./core-metrics-diversity/bray_curtis_distance_matrix.tsv \
            ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv \
            ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
    fi

    qiime diversity alpha-rarefaction --i-table $grouped_table \
        --m-metadata-file $grouped_metadata \
        --o-visualization alpha-rarefaction-curves.qzv \
        --p-min-depth 10 --p-max-depth $alpha_depth
    """
}
