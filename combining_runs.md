# Combining sequencing runs of the same sample

This pipeline was written for one library per sample. When the same DNA extract
is prepared and sequenced more than once (for example two library preps six
months apart), the runs should be **denoised separately and combined on the ASV
table**, not concatenated at the FASTQ level.

This document describes the supported route through the pipeline, and how to
merge outputs from two runs that were processed separately.

## Recommended design

1. Keep one sample ID per **library** (`sample1_PB334`, `sample1_PB479`). This
   preserves the ability to inspect each run and to detect a batch effect.
2. Denoise each sequencing run separately by setting the existing `pool` column
   in the metadata to the run name.
3. Add a column that identifies the **biological sample** (for example
   `subject`).
4. Run the pipeline with `--group_samples_by subject`. After denoising and ASV
   filtering, the merged table is summed per subject and the combined
   diversity, taxonomy and barplot outputs are written to `results/combined`.

Separate denoising is recommended because DADA2 learns an error model per
sequencing run. The combined table sums read counts for the same sample; ASVs
are identified by `md5(sequence)` in QIIME 2 2024.10, so an identical ASV found
in both runs collapses to one feature instead of being counted twice.

The libraries remain visible in the main report and in `results/dada2`, which
is useful to check whether one run contributed disproportionately many reads or
ASVs.

## Input files

`samples.tsv` maps every library to its FASTQ:

```tsv
sample-id	absolute-filepath
sample1_PB334	/data/PB334/sample1.fastq.gz
sample1_PB479	/data/PB479/sample1.fastq.gz
sample2_PB334	/data/PB334/sample2.fastq.gz
sample2_PB479	/data/PB479/sample2.fastq.gz
```

`metadata.tsv` has one row per library. `sample_name` must match `sample-id`,
`pool` groups libraries that are denoised together (here: one pool per run),
and `subject` is the biological sample that libraries are combined into:

```tsv
sample_name	condition	subject	pool
sample1_PB334	groupA	sample1	PB334
sample1_PB479	groupA	sample1	PB479
sample2_PB334	groupA	sample2	PB334
sample2_PB479	groupA	sample2	PB479
```

Columns that differ between libraries of the same subject (run name, run date,
library concentration, ...) are dropped from the combined metadata
automatically and listed in
`results/combined/combining_summary.txt`. Columns that are constant within a
subject (`condition`, ...) are kept.

## Running the pipeline

```bash
nextflow run main.nf \
    --input samples.tsv \
    --metadata metadata.tsv \
    --group_samples_by subject \
    -profile singularity
```

For runs that were already processed separately, running the pipeline once on
all libraries (with `-resume` if the work directory is still available) is
preferable, because ASV filtering is then applied once to the combined cohort.

## Outputs

The normal per-library outputs are unchanged. In addition:

- `results/combined/`
  - `dada2-ccs_table_combined.qza`: ASV counts summed per `subject`.
  - `metadata_combined.tsv`: metadata matching the combined table.
  - `vsearch_merged_freq_tax_combined.tsv`: ASV sequence, taxonomy and combined
    counts per sample.
  - `feature-table-tax_vsearch_combined.biom`: combined table with taxonomy,
    suitable for `phyloseq`/`qiime2R`.
  - `physeq_combined.rds`: `phyloseq` object with the combined ASV table,
    VSEARCH taxonomy, rooted tree and combined metadata. Load it in R with
    `readRDS()`.
  - `merged_freq_tax_combined.qzv`: interactive table.
  - `taxonomy_barplot_vsearch_combined.qzv`: taxonomy barplot for the combined
    samples.
  - `combining_summary.txt`: metadata columns dropped because they differed
    between libraries of the same subject.
  - `phylogeny_diversity/`: core diversity metrics and
    `alpha-rarefaction-curves.qzv` computed on the combined table.
    The phylogenetic tree is shared with the per-library analysis.

The HTML report (`results/visualize_biom.html`) intentionally remains
library-level. Use it to check per-run yield and ASV recovery, then use
`results/combined` for biological interpretation.

A `phyloseq` object for the per-library table is also written to
`results/physeq_vsearch.rds`, using the same VSEARCH taxonomy and rooted tree.

## Caveats

- **ASV filtering happens before combining.** `--min_asv_sample` counts
  libraries, not subjects. An ASV seen in both libraries of one subject counts
  as present in two samples. If you want "present in at least N subjects"
  semantics, set `--min_asv_sample 0` and filter the combined
  `dada2-ccs_table_combined.qza` afterwards.
- **Rarefaction depth is computed from the per-library table.** It is
  conservative for the combined table because summed counts are only larger.
  Use `--rarefaction_depth` if you want a fixed value for both analyses.
- **Do not combine raw FASTQ files.** Concatenating reads from different
  library preps forces DADA2 to fit one error model to mixed quality profiles
  and hides run-specific failure modes.
- **Same sample ID in multiple pools.** The pipeline's sample IDs are library
  IDs, so this should not happen. If you do reuse a sample ID across pools,
  `mergeASV` needs `--merge_overlap_method sum`; the default
  `error_on_overlapping_sample` stops the run instead of silently summing.

## Merging two already-finished runs manually

If rerunning is not possible, merge the **unfiltered** DADA2 outputs from each
run (`dada2-ccs_table_merged.qza` and `dada2-ccs_rep_merged.qza`), then apply
the ASV filter once to the merged table. Merging tables that were already
filtered with `--min_asv_sample`/`--min_asv_totalfreq` changes the filter
semantics and is not equivalent to a joint run.

```bash
# 1. Merge ASV tables (use --p-overlap-method sum if sample IDs overlap)
qiime feature-table merge \
    --i-tables run_PB334/dada2/dada2-ccs_table_merged.qza \
    --i-tables run_PB479/dada2/dada2-ccs_table_merged.qza \
    --p-overlap-method error_on_overlapping_sample \
    --o-merged-table merged/dada2-ccs_table_merged.qza

# 2. Merge representative sequences (identical ASVs collapse to one feature)
qiime feature-table merge-seqs \
    --i-data run_PB334/dada2/dada2-ccs_rep_merged.qza \
    --i-data run_PB479/dada2/dada2-ccs_rep_merged.qza \
    --o-merged-data merged/dada2-ccs_rep_merged.qza

# 3. Filter ASVs once, on the combined cohort
qiime feature-table filter-features \
    --i-table merged/dada2-ccs_table_merged.qza \
    --p-min-frequency 5 \
    --p-min-samples 1 \
    --p-filter-empty-samples \
    --o-filtered-table merged/dada2-ccs_table_filtered.qza

qiime feature-table filter-seqs \
    --i-data merged/dada2-ccs_rep_merged.qza \
    --i-table merged/dada2-ccs_table_filtered.qza \
    --o-filtered-data merged/dada2-ccs_rep_filtered.qza
```

Then sum libraries per subject and build the matching metadata:

```bash
qiime feature-table group \
    --i-table merged/dada2-ccs_table_filtered.qza \
    --p-axis sample \
    --p-mode sum \
    --m-metadata-file metadata.tsv \
    --m-metadata-column subject \
    --o-grouped-table merged/dada2-ccs_table_combined.qza

python3 scripts/combine_metadata.py \
    --input metadata.tsv \
    --output merged/metadata_combined.tsv \
    --column subject
```

`--p-overlap-method` only matters when a sample ID occurs in both inputs. With
distinct library IDs the default (`error_on_overlapping_sample`) is correct and
will flag accidental duplicates.

## Validation checklist

Before interpreting the combined results:

- Sum of reads in `dada2-ccs_table_combined.qza` equals the sum of the
  per-library filtered table; no reads should be gained or lost by grouping.
- ASV count in the combined table is at most the sum of the per-run ASV counts;
  it is lower when the same sequence was recovered in both runs.
- The number of samples in `metadata_combined.tsv` equals the number of unique
  `subject` values.
- Check `combining_summary.txt` to confirm that run-specific columns were
  dropped and that nothing biologically meaningful was lost.
- Inspect the per-library report for a run effect (reads, ASV richness) before
  treating the combined sample as one observation.
