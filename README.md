# Mouse brain Visium analysis in Duchenne muscular dystrophy

This repository accompanies the study *Spatial transcriptomics maps region and isoform specific brain alterations in mdx and mdx4cv models of Duchenne muscular dystrophy*. The study used 10x Visium to compare regional gene expression in the developing brains of wild-type, mdx, and mdx4cv mice.

## Where to start

| Folder | What it contains |
| --- | --- |
| [R/](R/README.md) | Seurat preprocessing and spot QC, summaries of regional differential expression, pathway and transcription-factor analyses, and figure generation. |
| [notebooks/](notebooks/README.md) | Python workflow for scVI integration and scANVI transfer of brain region annotations from reference sections to unannotated sections. |

- **To regenerate figures from existing results**, start with [the R script index](R/README.md#regional-mast-results-and-figures). Most of these scripts consume precomputed MAST tables and can be run independently of the Python workflow.
- **To run spot QC on annotated Seurat files**, follow [the preprocessing and QC steps](R/README.md#seurat-preprocessing-and-spot-qc).
- **To rerun label transfer**, follow [the Python instructions](notebooks/README.md). The per-sample CSVs contain the predicted label, sample ID, and original Visium barcode. When adding these annotations to individual Seurat objects, match by sample and barcode.

Run the examples from the repository root. By default, the scripts look for inputs under `data/` and write generated files under `results/`. The folder READMEs describe how to override these locations.

## Data and scope

Visium study data are indexed in the NCBI Gene Expression Omnibus under [GSE349849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349849). Consult the study's Data Availability statement and the GEO record for acquisition details.

The scripts expect specific prepared inputs, including annotated reference H5AD files, per-sample Seurat RDS files, and precomputed differential expression results. The GEO accession is a source for the study data, not a guarantee that it contains those exact intermediate files. In particular, the figure scripts that read `ALL_MAST_results.RData` do not fit the MAST models, and the Python label-transfer script does not create Seurat RDS files.

The integration and QC workflows include mdx52 sections. The main genotype comparisons described in the manuscript focus on wild-type, mdx (labeled `mdxbl6` in the code), and mdx4cv.

## Citation

Mao Q, Heezen LGM, de Vries S, Verwoerd N, Aartsma-Rus A, van Putten M, Spitali P, Mahfouz A. *Spatial transcriptomics maps region and isoform specific brain alterations in mdx and mdx4cv models of Duchenne muscular dystrophy*. **Neurobiology of Disease** (2026), 107677. [doi:10.1016/j.nbd.2026.107677](https://doi.org/10.1016/j.nbd.2026.107677).
