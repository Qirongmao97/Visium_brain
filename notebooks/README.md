# Visium annotation transfer in Python

[visium_scanvi.py](visium_scanvi.py) is the standalone version of the Visium scANVI notebook. It runs the label-transfer workflow without a Jupyter kernel, exploratory visualizations, or machine-specific paths. See the [repository overview](../README.md) for the publication and dataset.

The workflow uses four annotated sections from slide 1 as a reference and transfers labels to eight unannotated sections from slides 2 and 3. It intersects genes across the 12 sections, selects up to 2,000 highly variable genes by slide, integrates the sections with scVI using processing batch, and uses scANVI to predict brain region labels. Slide 1 is batch 1. Slides 2 and 3 are batch 2.

## Python environment

Use an isolated Python environment with `anndata`, `pandas`, `scanpy`, `scvi-tools`, `scikit-misc`, `igraph`, and `leidenalg`. The script uses `scanpy.read_visium` to read Space Ranger outputs, so use a Scanpy version that provides it. Jupyter is not required.

~~~bash
python -m pip install anndata pandas "scanpy<2" scvi-tools scikit-misc igraph leidenalg
~~~

## Inputs

The default input directory is `data/visium_scanvi/` relative to the repository root. Override it with `--input-dir` if the files are elsewhere.

| Slide | Expected input paths relative to the input directory |
| --- | --- |
| 1 (reference) | `reference/wildtype_noCB.h5ad`, `reference/mdxbl6_noCB.h5ad`, `reference/mdx4cv_noCB.h5ad`, and `reference/mdx52_noCB.h5ad` |
| 2 (query) | `query/<sample>/outs/` for `wildtype_2`, `mdxbl6_2`, `mdx4cv_2`, and `mdx52_2` |
| 3 (query) | `query/<sample>/outs/` for `wildtype_3`, `mdxbl6_3`, `mdx4cv_3`, and `wildtype_4` |

Each reference H5AD must contain raw counts in `X` and the reference region labels in `obs["annotation"]`. Use the prepared `noCB` reference files. The script does not apply an additional reference barcode exclusion list. Each query `outs/` directory must contain `filtered_feature_bc_matrix.h5` and the `spatial/` directory.

## Run

From the repository root:

~~~bash
python notebooks/visium_scanvi.py \
  --input-dir data/visium_scanvi \
  --output-dir results/visium_scanvi
~~~

Use `--help` to list the options. The script defaults to 191 scVI epochs and 10 scANVI epochs, matching the completed training runs in the source notebook. To change these settings, use `--scvi-epochs` and `--scanvi-epochs`. The script sets a random seed, but different package versions or hardware can still change model results.

## Outputs

| File or directory | Contents |
| --- | --- |
| `Visium_merged.h5ad` | Merged dataset with unique spot IDs, `X_scVI` and `X_scANVI` embeddings, `leiden_scVI` clusters, and `C_scANVI` predicted region labels. |
| `Visium_scvi/` | Trained scVI model checkpoint. |
| `sample_metadata/<sample>.csv` | Spot metadata and predictions for one section. `Barcode` holds the original Visium barcode and `spot_id` identifies the same spot in the merged H5AD. |
| `frequency_table.csv` | Counts of predicted scANVI labels within scVI Leiden clusters. |

For downstream work in R, use the per-sample CSV and match `C_scANVI` to Seurat spots using the sample and original `Barcode`. Do not match by barcode across all samples without also checking sample identity. Continue with the [R preprocessing and QC instructions](../R/README.md#seurat-preprocessing-and-spot-qc) once the annotated Seurat RDS files are prepared.
