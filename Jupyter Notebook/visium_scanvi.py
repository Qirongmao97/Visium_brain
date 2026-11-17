"""Transfer mouse brain region annotations across Visium slides with scVI/scANVI.

Run from the repository root:
    python visium_scanvi.py --input-dir data/visium_scanvi

Inputs relative to --input-dir:
    reference/{wildtype,mdxbl6,mdx4cv,mdx52}_noCB.h5ad
        Four annotated slide 1 sections with raw counts in X and an
        annotation column in obs. Use the prepared noCB files from the study.
    query/<sample>/outs/
        Space Ranger outputs for wildtype_2, mdxbl6_2, mdx4cv_2, mdx52_2,
        wildtype_3, mdxbl6_3, mdx4cv_3, and wildtype_4. Each outs directory
        contains filtered_feature_bc_matrix.h5 and the spatial directory.

The script merges on shared genes, selects up to 2,000 variable genes by
Slide, fits scVI with Batch as the covariate, then trains scANVI with the
reference annotations. Slide 1 is batch 1 and slides 2–3 are batch 2.

Outputs under --output-dir (default: results/visium_scanvi):
    Visium_merged.h5ad    Integrated data, predictions, and embeddings.
    Visium_scvi/         Trained scVI model.
    sample_metadata/    One CSV per sample, indexed by original barcode.
    frequency_table.csv  Cross-tabulation of scVI clusters and scANVI labels.

Requires anndata, pandas, Scanpy (with read_visium), scvi-tools, scikit-misc
(for Seurat v3 feature selection), igraph, and leidenalg. The training
epoch defaults match the completed runs recorded in the source notebook.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


SLIDES = {
    1: ("wildtype", "mdxbl6", "mdx4cv", "mdx52"),
    2: ("wildtype_2", "mdxbl6_2", "mdx4cv_2", "mdx52_2"),
    3: ("wildtype_3", "mdxbl6_3", "mdx4cv_3", "wildtype_4"),
}
UNLABELED = "Unpredicted"


def load_visium(input_dir: Path) -> AnnData:
    """Read the 12 sections and preserve the input barcode before merging."""
    import anndata as ad
    import pandas as pd
    import scanpy as sc

    samples = []
    for slide, names in SLIDES.items():
        for name in names:
            if slide == 1:
                path = input_dir / "reference" / f"{name}_noCB.h5ad"
                if not path.is_file():
                    raise FileNotFoundError(f"Missing reference H5AD: {path}")
                sample = sc.read_h5ad(path)
                if "annotation" not in sample.obs:
                    raise ValueError(f"{name} is missing obs['annotation'].")
            else:
                path = input_dir / "query" / name / "outs"
                if not (path / "filtered_feature_bc_matrix.h5").is_file():
                    raise FileNotFoundError(f"Missing Space Ranger count file under {path}")
                sample = sc.read_visium(path)

            if sample.n_obs == 0 or not sample.obs_names.is_unique:
                raise ValueError(f"{name} must have nonempty, unique spot barcodes.")
            sample.var_names_make_unique()
            # Keep the QC columns exported by the notebook, without adding filters.
            percent_top = [n for n in (50, 100, 200, 500) if n <= sample.n_vars]
            sc.pp.calculate_qc_metrics(sample, percent_top=percent_top, inplace=True)
            sample.obs["Barcode"] = sample.obs_names.to_numpy()
            sample.obs["Sample"] = name
            sample.obs["Slide"] = slide
            sample.obs["Batch"] = 1 if slide == 1 else 2
            samples.append(sample)

    # A sample index is appended to the merged IDs, as in the original
    # concatenation. The same raw barcode can occur on multiple slides.
    adata = ad.concat(samples, join="inner", label="batch", index_unique="-")
    if adata.n_vars == 0:
        raise ValueError("The 12 samples have no shared gene names.")
    if not adata.obs_names.is_unique:
        raise ValueError("The merged spot IDs are not unique.")
    # Intersect genes, not metadata columns. An inner AnnData join can drop
    # reference-only columns such as annotation. Retain the original union
    # of spot metadata, with missing annotations for the query sections.
    spot_metadata = pd.concat([sample.obs for sample in samples], ignore_index=True)
    spot_metadata.index = adata.obs_names
    spot_metadata["batch"] = adata.obs["batch"]
    adata.obs = spot_metadata
    adata.obs["Slide"] = adata.obs["Slide"].astype("category")
    return adata


def fit_models(adata: AnnData, scvi_epochs: int, scanvi_epochs: int, seed: int):
    """Prepare raw counts, fit both models, and store spot-level predictions."""
    import pandas as pd
    import scanpy as sc
    import scvi

    scvi.settings.seed = seed
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=min(2000, adata.n_vars),
        layer="counts",
        batch_key="Slide",
        subset=True,
    )

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Batch")
    scvi_model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    scvi_model.train(max_epochs=scvi_epochs)
    adata.obsm["X_scVI"] = scvi_model.get_latent_representation()

    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata, min_dist=0.3)
    sc.tl.leiden(adata, key_added="leiden_scVI", flavor="leidenalg")

    labels = adata.obs["annotation"].astype("category")
    if UNLABELED not in labels.cat.categories:
        labels = labels.cat.add_categories([UNLABELED])
    labels = labels.fillna(UNLABELED)
    if not labels.ne(UNLABELED).any() or not labels.eq(UNLABELED).any():
        raise ValueError("scANVI transfer requires both annotated and unannotated spots.")
    adata.obs["annotation"] = labels

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata,
        unlabeled_category=UNLABELED,
        labels_key="annotation",
    )
    scanvi_model.train(max_epochs=scanvi_epochs)
    adata.obsm["X_scANVI"] = scanvi_model.get_latent_representation(adata)
    adata.obs["C_scANVI"] = pd.Categorical(
        scanvi_model.predict(adata), categories=labels.cat.categories
    )
    return scvi_model


def write_outputs(adata: AnnData, scvi_model, output_dir: Path) -> None:
    """Save the integrated data, scVI checkpoint, and annotation tables."""
    import pandas as pd

    metadata_dir = output_dir / "sample_metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    for name, sample_obs in adata.obs.groupby("Sample", sort=False, observed=True):
        # Use the original barcode as the join key for the individual sample.
        # Preserve the unique merged spot ID in a separate column.
        sample_obs = sample_obs.assign(spot_id=sample_obs.index).set_index("Barcode")
        sample_obs.to_csv(metadata_dir / f"{name}.csv")

    pd.crosstab(adata.obs["leiden_scVI"], adata.obs["C_scANVI"]).to_csv(
        output_dir / "frequency_table.csv"
    )
    scvi_model.save(str(output_dir / "Visium_scvi"), overwrite=True)
    # The notebook retained array-grid columns in the CSVs, not in the H5AD.
    if "in_tissue" in adata.obs:
        adata.obs["in_tissue"] = adata.obs["in_tissue"].astype(str)
    adata.obs = adata.obs.drop(columns=["array_row", "array_col"], errors="ignore")
    adata.write_h5ad(output_dir / "Visium_merged.h5ad")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input-dir", type=Path, default=Path("data/visium_scanvi"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/visium_scanvi"))
    parser.add_argument("--scvi-epochs", type=int, default=191)
    parser.add_argument("--scanvi-epochs", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    if min(args.scvi_epochs, args.scanvi_epochs) < 1:
        parser.error("Training epoch counts must be positive.")

    adata = load_visium(args.input_dir)
    scvi_model = fit_models(adata, args.scvi_epochs, args.scanvi_epochs, args.seed)
    write_outputs(adata, scvi_model, args.output_dir)
    print(f"Saved Visium label-transfer results to {args.output_dir}")


if __name__ == "__main__":
    main()
