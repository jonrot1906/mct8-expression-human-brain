# Multiscale characterization and distribution of the thyroid hormone transporter MCT8 in the adult human brain

This is the data repository for the study "Multiscale characterization and distribution of the thyroid hormone transporter MCT8 in the adult and aged human brain".

### Authors:

Jonas Rotter (1), Toni Kirmann (2), Jana Brendler (1), Larissa Anthofer (1,3), Nina Maria Wilpert (4,5,6), Stefan Hallermann (2), Heike Biebermann (3), Ingo Bechmann (1,7)

1 Institute of Anatomy, Leipzig University, Leipzig, Germany.
2 Carl-Ludwig-Institute of Physiology, Leipzig University, Leipzig, Germany.
3 Charité Universitätsmedizin Berlin, corporate member of Freie Universität Berlin and Humboldt-Universität zu Berlin, Germany, Institute for Experimental Paediatric Endocrinology.
4 Department of Paediatric Neurology, Charité-Universitätsmedizin Berlin, Corporate member of Freie Universität Berlin, Humboldt-Universität zu Berlin, and Berlin Institute of Health (BIH), Berlin, Germany. 
5 NeuroCure Cluster of Excellence, Charité-Universitätsmedizin Berlin, Corporate member of Freie Universität Berlin, Humboldt-Universität zu Berlin, and Berlin Institute of Health (BIH), Berlin, Germany. 
6 Berlin Institute of Health at Charité - Universitätsmedizin Berlin, BIH Biomedical Innovation Academy, BIH Charité Junior Clinician Scientist Program, Berlin, Germany. 
7 LeiCeM – Leipzig Center of Metabolism, Leipzig University, Leipzig, Germany.

### Contents

```
Automatic MCT8 detection/
  MCT8_segmentation_with_own_lt.groovy    QuPath: StarDist nucleus detection + MCT8 intensity thresholding
  heatmap_processing.ijm                  Fiji/ImageJ: batch heatmap processing of an image directory tree
SLC16A2 expression in single-nucleus RNA sequencing data/
  slc16a2_umap.ipynb                      former UMAP and per-cell-type SLC16A2 positivity (initial submission)
  slc16a2_heatmap.py                      former cell type x region positivity heatmap (initial submission)
  revision/slc16a2_expression_fig7.py     main-text figure and robustness checks (ABC Atlas)
```

### Running the code

**Transcriptomics** (Python 3.9+):

```
pip install scanpy pandas numpy seaborn matplotlib abc_atlas_access
```

* `slc16a2_umap.ipynb` and `slc16a2_heatmap.py` read CELLxGENE `.h5ad` files from a `data/` folder next to the script; the dataset IDs are hard-coded at the top of each file.
* `revision/slc16a2_expression_fig7.py`: set `DOWNLOAD_BASE` to a local ABC Atlas cache directory, then run `python slc16a2_expression_fig7.py`. The pinned atlas release is fetched on first run (large download; the cache is reused afterwards). Figures and source-data CSVs are written to `revision/main_figure/`, `revision/robustness_threshold/` and `revision/robustness_pseudobulk/`.

Throughout, a nucleus counts as SLC16A2-positive at a raw count > 1, and groups with fewer than 20 nuclei are masked in the figures.

**Image analysis:**

* `MCT8_segmentation_with_own_lt.groovy`: run from the QuPath script editor with the StarDist extension installed; the `dsb2018_heavy_augment.pb` model must be available to the project.
* `heatmap_processing.ijm`: run in Fiji/ImageJ; it prompts for a source directory (containing per-image subdirectories) and a destination directory.

### License

MIT (see `LICENSE`).