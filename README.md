# scATAC

This repository contains general scripts for processing and analysis of single-cell ATAC-seq data.

## cellranger-atac
Running cellranger-atac is the first step in scATAC-seq analysis. This requires **fastq files**, with a specific naming convention ("{sample_id}_S1_{laneid}_{read_type}_001.fastq.gz"). 

If fastq files are not already named accordingly, the _Prepare_samples_for_cell_ranger.R_ script can help to modify them.

To run cellranger, create a short **bash script** containing something like the following, adjusting the paths to where your fastq files are located, and changing the sample ID accordingly. This can be adjusted to loop through several samples.
<pre> cellranger-atac count \
     --id=yoursamplename \
     --fastqs='/path/to/your/fastq/files' \
     --sample=yoursamplename \
     --reference=/path/to/your/reference/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/
</pre>

The cellranger **reference genome** can be downloaded from the 10X genomics website: https://www.10xgenomics.com/support/cn/software/cell-ranger-arc/latest/release-notes/reference-release-notes.

Running cellranger usually takes a couple of days for each sample. In the end the important outputs are the **fragments.tsv.gz** files, which are used as the first input to ArchR for downstream analysis. Before moving on you should also inspect the **web_summary.html**. This will give you an overview of the data content and indicate any quality concerns.

## ArchR

The main steps in processing and analysis of scATAC-seq data in ArchR are listed below and can be tested using the _ArchR_example.R_ script.

* loading the fragments files to create so-called **arrow files** and build an **ArchR project**
* detecting and removing **low quality cells** and **doublets**
* **dimensionality reduction** using iterative LSI
* **clustering**
* **peak calling** (requires installing/loading masc2 software)

After peak calling we obtain a **"peaks by cells"** matrix which can be used for downstream analyses. This is a sparse matrix with cells represented by columns and peaks/regions represented by rows, where 1's represent **open/accessible chromatin** (typically active promoters/enhancers) and 0's represent **closed/innaccessible chromatin** (or dropouts/missing vales).

In addition to the peaks matrix, ArchR also creates **"genes by cells"** and **"tiles by cells"** matrices.

For more details about these steps and further downstream analyses (e.g. identifying marker genes/peaks, motif and feature enrichment), the ArchR developers provide a very comprehensive **ArchR manual**, from which the ArchR example script here was adapted: https://www.archrproject.com/bookdown/index.html
