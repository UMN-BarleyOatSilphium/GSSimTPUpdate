
<!-- README.md is generated from README.Rmd. Please edit that file -->
GSSimTPUpdate
=============

### Project Title:

Evaluating Methods of Updating Training Data in Long-Term Genomewide Selection

### Description:

The objective of this project is to determine the dynamics of updating a training population and the long-term impact on a closed breeding program. We aimed to compare several methods of training population updating in terms of prediction accuracy, response to selection, and decay of genetic variance. This R package provide the code and data to run the simulation and create the figures used in the paper.

### Paper

Please read our [pre-print on *biorXiv*](http://biorxiv.org/content/early/2016/11/10/087163).

Our paper has been accepted at [G3: Genes|Genomes|Genetics](http://www.g3journal.org/)! You can read the [early online paper](http://www.g3journal.org/content/early/2017/03/15/g3.117.040550), and the full version should be out in May.

### Installation

To install, use the `devtools` package like so:

    devtools::install_github("UMN-BarleyOatSilphium/GSSimTPUpdate")

### Obtaining Data

We used data from The Triticeae Toolbox data base in our simulations. A step-by- step set of instructions for obtaining and downloading this data is provided in a package vignette. Use the following code to access the vignette:

    vignette("t3_data_access")

### Running the Simulations

The simulations were run on a supercomputing platform that accepts batch submissions. The file [`inst/scripts/run_simulations_use.job`](https://github.com/UMN-BarleyOatSilphium/GSSimTPUpdate/tree/master/inst/scripts/run_simulations_use.job) is a `Bash` script that calls upon `R` and the script [`inst/scripts/run_experiment_use.R`](https://github.com/UMN-BarleyOatSilphium/GSSimTPUpdate/tree/master/inst/scripts/run_experiment_use.R) to run the simulations. Both of these scripts need to be edited (i.e. with different directory paths) before they can be used.

### Replicating the Tables/Figures

Tables and figures can be replicated by using the procedure outlined in the Markdown file [`inst/scripts/plot_results_use.Rmd`](https://github.com/UMN-BarleyOatSilphium/GSSimTPUpdate/tree/master/inst/scripts/plot_results_use.Rmd)

### Support

Please [open an issue](https://github.com/UMN-BarleyOatSilphium/GSSimTPUpdate/issues/new) for support with the package or to add comments.

### Contact:

[Jeff Neyhart](neyha001@umn.edu)

### Citation

To cite the pre-print, please use the following (formatting for *Genetics*):

Neyhart, J. L., T. Tiede, A. J. Lorenz, and K. P. Smith. 2016 Evaluating Methods of Updating Training Data in Long-Term Genomewide Selection. bioRxiv. doi: <http://dx.doi.org/10.1101/087163>.
