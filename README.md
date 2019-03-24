# Human-Genotype-Visualization-GSoC-19-CCCG

This Repository Contains Code , Output Files and Data of the Selection Test by Canadian Centre for Computaional Genomics GSoC'19.

The Repository has 3 folders each containing LOG_Files , Images and Reduced 2D coordinates from Dimensionality Reduction.
#The Problem statement -

Selection tests
Write a program in python that does the following:

Imports the genotype data in VCF and population description files (suggested package: scikit-allel).
Dimensionally reduces the genotype data using principal components analysis and uniform manifold approximation and projection (there are libraries for this)
Outputs a file containing coordinates of the low dimensional projections
Outputs a 2D plot, which is coloured by the variables in the panel file and labelled using full variable names in the tsv file
Outputs a log containing the random seeds used to generate coordinates, the program's run-time, versions of all packages used, verbose output of methods used, time stamps, and all parameters specified (results should be reproducible)
All outputs should have a unique ID connecting them to the log file

#Data is located at:

#http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/
#http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/

Umap tutorial
https://umap-learn.readthedocs.io/en/latest/how_umap_works.html

Files are:

#ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf.gz
#affy_samples.20141118.panel
#20131219.populations.tsv

The below image displays the fact that Genotype Variations can be linked to Human Evolution and Genetic Variation can be classified by origin. This image file is one of the outputs acquired from the data.

![Image](https://github.com/DataGasmic/Human-Genotype-Visualization-GSoC-19-CCCG/blob/master/2D-Image-Plots/Image_25_15_0.5.pdf)
