# Dynamic-Strength-and-Fragmentation

This repository contains all data and plotting scripts used to generate the data presented in "Dynamic Compressive Strength and Fragmentation in Sedimentary and Metasedimentary Rocks" by Auriol S. P. Rae, T. Kenkmann, V. Padmanabha, M. H. Poelchau, F. Schäfer, M Dörffler, and L. Müller. This repository also contains data from an earlier paper "Dynamic Strength and Fragmentation and Fragmentation in Felsic Crystalline Rocks" by Auriol S. P. Rae, T. Kenkmann, V. Padmanabha, M. H. Poelchau, and F. Schäfer.

This repository contains three subdirectories:

1) Mechanical Testing and Analysis
  - PetrophysicalProperties.xlsx - Spreadsheet containing a summary of material properties determined in this study or derived from the literature, includes results from He-pycnometry.
  - Rae2021 - Directory containing spreadsheets with summaries of mechanical experiments in this study (Rae et al., 2021(?)) organised by lithology (SeeSst - Seeberger Sandstone, TaQu - Taunus Quartzite, SaLi - Savonnieres Limestone, CaMa - Carrara Marble). 
  - Rae2020 - Directory containing spreadsheets with summaries of mechanical experiments in the study of Rae et al. (2020) organised by experimental method. Also available at http://www.doi.org/10.5281/zenodo.3987223
  - SHPBdata - Directory containing all raw SHPB data. Data is organised into five columns (time in seconds; voltage from 1st strain gauge on incident bar; voltage from 2nd strain gauge on incident bar; voltage from 1st strain gauge on transmission bar; voltage from 2nd strain gauge on transmission bar). Directory contains the python script used to process the data (compatible with Python 2.7).

2) Fragment Size
  - AllSieveData.xlsx - Spreadsheet containing sieve results for all lithologies in this study and in the study of Rae et al. (2020; Also available at http://www.doi.org/10.5281/zenodo.3987223)

3) Fragment Shape
  - ImageAnalysis - Contains a list of the samples where fragment shapes were analysed, all original and final thresholded images for fragments > 2 mm (LargeFragments) and for fragments between 0.5 and 2 mm (SmallFragments), the output image analysis of each image from Fiji, and python scripts used to analyse the results.
  - GeometricFragmentationAlgorithms - Three python scripts, one for each geometric algorithm analysed in this study. Each script generates a set images with increasing numbers of fragments using the specified algorithm, the script then analyses those generated images to determine the fragment shape distributions for each of those images.

