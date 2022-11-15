# proteiNorm

The app has been updated to work with R 4.2.2. Data from both mass spectrometry quantitation methods, Tandem Mass Tag (TMT) and label-free mass spectrometry, are supported.

## Data Requirements
# TMT
As input, proteiNorm expects tab-separated peptide (optional) and protein data (not on logarithmic scale) as produced by software such as MaxQuant, where each row represents a peptide or protein and the column names of the measured intensities (samples) beginning with “Reporter intensity corrected” followed by an integer and an optional label (e.g. “Reporter intensity corrected 5 TMT2”) for TMT experiments. 
The TMT reporter intensities should be corrected using the error correction factors provided by Thermo Fisher for the specific TMT lot that was used during the sample preparation. These correction factors are imported into MaxQuant when setting up the database search parameters. The MaxQuant “proteinGroups.txt” output file will then contain column names for each TMT reporter ion labeled as “Reporter intensity corrected” and “Reporter intensity”.

# label free
The column names for samples in label free experiments should begin with “Intensity” followed by an integer (“Intensity 01”). 

# example data
In the following examples, we are using the “Reporter intensity corrected” columns from MaxQuant since these intensities apply the TMT correction factors. The reporter ion isotopic distributions (-2, -1, +1, +2) are primarily for carbon isotopes with reporter ion interference for each mass tag. ProteiNorm automatically detects these column names, determines if the data is TMT or label-free data type. ProteiNorm also removes proteins flagged as common contaminates and reverse sequences in the MaxQuant output using the column names “Potential contaminant” and “Reverse” provided in the proteinGroups.txt MaxQuant output file. If you prefer to load in a different matrix of data, then the column labels for the samples will need to be modified to match the naming structure provided here. An example of the data input files are hosted on Github for both TMT and label-free data sets. 

Due to the detection limit of mass spectrometry instruments, many measurements of peptides or proteins result in intensity levels of zero. These values will be considered as missing values (denoted as NA) and can be imputed with precaution. A modified heatmap of missing values (Figure 2C) from the DEP Bioconductor package 10 helps to determine if data is missing at random (MAR) or missing not at random (MNAR) and the MSnbase vignette describes different imputation methods used for different types of missing data 
