# INBIA

Inference Network Based on iRefIndex Analysis (INBIA) is a methodology to accurately correlate proteomic inferred relations to protein-protein interaction (PPI) networks. The release of The Cancer Proteome Atlas (TCPA) has provided proteomic expression data for 190 proteins in 16 cancer types using reverse-phase protein arrays (RPPA) technology. INBIA makes use of 14 network inference methods on protein expressions related to 16 cancer types and uses  the iRefIndex human PPI network as reference model:

| Category | Method | name |
| ------------------- |:------------------------------------:| ------------:|
| Correlation         | Spearman correlation                 | [SPEARMAN]   |
|	                    | Pearson correlation                  | [PEARSON]    |
|	                    | TOM Similarity                       | [WGCNA]      |
| Partial Correlation | Simple partial correlation           | [SPC]        |
|	 	                  | GeneNet shrunken                     | [GENENET]    |
|	 	                  | Graphical lasso                      | [GLASSO]     |
| Regression          | Partial least squares regression     | [PLS]        |
|	 	                  | Ridge regression                     | [RIDGE]      |
|	 	                  | Lasso regression                     | [LASSO]      |
|	 	                  | Elastic net regression               | [ELASTICNET] |
|Mutual Information   | ARACNE additive                      | [ARACNEA]    |
|		                  | ARACNE multiplicative                | [ARACNEM]    |
|	                 	  | Context likelihood of relatedness    | [CLR]        |
|	                 	  | Maximum relevance minimum redundancy | [MRNET]      |

## Load dataset
To load RPPA datasets used to train INBIA use:
```R
fs      <- list.files(path = "TCPA_2016-05-11") # list all files
tcpa.fs <- fs[-grep(pattern = "META", x = fs)] # get only files with expression, delete metadata

# use load.dataset function to load data from TCPA.
```

or download the latest from https://tcpaportal.org/tcpa/.

## Run INBIA

To run INBIA after loading datasets from TCPA use: 
```R
execute.all(dataset, ds.name)
```
dataset is retrieved from TCPA and loaded by  `load.dataset ` function.
