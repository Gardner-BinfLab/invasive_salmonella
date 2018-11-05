# invasive_salmonella
This repository contains data and code to accompany our paper on identifying patterns of gene degradation associated with invasiveness in Salmonella. 

The method relies on the DeltaBS metric (https://github.com/Gardner-BinfLab/deltaBS) for identifying deleterious mutations in protein-coding genes. The DeltaBS approach uses profile hidden Markov models to take a diverse collection of sequences of the same gene from different organisms and capture patterns of sequence variation that occur commonly in nature (in this case, within the Gammaproteobacteria). These models are then used to score the protein coding genes of these strains to produce bitscores, which indicate how well each protein conforms to modelled sequence constraints on that gene. The distribution of these bitscores can then be used by the random forest to identify genes which show a difference in adherence to sequence constraints between strains from different niches. A difference in bitscore distribution is likely to be caused by a change in selective pressures on a gene in one niche, leading to increased accumulation of mutations not commonly observed in nature (called rare mutations in this manuscript, for brevity). This accumulation of rare mutations is likely to result in the degradation of the gene, however these mutations could also result in a change to the function or structural stability of the protein. For a given protein family, the difference in bitscore for one strain and the median bitscore for the strain collection (DeltaBS) can be used as an indication of the likelihood that this accumulation of mutations is deleterious. 

The repository is separated into three directories:
* One for building your own model for identifying genes associated with a phenotype.
* One containing code for the original analysis.
* One containing code to analyse your own Salmonella strains using our model.

See the Wiki page for instructions on building your own random forest model or testing your own Salmonella isolates using our invasiveness index:
https://github.com/UCanCompBio/invasive_salmonella/wiki
