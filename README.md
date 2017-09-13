# invasive_salmonella
This repository contains data and code to accompany our paper on identifying patterns of gene degradation associated with invasiveness in Salmonella. 

The code for producing the random forest model is contained in the R notebook rf_building.Rmd. The notebook is designed to be adapted for new analyses, so has commentary and options for changing the analysis, as well as optimisation steps for assessing the ideal parameters for building your random forest. 

To use the model to determine the invasiveness index of your own samples, you can load the Rdata file finalmodel.Rdata. This contains the original training data for the model (in case you wish to compare to your samples) as well as the model itself. We recommend following the workflow outlined in the R notebook for building the model to generate a test dataset, then rather than training a model, using the predict() function in the R package randomForest.
