# glaciers
Combine projected thermal forcing and subglacial runoff into projected submarine melt rate for each glacier in Greenland. Note that at present not all of the datasets required to run the prepare_past_datasets.m script have been uploaded to the repository, because these datasets were personally requested and it has not been checked if the owners would be happy to make these datasets fully public. But it is possible to run add_future_data.m in isolation, for example if it was wanted to add new CMIP models

## process_terminus_position.mlx (note script has now been replaced by the combination of prepare_past_data.m & add_future_data.m)
Matlab live script that creates the glaciers.mat structure used to do the projections. The structure contains lots of information about each glacier in Greenland, including e.g. past ocean conditions, subglacial runoff, submarine melt rate, terminus position, ice flux etc. The most relevant pieces here are that for each CMIP model and scenario, this script associates projected thermal forcing and subglacial runoff to each glacier, bias-corrects these quantities based on the present day, and then combines these into projected submarine melt rate, which is subsequently used to project retreat

## prepare_past_data.m
Create a glaciers structure that contains past terminus positions, runoff, ocean thermal forcing and submarine melt rate. This past data is used in the calculation of the distribution for kappa. Note that at present not all of the datasets required to run this script have been uploaded to the repository, because these datasets were personally requested and it has not been checked if the owners would be happy to make these datasets fully public

## add_future_data.m
Add CMIP5/6 projected ocean temperature and runoff to the glaciers structure, do bias corrections and calculate future submarine melt. Saves as glaciers.mat, which is subsequently used to do the projections

## EN4_ISMIP6.mat
EN4 ocean temperature data used to bias-correct the ocean projections

## glaciers.mat
The glaciers structure containing all information about the glaciers
