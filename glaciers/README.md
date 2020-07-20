# glaciers
Combine projected thermal forcing and subglacial runoff into projected submarine melt rate for each glacier in Greenland

## process_terminus_position.mlx
Matlab live script that creates the glaciers.mat structure used to do the projections. The structure contains lots of information about each glacier in Greenland, including e.g. past ocean conditions, subglacial runoff, submarine melt rate, terminus position, ice flux etc. The most relevant pieces here are that for each CMIP model and scenario, this script associates projected thermal forcing and subglacial runoff to each glacier, bias-corrects these quantities based on the present day, and then combines these into projected submarine melt rate, which is subsequently used to project retreat

## glaciers.mat
The glaciers structure containing all information about the glaciers
