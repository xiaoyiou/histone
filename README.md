# histone
Histone Modification Experiment (code only)
This repository contains the code for running epigenetics analysis using various data mining algorithms. 
## Loading data
Loadfiles.py performs the preprocessing of the loaded track data and overlay them with gene annotations from TAIR database

## Searching for patterns
mdata.Madata is the class which encapsulate "experiements" objectcs. 

## Find patterns 
In order to find the function-specfic patterns, we need to addlabels first and then call the calcPScores

## Function specific scores 
The final result will be stored in obj.PScores 


----
## Shape-based method
We apply the wavelete transformation for series like data and cluster using Affinity Propagation. Each centroid is used as the 
"signature" shape for histone modfiication series.

