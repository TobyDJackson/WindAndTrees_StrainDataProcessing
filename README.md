# WindAndTrees_StrainDataProcessing
Complete strain data collected in Wytham Woods is available here: https://catalogue.ceh.ac.uk/documents/533d87d3-48c1-4c6e-9f2f-fda273ab45bc

Complete strain data collected in Danum Valley is available here: https://doi.org/10.5285/657f420e-f956-4c33-b7d6-98c7a18aa07a

Sample data from 4 pairs of strain gauges at 1.3m on 4 Birch trees in Wytham hosted here https://goo.gl/jr5jMT. I would suggest starting with this if you are unfamiliar with strain data. 

This data is in matlab .mat format and contains a datenum column followed by 8 columns of strain data. The code 'strain_process' takes this data through the steps to estimating critical wind speeds, if you run it section by section it should be simple enough - the other codes are dependencies and there is a file of wind data. Choices along the way are explored in the issues.
