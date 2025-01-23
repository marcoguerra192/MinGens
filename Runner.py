## RUN THE COMPUTATION OF MIN GENS AND DEATH POINTS
import os
from ComputeMinGensDP import computeMinGensAndDeathPoints


data_folder = os.path.join(os.getcwd(), 'data')
output_folder = os.path.join(os.getcwd(), 'results')
picture_folder = os.path.join(os.getcwd(), 'results')

Thresh = 3.7
plot = True
deathpoints=True
verbose=True
square=False


computeMinGensAndDeathPoints(Thresh = Thresh, 
                             data_folder = data_folder,
                             output_folder = output_folder,
                             plot=plot,
                             picture_folder=picture_folder, 
                             deathpoints=deathpoints,
                             verbose=verbose,
                             square=square)

