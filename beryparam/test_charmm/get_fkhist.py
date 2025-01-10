

import matplotlib.pylab as plt
import seaborn as sns
import numpy as np 
import pandas as pd
import sys
import numpy.ma as ma

from openmmtools.multistate import MultiStateReporter
from netCDF4 import Dataset

fnum=sys.argv[1]

lambda_electrostatics = [1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50, 0.45, 0.40, 0.35, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00]
lambda_sterics =  [1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.00]

for i in np.arange(1,3):
    ncfile = Dataset('explicit_'+str(fnum)+'/experiments/solvent'+str(i)+'.nc', 'r')
    
    fk_hist = ncfile.groups['online_analysis'].variables['f_k_history']
    n_iterations, n_states = fk_hist.shape
    
    #lastLine1 = ma.getdata(fk_hist)
    #df1 = pd.DataFrame({'lambE' : lambda_electrostatics, 'lambS' : lambda_sterics, 'val': lastLine1})
    #df1.to_csv('lambda_'+str(fnum)+'_solv'+str(i)+'_all.csv')

    # Plot the free energy history
    
    fig1 = plt.figure(figsize=(10,5));
    plt.plot(fk_hist, marker='.')

    plt.xlabel('iteration');
    plt.ylabel('Free energy (kT)');
    plt.savefig('fk_hist_'+str(fnum)+'_solv'+str(i)+'.png', dpi=300)
    
    lastLine2 = ma.getdata(fk_hist[-2])
    df2 = pd.DataFrame({'lambE' : lambda_electrostatics, 'lambS' : lambda_sterics, 'val': lastLine2})
    df2.to_csv('lambda_'+str(fnum)+'_solv'+str(i)+'_last.csv')

    fig3 = plt.figure(figsize=(10,5)) 
    plt.plot(lambda_electrostatics, lastLine2, marker='.')
    plt.savefig('lambdaE_'+str(fnum)+'_solv'+str(i)+'_last.png', dpi=300)
