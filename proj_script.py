'''
Author: Kendra Lyons
Date: 12/1/2020
Class: ISTA 131
Section Leader: Jocelyn Connors
Final Project Script

Description:
Script for analysis and creation of figures for final project. 
'''
import pandas as pd
import numpy  as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
# import * from proj_script

def read_black_holes():
    black_holes = pd.read_csv('https://raw.githubusercontent.com/eventhorizontelescope/2019-D01-01/master/csv/SR1_M87_2017_095_hi_hops_netcal_StokesI.csv', header=1)
    return black_holes

def read_exoplanets():
    '''
    Reads the relevant data from 'exoplanet_cumulative.csv' into a new DataFrame and returns it. 
    '''
    epdf = pd.read_csv('~/Documents/UA_classes_FA20/ISTA131/Final Project/exoplanet_cumulative.csv', 
        index_col = 'kepoi_name', usecols = ('kepoi_name','koi_pdisposition','koi_score', 'koi_period', 'koi_prad'))
    return epdf
    
def get_candidates(epdf):
    '''
    Makes a DataFrame with data for all KOI Candidates
    '''
    cds = epdf.loc[epdf.koi_pdisposition == 'CANDIDATE']
    return cds
    
def get_false_positives(epdf):
    '''
    Makes a DataFrame with data for all KOI False Positives
    '''
    fps = epdf.loc[epdf.koi_pdisposition == 'FALSE POSITIVE']
    return fps

def koi_score_frequency(cds, fps):
    '''
    '''
    ax1 = plt.hist(cds.koi_score, bins = 10, alpha = .5, ec = 'green', color = 'lightgreen', label = 'Candidates')
    ax2 = plt.hist(fps.koi_score, bins = 10, alpha = .5, ec = 'orange', color = 'gold', label = 'False Positives')
    # ax3 = ax1.twinx()
    plt.xlabel('KOI Score', labelpad = 5, fontsize = 20)
    plt.ylabel('Frequency', labelpad = 5, fontsize = 20)
    plt.title('KOI Score Frequency Distribution', pad = 10, fontsize = 24)
    plt.legend(loc = 'center', title = 'KOI Disposition')

def koi_score_v_period(cds, fps):
    '''
    Makes a scatterplot of the KOI Score by planetary orbital period in days. 
    '''
    cds_period = cds.loc[:, 'koi_period']
    cds_score = cds.loc[:, 'koi_score']
    fps_period = fps.loc[:, 'koi_period']
    fps_score = fps.loc[:, 'koi_score']
    plt.scatter(cds_period, cds_score, marker = 'o', color = 'green', alpha = .2, label = 'Candidates')
    plt.scatter(fps_period, fps_score, marker = 'o', color = 'orange', alpha = .2, label = 'False Positives')
    plt.xlim(0,700)
    plt.xlabel('Orbital Period in Days', labelpad = 5, fontsize = 20)
    plt.ylabel('KOI Score', labelpad = 5, fontsize = 20)
    plt.title('KOI Score by Planetary Orbital Period', pad = 10, fontsize = 24)
    plt.legend(title = 'KOI Disposition')
    df = pd.concat([cds_period, cds_score], axis = 1)
    df = df.dropna()
    #candidates linear regression line
    cds_df = pd.concat([cds_period, cds_score], axis = 1)
    cds_df = cds_df.dropna()
    cds_period = sm.add_constant(cds_df['koi_period'].values)
    model = sm.OLS(cds_df['koi_score'], cds_period)
    results = model.fit()
    xs = np.arange(0, 701)
    ys = results.params.loc['x1'] * xs + results.params.loc['const'] 
    plt.plot(xs, ys, linewidth = 3, color = 'lightgreen')
    #false positives linear regression line 
    fps_df = pd.concat([fps_period, fps_score], axis = 1)
    fps_df = fps_df.dropna()
    fps_period = sm.add_constant(fps_df['koi_period'].values)
    model = sm.OLS(fps_df['koi_score'], fps_period)
    results = model.fit()
    xs = np.arange(0, 701)
    ys = results.params.loc['x1'] * xs + results.params.loc['const'] 
    plt.plot(xs, ys, linewidth = 3, color = 'gold')

def koi_score_v_radius(cds, fps):
    '''
    Makes a scatterplot of the KOI Score by planetary radius in Earth radii. 
    '''
    cds_radius = cds.loc[:, 'koi_prad']
    cds_score = cds.loc[:, 'koi_score']
    fps_radius = fps.loc[:, 'koi_prad']
    fps_score = fps.loc[:, 'koi_score']
    plt.scatter(cds_radius, cds_score, marker = 'o', color = 'green', alpha = .2, label = 'Candidates')
    plt.scatter(fps_radius, fps_score, marker = 'o', color = 'orange', alpha = .2, label = 'False Positives')
    plt.xlim(0, 500)
    plt.xlabel('Planetary Radius in Earth Radii', labelpad = 5, fontsize = 20)
    plt.ylabel('KOI Score', labelpad = 5, fontsize = 20)
    plt.title('KOI Score by Planetary Radius', pad = 10, fontsize = 24)
    plt.legend(title = 'KOI Disposition')
    #candidates linear regression line
    cds_df2 = pd.concat([cds_radius, cds_score], axis = 1)
    cds_df2 = cds_df2.dropna()
    cds_df2.head()
    cds_radius = sm.add_constant(cds_df2['koi_prad'].values)
    model = sm.OLS(cds_df2['koi_score'], cds_radius)
    results = model.fit()
    xs = np.arange(0, 501)
    ys = results.params.loc['x1'] * xs + results.params.loc['const'] 
    plt.plot(xs, ys, linewidth = 3, color = 'lightgreen')
    #false positives linear regression line
    fps_df2 = pd.concat([fps_radius, fps_score], axis = 1)
    fps_df2 = fps_df2.dropna()
    fps_df2.head()
    fps_radius = sm.add_constant(fps_df2['koi_prad'].values)
    model = sm.OLS(fps_df2['koi_score'], fps_radius)
    results = model.fit()
    xs = np.arange(0, 501)
    ys = results.params.loc['x1'] * xs + results.params.loc['const'] 
    plt.plot(xs, ys, linewidth = 3, color = 'gold')

#==========================================================
def main():
    '''
    Write a description of what happens when you run
    this file here.
    '''
    expl = read_exoplanets()
    candidates = get_candidates(expl)
    false_positives = get_false_positives(expl)
    koi_score_frequency(candidates, false_positives)
    plt.show()
    plt.figure()
    koi_score_v_period(candidates, false_positives)
    plt.show()
    plt.figure()
    koi_score_v_radius(candidates, false_positives)
    plt.show()
  

if __name__ == '__main__':
    main()
