# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 13:06:55 2022

@author: Arnaud Reboud, IGE
"""

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.dates as mdate
import seaborn as sn
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from copy import copy, deepcopy
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import func # func.py: where the internal functions are located
import pydsd

path_save = './Plots/'
station_name = 'Chamrousse' #2020-12
Parsi_filename='./data/OHMCV_fr_Chamrousse_Disdro_Parsivel_202012.dat'
Spect_filename='./data/OHMCV_fr_Chamrousse_Disdro_Parsivel_Spectre_202012.dat'


#load the laws fall velocity - diameter
laws_vd = pd.read_csv('./utils/Laws_V_vs_D.csv',sep=';')
#plot the laws
func.plot_laws_vd(laws_vd)

#load the class of sizes and velocities of the Parsivel
class_size_vd = pd.read_csv('./utils/Parsivel_vd_class.csv',sep=';')

matrix_vd = func.build_matrix_precip_type(Parsivel_class=class_size_vd)

#compute the DSD from the Parsivel raw records
dsd = pydsd.read_parsivel_Campbell(Parsi_filename,
                                   Spect_filename,
                                      # start='2020-12-11 06:00:00',
                                      # stop='2020-12-13 06:00:00',
                                     resampling='1min'
)


#define the starting end ending time of the precipitation event
start_time = '2020-12-11 06:00:00'
end_time = '2020-12-13 06:00:00'

func.plot_dsd_vs_matrix(dsd, matrix_vd, spec_time = '2020-12-11 11:10:00', step_time=1,
                   classif= 'Type', xlim=(0,10), ylim=(0,12))

#Examples of use with different snowflake densities 
#Using the Parsivel spectrum (particle size and velocity distribution)
df_phase_50 = func.compute_drop_by_step(dsd, matrix_vd, 
                             start_time ='2020-12-11 06:00:00', 
                             end_time = '2020-12-13 06:00:00',
                             step_time = 1, rho_snow=50)
df_phase_Rees = func.compute_drop_by_step(dsd, matrix_vd, 
                             start_time ='2020-12-11 06:00:00', 
                             end_time = '2020-12-13 06:00:00',
                             step_time = 1, rho_snow=None)

#Using the Parsivel Code Synop information
df_Parsi = func.read_codesynop(filename=Parsi_filename, 
                             start_time ='2020-12-11 06:00:00', 
                             end_time = '2020-12-13 06:00:00')

#merge the data from the experimental process with the manufacturer-based data
df_merged_data = func.get_phase_mix(df_predicted = df_phase_50,
                               df_Parsivel = df_Parsi)
#Plot the cumulative density of the major phase contribution based on the experimental classification
func.plot_major_phase_contribution(df_phase_Rees)

#compute the confusion matrix to compare the manufacturer-based phase classification with the experimental classification
conf_mat, scores = func.get_confusion_matrix(df_compare = df_merged_data, var='Masspart')

#Plot the confusion matrix
func.plot_confusion_matrix(conf_mat)
