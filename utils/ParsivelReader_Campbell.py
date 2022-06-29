# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 17:39:23 2022

@author: fred
"""

import numpy as np
import pandas as pd
from math import log10
from datetime import datetime, timedelta
import os.path
import re
import pydsd

import sys
  
# setting path
sys.path.append('../pydsd')
from pydsd.DropSizeDistribution import DropSizeDistribution
sys.path.append('../pydsd/io')
from pydsd.io import common

def lin(x) : return np.power(10.,x)
def lin10(x) : return np.power(10.,x/10)
def toDBz(x) : return 10*log10(x)

def remove_duplicate(df):
    """  remove duplicate lines in a data frame where index is timestamp
    Returns: df

    """
    df=df.reset_index()
    dup = df[df.duplicated()]
    if len(dup) > 0 : 
        print("\nAll the following rows are duplicated and will be removed : ")
        print(dup)
        df=df.drop_duplicates(['TIMESTAMP'],keep='first')
    df=df.set_index('TIMESTAMP')
    return df

def remove_nan(df):
    """  remove  lines with NAN values in a data frame where index is timestamp
    Returns: df

    """
    df=df.reset_index()
    row_nan = df[df.isna().any(axis=1)]
    if len(row_nan) > 0 : 
        print("\nAll the following rows contain NaN values and will be removed : ")
        print(row_nan)    
        df=df.drop(df[df.isna().any(axis=1)].index)
    df=df.set_index('TIMESTAMP')
    return df

def remove_spectrum_Nd_zero(df):
    """  remove timestamp where the number of drops in the spectrum is zero
    Returns: df

    """  
    df=df.reset_index()
    num_raw = df.columns[df.columns.str.contains('var93split', regex=False)]
    Sp_Nd_zero = df.loc[df[num_raw].loc[:].sum(axis=1) == 0]
    if len(Sp_Nd_zero ) > 0 : 
        print("\nAll the following Timestamp with 0 drops in spectrum will be removed : ")
        print(Sp_Nd_zero)
        df=df.loc[df[num_raw].loc[:].sum(axis=1) > 0]
    df=df.set_index('TIMESTAMP')
    return df   
        

def read_parsivel_Campbell(Parsi_filename,Spect_filename,start=None,stop=None,v_filter=None,resampling=None):
    """
    Takes 2 filenames pointing to a parsivel raw file and a parsivel raw spectrum  
    and returns a drop size distribution object.
    start : define the starting date of sampling extraction - the nearest will be used
    if start is not set, extraction will start at the begining of the sampling 
    stop  : define the stoping date of the sampling extraction - the nearest will be used
    if stop is not set, extraction will end at the ending of the sampling 
    if start and stop are not set, all the sampling will be passed to the dsd object
    start and stop should be '%Y-%m-%d %H:%M:%S'
    
    if v_filter is set to True, the raw spectrum is multiplyed by de +or- 50% speed matrix
    
    resampling : resampling time of the dataset (in minutes)

    Usage:
    dsd = read_parsivel_Campbell(Parsi_filename,Spect_filename,start=None,stop=None,v_filter=None,resampling=None)

    Returns:
    DropSizeDistrometer object

    """   
    if os.path.exists(Parsi_filename) != True :
        print ("Parsi_filename doesn't exist !")
        return None
    
    if os.path.exists(Spect_filename) != True :
        print ("Spect_filename doesn't exist !")
        return None

    if start != None :
        try :
            start = pd.to_datetime(start, utc= True)
        except :
                print ("Start date format is not correct, try with %Y-%m-%d %H:%M:%S")
                return None
    
    if stop != None :
        try :
            stop = pd.to_datetime(stop, utc= True)
        except :
            print ("Stop date format is not correct, try with %Y-%m-%d %H:%M:%S")
            return None
    
    reader = ParsivelReaderCampbell(Parsi_filename,Spect_filename,start,stop,v_filter,resampling)
  
    dsd = DropSizeDistribution(reader)
    return dsd  


class ParsivelReaderCampbell(object):

    """
    ParsivelReader class takes  two filenames as only arguments(for now).
    The first one should be a parsivel raw datafile at Campbell format (output from the Campbell).
    The second one should be a parsivel raw spectrum datafile at Campbell format (output from the Campbell)
    """

    def __init__(self, Parsi_filename,Spect_filename,start,stop,v_filter,resampling):
        self.Parsi_filename = Parsi_filename
        self.Spect_filename = Spect_filename
        self.resampling=resampling
        self.rain_rate = []
        self.Z = []
        self.num_particles = []
        self._base_time = []

        self.nd = []
        self.vd = []
        self.raw = []
        self.code = []
        self.time = []

        self.ndt = []
        self.spect_index = []
        self.smp_interval = []

        self.pcm = np.reshape(self.pcm_matrix, (32, 32))
        


        self._read_file_Spectre(start,stop)
        self._read_file_Parsivel()
        
        if np.any(self.smp_interval == None) :
            print('test')
            new_sampling_sec = np.nanmin(np.diff(self.time['data']))
            self.smp_interval = np.full(len(self.spect_index), new_sampling_sec)     
            
        if self.resampling != None :
            try : 
                x = int(re.findall(r'\d+', self.resampling)[0])
            except : 
                x = 1
            if len(re.findall('min', self.resampling)) == 1 :
                new_sampling_sec = x * 60
            else :
                if len(re.findall('H', self.resampling)) == 1 :
                    new_sampling_sec = x * 3600
                else : 
                    if len(re.findall('D', self.resampling)) == 1 :
                        new_sampling_sec = x * 86400   
            self.div_sampling = np.nanmin(self.smp_interval) / new_sampling_sec
            self.smp_interval = new_sampling_sec*np.ones(len(self.spect_index))
        else:
            self.div_sampling  = 1
            new_sampling_sec = 0
            


        print("Resampling = ",new_sampling_sec,"s") 
        print("div_sampling = ",self.div_sampling) 
        
        if v_filter == True :
            self._apply_pcm_matrix()
            self.raw = self.filtered_raw_matrix 
            del self.filtered_raw_matrix 
        else :
            self.raw=np.reshape(self.raw, (len(self.spect_index),32, 32))

        self._prep_data()

        self.bin_edges = np.hstack(
             (0, self.diameter["data"] + np.array(self.spread["data"]) / 2)
        )

        self.bin_edges = common.var_to_dict(
             "bin_edges", self.bin_edges, "mm", "Bin Edges"
        )
        
        
        self.spectrum_fall_velocity=self.fields["spectrum_fall_velocity"]
        
                        
        print("\nDSD sampling available from "
              +str(pd.to_datetime(self.spect_index[0]))
              +" to "
              +str(pd.to_datetime(self.spect_index[-1])))
        
    def _read_file_Parsivel(self):
        """  Read the Parsivel Data file and store it in internal structure.
        Returns: None

        """
        print("Parsivel2 file in process : ")
        
        df = pd.read_csv(self.Parsi_filename, skiprows=[0,2,3], index_col=0, parse_dates=True,
                         na_values=["NAN", "INF", "-INF"], low_memory=False)

        df = remove_duplicate(df)
        
        df['Radar_Reflectivity']=df['Radar_Reflectivity'].apply(lin10)
        
        if self.resampling != None :
            saved_df = df.copy(deep=True)
            pd.options.mode.chained_assignment = None
            agg_dict = {'Intensity' : lambda x: x.sum(min_count=1),
                        'Radar_Reflectivity' : lambda x: x.sum(min_count=1),
                        'Number_of_particles' : lambda x: x.sum(min_count=1),
                        'Smp_interval' : 'mean'}
            try:
                df=df.resample(self.resampling).agg(agg_dict)
                del saved_df
            except:
                df = saved_df.copy(deep=True)
                del saved_df
        
        
        df=df.reindex(self.spect_index,copy=True)
        
        self.rain_rate = df.Intensity[self.spect_index]
        self.Z = df.Radar_Reflectivity[self.spect_index]
        self.num_particles = df.Number_of_particles.loc[self.spect_index].values
        
        try:
            self.smp_interval = df.Smp_interval.loc[self.spect_index].values
        except:
            self.smp_interval = None
        
    def _read_file_Spectre(self,start,stop):
        """  Read the Parsivel Data file and store it in internal structure.
        Returns: None

        """
        cols_Nd = [u'mm062', u'mm187', u'mm312', u'mm437', u'mm562', u'mm687', u'mm812',
               u'mm937', u'mm1062', u'mm1187', u'mm1375', u'mm1625', u'mm1875',
               u'mm2125', u'mm2375', u'mm2750', u'mm3250', u'mm3750', u'mm4250',
               u'mm4750', u'mm5500', u'mm6500', u'mm7500', u'mm8500', u'mm9500',
               u'mm11000', u'mm13000', u'mm15000', u'mm17000', u'mm19000', u'mm21500',
               u'mm24500']

        cols_Vd= ["ms250","ms350","ms450","ms550","ms650","ms750","ms850","ms950","ms1100",
                  "ms1300","ms1500","ms1700","ms1900","ms2200","ms2600","ms3000","ms3400",
                  "ms3800","ms4400","ms5200","ms6000","ms6800","ms7600","ms8800","ms10400",
                  "ms12000","ms13600","ms15200","ms17600","ms20800"]
        
        bin_edges = [0,0.1245,0.2495,0.4995,0.6245,0.7495,0.8745,0.9995,1.1245,
                   1.50,1.75,2.00,2.25,2.50,3.00,3.50,4.00,4.50,5.00,6.00,
                   7.00,8.00,9.00,10.0,12.0,14.0,16.0,18.0,20.0,23.0,26.0]
        
        print("Parsivel2 spectrum raw file in process : ")
        
        df = pd.read_csv(self.Spect_filename, skiprows=[0,2,3], index_col=0, parse_dates=True,
                         na_values=["NAN", "INF", "-INF"], low_memory=False)
        
        df = remove_spectrum_Nd_zero(df)
        df = remove_duplicate(df)
        df = remove_nan(df)

        num_raw = df.columns[df.columns.str.contains('var93split', regex=False)]
        
        if self.resampling != None : 
            agg_dict = {name : 'sum' for name in cols_Nd}
            agg_dict |= {name : 'sum' for name in cols_Vd}
            agg_dict |= {name : 'sum' for name in num_raw}
            try:
                print ("Dataset resampled over "+self.resampling+" period")
                df[cols_Nd]=df[cols_Nd].replace(-10.,np.nan).apply(lin)
                df=df.resample(self.resampling).agg(agg_dict)
                df=df.fillna(np.power(10, -10.))
                # to clean empty spectrum after resampling operation 
                df = remove_spectrum_Nd_zero(df)

            except:
                print ("Resampling format is probably not correct.")
                print ("Should be like \'5min\’ or \'H\' or \'D\' ..")
                print ("see Offset aliases in df.resample Python documentation.")
                print ("Resamplig aborted.")
                self.resampling = None
                df=df.fillna(np.power(10, -10.))
        else:
            df[cols_Nd]=df[cols_Nd].apply(lin)

        index=df.index
        
        if start != None :
           index=df.index[np.where(pd.to_datetime(df.index, utc = True) >= start)]
    
        if stop != None :
           index=index[np.where(pd.to_datetime(index, utc = True) <= stop)]
           
        self.spect_index = index
        
        self.nd = np.reshape(df[cols_Nd].loc[index].values,(len(index),32)) #Nd
        self.vd = np.reshape(df[cols_Vd].loc[index].values,(len(index),30)) # Vd
        self.raw = df[num_raw].loc[index].values

        time_epoch = (index - datetime.utcfromtimestamp(0)) / np.timedelta64(1, 's')
        self.time = {
            "data": time_epoch.values,
            "units": common.EPOCH_UNITS,
            "title": "Time",
            "long_name": "time",
        }
        
    def _apply_pcm_matrix(self):
        """ Apply Data Quality matrix from Ali Tokay
        Returns: None

        """
        print("PCM Matrix applied")
        self.filtered_raw_matrix = np.ndarray(
            shape=(len(self.raw), 32, 32), dtype=float
        )
        for i in range(len(self.raw)):
            self.filtered_raw_matrix[i] = np.multiply(
                self.pcm, np.reshape(self.raw[i], (32, 32))
            )

    def _prep_data(self):
        self.fields = {}

        self.fields["rain_rate"] = common.var_to_dict(
            "Rain rate", np.ma.array(self.rain_rate*self.div_sampling), "mm/h", "Rain rate"
        )
        mask=self.Z.isna()
        self.fields["reflectivity"] = common.var_to_dict(
            "Reflectivity",
            np.ma.array((self.Z*self.div_sampling).apply(toDBz),mask=mask),
            "dBZ",
            "Equivalent reflectivity factor",
        )
        Nd = np.ma.masked_equal(self.nd, np.power(10, -10.))
        Nd= Nd.__mul__(self.div_sampling)
        self.fields["Nd"] = common.var_to_dict(
            "Nd",
            Nd,
            "m^-3 mm^-1",
            "Liquid water particle concentration",
        )
        self.fields["Nd"]["data"].set_fill_value(0)
        
        mask=np.isnan(self.num_particles)
        self.fields["num_particles"] = common.var_to_dict(
            "Number of Particles",
            np.ma.array(self.num_particles,mask=mask),
            "",
            "Number of particles",
        )
        self.fields["drop_spectrum"] = common.var_to_dict(
            "drop spectrum",
            np.array(self.raw),  
            "",
            "Number of drops in a 32 x 32 size versus fall velocity matrix",
        )
        self.smp_interval = {"data": np.ma.masked_equal(self.smp_interval,np.nan)} 
        
 # # # A modifier 
  #       self.fields["spectrum_fall_velocity"] = common.var_to_dict(
  #           "Terminal Fall Velocity",
  #           np.array(
  # #              self.vd[0]
  #               self.vd   # if timeserie is needed
  #           ),  # Should we do something different here? Don't think we want the time series.
  #           "m/s",
  #           "Terminal fall velocity for each bin",
  #       )
  
        # try:  # check if filtered_raw_matrix exits 
        #     raw_matrix = self.filtered_raw_matrix
        # except:
        #     raw_matrix  = np.reshape(self.raw, (len(self.spect_index),32, 32))
            
     
        
        diameter= np.array([ 0.062,  0.187,  0.312,  0.437,  0.562,  0.687,  0.812,  0.937,
                1.062,  1.187,  1.375,  1.625,  1.875,  2.125,  2.375,  2.75 ,
                3.25 ,  3.75 ,  4.25 ,  4.75 ,  5.5  ,  6.5  ,  7.5  ,  8.5  ,
                9.5  , 11.   , 13.   , 15.   , 17.   , 19.   , 21.5  , 24.5  ])
                
        self.diameter = common.var_to_dict(
            "diameter", diameter, "mm", "Particle diameter of bins",
        )
        spread = np.array([ 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
               0.125, 0.25 , 0.25 , 0.25 , 0.25 , 0.25 , 0.5  , 0.5  , 0.5  ,
               0.5  , 0.5  , 1.   , 1.   , 1.   , 1.   , 1.   , 2.   , 2.   ,
               2.   , 2.   , 2.   , 3.   , 3.   ])
        self.spread = common.var_to_dict(
            "spread", spread, "mm", "Bin size spread of bins",
        )        
        # bin_edges = [  124.25,   249.5 ,   374.5 ,   499.5 ,   624.5 ,   749.5 ,
        #          874.5 ,   999.5 ,  1124.5 ,  1249.75,  1500.  ,  1750.  ,
        #         2000.  ,  2250.  ,  2500.  ,  3000.  ,  3500.  ,  4000.  ,
        #         4500.  ,  5000.  ,  6000.  ,  7000.  ,  8000.  ,  9000.  ,
        #        10000.  , 12000.  , 14000.  , 16000.  , 18000.  , 20000.  ,
        #        23000.  , 26000.  ]
        self.bin_edges = common.var_to_dict(
            "bin_edges",
            np.hstack((0, self.diameter["data"] + np.array(self.spread["data"]) / 2)),
            "mm",
            "Boundaries of bin sizes",
        )
        
        velocity = {'data': np.array([ 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 
                                   0.75, 0.85, 0.95, 1.1 , 1.3 , 1.5 , 1.7 , 
                                   1.9 , 2.2 , 2.6 , 3. , 3.4 , 3.8 , 4.4 , 
                                   5.2 , 6. , 6.8 , 7.6 , 8.8 , 10.4 , 12. , 
                                   13.6 , 15.2 , 17.6 , 20.8 ]), 
                    'long_name': 'Terminal fall velocity for each bin', 
                    'standard_name': 'velocity', 'units': 'm s^-1'}
        
        self.fields["spectrum_fall_velocity"] = velocity

        
        v_spread = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 
                    0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.4, 0.8, 0.8, 
                    0.8, 0.8, 0.8, 1.6, 1.6, 1.6, 1.6, 1.6, 3.2, 3.2]
        
 

    diameter = common.var_to_dict(
        "diameter",
        np.array(
            [
                0.06,
                0.19,
                0.32,
                0.45,
                0.58,
                0.71,
                0.84,
                0.96,
                1.09,
                1.22,
                1.42,
                1.67,
                1.93,
                2.19,
                2.45,
                2.83,
                3.35,
                3.86,
                4.38,
                4.89,
                5.66,
                6.7,
                7.72,
                8.76,
                9.78,
                11.33,
                13.39,
                15.45,
                17.51,
                19.57,
                22.15,
                25.24,
            ]
        ),
        "mm",
        "Particle diameter of bins",
    )
 
    pcm_matrix =   (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

