# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 23:13:11 2022

@author: rebouda
"""

def fall_velocity(d):
    """
    Return the terminal fall speed of a hydrometeor according its diameter

    Parameters
    ----------
    d : float
        Diameter in mm of the hydrometeor.

    Returns
    -------
    float
        Vertical velocity of the falling particles in m/s.

    """
    import numpy as np
    #d= diameter in mm, return v=velocity in m/s
    #Empirical formula for rain: v= 9.65 - 10.3*np.exp(-0.6*d) #fall velocity for Rain part. from Atlas(1973)
    if d<1.34 :
        return 2.158*d + 0.095
    elif d<5 :
        return 0.0919*np.power(d,3) - 1.2145*np.power(d,2) + 5.8075*d - 2.8398
    else : 
        return 0.0389*d + 7.2293

def plot_laws_vd(laws_vd):
    """
    Plot the laws velocity (diameters) from the file provided.

    Parameters
    ----------
    laws_vd : Dataframe
        Some v(d) laws for different types of precipitation. V is in m/s and
        d in mm.

    Returns
    -------
    None.

    """
    import matplotlib.pyplot as plt
    import numpy as np
    
    cmap = plt.get_cmap('tab20c')
    for laws in laws_vd.index:
        if str(laws_vd.loc[laws]['Type']) != 'nan' :
            alpha = laws_vd.loc[laws]['alpha']
            lamb = laws_vd.loc[laws]['lambda']
            dmin = laws_vd.loc[laws]['Dmin']
            dmax = laws_vd.loc[laws]['Dmax']
            d = np.arange(dmin,dmax,0.1)
            if laws_vd.loc[laws]['Law'] == 'Exponential' :
                C = laws_vd.loc[laws]['C']
                v = C + alpha * np.exp(lamb*d)
            if laws_vd.loc[laws]['Law'] == 'Power' :
                v = alpha * np.power(d,lamb)
            name = str(laws_vd.loc[laws]['Name'])
            plt.plot(d, v, label=name, c= cmap(laws), ls='--')

    plt.xlabel('Diameter [mm]')
    plt.ylabel('Velocity [m/s]')
    plt.xlim(0,10)
    plt.ylim(0,12)
    plt.grid(alpha=0.5, dashes = (5,5))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol = 3)
    plt.tight_layout()
    plt.show()
    #plt.savefig(path_save+'matrix_vs_laws.png',format='png',dpi=100)
    
def build_matrix_precip_type(Parsivel_class, v_g=10, v_s=3, d_h=5):
    """
    Return the matrix of classification of Precipitation types according the
    Parsivel classes of particle sizes and velocities.

    Parameters
    ----------
    Parsivel_class : Dataframe
        Data list of the class of velocity and particle size of Parsivel.
    v_g : float, optional
        Graupel/Rain velocity limit in m/s. The default is 10.
    v_s : float, optional
        Snow/rain velocity limit in m/s. The default is 3.
    d_h : float, optional
        Hail/Rain limit of the particle diameters in mm. The default is 5.

    Returns
    -------
    matrix_vd : DataFrame
        Classification matrx of the 6 precipitation types.

    """
    import pandas as pd
    import numpy as np
    
    matrix_vd=pd.DataFrame(index=Parsivel_class.Size_Class_average, 
                           columns=Parsivel_class.Velocity_Class_average)
    for diam in matrix_vd.index:
        for vel in matrix_vd.columns:
            if vel > fall_velocity(diam):
                if diam > d_h: #for Hail
                    matrix_vd.loc[diam][vel]=1
                else:
                    if vel > v_g: #for Graupel
                        matrix_vd.loc[diam][vel]=2
                    else:
                        if vel < v_s: #for Drizzle
                            matrix_vd.loc[diam][vel]=4
                        else: #for Rain
                            matrix_vd.loc[diam][vel]=3
            else:
                if vel > v_s: #for Mixture
                    matrix_vd.loc[diam][vel]=5
                else: #for Snow
                    matrix_vd.loc[diam][vel]=6
    matrix_vd = matrix_vd.astype('int')
    matrix_vd.iloc[0:2,:] = np.nan #first 2 classes of Parsivel are below the sensor sensitivity
    return matrix_vd

def Synop_to_Type(Code_4680):
    """
    Convert the Synop Code from the manufacturer Parsivel output to a 
    corresponding precipitation type.

    Parameters
    ----------
    Code_4680 : Int
        Synop Code from Parsivel.

    Returns
    -------
    str
        Name of the precipitation type.

    """
    if Code_4680 == 0:
        return 'No precipitation'
    elif (Code_4680 >= 50) & (Code_4680 < 55):
        return 'Drizzle'
    elif (Code_4680 >= 55) & (Code_4680 < 60):
        return 'Drizzle with rain'
    elif (Code_4680 >= 60) & (Code_4680 < 65):
        return 'Rain'
    elif (Code_4680 >= 65) & (Code_4680 < 70):
        return 'Mixture'
    elif (Code_4680 >= 70) & (Code_4680 < 80):
        return 'Snow'
    elif (Code_4680 == 87) | (Code_4680 == 88):
        return 'Freezing Rain'
    elif (Code_4680 == 89) | (Code_4680 == 90):
        return 'Hail'
    else:
        return 'error: wrong Code Synop'

def Synop_to_Phase(Code_4680):
    """
    Convert the Synop Code from the manufacturer Parsivel output to a 
    corresponding precipitation phase.

    Parameters
    ----------
    Code_4680 : Int
        Synop Code from Parsivel.

    Returns
    -------
    str
        Name of the precipitation phase.

    """
    if Code_4680 == 0:
        return 'No precipitation'
    elif (Code_4680 < 66):
        return 'Liquid'
    elif (Code_4680 >= 66) & (Code_4680 < 70):
        return 'Mixture'
    elif (Code_4680 == 87) | (Code_4680 == 88):
        return 'Liquid'
    elif (Code_4680 >= 70) :
        return 'Solid'
    else:
        return 'error: wrong Code Synop'

def read_codesynop(filename, start_time, end_time, step_time=1):
    """
    Return the precipitation type and phase at each time step of a Parsivel record,
    based on the Synop Code returned by the Parsivel.

    Parameters
    ----------
    filename : DataFrame
        Parsivel file.
    start_time : str
        Time at format %Y-%m-d %h:%M:%s of the desired starting time.
    end_time : str
        Time at format %Y-%m-d %h:%M:%s of the desired ending time.
    step_time : Int, optional
        Step time of the desired integration in min. The default is 1.

    Returns
    -------
    data_Parsi : DataFrame
        Same Dataframe with added information of precipitation Type and Phase.

    """
    import pandas as pd
    import numpy as np
    
    data_Parsi = pd.read_csv(filename, skiprows=[2, 3], index_col=0, header=1,
                           parse_dates=[0], infer_datetime_format=True)
    data_Parsi.index.name = 'Time'
    data_Parsi.Code_4680 = data_Parsi.Code_4680.astype('float')
    
    istart_time = (data_Parsi.index.get_loc(start_time, method='nearest'))
    iend_time = data_Parsi.index.get_loc(end_time, method='nearest')
    if np.isscalar(istart_time) == False :
        istart_time = istart_time.item()
        iend_time = iend_time.item()
    print('start: '+ str(data_Parsi.index[istart_time]))
    #istart_time = np.argwhere(pd.to_datetime(data_Parsi.index, unit='s')==start_time)[0][0]
    print('end :'+ str(data_Parsi.index[iend_time]))
    #iend_time = np.argwhere(pd.to_datetime(data_Parsi.index, unit='s')==end_time)[0][0]
    
    data_Parsi = pd.DataFrame(data_Parsi[istart_time:iend_time].Code_4680)
    data_Parsi = data_Parsi.dropna(axis=0)
    data_Parsi['Precip_type'] = (np.vectorize(Synop_to_Type))(data_Parsi.Code_4680)
    data_Parsi['Precip_phase'] = (np.vectorize(Synop_to_Phase))(data_Parsi.Code_4680)
    return data_Parsi

def compute_drop_by_step(dsd, matrix_vd,
                        start_time, end_time, step_time= 1, rho_snow = None,
                        rho_solid = 900, rho_liquid = 1000, rho_mixture = 300):
    """
    Calculate the number and water equivalent (in kg) of each phase 
    at each desired time step (default= 1 min). Also returns the major phase
    according the number of particles and the mass of water. Also computes the
    proportion of each of those 3 phases (liquid, solid, mixture).
    Based on the raw values measured by the Parsivel: the distribution 
    of the particule size and velocity.

    Parameters
    ----------
    dsd : DSD
        DSD computed from the Pasivel spectrum file.
    matrix_vd : Dataframe
        Classification matricx of the Precipitation type.
    start_time : str
        Time at format %Y-%m-d %h:%M:%s of the desired starting time.
    end_time : str
        Time at format %Y-%m-d %h:%M:%s of the desired ending time.
    step_time : int, optional
        Step time of the desired integration in min. The default is 1.
    rho_snow : float, optional
        Snowflake density in kg/m3. If None, the Rees(2021) density 
        snowflake formula is used. The default is None.
    rho_solid : float, optional
        Ice (hail) density in kg/m3. The default is 900.
    rho_liquid : float, optional
        Liquid water (drizzle, rain) density in kg/m3. The default is 1000.
    rho_mixture : float, optional
        Mixture phase density in kg/m3. The default is 300.

    Returns
    -------
    df_matrix : DataFrame
        Dataframe with the number and the mass (in kg) of particles 
        for each phases at every time step (1 min by default). Also the
        major phase observed during the time step. Based on the raw values
        measured by the Parsivel: the distribution of the particule size and
        velocity.

    """
    import datetime
    from copy import copy, deepcopy
    import pandas as pd
    import numpy as np
    
    dsd_local = deepcopy(dsd)
    
    if pd.to_datetime(arg=start_time) < pd.to_datetime(arg=dsd_local.time["data"][0], unit='s'):
        start_time = str(pd.to_datetime(arg=dsd_local.time["data"][0], unit='s'))
    if pd.to_datetime(arg=end_time) > pd.to_datetime(arg=dsd_local.time["data"][-1], unit='s'):
        end_time = str(pd.to_datetime(arg=dsd_local.time["data"][-1], unit='s'))
        
    run_time = datetime.datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
    step_time = datetime.timedelta(minutes= step_time)
    print('start= '+str(start_time))
    print('end= '+str(end_time))
    it_time = 0
    while run_time < end_time:
        print(run_time)
        
        try : 
            irun_time = np.argwhere(pd.to_datetime(arg=dsd_local.time["data"][:], unit='s')==run_time)[0][0]
            #inext_time = np.argwhere(pd.to_datetime(arg=dsd_local.time["data"][:], unit='s')==run_time+step_time)[0][0]
            it_time += 1
            
            Nd_matrix = dsd_local.fields["drop_spectrum"]["data"][irun_time:irun_time+1,:,:].sum(axis=0)
            Nd_tot = Nd_matrix.sum()

            Nd_hail = np.ma.masked_where(matrix_vd.T != 1, Nd_matrix).sum()
            Nd_graupel = np.ma.masked_where(matrix_vd.T != 2, Nd_matrix).sum()
            Nd_rain = np.ma.masked_where(matrix_vd.T != 3, Nd_matrix).sum()
            Nd_drizzle = np.ma.masked_where(matrix_vd.T != 4, Nd_matrix).sum()
            Nd_mix = np.ma.masked_where(matrix_vd.T != 5, Nd_matrix).sum()
            Nd_snow = np.ma.masked_where(matrix_vd.T != 6, Nd_matrix).sum()
            
            # find the major phase according to drop numbers
            if np.argmax([Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mix,Nd_snow]) == 0 :
                Nd_majo_loffler = 'Hail'
            elif np.argmax([Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mix,Nd_snow]) == 1 :
                Nd_majo_loffler = 'Graupel'
            elif np.argmax([Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mix,Nd_snow]) == 2 :
                Nd_majo_loffler = 'Rain'
            elif np.argmax([Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mix,Nd_snow]) == 3 :
                Nd_majo_loffler = 'Drizzle'
            elif np.argmax([Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mix,Nd_snow]) == 4 :
                Nd_majo_loffler = 'Mixture'
            else :
                Nd_majo_loffler = 'Snow'
          
            #mask the non-solid phase (3,4,5), and sum the number of drops
            mask_solid = (matrix_vd.T == 3) |(matrix_vd.T == 4)| (matrix_vd.T == 5)
            mask_hail = (matrix_vd.T != 1)
            mask_snow = (matrix_vd.T != 6) & (matrix_vd.T != 2)
            Nd_solid = np.ma.masked_where(mask_solid, Nd_matrix).sum()
            #mask the non-liquid phase (1,2,6,5), and sum the nuber of drops
            mask_liquid = (matrix_vd.T != 3) & (matrix_vd.T != 4)
            Nd_liquid = np.ma.masked_where(mask_liquid, Nd_matrix).sum()
            #mask the non-mix phase (1,2,3,4,6), and sum the nuber of drops
            mask_mixture = (matrix_vd.T != 5)
            Nd_mixture = np.ma.masked_where(mask_mixture, Nd_matrix).sum()
            
            # find the major phase according to drop numbers
            if np.argmax([Nd_solid,Nd_liquid,Nd_mixture]) == 0 :
                Nd_majo = 'Solid'
            elif np.argmax([Nd_solid,Nd_liquid,Nd_mixture]) == 1 :
                Nd_majo = 'Liquid'
            else :
                Nd_majo = 'Mixture'
                
            # compute the proportion of drops of the major phase
            if (Nd_solid+Nd_liquid+Nd_mixture) != 0:
                Frac_Nd_majo = np.max([Nd_solid,Nd_liquid,Nd_mixture])/(Nd_solid+Nd_liquid+Nd_mixture)
            else :
                Frac_Nd_majo = np.nan
            
            #Volume of all drops in a phase: Number of drops in each class of diameter in the phase * diameter^3 * pi /6
            Vol_solid = (np.ma.masked_where(mask_solid, Nd_matrix).sum(axis=0)*(np.power(matrix_vd.T.columns.values,3))*np.pi/6).sum() * 1e-9 #in m^3
            Vol_snow = (np.ma.masked_where(mask_snow, Nd_matrix).sum(axis=0)*(np.power(matrix_vd.T.columns.values,3))*np.pi/6).sum() * 1e-9 #in m^3
            Vol_hail = (np.ma.masked_where(mask_hail, Nd_matrix).sum(axis=0)*(np.power(matrix_vd.T.columns.values,3))*np.pi/6).sum() * 1e-9 #in m^3
            Vol_liquid = (np.ma.masked_where(mask_liquid, Nd_matrix).sum(axis=0)*(np.power(matrix_vd.T.columns.values,3))*np.pi/6).sum() * 1e-9 #in m^3
            Vol_mixture = (np.ma.masked_where(mask_mixture, Nd_matrix).sum(axis=0)*(np.power(matrix_vd.T.columns.values,3))*np.pi/6).sum() * 1e-9 #in m^3
            
            Mass_liquid = rho_liquid * Vol_liquid
            Mass_hail = rho_solid * Vol_hail
            Mass_mixture = rho_mixture * Vol_mixture
            if rho_snow == None : #if no snow density provided, compute the snowflakes mass acoording to Rees (2021): M_snowflake(mg)=0.018*D(mm)^2.75
                Mass_snow = (0.018*np.ma.masked_where(mask_snow, Nd_matrix).sum(axis=0)*np.power(matrix_vd.T.columns.values,2.52)).sum() * 1e-6 #in kg
            else:
                Mass_snow = rho_snow * Vol_snow
            Mass_solid = Mass_snow + Mass_hail
            Mass_total = Mass_solid + Mass_liquid + Mass_mixture
            
            # find the major phase according to mass of drops
            if np.argmax([Mass_solid,Mass_liquid,Mass_mixture]) == 0 :
                Mass_majo = 'Solid'
            elif np.argmax([Mass_solid,Mass_liquid,Mass_mixture]) == 1 :
                Mass_majo = 'Liquid'
            else :
                Mass_majo = 'Mixture'
                
            # compute the proportion of drop mass of the major phase
            if Mass_total != 0:
                Frac_Mass_majo = np.max([Mass_solid,Mass_liquid,Mass_mixture])/Mass_total
            else :
                Frac_Mass_majo = np.nan

            df_matrix_temp = pd.DataFrame([[Nd_solid,Nd_liquid,Nd_mixture,
                                        Vol_solid,Vol_liquid,Vol_mixture,
                                        Mass_solid,Mass_liquid,Mass_mixture,Mass_total,
                                        Nd_majo, Nd_majo_loffler, Frac_Nd_majo,Mass_majo,Frac_Mass_majo]],
                                 columns=['Solid','Liquid','Mixture',
                                          'Vol_solid','Vol_liquid','Vol_mixture',
                                          'Mass_solid','Mass_liquid','Mass_mixture','Mass_total',
                                          'Nd_majo', 'Nd_majo_Loffler','Frac_Nd_majo','Mass_majo','Frac_Mass_majo'],
                                 index=[run_time])
            
            if it_time == 1 : #for the first good iteration, create the dataframe
                df_matrix = df_matrix_temp
            else :
                df_matrix = df_matrix.append(df_matrix_temp)
        except : 
            print(str(run_time)+' does not exist in specrtum')
            
        run_time += step_time
        
    df_matrix.index.name = 'Time'
    
    return df_matrix

def plot_dsd_vs_matrix(dsd, matrix_vd, 
                        spec_time, step_time= 1, classif = 'Type',
                        log_drop = False, xlim=(0,6), ylim=(0,15)):
    """
    Plot the raw distribution of the particles size and velocity at a specific
    time during a specified duration step_time in min(default 1 min).
    It can classify the precipitation by considering the types or the phase 
    (default is 'Type') of the particles.

    Parameters
    ----------
    dsd : DSD
        DSD computed from the Pasivel spectrum file..
    matrix_vd : DataFrame
        Classification matricx of the Precipitation type.
    start_time : str
        Time at format %Y-%m-d %h:%M:%s of the desired time to investigate.
    step_time : int, optional
        Time interval to consider to integrate the data in min. The default is 1.
    classif : str, optional
        'Type' or 'Phase', depending on the classification you want. 
        The default is 'Type'.
    log_drop : bool, optional
        Wheter to use or not a logarithmic scale for particle number. Useful
        if several minutes of precipitation are observed.
        The default is False.
    xlim : tuple, optional
        Set the particle size limits of the plot in mm. The default is (0,6).
    ylim : tuple, optional
        Set the particle velocity limits of the plot in m/s. The default is (0,15).

    Returns
    -------
    Plot the distribution of the particule size and velocity above the classification
    matrix (in type or phase) and plot the barplot of the contribution in each type/phase.

    """
    
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    from matplotlib import cm
    
    ispec_time = np.argwhere(pd.to_datetime(arg=dsd.time["data"][:], unit='s')==spec_time)[0][0]
    
    #compute drop numbers for each precip types
    Nd_matrix = dsd.fields["drop_spectrum"]["data"][ispec_time:ispec_time+step_time,:,:].sum(axis=0)
    Nd_tot = Nd_matrix.sum()
    Nd_hail = np.ma.masked_where(matrix_vd.T != 1, Nd_matrix).sum()
    Nd_graupel = np.ma.masked_where(matrix_vd.T != 2, Nd_matrix).sum()
    Nd_rain = np.ma.masked_where(matrix_vd.T != 3, Nd_matrix).sum()
    Nd_drizzle = np.ma.masked_where(matrix_vd.T != 4, Nd_matrix).sum()
    Nd_mixture = np.ma.masked_where(matrix_vd.T != 5, Nd_matrix).sum()
    Nd_snow = np.ma.masked_where(matrix_vd.T != 6, Nd_matrix).sum()
    
    Nd_liquid = Nd_drizzle + Nd_rain
    Nd_solid = Nd_snow + Nd_graupel + Nd_hail
    df_Nd = pd.DataFrame([[Nd_hail,Nd_graupel,Nd_rain,Nd_drizzle,Nd_mixture,Nd_snow]],
                         columns=['Hail','Graupel','Rain','Drizzle','Mixture','Snow'])
    df_phase = pd.DataFrame([[Nd_liquid,Nd_solid,Nd_mixture]],
                         columns=['Liquid','Solid','Mixture'])
    
    #plot1
    if classif == 'Type':
        plt.bar(df_Nd.columns,height=df_Nd.loc[0])
        plt.grid()
        plt.title('Number of drops by type\nfrom '+ spec_time+
                  ' during '+str(step_time)+'min'+'\n#drops = '+str(Nd_tot))
        plt.show()
    if classif == 'Phase' : 
        plt.figure(figsize=(4.5,4))
        plt.bar(df_phase.columns,height=df_phase.loc[0],color=['tab:red', 'tab:blue', 'grey'],
                alpha=0.3)
        plt.grid(ls='--', alpha=0.5)
        plt.ylabel('Particles number')
        plt.title('Number of drops by phase\nfrom '+ spec_time+
                  ' during '+str(step_time)+'min'+'\n#drops = '+str(Nd_tot))
        plt.tight_layout()
        # plt.savefig(path_save+'frac_phase_'+start_time.replace(':','')+'.png', 
        #             format='png',dpi=300)
        plt.show()
    
    #plot2
    if classif == 'Type':
        plt.figure(figsize=(8,5))
        pal = plt.cm.get_cmap("viridis_r").copy()
        pal.set_bad(color='white', alpha=0) #set nan or 0 colors as transparent
        #pal.colors[0]=[1,1,1]
        fig, ax = plt.subplots()
        
        precip_type = plt.pcolormesh(matrix_vd.index,matrix_vd.columns,matrix_vd.T, 
                       shading='auto',cmap=cm.get_cmap('Set2',6), alpha=0.4)
        cbar_type = plt.colorbar(precip_type, label='Precipitation type',
                                 ticks=np.arange(1.5-1/12,6,5/6), orientation='horizontal')
        cbar_type.ax.set_xticklabels(['Hail', 'Graupel', 'Rain', 'Drizzle', 'Mixture','Snow']) 
        
    if classif == 'Phase':
        matrix_phase = matrix_vd.replace(to_replace = [[3,4],[1,2,6],5], value = [0,1,2])
        pal = plt.cm.get_cmap("turbo").copy()
        pal.set_bad(color='white', alpha=0) #set nan or 0 colors as transparent
        #pal.colors[0]=[1,1,1]
        fig, ax = plt.subplots(facecolor= 'white', dpi=300, figsize=(5,4))
        
        precip_phase = ax.pcolormesh(matrix_phase.index,matrix_phase.columns,matrix_phase.T, 
                        shading= 'auto', alpha= 0.3,
                        cmap= mpl.colors.ListedColormap(['tab:red', 'tab:blue', 'grey']))
        cbar_type = plt.colorbar(precip_phase, label='Precipitation phase',
                                  ticks= np.arange(1/3, 2, 2/3), 
                                  orientation='horizontal', 
                                  pad=0.2, fraction=0.08
                                  )
        cbar_type.ax.set_xticklabels(['Liquid', 'Solid', 'Mixture']) 
        
    if log_drop == False: 
        im=ax.pcolormesh(np.insert(dsd.diameter["data"],0,0),
                          np.insert(dsd.spectrum_fall_velocity["data"],0,0),
                          np.ma.masked_values(dsd.fields["drop_spectrum"]["data"][ispec_time:ispec_time+step_time,:,:].sum(axis=0),0),
                          shading='auto', cmap=pal)
        color_label = 'Number of particles'
    else:
        im=ax.pcolormesh(np.insert(dsd.diameter["data"],0,0),
                          np.insert(dsd.spectrum_fall_velocity["data"],0,0),
                          np.log10(dsd.fields["drop_spectrum"]["data"][ispec_time:ispec_time+step_time,:,:].sum(axis=0)),
                          shading='auto', cmap=pal)
        color_label = 'log10(m$^{-3}$ mm$^{-1}$)'
    if xlim != None :
        ax.set_xlim(xlim)
    if ylim != None :
        ax.set_ylim(ylim)
    plt.xlabel('Particle size class [mm]')
    plt.ylabel('Fall velocity [m/s]')
    plt.title('Drop size-velocity distribution \nfrom '+ spec_time+' during '+str(step_time)+'min')
    plt.colorbar(im,label=color_label) 
        
    #ax.plot(dsd.diameter["data"],dsd.velocity["data"],'r-..',label='Velocity model (Atlas)')
    #plt.legend(loc='best')

    plt.tight_layout()
    # plt.savefig(path_save+'PVSD_big_'+start_time.replace(':','')+'.png', 
    #             format='png',dpi=300)
    plt.show()

def get_phase_mix(df_predicted, df_Parsivel,minMass=4e-6):
    """
    Merge the data of the major phase observed by the experimental classification
    matrix df_predicted with the data of the major phase observed by the Parsivel
    according the Synop Code.

    Parameters
    ----------
    df_predicted : DataFrame
        Data computed from the Parsivel DSD and the experimental phase 
        classification.
    df_Parsivel : DataFrame
        Data computed from the Parsivel Code Synop.
    minMass : float, optional
        Mass threshold in kg to define the absence of precipitation when calculating
        the phase from the spectrum data of Parsivel. e.g. below 4e-6 kg of water
        precipited during 1 min, it is considered that 'No precipitaion' is detected.
        The default is 4e-6.

    Returns
    -------
    df_compare : DataFrame
        Merged Dataframe with both data from the experimental classification
        and the manufacturer-based determination of the particle's phase.

    """
    
    import numpy as np
    import pandas as pd
    
    #Min Mass: filter the precip intensity lower than a threshold in kg/min

    #create a local dataframe, not mutable
    df_phase_predict = df_predicted.copy(deep=True)
    df_phase_predict[df_phase_predict['Mass_total'] < minMass] = np.nan
    df_compare = pd.merge(df_phase_predict, df_Parsivel, how='outer', on='Time', sort=True)
    df_compare.fillna({'Mass_majo':'No precipitation','Nd_majo':'No precipitation'}, inplace=True)
    df_compare.dropna(subset=['Precip_phase'],inplace=True)
    return df_compare


def plot_major_phase_contribution(df_compare, density=True):
    """
    Plot the cumulative density of the major pahse contribution according the DSD
    from the spectrum file of the Parsivel and the experimental matrix of phase
    classification.

    Parameters
    ----------
    df_compare : DataFrame
        Merged Dataframe with both data from the experimental classification
        and the manufacturer-based determination of the particle's phase.
    density : bool, optional
        Wether to plot the results in a cumulative density or not. The default is True.

    Returns
    -------
    Plots of the cumulative density of the major phase contribution according the 
    Number of particles or the mass equivalent.

    """
    
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sn
    
    minMass = 4e-6 #kg/1min to be modified if other value has been used for the thereshold

    plt.hist(df_compare[['Frac_Nd_majo','Frac_Mass_majo']], bins=np.linspace(0.3,1,15), 
             alpha=0.4,density=False, rwidth=0.9, cumulative=False, 
             label=['Particle number','Particle mass'])
    plt.title('Major phase contribution'+'\nMass intensity threshold = '+str(minMass)+' kg/min')
    plt.ylabel('# records')
    plt.xlabel('Major phase fraction')
    plt.legend()
    plt.show()
    
    if density:
        plt.hist(df_compare['Frac_Mass_majo'],bins=np.linspace(0.3,1,15), 
                 alpha=0.5,density=True, rwidth=0.9, cumulative=True, 
                 histtype='bar', label='Particle mass')
        plt.hist(df_compare['Frac_Nd_majo'], bins=np.linspace(0.3,1,15), 
                 alpha=0.5,density=True, rwidth=0.9, cumulative=True, 
                 histtype='bar', label='Particle number')
        plt.title('Cumulative distribution of major phase contribution'+'\nMass intensity threshold = '+str(minMass)+' kg/min')
        plt.ylabel('density')
    else :
        plt.hist(df_compare['Frac_Nd_majo'], bins=np.linspace(0.3,1,15), 
                 alpha=0.4,density=False, rwidth=0.9, cumulative=False, label='Particle number')
        plt.hist(df_compare['Frac_Mass_majo'],bins=np.linspace(0.3,1,15), 
                 alpha=0.4,density=False, rwidth=0.9, cumulative=False, label='Particle mass')
        plt.title('Major phase contribution'+'\nMass intensity threshold = '+str(minMass)+' kg/min')
        plt.ylabel('# records')
    plt.xlabel('Major phase fraction')
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(5,4))
    sn.ecdfplot(df_compare['Frac_Mass_majo'], label= 'Particle mass',stat='count')
    sn.ecdfplot(df_compare['Frac_Nd_majo'], label='Particle number',stat='count')
    plt.title('Cumulative distribution of major phase contribution'+'\nMass intensity threshold = '+str(minMass)+' kg/min')
    plt.xlabel('Major phase contribution')
    plt.ylabel('Cumulative counts of records')
    plt.text(0.7,0.9,'#records= '+str(df_compare['Frac_Mass_majo'].size))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(5,4))
    sn.ecdfplot(df_compare['Frac_Mass_majo'], label= 'Particle mass',stat='proportion')
    sn.ecdfplot(df_compare['Frac_Nd_majo'], label='Particle number',stat='proportion')
    plt.title('Cumulative distribution of major phase contribution'+'\nMass intensity threshold = '+str(minMass)+' kg/min')
    plt.xlabel('Major phase contribution')
    plt.ylabel('Cumulative density of records')
    plt.text(0.7,0.9,'#records= '+str(df_compare['Frac_Mass_majo'].size))
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()

def get_confusion_matrix(df_compare, var='Npart'):
    """
    Min Mass = filter the precip intensity lower than a threshold in kg/min
    PrecipType = "solid",'Liqud','Mixture' --> the phase wanted to be compared between the Parsivel output and the experimental matrix 
    """
    
    import numpy as np
    import pandas as pd
    from sklearn import metrics
    from sklearn.metrics import confusion_matrix
    
    if var == 'Npart' :

        conf_mat = confusion_matrix(df_compare['Precip_phase'],df_compare['Nd_majo'],
                                    labels=['Liquid','Solid','Mixture','No precipitation'])
        print(conf_mat)
        
        Accuracy = np.round(metrics.accuracy_score(df_compare['Precip_phase'], 
                                      df_compare['Nd_majo']),3)
        Kappa = np.round(metrics.cohen_kappa_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo']),3)
        Precision_w = np.round(metrics.precision_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average='weighted', 
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        Recall_w = np.round(metrics.recall_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average='weighted',
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        F1_score_w = np.round(metrics.f1_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average='weighted',
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        Precision_l, Precision_s = np.round(metrics.precision_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        Recall_l, Recall_s = np.round(metrics.recall_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        F1_score_l, F1_score_s = np.round(metrics.f1_score(df_compare['Precip_phase'], 
                                     df_compare['Nd_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        Test_score = (conf_mat[0,0]+conf_mat[1,1]) / conf_mat.sum()
        vec_score = pd.DataFrame(data=[[Accuracy,Kappa,Precision_w,Recall_w,F1_score_w,
                                        Precision_l, Precision_s,Recall_l,Recall_s,
                                        F1_score_l, F1_score_s,Test_score]]
                                 ,columns= ['Accuracy','Kappa','Precision_w','Recall_w','F1_score_w',
                                            'Precision_l', 'Precision_s','Recall_l','Recall_s',
                                            'F1_score_l', 'F1_score_s','Test_score'])
    
    if var == 'Masspart' :
        conf_mat = confusion_matrix(df_compare['Precip_phase'],df_compare['Mass_majo'],
                                    labels=['Liquid','Solid','Mixture','No precipitation'])
        print(conf_mat)
        Accuracy = np.round(metrics.accuracy_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo']),3)
        Kappa = np.round(metrics.cohen_kappa_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo']),3)
        Precision_w = np.round(metrics.precision_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average='weighted', 
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        Recall_w = np.round(metrics.recall_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average='weighted',
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        F1_score_w = np.round(metrics.f1_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average='weighted',
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)
        Precision_l, Precision_s = np.round(metrics.precision_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        Recall_l, Recall_s = np.round(metrics.recall_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        F1_score_l, F1_score_s = np.round(metrics.f1_score(df_compare['Precip_phase'], 
                                     df_compare['Mass_majo'], average=None,
                                     labels=['Liquid','Solid','Mixture','No precipitation']),3)[0:2]
        Test_score = (conf_mat[0,0]+conf_mat[1,1]) / conf_mat.sum()
        print(conf_mat[0,0])
        print(conf_mat[1,1])
        vec_score = pd.DataFrame(data=[[Accuracy,Kappa,Precision_w,Recall_w,F1_score_w,
                                        Precision_l, Precision_s,Recall_l,Recall_s,
                                        F1_score_l, F1_score_s,Test_score]]
                                 ,columns= ['Accuracy','Kappa','Precision_w','Recall_w','F1_score_w',
                                            'Precision_l', 'Precision_s','Recall_l','Recall_s',
                                            'F1_score_l', 'F1_score_s','Test_score'])
        
    return conf_mat, vec_score

def plot_confusion_matrix(mat_confusion):
    
    import matplotlib.pyplot as plt
    import seaborn as sn
    from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
    
    fig, ax = plt.subplots(facecolor='white', dpi=300, figsize = (6,5))
    sn.heatmap(mat_confusion, annot=True, cmap='Blues', fmt='.4g',square=True, robust=True,
               xticklabels=['Liquid','Solid','Mixture','No precipitation'],
               yticklabels=['Liquid','Solid','Mixture','No precipitation'],
               cbar_kws={'label':'Particle number'
               #, 'shrink':0.75
               #, 'extend':'max', 'format':'%.2g',
               }, 
               ax=ax
               #, vmax = 5e4
               ) 
    # plt.title(label='Confusion matrix - '+station_name+' '+event_date+
    #           '\nMass intensity threshold = '+str(mass_thrshld)+' kg/min - '+
    #           r'$\rho_{snow}$='+str('Rees'))
    # plt.title(label='Confusion matrix - '+station_name+' '+event_date+
    #           '\nMass intensity threshold = '+str(mass_thrshld)+' kg/min - '+
    #           r'$\rho_{snow}$='+str(rho_snow)+' kg/m3')
    plt.xlabel('Experimental-based\n predictions')
    plt.ylabel('Manufacturer-based\n observations')
    plt.xticks(rotation=90)
    #ax.xaxis.grid(False, which='minor')
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
    ax.tick_params(which='major')
    plt.yticks(rotation=0)
    plt.tight_layout()
    #plt.savefig(path_save+'conf_mat_Cham_snow50.png',format='png', dpi=300)
    plt.show()

