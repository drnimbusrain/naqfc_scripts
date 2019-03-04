#!/data/aqf2/barryb/anaconda2/envs/patrick_monet/bin/python

__author__  = 'Patrick Campbell'
__email__   = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'



#Simple MONET utility to calculate statistics from paired hdf file

import os
from glob import glob
import sys
sys.path.append('/data/aqf/patrickc/MONET/')

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave
from monet.util.mystats import NO,NP,NOP,MO,MP,MdnO,MdnP,STDO,STDP,MB,NMB,NME_m,RMSE,IOA_m,R2 
import pandas as pd
import numpy as np
from numpy import sqrt

#Define all statistics and statistical plots desired

def  calc_r2(obs,mod):
     """ Coefficient of Determination (unit squared) """
     return R2(obs,mod,axis=None)

def  calc_ioa(obs,mod):
     """ Index of Agreement """
     return IOA_m(obs,mod,axis=0)

def  calc_rmse(obs,mod):
     """ Root  Mean Square Error """
     return RMSE(obs,mod,axis=0)

def  calc_nme(obs,mod):
     """ Normalized Mean Error (%)"""
     return NME_m(obs,mod,axis=0)

def  calc_nmb(obs,mod):
     """ Normalized Mean Bias (%)"""
     return NMB(obs,mod,axis=0)

def  calc_mb(obs,mod):
     """ Mean Bias """
     return MB(obs,mod,axis=0)

def  calc_stdp(obs,mod):
     """ Standard deviation of Predictions """
     return STDP(obs,mod,axis=0)

def  calc_stdo(obs,mod):
     """ Standard deviation of Observations """
     return STDO(obs,mod,axis=0)

def  calc_MdnP(obs,mod):
     """ Median Predictions (obs unit) """
     return MdnP(obs,mod,axis=0)

def  calc_MdnO(obs,mod):
     """ Median Observations (model unit) """
     return MdnO(obs,mod,axis=0)

def  calc_MP(obs,mod):
     """ Mean Predictions (model unit) """
     return MP(obs,mod,axis=0)

def  calc_MO(obs,mod):
     """ Mean Observations (obs unit) """
     return MO(obs,mod,axis=0)

def  calc_NOP(obs,mod):
     """ N Observations/Prediction Pairs (#)"""
     return NOP(obs,mod,axis=0)

def  calc_NP(obs,mod):
     """ N Predictions (#) """
     return NP(obs,mod,axis=0)

def  calc_NO(obs,mod):
     """ N Observations (#) """
     return NO(obs,mod,axis=0)

def  make_24hr_regulatory(df,col=None):
     """ Make 24-hour averages """
     return calc_24hr_ave(df,col)

def  make_8hr_regulatory(df,col=None):
     """ Make 8-hour rolling average daily """
     return calc_8hr_rolling_max(df,col,window=8)


if __name__ == '__main__':

    parser = ArgumentParser(description='calculates statistics from paired file', formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-f',   '--files',       help='string input paired file directory/names', type=str, required=True)
    parser.add_argument('-r',   '--regulatory',  help='boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs', type=bool, required=False, default=False)
    parser.add_argument('-sd',  '--startdate',   help='string start date to isolate periods for statistics YYYY-MM-DD HH:MM:SS', type=str, required=False, default=None)
    parser.add_argument('-ed',  '--enddate',     help='string end date to isolate periods for statistics YYYY-MM-DD HH:MM:SS', type=str, required=False, default=None)
    parser.add_argument('-n',   '--networks',    help='string/list input for networks to do stats',type=str,nargs='+', required=False, default={'airnow'})
    parser.add_argument('-m',   '--models',      help='string/list input models: cmaq, fv3, hysplit (not-ready), or camx (not-ready)', type=str,nargs='+', required=False, default={'cmaq'}) 
    parser.add_argument('-s',   '--species',     help='string/list input for obs species-variables to create stats',type=str,nargs='+', required=False, default={'OZONE','PM2.5'})
    parser.add_argument('-b',   '--subset_epa',  help='boolean set to True for subsetting by U.S. EPA region', type=bool, required=False, default=False)
    parser.add_argument('-e',   '--epa_regions', help='string/list input for set U.S. EPA regions',type=str,nargs='+', required=False, default={'R1'})
    parser.add_argument('-v',   '--verbose',     help='print debugging information', action='store_true', required=False)
    args = parser.parse_args()

    finput       = args.files
    reg          = args.regulatory
    networks     = args.networks
    startdate    = args.startdate
    enddate      = args.enddate
    models       = args.models
    species      = args.species
    subset_epa   = args.subset_epa
    epa_regions  = args.epa_regions
    verbose      = args.verbose

    for ee in epa_regions:
     for xx in models:
#reads paired network
      for ii in networks:
       if ii == 'airnow' and xx == 'cmaq':
          df = pd.read_hdf(finput)
          mapping_table = {'OZONE':'O3', 'PM2.5':'PM25_TOT', 'PM10':'PMC_TOT', 'CO':'CO', 'NO':'NO', 'NO2':'NO2', 'SO2':'SO2','NOX':'NOX','NO2Y':'NOY'}
          sub_map = {i: mapping_table[i] for i in species if i in mapping_table}   
#subsetting data for dates, regulatory calc, and/or epa regions    
       if startdate != None and enddate != None:
          mask = (df['time'] >= startdate) & (df['time'] <= enddate)
          df =df.loc[mask]
          import datetime
          startdatename_obj = datetime.datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S')
          enddatename_obj   = datetime.datetime.strptime(enddate, '%Y-%m-%d %H:%M:%S')
          startdatename = str(datetime.datetime.strftime(startdatename_obj,'%Y-%m-%d_%H'))
          enddatename = str(datetime.datetime.strftime(enddatename_obj,'%Y-%m-%d_%H'))
       else:
          startdatename='Entire'
          enddatename  ='Period'

       if subset_epa is True:
          df.query('epa_region == '+'"'+ee+'"',inplace=True)
       if reg is True and subset_epa is False:
          stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_domain.txt','w')
       elif reg is True and subset_epa is True:
          stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_'+ee+'.txt','w')
       elif reg is False and subset_epa is True:
          stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_'+ee+'.txt','w')
       else:
          stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_domain.txt','w')          
#Converts OZONE, PM10, or PM2.5 dataframe to NAAQS regulatory values
       for jj in species: 
        df_replace = df.replace(0.0,np.nan) #Replace all values with exactly 0.0 (non-physical)
        df_drop=df_replace.dropna(subset=[jj,sub_map.get(jj)]) #Drops all rows with obs species = NaN        
       
        if jj == 'OZONE' and reg is True:
       	 df2 = make_8hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)}) 
        elif jj == 'PM2.5' and reg is True:
       	 df2 = make_24hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
        elif jj == 'PM10' and reg is True:
       	 df2 = make_24hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
        else:
         df2=df_drop  
#Calculates average statistics over entire file time
        if reg is True and subset_epa is False:
       	 stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_domain.txt','a')
        elif reg is True and subset_epa is True:
         stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_'+ee+'.txt','a')
        elif reg is False and subset_epa is True:
         stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_'+ee+'.txt','a')
        else:
       	 stats=open(finput.replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_domain.txt','a')
         
        print('---------------------------------')
        print('Statistics of ',jj,' and ',sub_map.get(jj),' pair over ', startdatename,' to ', enddatename)
        stats.write('Statistics of '+jj+' and '+sub_map.get(jj)+' pair over file period '+startdatename+ ' to '+enddatename+'\n')
        stats.write('---------------------------------'+'\n')

        obs_stats = df2[jj]
        mod_stats = df2[sub_map.get(jj)]
        
        no_stats=calc_NO(obs_stats,mod_stats)
        print('Number of',jj,' Observations = ', no_stats)
        stats.write('Number of '+jj+' Observations = '+str(no_stats)+'\n')
        
        np_stats=calc_NP(obs_stats,mod_stats)
        print('Number of',sub_map.get(jj),' Predictions = ', np_stats)
        stats.write('Number of '+sub_map.get(jj)+' Predictions = '+str(np_stats)+'\n')        

        nop_stats=calc_NOP(obs_stats,mod_stats)
        print('Number of',jj,'/',sub_map.get(jj),' Observations/Prediction Pairs (#) = ', nop_stats)
        stats.write('Number of '+jj+'/'+sub_map.get(jj)+' Observations/Prediction Pairs (#) = '+str(nop_stats)+'\n')        

        mo_stats=calc_MO(obs_stats,mod_stats)
        print('Mean of',jj,' Observations = ', "{:8.2f}".format(mo_stats))
        stats.write('Mean of '+jj+' Observations = '+str("{:8.2f}".format(mo_stats))+'\n') 
       
        mp_stats=calc_MP(obs_stats,mod_stats)
        print('Mean of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(mp_stats))
        stats.write('Mean of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(mp_stats))+'\n')       
 
        mdno_stats=calc_MdnO(obs_stats,mod_stats)
        print('Median of',jj,' Observations = ', "{:8.2f}".format(mdno_stats))
        stats.write('Median of '+jj+' Observations = '+str("{:8.2f}".format(mdno_stats))+'\n')        

        mdnp_stats=calc_MdnP(obs_stats,mod_stats)
        print('Median of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(mdnp_stats))
        stats.write('Median of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(mdnp_stats))+'\n')

        stdo_stats=calc_stdo(obs_stats,mod_stats)
        print('Standard deviation of',jj,' Observations = ', "{:8.2f}".format(stdo_stats))
        stats.write('Standard deviation of '+jj+' Observations = '+str("{:8.2f}".format(stdo_stats))+'\n')
       
        stdp_stats=calc_stdp(obs_stats,mod_stats)
        print('Standard deviation of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(stdp_stats))
        stats.write('Standard deviation of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(stdp_stats))+'\n')       

        mb_stats=calc_mb(obs_stats,mod_stats)
        print('Mean Bias of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(mb_stats))
        stats.write('Mean Bias of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(mb_stats))+'\n')        

        nmb_stats=calc_nmb(obs_stats,mod_stats)
        print('Normalized Mean Bias (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nmb_stats))        
        stats.write('Normalized Mean Bias (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nmb_stats))+'\n')

        nme_stats=calc_nme(obs_stats,mod_stats)
        print('Normalized Mean Error (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nme_stats))
        stats.write('Normalized Mean Error (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nme_stats))+'\n')

        rmse_stats=calc_rmse(obs_stats,mod_stats)
        print('Root Mean Square Error of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(rmse_stats))
        stats.write('Root Mean Square Error of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(rmse_stats))+'\n')

        ioa_stats=calc_ioa(obs_stats,mod_stats)
        print('Index of Agreement of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(ioa_stats))
        stats.write('Index of Agreement of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(ioa_stats))+'\n')

        r_stats=sqrt(calc_r2(obs_stats,mod_stats))
        print('Pearsons Correlation Coefficient of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(r_stats))
        stats.write('Pearsons Correlation Coefficient of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(r_stats))+'\n')

        print('Statistics done!')
        print('---------------------------------')
        stats.write('---------------------------------'+'\n')
        stats.close()







    sys.exit(0)
    


