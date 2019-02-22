#!/data/aqf2/barryb/anaconda2/envs/patrick_monet/bin/python

__author__  = 'Patrick Campbell'
__email__   = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'



#Simple MONET utility to calculate statistics from paired hdf file

import os
from glob import glob
import sys
sys.path.append('/data/aqf/patrickc/MONET/')
#os.chdir('/data/aqf/patrickc/MONET/scripts/')

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave
from monet.util.mystats import NO,NP,NOP,MO,MP,MdnO,MdnP,STDO,STDP,MB,NMB,NME,RMSE,IOA,R2 
#from monet.plots import TaylorDiagram
import pandas as pd
import numpy as np
from numpy import sqrt

#Define all statistics and statistical plots desired

#def  plot_taylor(finput,verbose=False):
#     dset=monet.models.cmaq.open_mfdataset(finput)
#     return dset

def  calc_r2(obs,mod):
     """ Coefficient of Determination (unit squared) """
     return R2(obs,mod,axis=None)

def  calc_ioa(obs,mod):
     """ Index of Agreement """
     return IOA(obs,mod,axis=0)

def  calc_rmse(obs,mod):
     """ Root  Mean Square Error """
     return RMSE(obs,mod,axis=0)

def  calc_nme(obs,mod):
     """ Normalized Mean Error (%)"""
     return NME(obs,mod,axis=0)

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
    
    parser.add_argument('-f', '--files',       help='string input paired file directory/names', type=str, required=True)
    parser.add_argument('-r', '--regulatory',  help='boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs', type=bool, required=False, default=False)
    parser.add_argument('-n', '--networks',    help='string input for networks to do stats',type=str,nargs='+', required=False, default=['airnow'])
    parser.add_argument('-m', '--models',      help='string input models: cmaq, fv3, hysplit (not-ready), or camx (not-ready)', type=str,nargs='+', required=False, default=['cmaq']) 
    parser.add_argument('-s', '--species',     help='string input for obs species-variables to create stats',type=str,nargs='+', required=False, default=['OZONE','PM2.5'])
    parser.add_argument('-v', '--verbose',     help='print debugging information', action='store_true', required=False)
    args = parser.parse_args()

    finput       = args.files
    reg          = args.regulatory
    networks     = args.networks
    models       = args.models
    species      = args.species
    verbose      = args.verbose


    for xx in models:
#reads paired network
     for ii in networks:
      if ii == 'airnow' and xx == 'cmaq':
         df = pd.read_hdf(finput+'.hdf')
         mapping_table = {'OZONE':'O3', 'PM2.5':'PM25_TOT', 'PM10':'PMC_TOT', 'CO':'CO', 'NO':'NO', 'NO2':'NO2', 'SO2':'SO2','NOX':'NOX','NO2Y':'NOY'}
         sub_map = {i: mapping_table[i] for i in species if i in mapping_table}   
         
      if reg is True:
         stats=open(finput+'_reg_stats_domain.txt','w')
      else:
         stats=open(finput+'_stats_domain.txt','w')
          
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
        
#Calculates domain-wide average statistics over entire file
       if reg is True:
       	stats=open(finput+'_reg_stats_domain.txt','a')
       else:
       	stats=open(finput+'_stats_domain.txt','a')
         
       print('---------------------------------')
       print('Domain-wide statistics of ',jj,' and ',sub_map.get(jj),' pair over file period')
       stats.write('Domain-wide statistics of '+jj+' and '+sub_map.get(jj)+' pair over file period'+'\n')
       stats.write('---------------------------------'+'\n')

       obs_domain = df2[jj]
       mod_domain = df2[sub_map.get(jj)]
        
       no_domain=calc_NO(obs_domain,mod_domain)
       print('Number of',jj,' Observations = ', no_domain)
       stats.write('Number of '+jj+' Observations = '+str(no_domain)+'\n')
        
       np_domain=calc_NP(obs_domain,mod_domain)
       print('Number of',sub_map.get(jj),' Predictions = ', np_domain)
       stats.write('Number of '+sub_map.get(jj)+' Predictions = '+str(np_domain)+'\n')        

       nop_domain=calc_NOP(obs_domain,mod_domain)
       print('Number of',jj,'/',sub_map.get(jj),' Observations/Prediction Pairs (#) = ', nop_domain)
       stats.write('Number of '+jj+'/'+sub_map.get(jj)+' Observations/Prediction Pairs (#) = '+str(nop_domain)+'\n')        

       mo_domain=calc_MO(obs_domain,mod_domain)
       print('Mean of',jj,' Observations = ', mo_domain)
       stats.write('Mean of '+jj+' Observations = '+str(mo_domain)+'\n') 
       
       mp_domain=calc_MP(obs_domain,mod_domain)
       print('Mean of',sub_map.get(jj),' Predictions = ', mp_domain)
       stats.write('Mean of '+sub_map.get(jj)+' Predictions = '+str(mp_domain)+'\n')       
 
       mdno_domain=calc_MdnO(obs_domain,mod_domain)
       print('Median of',jj,' Observations = ', mdno_domain)
       stats.write('Median of '+jj+' Observations = '+str(mdno_domain)+'\n')        

       mdnp_domain=calc_MdnP(obs_domain,mod_domain)
       print('Median of',sub_map.get(jj),' Predictions = ', mdnp_domain)
       stats.write('Median of '+sub_map.get(jj)+' Predictions = '+str(mdnp_domain)+'\n')

       stdo_domain=calc_stdo(obs_domain,mod_domain)
       print('Standard deviation of',jj,' Observations = ', stdo_domain)
       stats.write('Standard deviation of '+jj+' Observations = '+str(stdo_domain)+'\n')
       
       stdp_domain=calc_stdp(obs_domain,mod_domain)
       print('Standard deviation of',sub_map.get(jj),' Predictions = ', stdp_domain)
       stats.write('Standard deviation of '+sub_map.get(jj)+' Predictions = '+str(stdp_domain)+'\n')       

       mb_domain=calc_mb(obs_domain,mod_domain)
       print('Mean Bias of ',sub_map.get(jj),'-',jj,' =  ', mb_domain)
       stats.write('Mean Bias of '+sub_map.get(jj)+'-'+jj+' =  '+str(mb_domain)+'\n')        

       nmb_domain=calc_nmb(obs_domain,mod_domain)
       print('Normalized Mean Bias (%) of ',sub_map.get(jj),'-',jj,' =  ', nmb_domain)        
       stats.write('Normalized Mean Bias (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str(nmb_domain)+'\n')

       nme_domain=calc_nme(obs_domain,mod_domain)
       print('Normalized Mean Error (%) of ',sub_map.get(jj),'-',jj,' =  ', nme_domain)
       stats.write('Normalized Mean Error (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str(nme_domain)+'\n')

       rmse_domain=calc_rmse(obs_domain,mod_domain)
       print('Root Mean Square Error of ',sub_map.get(jj),'-',jj,' =  ', rmse_domain)
       stats.write('Root Mean Square Error of '+sub_map.get(jj)+'-'+jj+' =  '+str(rmse_domain)+'\n')

       ioa_domain=calc_ioa(obs_domain,mod_domain)
       print('Index of Agreement of ',sub_map.get(jj),'-',jj,' =  ', ioa_domain)
       stats.write('Index of Agreement of '+sub_map.get(jj)+'-'+jj+' =  '+str(ioa_domain)+'\n')

       r_domain=sqrt(calc_r2(obs_domain,mod_domain))
       print('Pearsons Correlation Coefficient of ',sub_map.get(jj),'-',jj,' =  ', r_domain)
       stats.write('Pearsons Correlation Coefficient of '+sub_map.get(jj)+'-'+jj+' =  '+str(r_domain)+'\n')

       print('Statistics done!')
       print('---------------------------------')
       stats.write('---------------------------------'+'\n')
       stats.close()







    sys.exit(0)
    


