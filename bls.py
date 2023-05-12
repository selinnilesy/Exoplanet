import lightkurve as lk
import astropy

from astropy.io import fits
import csv
import sys
from timeit import default_timer as timer
from sklearn.preprocessing import normalize

import numpy as np

def weightSoma(fluxErr):
    return np.power(np.sum(np.power(fluxErr,-2)),-1)

def rValue(i1,i2,weight):
    return np.sum(np.fromiter((weight[i] for i in range(i1,i2)), dtype=float))

def sValue(i1,i2,weight,fluxo):
    return np.sum(np.fromiter((weight[i]*fluxo[i] for i in range(i1,i2)), dtype=float))

def dValue(lisWeight,fluxo,r,s):
    soma = np.sum(np.fromiter((lisWeight[i]*fluxo[i]**2 for i in range(len(lisWeight))), dtype=float))
    return (soma - (s**2)/(r*(1-r)))

def myBls(flux, fluxerr):
  lisWeight = []
  somaW = weightSoma(fluxerr)

  for i1 in range(len(flux)):
    wi = somaW*(fluxerr[i1]**(-2))
    lisWeight.append(wi)

  min_d = sys.float_info.max
  min_i1 = -1
  min_i2 = -1
  for i1 in range(0,len(flux)):
    for i2 in range(i1+1,len(flux)):
      r = rValue(i1,i2,lisWeight)
      s = sValue(i1,i2,lisWeight,flux)
      d = dValue(lisWeight,flux,r,s)
      if d<min_d :
        min_d=d
        min_i1=i1
        min_i2=i2
        #print(f"min found: s {s}, r = {r} , d = {d}")
  print(f"FINAL D: {d}, with i1: {min_i1} , i2: {min_i2} ")
  sys.stdout.flush()
  return min_i1,min_i2

def main():
  search_result = lk.search_lightcurve('Kepler-37', author='Kepler', cadence='long')
  # Download all available Kepler light curves
  lc_collection = search_result.download_all()

  for x in range(len(lc_collection)) :
    dataf = lc_collection[x]
    #dataf.info()
    #cols = dataf[1].columns
    print(f"Number of Lines in File: {len(dataf)}")
    sys.stdout.flush()
    #cols.info()
    time = np.asarray(dataf['time'])
    flux = np.asarray(dataf['sap_flux'])
    fluxerr = np.asarray(dataf['sap_flux_err'])
    #print(f"datafile: {d}, with i1: {min_i1} , i2: {min_i2} ")

    #clean nans
    flux = flux[np.logical_not(np.isnan(flux))]
    fluxerr = fluxerr[np.logical_not(np.isnan(fluxerr))]
    #normalize
    flux = flux / np.linalg.norm(flux)
    fluxerr = fluxerr / np.linalg.norm(fluxerr)
    #clean 0.0s
    #flux = flux[np.logical_not(flux == 10^-13)]
    #fluxerr = fluxerr[np.logical_not(fluxerr == 10^-13)]
    print(f"Remaining Lines in File after Cleansing: {len(flux)}")
    sys.stdout.flush()

    start = timer()
    i1,i2 = myBls(flux,fluxerr )
    end = timer()
    print(f"Execution Time: {end-start} with Period: {time[i2]-time[i1]}")
    sys.stdout.flush()
    #break
main()
