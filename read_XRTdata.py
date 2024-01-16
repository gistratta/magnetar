import pandas as pd
import numpy as np
import os

"""

Reads 0.3 - 10 keV unabsorbed flux XRT light curves and photon indexes
in the Swift XRT Repository format

INPUT: GRB name ID

OUTPUT:
tx = time
dtxp,dtxm = time bin
fx = 0.3-10 keV integrated flux
dfxp,dfxm = integrated flux errors
FnuxmJy = flux density in mJyansky
dFnuxp_mJy,dFnuxm_mJy = flux density errors
gammaxrt = Photon Index
dgammaxrtp,dgammaxrtm = Photon index errors

"""


def read_data(grb):

    # --- Directory where the light curve data files are stored
    outpath = './data/'

    # --- Check if GRB light curve is present in the directory
    if grb[-1] == 'A':
        print('GRB',grb,'ends with A: check if it is quoted without A....')
        if not os.path.isfile(outpath+grb+'.txt'):
            if not os.path.isfile(outpath+grb[:-1]+'.txt'):
                print('GRB', grb,'is not present in database')
            else:
                grb=grb[:-1]
        else:
            print('GRB',grb,' is quoted correctly')

    # --- read the file
    db=pd.read_csv(outpath+grb+'.txt',skiprows=[0,1],header=None,sep='\s+')

    # --- Find index where data block starts after separation line (contining NO)
    idx = db[db[0]=='NO'].index.values.astype(int)
    idx_name = db[db[0]=='!'].index

    #--- data up to the first block index (wt or pc)
    if 'flux' in db.iloc[idx[0]+1][1]:
        db01gamma=db[0:idx[0]]
    else:
        db01=db[0:idx[0]]
    if len(idx) > 1:
        if 'flux' in db.iloc[idx[0]+1][1]:
            db01=db[idx[0]+2:idx[1]]
        else:
            db01gamma=db[idx[0]+2:idx[1]]
    else:
        if 'flux' in db.iloc[idx[0]+1][1]:
            db01=db[idx[0]+2:]
        else:
            db01gamma=db[idx[0]+2:]
    db0=db01
    db0gamma=db01gamma
    if len(idx_name) > 1:
        if 'wt' in db[1].iloc[idx_name[1]]:
            print('reading wt data')
            # if wt then first block is slew
            db02=db[idx[1]+2:idx[2]]
            db02gamma=db[idx[2]+2:idx[3]]
            db03=db[idx[3]+2:idx[4]]
            db03gamma=db[idx[4]+2:]
            db0=pd.concat([db01,db02,db03])
            db0gamma=pd.concat([db01gamma,db02gamma,db03gamma])
        elif 'xrtpcflux' in db[1].iloc[idx_name[1]]:
            # if pcflux then first block is wt
            db02=db[idx[1]+2:idx[2]]
            db02gamma=db[idx[2]+2:]
            db0=pd.concat([db01,db02])
            db0gamma=pd.concat([db01gamma,db02gamma])
    else:
        print('only pc data')
    tx=np.array([float(x) for x in db0[0]])
    dtxp=np.array([float(x) for x in db0[1]])
    dtxm=np.array([float(x) for x in db0[2]])
    fx=np.array([float(x) for x in db0[3]])
    dfxp=np.array([float(x) for x in db0[4]])
    dfxm=np.array([float(x) for x in db0[5]])
    gammaxrt=np.array([float(x) for x in db0gamma[3]])
    dgammaxrtp=np.array([float(x) for x in db0gamma[4]])
    dgammaxrtm=np.array([float(x) for x in db0gamma[5]])

    # --- Compute flux density in mJy at 1 keV by assuming a power law spectrum with index beta
    beta = gammaxrt-1
    dbetap = dgammaxrtp-1
    dbetam = dgammaxrtm-1
    keV2Hz = 2.4180e17
    nu1 = 0.3*keV2Hz
    nu2 = 10.0*keV2Hz
    nu0 = 1*keV2Hz
    K = fx*(nu0**(-beta))/((nu2**(1-beta) - nu1**(1-beta))/(1-beta) )
    FnuxmJy = 1e26*K
    Kp = (fx+dfxp)/( (nu2**(1-beta) - nu1**(1-beta))/(1-beta) )
    Km = (fx+dfxm)/( (nu2**(1-beta) - nu1**(1-beta))/(1-beta) )
    Fnuxp = 1e26*Kp*(nu0)**(-beta)
    Fnuxm = 1e26*Km*(nu0)**(-beta)
    dFnuxp_mJy = Fnuxp - FnuxmJy
    dFnuxm_mJy = FnuxmJy - Fnuxm

    return tx,dtxp,dtxm,fx,dfxp,dfxm,FnuxmJy,dFnuxp_mJy,dFnuxm_mJy,gammaxrt,dgammaxrtp,dgammaxrtm
