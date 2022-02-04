import pandas as pd
import numpy as np

def read_data(grb):
    outpath='./data/'
    db=pd.read_csv(outpath+grb+'.txt',skiprows=[0,1],header=None,sep='\s+')
    # find index where each block starts
    idx = db[db[0]=='NO'].index.values.astype(int)
    idx_name = db[db[0]=='!'].index
    # NOTE: 140903 has inverted blocks
    if grb == '140903A':
        db01gamma=db[0:idx[0]]
    else:
        db01=db[0:idx[0]]
    if len(idx) > 1:
        if grb == '140903A':
            db01=db[idx[0]+2:idx[1]]
        else:
            db01gamma=db[idx[0]+2:idx[1]]
    else:
        if grb == '140903A':
            db01=db[idx[0]+2:]
        else:
            db01gamma=db[idx[0]+2:]
    db0=db01
    db0gamma=db01gamma
    if len(idx_name) > 1:
        if 'wt' in db[1].iloc[idx_name[1]]:
            # if second block is wt then first block is slew
            db02=db[idx[1]+2:idx[2]]
            db02gamma=db[idx[2]+2:idx[3]]
            db03=db[idx[3]+2:idx[4]]
            db03gamma=db[idx[4]+2:]
            db0=pd.concat([db01,db02,db03])
            db0gamma=pd.concat([db01gamma,db02gamma,db03gamma])
        elif 'xrtpcflux' in db[1].iloc[idx_name[1]]:
            # if second block is pc then first block is wt
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
    return tx,dtxp,dtxm,fx,dfxp,dfxm,gammaxrt,dgammaxrtp,dgammaxrtm
