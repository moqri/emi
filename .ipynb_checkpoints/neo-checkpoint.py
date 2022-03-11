import pandas as pd
import numpy as np
import time
source='data/source/'
neo_=source+'pbmc_neo/GSE140730-GPL20795_series_matrix.txt.gz'
n=8
neos=[]
t=time.time()
neo=pd.read_csv(neo_,sep='\t',skiprows=31).T
neo=neo[neo[9]=='diagnosis: TD']
gsms=neo[43].index
s3s=[]
s4s=[]
for i in range(len(gsms)):
    gsm=gsms[i]           
    print(gsm,end=',')
    neo=pd.read_table('data/source/pbmc_neo/'+gsm+'.txt.gz',header=None,usecols=[0,1,2,3,4],nrows=10**n)
    neo=neo[~neo[0].isin(['chrX','chrY','chrM'])]
    neo=neo[~neo[0].str.startswith('chrU')]
    neo=neo[~neo[0].str.endswith('m')]
    neo=neo[~neo[0].str.endswith('t')]
    neo.index=neo[0].str[3:]+'_'+(neo[1]-1).astype(str)    
    neo=neo.drop([0,1],axis=1)
    if i==0:
        ind=neo.index
        s3=np.array(neo.loc[:,3])
        s4=np.array(neo.loc[:,4])
    else:
        s3=np.array(neo.loc[ind,3])
        s4=np.array(neo.loc[ind,4])
    s3s.append(s3)
    s4s.append(s4)
    print(round(time.time()-t))
s3ss=np.array(s3s).sum(0)
s4ss=np.array(s4s).sum(0)
df=pd.concat([pd.DataFrame(s3ss),pd.DataFrame(s4ss)],1)
df.index=ind
df.to_csv('data/pmbc_neo_signals.csv')  
print(round(time.time()-t))
