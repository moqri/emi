import pandas as pd
import numpy as np
import pyBigWig
from bisect import bisect_left
import seaborn as sns

def take_closest(myList, myNumber):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before
    
def bw2df(bw_url,tss):
    bw=pyBigWig.open(bw_url)
    bws=[]
    count=0
    for i,r in tss.iterrows():
        count+=1
        if (count % 1000) == 0:        
            print(r['chrom'],end='.')
        bwi=bw.intervals("chr"+str(r['chrom']),r['chromStart'],r['chromEnd'])
        if bwi is not None:
            for b in bwi:
                bws.append([r['chrom'],b[0],b[2]])
    df=pd.DataFrame(bws)
    df.drop_duplicates(inplace=True)
    df.index=df[0].astype(str)+'_'+df[1].astype(str)
    df.drop([0,1],axis=1,inplace=True)
    df.columns=['beta']
    print(df.shape)
    return df
def palette(values, palette_name):
    # normalize the values to range [0, 1]
    normalized = (values - min(values)) / (max(values) - min(values))
    # convert to indices
    indices = np.round(normalized * (len(values) - 1)).astype(np.int32)
    # use the indices to get the colors
    palette = sns.color_palette(palette_name, len(values))
    return np.array(palette).take(indices, axis=0)    

    
