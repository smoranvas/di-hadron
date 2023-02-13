import numpy as np, pandas as pd, root_pandas,sys


nbins = 100
df1 = root_pandas.read_root(sys.argv[1])
df2 = root_pandas.read_root(sys.argv[2])

y1,x = np.histogram(df1.diff_phi_cm,range=(-np.pi,np.pi),bins=nbins)
y2,x = np.histogram(df2.diff_phi_cm,range=(-np.pi,np.pi),bins=nbins)

y1 = pd.Series(y1)
y2 = pd.Series(y2)

sumy1 = np.sum(y1)
sumy2 = np.sum(y2)

dy1 = np.sqrt(y1)/sumy1
dy2 = np.sqrt(y2)/sumy2

y1 = y1/sumy1
y2 = y2/sumy2

from scipy.stats import chi2 as chi2lib
chi2 = sum((y1-y2)**2/(dy1**2+dy2**2))
print("chi2=",chi2, ",  ndof = ", nbins-1, ",  p=",chi2lib.sf(chi2,nbins-1))


