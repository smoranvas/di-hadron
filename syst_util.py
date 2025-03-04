import pandas as pd, numpy as np
from scipy.stats import chi2 as chi2lib
from scipy.stats import norm as normlib
from scipy.stats import binom as binomlib

#tightVsLoose=false (default):  delta(y1-y2)^2=(delta(y1)^2+delta(y2)^2
#  tightLoose=true  (when you're comparing results with a tight cut and a loose cut, so the former is a subset of the latter, and thus the results are strongly correlated): delta(y1-y2)^2=(delta(y1)^2-delta(y2)^2
def check_for_equality(y1, dy1,y2,dy2=0, ndof=None, printInfo=True, tightVsLoose=False, assignSystToBins=None,threshold=0.05):
    if ndof is None:
        ndof = len(y1)
    if not tightVsLoose:
        normedResiduals = (y1-y2)/np.hypot(dy1,dy2)
    else:
        normedResiduals = (y1-y2)/np.sqrt(abs(dy1**2-dy2**2))
    #print(normedResiduals)
    chi2 = np.sum(normedResiduals**2)
    pval_chi2 = chi2lib.sf(chi2,ndof)
    if printInfo == 1:
        print(f"chi2 test\n: chi2={chi2:2f}\nndof={ndof}\npval={pval_chi2:3f}")
    
    if pval_chi2<threshold:
        # estimate the systematic error by searching for a chi2
        # value such that adding it to both halves would bring the pval to 0.5
        dy_syst = np.min(dy1)/10.
        change = 2
        if assignSystToBins == None:
            assignSystToBins=[1]*len(dy1)
        assignSystToBins=pd.Series(assignSystToBins)
        for i in range(100):
            if not tightVsLoose:
                normedResiduals_mod = (y1-y2)/np.sqrt(dy1**2+dy2**2+2*((y1+y2)/2*dy_syst*assignSystToBins)**2)
            else:
                normedResiduals_mod = (y1-y2)/np.sqrt(dy1**2-dy2**2+2*((y1+y2)/2*dy_syst*assignSystToBins)**2)
            pval_mod = chi2lib.sf(sum(normedResiduals_mod**2),ndof)
            #print(pval_mod)
            if pval_mod >=threshold*.8 and pval_mod<=threshold*1.2:
                break
            if pval_mod >threshold*1.2:
                change = 1.05
                dy_syst/=change
            else :
                dy_syst*=change
        print(f"Failed chi2 test (p<{threshold}).  estimated systematic error:", dy_syst)
    #now do look-elsewhere-effect test on the largest residual
    largestNormedResidual = 0
    i = 0
    for ii,zi in enumerate(normedResiduals):
        if abs(zi) > abs(largestNormedResidual):
            largestNormedResidual=abs(zi)
            i = ii
    if printInfo == 1:
        print(f"\nLEE test:\n largest normed residual ={largestNormedResidual:2f} at bin {i}")
    pval_local = 2*normlib.sf(abs(largestNormedResidual))
    if printInfo == 1:
        print(f"local pval: {pval_local:3f}")
    pval_global = 1-(1-pval_local)**ndof
    if printInfo == 1:
        print(f"global pval: {pval_global:3f}")
    
    #now check if there is a systematic trend in the signs of the residuals
    n_positives = sum(normedResiduals>0)
    
    pval_signs = binomlib.sf(n_positives,len(y1), 0.5)
    if pval_signs < 0.5:
        pval_signs = 2*pval_signs
    if pval_signs > 0.5:
        pval_signs = 2*(1-pval_signs)
    if printInfo == 1:
        print(f"\nsigns test:\npval:  {pval_signs}")
    #more compact level of information printing
    if printInfo == 2:
        print(f"p_chi2 ={pval_chi2:3f}; p_lee ={pval_global:3f}; p_signs ={pval_signs:3f}")
    return pval_chi2, pval_global, pval_signs

def check_for_equality_of_dfs(df1,df2,**arg):
    y1 = df1.y
    dy1 = df1.dy
    y2 = df2.y
    dy2 = df2.dy
    
    return check_for_equality(y1, dy1,y2,dy2, **arg)

def check_for_equality_of_dfs_multi(dfs,**arg):
    y1 = pd.concat([df.y for df in dfs])
    dy1 = pd.concat([df.dy for df in dfs])
    
    y2 = pd.Series([0]*len(dfs[0]))
    dy2 = pd.Series([0]*len(dfs[0]))
    wsum = 0
    for df in dfs:
        w = 1/df.dy**2
        y2+= df.y*w
        dy2+= w
        wsum +=w
    y2/=wsum
    dy2 =np.sqrt(1/dy2)
    y2 = pd.concat([y2 for df in dfs])
    dy2 = pd.concat([dy2 for df in dfs])
    ndof = len(dfs[0])*(len(dfs)-1)
    #print(y2,dy2)
    return check_for_equality(y1, dy1,y2,dy2,ndof=ndof, **arg)
