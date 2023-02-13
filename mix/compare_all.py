import numpy as np, pandas as pd, root_pandas,sys,matplotlib , matplotlib.pyplot as plt

var = "diff_phi_cm"
#"h1_ph-h2_ph+2*3.14159*(h1_ph-h2_ph<-3.14159)-2*3.14159*(h1_ph-h2_ph>3.14159)"

nbins = 50
dys = []
ys = []
rang = None
for fname in sys.argv[1:]:
    df = root_pandas.read_root(fname)
    a = df.eval(var)
    if rang == None:
        rang=(min(a),max(a))
    y,x = np.histogram(a,range=rang,bins=nbins)
    y = pd.Series(y)
    denom = np.sum(y)*(x[1]-x[0])/(rang[1]-rang[0])
    dy = np.sqrt(y)/denom
    y = y/denom
    
    dys.append(dy)
    ys.append(y)
    x = np.add(x[1:],x[:-1])/2
    plt.errorbar(x,y,dy,linestyle='',marker='o',label=fname)

from scipy.stats import chi2 as chi2lib

for i in range(len(sys.argv)-1):
    for j in range(i+1,len(sys.argv)-1):
        chi2 = sum((ys[i]-ys[j])**2/(dys[i]**2+dys[j]**2))
        print(sys.argv[1+i]," vs ",sys.argv[1+j])
        print("chi2=",chi2, ",  ndof = ", nbins-1, ",  p=",chi2lib.sf(chi2,nbins-1))

plt.legend()
plt.gca().set_ylim(0,2)
plt.gca().set_ylabel("$dN/d\\Delta\\phi_{\\mathrm{cm}}\\times 2\\pi/N$")
plt.gca().set_xlabel("$\\Delta\\phi_{\\mathrm{cm}}$")
plt.savefig("mixed_dphi.png")
plt.show()
