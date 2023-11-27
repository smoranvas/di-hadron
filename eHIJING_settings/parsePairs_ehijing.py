import numpy as np, pandas as pd


def get_df(filename,nevents=1e17):
    with open(filename, "r") as f:
        d={"Q2":[], "x":[], "pT1":[], "pT2":[], "y1":[], "y2":[], "z1":[], "z2":[], "phi1":[], "phi2":[],
          "p1":[], "p2":[], "thetalab1":[], "thetalab2":[],"philab1":[], "philab2":[],}
        d_lead={"Q2":[], "x":[], "pT1":[], "y1":[], "z1":[], "phi1":[],
          "p1":[],  "thetalab1":[],"philab1":[],}
        pT2s=[]
        phi2s=[]
        y2s=[]
        z2s=[]
        p2s=[]
        thetalab2s=[]
        philab22s=[]
        foundLeading=False
        events=0
        while True:
            line=f.readline()
            if line is not None:
                s=line.split()
            else : break
            if events>nevents:
                break
            if len(s)==0:
                break
            if s[0]=='#':
                nh2=len(pT2s)
                if foundLeading==True:
                    
                    d_lead['Q2'].append(Q2)
                    d_lead['x'].append(x)
                    d_lead['pT1'].append(pT1)
                    d_lead['y1'].append(y1)
                    d_lead['z1'].append(z1)
                    d_lead['p1'].append(p1)
                    d_lead['thetalab1'].append(thetalab1)
                    d_lead['philab1'].append(philab1)
                    d_lead['phi1'].append(phi1)
                    if nh2 >0:
                        d['Q2']+=[Q2]*nh2
                        d['x']+=[x]*nh2
                        d['pT1']+=[pT1]*nh2
                        d['y1']+=[y1]*nh2
                        d['z1']+=[z1]*nh2
                        d['p1']+=[p1]*nh2
                        d['thetalab1']+=[thetalab1]*nh2
                        d['philab1']+=[philab1]*nh2
                        d['phi1']+=[phi1]*nh2
                        d['pT2']+=pT2s
                        d['y2']+=y2s
                        d['z2']+=z2s
                        d['phi2']+=phi2s
                        d['p2']+=p2s
                        d['philab2']+=philab2s
                        d['thetalab2']+=thetalab2s
                    foundLeading=False

                Q2=float(s[1])
                x=float(s[2])
                
                #reset
                pT2s=[]
                phi2s=[]
                y2s=[]
                z2s=[]
                p2s=[]
                thetalab2s=[]
                philab2s=[]
                foundLeading=False
                events+=1
                
                    
            elif s[0]=="211" and float(s[1])>0.5: #leading pi+
                z1=float(s[1])
                pT1=float(s[2])
                phi1=float(s[3])
                y1=float(s[4])
                p1=float(s[5])
                thetalab1=np.pi-float(s[6])
                philab1=float(s[7])
                foundLeading=True
            elif s[0]=="-211" and float(s[1])<0.5: #subleading pi-
                z2s.append(float(s[1]))
                pT2s.append(float(s[2]))
                phi2s.append(float(s[3]))
                y2s.append(float(s[4]))
                p2s.append(float(s[5]))
                thetalab2s.append(np.pi-float(s[6]))
                philab2s.append(float(s[7]))
    return pd.DataFrame(d), pd.DataFrame(d_lead), events


if __name__=="__main__":
    import glob, sys
    infile=sys.argv[1]
    outfile=sys.argv[2]
    n=1e9
    if len(sys.argv)>3:
        n=int(float(sys.argv[3]))
    pairs=None
    leading=None
    for f in glob.glob(infile):
        print("reading file", f)
        pairs1, leading1,n1=get_df(f,n)


        if pairs is None:
            pairs=pairs1
            leading=leading1
        else:
            pairs=pd.concat([pairs, pairs1])
            leading=pd.concat([leading, leading1])
        n-=n1
        if n<0:
            break
    import pickle
    with open(outfile, "wb") as f:
        pickle.dump([pairs, leading], f)
