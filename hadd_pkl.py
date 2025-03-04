import pandas as pd, pickle as pkl, sys

a=sys.argv[1]
df_mc_pairs=[]
df_mc_trigger=[]
for i in range(1,9+1):
    
    fname=f"MC_pion_proton_pairs_{a}_{i}.pkl"
    import glob
    if len(glob.glob(fname))==0:
        continue
    with open(fname, "rb") as f:
        print(fname)
        d=pkl.load(f)
        df_mc_pairs.append(d[a])
        df_mc_trigger.append(d[a+"_trigger"])
        print(len(df_mc_pairs[-1]), len(df_mc_trigger[-1]))
df_mc={}
df_mc[a]=pd.concat(df_mc_pairs)
df_mc[a+"_trigger"]=pd.concat(df_mc_trigger)

with open(f"MC_pion_proton_pairs_{a}.pkl", "wb") as f:
    print(df_mc)
    pkl.dump(df_mc, f, pkl.HIGHEST_PROTOCOL)
