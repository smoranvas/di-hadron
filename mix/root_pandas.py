import uproot3,ROOT,pandas as pd

def to_root(df,filename, treename):
    with uproot3.recreate(filename) as f:
        f[treename] = uproot3.newtree({col:df[col].dtype for col in df.columns})
        f[treename].extend(dict(df))

pd.DataFrame.to_root = to_root

def read_root(filename,treename=None,N=None):
    with uproot3.open(filename) as f:
        if treename is None:
            if len(f.keys()) == 1:
                treename = f.keys()[0]
            elif len(f.keys()) == 0:
                raise Exception("error: no tree names in file " + filename)
            else:
                raise Exception("error: treename must be specified; more than one tree in " + filename)
            
        df = f[treename].pandas.df(entrystop=N)
    return df
if __name__=='__main__':
    df = pd.DataFrame({'a':[1.0,2.2,3.3],'b':[0,1,2]})
    #to_root(df,'out.root', 'tree')
    df.to_root('out.root', 'tree')  

    df = read_root('out.root')
    print("reading a one tree file without specifying tree name\n",df)
    df = read_root('out.root',"tree")
    print("reading a one tree file with specifying tree name\n",df)
