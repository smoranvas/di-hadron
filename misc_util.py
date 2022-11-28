import pandas as pd, numpy as np


#mode options are
# df:  return the slice dataframes (default)
# len: return the numbers of entries in each slice
class BinIterator:
    def __init__(self,df,xvar,min,max,bins,mode='df'):
        self._df = df
        self._xvar = xvar
        self._min = min
        self._max = max
        self._bins = bins
        self._i = 0
        self._mode = mode
    def __iter__(self):
        self._i = 0
        return self
    def __next__(self):
        if(self._i >= self._bins):
            raise StopIteration
        dx = (self._max-self._min)/self._bins
        #case:  single dataframe
        if type(self._df)==pd.DataFrame:
            ret = self._df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx))
            if self._mode == 'len':
                ret = len(ret)
        #case:  list of dataframes
        else: 
            ret = [df.query("%s >= %s and %s < %s"\
                       %(self._xvar,
                         self._min + self._i*dx,
                         self._xvar,
                         self._min +(self._i+1)*dx)) for df in self._df]
            if self._mode == 'len':
                ret = [len(r) for r in ret]
        self._i+=1
        return self._min+(self._i-1/2)*dx, ret


def query_or_all(df,q):
    if q != "":
        return df.query(q)
    else:
        return df

class BinnedVariable:
    def __init__(self,var,tex,bins,unit):
        self._var=var
        self._tex=tex
        self._bins=bins
        self._unit=unit
        
    #formats the variable to display a label for an axis
    def axis_label(self,):
        if self._unit is not None:
            return f"{self._tex} [{self._unit}]"
        else :
            return self._tex
    #formats the label for a slice
    def slice_label(self, i,f="%.2f"):
        ret = f"{f}<{self._tex}<{f}" %(self._bins[i],self._bins[i+1])
        if self._unit is not None:
            ret += " " + self._unit
        return ret
    #formats the query string for slicing
    def slice_query(self,i):
        return f"{self._bins[i]}<{self._var} and {self._var}<{self._bins[i+1]}"
    #returns a queried dataframe
    def query(self,df,i):
        return df.query(self.slice_query(i))
    def count_query(self,df,i):
        return len(self.query(self,df,i))
    def bin_center(self,i,shift_frac=0):
        return (self._bins[i]+self._bins[i+1])/2 + \
            shift_frac*(self._bins[i+1]-self._bins[i])
    def auto_bins(self,df, nbins,style="linear"):
        a = df.eval(self._var)
        if style=="linear":
            self._bins=np.linspace(np.min(a),np.max(a),nbins+1)
        if style=="quantile":
            self._bins=[a.quantile(q) for q in np.linspace(0,1,nbins+1)]
        return self._bins
    def __len__(self):
        return len(self._bins)-1

class RatioHistogram :
    def __init__(self, binvar,df_num, df_denom, df_denom2=None,df_num2=None,scale=1, ratioName=""):
        self._binvar = binvar
        y1,self._xedges = np.histogram(df_num.eval(binvar._var),bins=binvar._bins)
        y2,_ = np.histogram(df_denom.eval(binvar._var),bins=binvar._bins)
        y3,y4 = None,None
        if df_denom2 is not None:
            y3,_ = np.histogram(df_denom2.eval(binvar._var),bins=binvar._bins)
            y4,_ = np.histogram(df_num2.eval(binvar._var),bins=binvar._bins)

        self._xcenter = np.add(self._xedges[:-1],self._xedges[1:])/2.
        self._y = y1/y2*scale
        if df_denom2 is not None:
            self._y /= y3/y4
        self._dy = self._y*np.sqrt(1/y1+1/y2+(0 if df_denom2 is None else (1/y3+1/y4)))
        self._ratioName = ratioName
    #def to_tex_table(self, filename=None, syst=None):
    #    ret = ""
    #    ret += f"{self._binvar._tex} range & {ratioName}" 
                                   
    
    
    
    
    