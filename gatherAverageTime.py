import pandas as pa
import numpy as np

data = pa.read_csv('processTimes_Np-40.txt',  delim_whitespace=True, header=None)
data = data.values

for i in range(0,40):
  sumTot = 0.
  for j in range(0,len(data)):
    if (data[j,0] == i):
      sumTot = sumTot+data[j,1]
  avTot = sumTot/150.
  print(avTot)
