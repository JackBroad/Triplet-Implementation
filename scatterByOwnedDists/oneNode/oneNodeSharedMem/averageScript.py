import pandas as pa
import numpy as np

data = pa.read_csv('30-data.txt',  delim_whitespace=True, header=None)
data = data.values
data = data[1:,:]

for i in range(0,8):
  dataVec = data[:,i]
  sumData = np.sum(dataVec)
  avData = sumData/len(data)
  print(avData, '\n')
