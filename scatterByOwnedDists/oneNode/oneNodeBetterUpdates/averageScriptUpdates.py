import pandas as pa
import numpy as np

data = pa.read_csv('40-data.txt',  delim_whitespace=True, header=None)
data = data.values
data = data[1:,:]

dataVec = np.zeros(len(data))
updateVec = np.zeros(len(data))
noUpdateVec = np.zeros(len(data))
for i in range (0,len(data)):
  dataVec[i] = data[i,1]
  updateVec[i] = data[i,2]
  noUpdateVec[i] = data[i,3]

sumUpdateData = np.sum(updateVec)
nUpdates = len(data)/2
avUpdateData = sumUpdateData/nUpdates
print('average accept move update time', avUpdateData)

sumNoUpData = np.sum(noUpdateVec)
avNoUpData = sumNoUpData/nUpdates
print('average reject move update time', avNoUpData)

sumData = np.sum(dataVec)
avData = sumData/len(data)
print('average move time', avData)
