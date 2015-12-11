from os import remove, close
import argparse
import numpy
import math

#definitions
datapath = '../data/PixelXYZ'
deltapath = '../data/DeltaToMax'
meandeltapath = '../data/MeanDelta'
labelpath = '../data/Killeroo/Label'
rawdatapath = '../data/AllRawData'
svmdatapath = '../data/Killeroo/SVMData'

datacount = 0
label = 0
maxArray = []
labelArray = []
indexArray = []
total_count = 0
meanDelta = []

def cleanData(layers):
  global total_count
  total_count = 0
  global indexArray
  indexArray = [0]*datacount
  firstArray = 1
  for i in range (1,layers+1):
    print(str(i))
    with open(datapath + str(i)) as f: 
      if firstArray == 1:
        for line in f:
          s,a1, a2, a3 = [float(x) for x in line.split()]
          s = int(s)
          total_count = total_count +1
          indexArray[s]=1
        firstArray = 0
      else:
        buffArray = [0]*datacount
        for line in f:
          s, a1, a2, a3 = [float(x) for x in line.split()]
          s = int(s)
          buffArray[s]=1
        for j in range(0,datacount):
          if indexArray[j] == 1 and buffArray[j] == 0:
            indexArray[j] = 0
            total_count = total_count -1
  print(str(total_count))

def delta(x,y,z,x1,y1,z1):
  return ((x1-x)**2 + (y1-y)**2 + (z1-z)**2)**0.5

def getSampleCount():
  global datacount
  with open('../data/datacount') as f:  
    datacount = int(f.readline())

def generateMaxArray(layers):
  global maxArray 
  maxArray = [0]*(3*datacount)
  with open(datapath + str(layers)) as f:
    for line in f:
      s,x,y,z = [float(x) for x in line.split()]
      s = int(s)
      maxArray[3*s]=x
      maxArray[3*s+1]=y
      maxArray[3*s+2]=z

def func(x):
  t = pow(x,1.0/3)
  if t <= (6.0/29):
    t = 1.0/3 * (841.0/36) * x + 4.0/29
  return t

def Lab(x,y,z):
  nx = x/95.047
  ny = y/100
  nz = z/108.883
  L = 116 * func(ny) - 16
  a = 500*(func(nx)-func(ny))
  b = 200*(func(ny) - func(nz))
  return (L,a,b)

def deltaE(x,y,z,x1,y1,z1):
  [L,a,b] = Lab(x,y,z)
  [L1,a1,b1] = Lab(x1,y1,z1)
  return 0.5*(pow(L1-L,2) + pow(a1-a,2) + pow(b1-b,2))

def generateDistanceFromMax(layer):
  mdelta = 0
  with open(deltapath + str(layer),'w') as result:
    with open(datapath + str(layer)) as layer_file:
      for line in layer_file:
        s,x1,y1,z1 = [float(x) for x in line.split()]
        s = int(s)
        if indexArray[s]==1:
          d = delta(maxArray[3*s],maxArray[3*s+1],maxArray[3*s+2],x1,y1,z1)
          mdelta = mdelta + d
          line = str(s) + ' ' + str(d) + '\n'
          result.write(line)

  mdelta = mdelta/datacount
  meanDelta[layer-1] = mdelta

def generateLabelsAndData(layer):
  global labelArray
  labelArray = [0]*datacount
  with open(labelpath + str(layer),'w') as result:
    with open(deltapath + str(layer)) as layer_file:
      for line in layer_file:
        s,d = [float(x) for x in line.split()]  
        s = int(s)
        if indexArray[s] == 1:
          difference = -1
          if d > meanDelta[layer-1]:
            difference = 1
            labelArray[s] = difference
          result.write(str(s) + ' ' + str(difference) + '\n')

  with open(svmdatapath+'-' + str(label) +'-' + str(layer),'w') as result:
    with open(rawdatapath + str(layer)) as layer_file:
      for line in layer_file:
        #hacky way for features
        s,a1, a2, a3, a4, a5, a6, a7, a8 = [float(x) for x in line.split()]
        s = int(s)
        if indexArray[s]== 1: # and s % 40 == 0:
          difference = labelArray[s]
          result.write(str(difference) 
            + ' 1:' + str(a1)+ ' 2:' + str(a2)+ ' 3:' + str(a3)+ ' 4:' + str(a4)
            + ' 5:' + str(a5) + ' 6:' + str(a6) + ' 7:' + str(a7) + ' 8:' + str(a8) + '\n') 

def ReadMean():
  global meanDelta;
  with open(meandeltapath) as f:
    for line in f:
      layer, layerthreshold = line.split()
      meanDelta[layer-1] = layerthreshold



## function runfile
def runfile(s_label,s_layers):
  getSampleCount()
  layers = int(s_layers)
  global indexArray
  indexArray = [1]*datacount
  global meanDelta
  meanDelta = [0]*layers
  global label
  label = int(s_label)
  if label > 0:
    for i in range(0,layers-1):
      generateLabelsAndData(i+1);
  else:
    generateMaxArray(layers);
    for i in range(0,layers-1):
      generateDistanceFromMax(i+1);
    with open(meandeltapath,'w') as result:
      for i in range(0,layers-1):
        result.write(str(i+1)+' '+str(meanDelta[i]) + '\n')



#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'filelabel')
parser.add_argument( 'layers')
args= parser.parse_args()
runfile(args.filelabel,args.layers)
