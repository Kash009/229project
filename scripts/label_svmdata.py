from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call, Popen, PIPE
import argparse
import numpy

#definitions
datapath = '../data/PixelXYZ'
deltapath = '../data/DeltaEToMax'
labelpath = '../data/Label'
rawdatapath = '../data/AllRawData'
svmdatapath = '../model/SVMData'

datacount = 0
threshold = 0
maxArray = []
labelArray = []

def getSampleCount():
  global datacount
  with open('../data/datacount') as file:  
    datacount = int(file.readline())

def generateMaxArray(layers):
  global maxArray 
  maxArray = [0]*(3*datacount)
  with open(datapath + str(layers)) as f:
    for i in range(0,datacount):
      s,x,y,z = [float(x) for x in f.readline().split()]
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
  with open(deltapath + str(layer),'w') as result:
    with open(datapath + str(layer)) as layer_file:
      for i in range(0,datacount):
        s,x1,y1,z1 = [float(x) for x in layer_file.readline().split()]
        s = int(s)
        delta = deltaE(maxArray[3*s],maxArray[3*s+1],maxArray[3*s+2],x1,y1,z1)
        line = str(s) + ' ' + str(delta) + '\n'
        result.write(line)

def generateLabelsAndData(layer):
  global labelArray
  labelArray = [0]*datacount
  with open(labelpath + str(layer),'w') as result:
    with open(deltapath + str(layer)) as layer_file:
      for i in range(0,datacount):
        s,delta = [float(x) for x in layer_file.readline().split()]
        s = int(s)
        difference = 0
        if delta > threshold:
          difference = 1
          labelArray[s] = difference
        result.write(str(s) + ' ' + str(difference) + '\n') 
  with open(svmdatapath+'-' + str(threshold) +'-' + str(layer),'w') as result:
    with open(rawdatapath + str(layer)) as layer_file:
      for i in range(0,datacount):
        #hacky way for features
        s,a1, a2, a3, a4, a5, a6 = [float(x) for x in layer_file.readline().split()]
        s = int(s)
        difference = labelArray[s]
        result.write(str(difference) + ' 1:' + str(a1)+ ' 2:' + str(a2)+ ' 3:' + str(a3)+ ' 4:' + str(a4)+ ' 5:' + str(a5)+ ' 6:' + str(a6) + '\n') 

## function runfile
def runfile(s_threshold,s_layers):
  getSampleCount()
  layers = int(s_layers)
  global threshold
  threshold = float(s_threshold)
  if threshold > 0:
    for i in range(0,layers-1):
      generateLabelsAndData(i+1);
  else:
    generateMaxArray(layers);
    for i in range(0,layers-1):
      generateDistanceFromMax(i+1);


#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'threshold')
parser.add_argument( 'layers')
args= parser.parse_args()
runfile(args.threshold,args.layers)
