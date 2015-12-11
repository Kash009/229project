from os import remove, close
import argparse
import numpy
from svmutil import *
#definitions
svmmodelpath = '../model/Full/SVMModel'

maxArray = []
labelArray = []


## function runfile
def runfile(modelpath, datapath, layers):
  layers = int(layers)
  for i in range(1, layers):
    [y,x] = svm_read_problem(str(datapath) +'-'+str(i))
    model = svm_load_model(str(modelpath) +'-'+str(i))#str(str(model) + '-'+str(i)))
    svm_predict(y,x,model)

#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'modelpath')
parser.add_argument( 'datapath')
parser.add_argument( 'layers')
args= parser.parse_args()
runfile(args.modelpath,args.datapath, args.layers)
