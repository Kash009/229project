from os import remove, close
import argparse
import numpy
from svmutil import *
#definitions
svmdatapath = '../model/SVMData'
svmmodelpath = '../model/SVMModel'

maxArray = []
labelArray = []


## function runfile
def runfile(s_threshold,s_layers):
  layers = int(s_layers)
  for i in range(1, layers+1):
    [y,x] = svm_read_problem(svmdatapath + '-'+s_threshold + '-'+str(i))
    model = svm_train(y,x)
    svm_predict(y,x,model)
    svm_save_model(svmmodelpath + '-'+s_threshold + '-'+str(i),model)


#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'threshold')
parser.add_argument( 'layers')
args= parser.parse_args()
runfile(args.threshold,args.layers)
