from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call, Popen, PIPE
import argparse

##function modifysample
#edits pbrt file corresponding to values
def modifysample(path, sampler, value1, value2):
  #temp file
  fh, result = mkstemp()
  nline = 1;
  with open(result,'w') as new_file:
    with open(path) as old_file:
      for line in old_file:
        if nline == 1:
          line = 'Film "image" "string filename" ["../images/' +str(value1) +'-'+str(value2)+ '.exr"] "integer xresolution" ['+str(value2)+'] "integer yresolution" ['+str(value2)+']\n'
          first = False
          second = True
      	if nline == 2:
      	  line = str(sampler);
      	  second = False
        nline = nline + 1
	new_file.write(line)
  close(fh)
  remove(path)
  move(result, path)

## function runfile
#runs through all sample combinations
def runfile(path,output,resolution, samples, modelpath):
  eye = 1
  svm = "SVM"
  standard = "Standard"
  ld = "LD"
  svmsampler = 'Sampler "adaptive" "string datatype" ["test"] "string modelpath" ["' + modelpath + '"]\n'
  standardsampler = 'Sampler "random" "integer pixelsamples" ['+str(samples)+']\n'
  ldsampler = 'Sampler "lowdiscrepancy" "integer pixelsamples" ['+str(samples)+']\n'
  modifysample(path,svmsampler,output+svm,resolution)
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  modifysample(path,standardsampler,output+standard,resolution)
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  modifysample(path,ldsampler,output+ld,resolution)
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  call(['../pbrt-v2/src/bin/pbrt',path])
  #p = Popen(['exrdiff', stdout=PIPE, stderr=PIPE)
  #stdout, stderr = p.communicate()
  #print >>f, stdout


#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'path')
parser.add_argument( 'output')
parser.add_argument( 'resolution')
parser.add_argument( 'samples')
parser.add_argument( 'modelpath')
args= parser.parse_args()
runfile(args.path, args.output, args.resolution, args.samples, args.modelpath)
