from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call, Popen, PIPE
import argparse

##function modifysample
#edits pbrt file corresponding to values
def modifysample(path, value1, value2):
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
      	  line = 'Sampler "adaptive" "string datatype" ["train"]\n'
      	  second = False
        nline = nline + 1
	new_file.write(line)
  close(fh)
  remove(path)
  move(result, path)

## function runfile
#runs through all sample combinations
def runfile(path,output,resolution):
  eye = 1
  f = open(output,'w')
  modifysample(path,output,resolution)
  call(['../pbrt-v2/src/bin/pbrt',path])
  #p = Popen(['exrdiff', stdout=PIPE, stderr=PIPE)
  #stdout, stderr = p.communicate()
  #print >>f, stdout
  f.close()


#main script
parser = argparse.ArgumentParser()
parser.add_argument( 'path')
parser.add_argument( 'output')
parser.add_argument( 'resolution')
args= parser.parse_args()
runfile(args.path, args.output, args.resolution)
