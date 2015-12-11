#include <pngwriter.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <ctime>
#include <cstdlib>
#include <vector>
#include <algorithm>
string files[5]; 
string imagepath = "../images/";

void writeFeature(string feature,int layer, int res){
    string datafile = feature + to_string(layer);
    string filename = imagepath + datafile +".png";
    pngwriter im(res,res,0.0,filename.c_str()); 
    std::ifstream file("../data/Z-bunny/" + datafile);
      int id = 0;
      double value = 0;
      char formatting_tab;
      char formatting_line;
      int x,y;
      double a,c;
      if (file.is_open())
        {
          while((file >> id) && (file >> value)){
            x = id /res;
            y = id - res*x;
            if(feature == files[3])value /=5;//for yvariance
            im.plot_blend((x+1), res-y, 1, value, value, value); 
          }
        }else std::cout <<"cannot open " << feature;
      file.close();
   std::cout << " done. Writing to disk...";
   im.close();
   std::cout << " done." << std::endl;
}


void writeCount(){
  int res = 201;
    string filename = imagepath + "result.png";
    pngwriter im(res,res,0.0,filename.c_str()); 
    float value;
    for(int x = 0; x < res; x++){
      for(int y = 0; y < res; y++){
        std::ifstream file("../data/dataVS/bunny/p" + std::to_string(x*res+y));
        if (file.is_open()){
          file >> value;
          value /= 512;
          value = powf(value, 0.2);
            im.plot_blend((x+1), res-y, 1, value, value, value); 
        }else std::cout << "meep!";
      file.close();
      }
    }  
   std::cout << " done. Writing to disk...";
   im.close();
   std::cout << " done." << std::endl;
}

void writeRGB(int layer, int res){
    std::ifstream file("../data/PixelXYZ"+std::to_string(layer));
    string s = imagepath + "rgb"+std::to_string(layer)+".png";
    pngwriter im(res,res,0.0,s.c_str()); 
      int id,x,y;
      double r,g,b;
      if (file.is_open())
        {
          while((file >> id) && (file >> r)&& (file >> g)&& (file >> b)){
            x = id /res;
            y = id - res*x;
            im.plot_blend((x+1), res-y, 1, 100*r,100*g,100*b); 
          }
        }else std::cout <<"cannot open file...";
      file.close();
      im.close();
}


int main(){  
  writeCount();
  /*int maxlayer = 1;
  files[0] = "SingleShape";
  files[1] = "SinglePrimitive";
  files[2] = "DeltaToMax";
  files[3] = "VarY";
  files[4] = "DeltaY";
  int res=1;

  std::cout << "layer: ";
  std::cin >> maxlayer;
  std::cout << "res: ";
  std::cin >> res;
  res++;

  for(int layer = 1; layer <= maxlayer; layer++){
       writeRGB(layer,res);
      for(int i = 0; i < 5; i++){
         writeFeature(files[i],layer,res);
      }
    }*/

   std::cout << "\n\npngtest has finished. Take a look at the PNG images that have been created!\n";
   return 0;
}