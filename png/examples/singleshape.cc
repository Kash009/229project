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
string files[4]; 

int main()
{  
  int ver = 1;
  files[0] = "SingleShape";
  files[1] = "SinglePrimitive";
  files[2] = "DeltaY";
  files[3] = "VarY";
  files[4] = "Label";
  int res=1;
  ////////////////////////////////////////////////////////////////////////////////////
   /*  two.png
    *  This will be a 300x300 image with a black background, it will be called two.png.
    * Note that we are using 0.0 as the background colour and that we are using a string
    * type object as the filename, which can convert itself into a const char * with 
    * filename.c_str().
    * */
     std::cout << "layer: ";
     std::cin >> ver;
     std::cout << "res: ";
     std::cin >> res;
     res++;

  std::ifstream file("../../data/RGB");
  pngwriter im(res,res,0.0,"rgb.png"); 
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

  /*for(int i = 0; i < 5; i++){
     string datafile = files[i] + to_string(ver);
     string filename = "../../images/" + datafile+".png";
     pngwriter im(res,res,0.0,filename.c_str()); 
    
     std::ifstream file("../../data/" + datafile);
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
            if(i==3)value /=5;//for yvariance
            im.plot_blend((x+1), res-y, 1, value, value, value); 
          }
        }else std::cout <<"cannot open file...";
      file.close();
    
   std::cout << " done. Writing to disk...";
   im.close();
   std::cout << " done." << std::endl;
    }*/

   std::cout << "\n\npngtest has finished. Take a look at the PNG images that have been created!\n";
   return 0;
}
