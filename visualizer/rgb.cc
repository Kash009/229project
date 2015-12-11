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
  int res=1;
  float scale = 1;
     std::cout << "res: ";
     std::cin >> res;
    std::cout << "scale: ";
     std::cin >> scale;
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
            im.plot_blend((x+1), res-y, 1, scale*r,scale*g,scale*b); 
          }
        }else std::cout <<"cannot open file...";
      file.close();
      im.close();

   std::cout << "\n\npngtest has finished. Take a look at the PNG images that have been created!\n";
   return 0;
}
