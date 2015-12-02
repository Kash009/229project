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

int main()
{  
  ////////////////////////////////////////////////////////////////////////////////////
   /*  two.png
    *  This will be a 300x300 image with a black background, it will be called two.png.
    * Note that we are using 0.0 as the background colour and that we are using a string
    * type object as the filename, which can convert itself into a const char * with 
    * filename.c_str().
    * */
   std::cout << "Generating im.png...";
   string filename = "im.png";
   pngwriter im(101,101,0.0,filename.c_str()); 
  
   std::ifstream file("/Users/winnielin/Schoolwork/2015/Fall/CS229/Project/code/pbrt-v2/src/samplers/mldata/VarY3");
    int id = 0;
    double value = 0;
    char formatting_tab;
    char formatting_line;
    int x,y;
    double a,c;
    int debug_count = 1;
    int drawn[101][101];
    for(int i = 0; i < 101;i++){
      for(int j = 0; j < 101; j++){
        drawn[i][j] = 0;
      }
    }
    vector<int>vec;

    if (file.is_open())
      {
//        while((file >> id) && (file >> a) && (file >> value) && (file >> c)){
        while((file >> id) && (file >> value)){
          x = id /101;
          y = id - 101*x;
          value /=20;//for yvariance
          //std::cout<<x <<" " <<y<<"\n";
          drawn[x][y] +=1;
          debug_count++;
          im.plot_blend((x+1), 101-y, 1, value, value, value); 
          vec.push_back(id);
        }
      }else std::cout <<"cannot open file...";
    file.close();
    std::sort (vec.begin(), vec.end());
    int size = vec.size();
    for(int i = 0; i < size; i++){
      cout << vec.at(i) <<"\n";
    }
    for(int i = 0; i < 101;i++){
      for(int j = 0; j < 101; j++){
        std::cout<<drawn[i][j];
      }
      std::cout<<"\n";
    }
   /* Now, just as an example, we will use a lower compression on this image.
    * The default is 6 (from 0 to 9) and we will set it to 3. The lower the compression used
    * the faster the image will be close()d. Complex images will take longer to 
    * close() than simple ones.
    * */
   //two.setcompressionlevel(3);
   std::cout << " done. Writing to disk...";
   im.close();
   std::cout << " done." << std::endl;
    

   std::cout << "\n\npngtest has finished. Take a look at the PNG images that have been created!\n";
   return 0;
}
