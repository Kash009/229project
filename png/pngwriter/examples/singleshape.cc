#include <pngwriter.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
#include <ctime>
#include <cstdlib>

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
   string filename = "im2.png";
   pngwriter im(51,51,0.0,filename.c_str()); 
  
   std::ifstream file("../../pbrt-v2/src/svmdata/32singleshape");
    unsigned int id = 0;
    double value = 0;
    char formatting_tab;
    char formatting_line;
    if (file.is_open())
      {
        while((file >> id) && (file >> formatting_tab) && (file >> value) && (file >> formatting_line) ){
          int x = id/51;
          int y = id - x*51;
          im.plot_blend(x, y, 1, value, value, value); 

        }
      }else std::cout <<"cannot open file...";
    file.close();
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
