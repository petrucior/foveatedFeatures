#include <iostream>
#include "../foveatedHessianDetector.h"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

int main(int argc, char** argv){
  Mat image = imread(argv[1]);
  if ( image.empty() ) return -1;
  
  return 0;
}
