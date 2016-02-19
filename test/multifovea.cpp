#include <iostream>
#include "../multiFoveation.h"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

int main(int argc, char** argv){
  //MultiFoveation m(2);
  //m.intersection(0, 2);
  Mat image = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  if ( image.empty() ) return -1;
  std::vector<String> s;
  s.push_back(argv[2]);
  s.push_back(argv[3]);
  MultiFoveation foveas(2, image, s);
  //FoveatedHessianDetectorParams p = foveas.getParams(0);
  //std::cout << p.nOctaveLayers << std::endl;
  //std::cout << p.foveaModel.fx << std::endl;
  //Point pointIntersection = foveas.intersection(1, 1, image.size(), 0, 1);
  //std::cout << "Intersection is realized in: (" << pointIntersection.x << ", " << pointIntersection.y << ")" << std::endl;
  return 0;
}
