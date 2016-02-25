#include <iostream>
#include "../multiFoveation.h"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

using namespace cv;

static void helpExtract(){
  printf("\nThis program demonstrates using multi foveated features2d detector and descriptor extractor\n"
	 "Using the SURF desriptor:\n"
	 "\n"
	 "Usage:\n multi <image1> <fovea yml files>\n");
}

static void on_mouse(int event, int x, int y, int flags, void *_param){
  FoveatedHessianDetectorParams *params = (FoveatedHessianDetectorParams *) _param;
  params->foveaModel.setFovea(x, y);
  params->foveaModel.fixFovea();
}


int main(int argc, char** argv){
  //MultiFoveation m(2);
  //m.intersection(0, 2);
  //  Mat image = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  //  if ( image.empty() ) return -1;
  //  std::vector<String> s;
  //  s.push_back(argv[2]);
  //  s.push_back(argv[3]);
  //  MultiFoveation foveas(2, image, s);
  //FoveatedHessianDetectorParams p = foveas.getParams(0);
  //std::cout << p.nOctaveLayers << std::endl;
  //std::cout << p.foveaModel.fx << std::endl;
  //Point pointIntersection = foveas.intersection(1, 1, image.size(), 0, 1);
  //std::cout << "Intersection is realized in: (" << pointIntersection.x << ", " << pointIntersection.y << ")" << std::endl;

  if ( argc < 3 ){
    helpExtract();
    return -1;
  }
  
  Mat image = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  if ( image.empty() ){
    printf("Can't read one of the images\n");
    return -1;
  }

  std::vector<String> s;
  for (int i = 2; i < argc; i++)
    s.push_back(argv[i]);
  
  MultiFoveation foveas(argc-2, image, s);
  std::cout << argc-2 << " Fóveas criadas" << std::endl;
  namedWindow("keypoints", 1);
  
  while(true){
    vector<KeyPoint> keypointSave;
    keypointSave.clear();
    // Detecting keypoints
    for (int i = 0; i < argc-2; i++){
      vector<KeyPoint> keypoints;
      foveatedHessianDetector(image, Mat(), keypoints, foveas.getParams(i));
      for (unsigned int k = 0; k < keypoints.size(); k++)
	keypointSave.push_back(keypoints[k]);
    }
    
    // Computing descriptors
    SurfDescriptorExtractor extractor;
    Mat descriptors;
    extractor.compute(image, keypointSave, descriptors);
    
    // Drawing the results
    Mat outputImg;
    drawKeypoints(image, keypointSave, outputImg, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    for (int i = 0; i < argc-2; i++)
      (i == 0) ? drawMultiFoveatedLevels(outputImg, foveas.getParams(i), Scalar(255, 255, 255)) : drawMultiFoveatedLevels(outputImg, foveas.getParams(i), Scalar(0, 0, 255));
    
    imshow("keypoints", outputImg);
    
    char key = waitKey(33);
    if ( key == 'q' ) break;
    
  }
  
  return 0;
}
