#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
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

  std::vector<Scalar> colors;
  srand (time(NULL)); // Initialize random seed
  for (int i = 0; i < argc-2; i++){
    int r = rand() % 256; // 0 - 255
    int g = rand() % 256; // 0 - 255
    int b = rand() % 256; // 0 - 255
    colors.push_back(Scalar(b, g, r));
  }

  // Eliminar a extração das features do código, pois estou testando o algoritmo
  // de regiões
  //return -1;
  
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
      drawMultiFoveatedLevels(outputImg, foveas.getParams(i), i, colors[i]);
    
    imshow("keypoints", outputImg);
    
    char key = waitKey(33);
    if ( key == 'q' ) break;
    
  }
  
  return 0;
}
