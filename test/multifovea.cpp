#include <stdio.h>
#include "../multiFoveation.h"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

using namespace cv;

struct strutureUnion{
  MultiFoveation* m;
  int* positionParam;
  void setStrutureUnion(MultiFoveation *_m, int *_positionParam){
    m = _m;
    positionParam = _positionParam;
  }
  MultiFoveation* getStrutureMultiFoveation(){ return m; }
  int* getStrutureParam(){ return positionParam; }
  void setStrutureMultiFoveation(MultiFoveation *_m){ m = _m; }
  void setStrutureParam(int *_positionParam){ positionParam = _positionParam; }
};

static void helpExtract(){
  printf("\nThis program demonstrates using multi foveated features2d detector and descriptor extractor\n"
	 "Using the SURF desriptor:\n"
	 "\n"
	 "Usage:\n multi <image1> <fovea yml files>\n");
}

static void on_mouse(int event, int x, int y, int flags, void *_param){
  strutureUnion* params = (strutureUnion *) _param;
  
  int* p = (int *) params->getStrutureParam();
  // Update fovea position
  params->getStrutureMultiFoveation()->params[*p].foveaModel.setFovea(x, y);
  // Update intersections
  params->getStrutureMultiFoveation()->updateParams(*p);
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

  // Foveae params
  int indice = 0;
  
  // Struture with multifoveation and parameter address
  strutureUnion su;
  su.setStrutureUnion(&foveas, &indice);

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
  
  int controlFoveaMove = 0;
  
  vector<KeyPoint> keypointSave;
  while(true){
    keypointSave.clear();
    int64 t = cv::getTickCount();
    foveas.extractKeypoints(image, keypointSave);
    t = cv::getTickCount() - t;
    std::cout << "Feature extraction = " << t*1000/cv::getTickFrequency() << " milliseconds" << std::endl;
    
    // Computing descriptors
    SurfDescriptorExtractor extractor;
    Mat descriptors;
    t = cv::getTickCount();
    extractor.compute(image, keypointSave, descriptors);
    t = cv::getTickCount() - t;
    std::cout << "Feature description = " << t*1000/cv::getTickFrequency() << " milliseconds" << std::endl;
    
    // Drawing the results
    Mat outputImg;
    drawKeypoints(image, keypointSave, outputImg, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    foveas.drawLevels(outputImg, colors);
    
    imshow("keypoints", outputImg);
    
    char key = waitKey(33);
    if ( key == 'q' ) break;
    if ( key == 'm' ){ // Control foveae
      ( controlFoveaMove + 1 < (int)foveas.params.size() ) ? controlFoveaMove++ : controlFoveaMove = 0;
      su.setStrutureParam(&controlFoveaMove);
    }
    if ( key == 'a' ){ // Add foveae
      String file;
      std::cin >> file;
      foveas.addFovea(image, file);
    }
    if ( key == 'r' ){ // Remove foveae
      foveas.removeFovea(controlFoveaMove);
    }
    
    cvSetMouseCallback("keypoints", &on_mouse, &su);
    
  }
  
  return 0;
}
