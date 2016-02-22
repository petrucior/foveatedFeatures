#include <iostream>
#include "../foveatedHessianDetector.h"
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

int main(int argc, char** argv){
  // Image e template
  Mat image = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
  Mat temp = imread(argv[2], CV_LOAD_IMAGE_GRAYSCALE);
  if ( ( image.empty() ) || ( temp.empty() ) ) return -1;
  
  // detecting keypoints
  SurfFeatureDetector detector(400);
  vector<KeyPoint> keypoints1, keypoints2;
  vector<Point2f> modelPoints, imgPoints;
  
  FoveatedHessianDetectorParams params(image.cols, image.rows, "fovea3.yml");
  
  unsigned int maximo = 0;
  vector<unsigned int> maximos;
  vector<int> pointx, pointy;
  for ( int i = 0; i < image.cols; i++ ){
    for ( int j = 0; j < image.rows; j++ ){
      pointx.push_back(i);
      pointy.push_back(j);
      
      // Setando a posição da fóvea
      params.foveaModel.setFovea(j, i);
      params.foveaModel.fixFovea();
      
      foveatedHessianDetector(image, Mat(), keypoints1, params);
      //foveatedHessianDetector(temp, Mat(), keypoints2, params);
      detector.detect(temp, keypoints2);
      
      // computing descriptors
      SurfDescriptorExtractor extractor;
      Mat descriptors1, descriptors2;
      extractor.compute(image, keypoints1, descriptors1);
      extractor.compute(temp, keypoints2, descriptors2);

      // matching descriptors
      BFMatcher matcher(NORM_L2);
      vector<DMatch> matches;
      vector<DMatch> good_matches;
      if ( !descriptors2.empty() )
	matcher.match(descriptors1, descriptors2, matches);
      
      if ( matches.size() > 4 ){
	modelPoints.clear();
	imgPoints.clear();
	for(unsigned int k = 0; k < matches.size(); k++){
	  DMatch m = matches[k];
	  imgPoints.push_back(keypoints1[m.trainIdx].pt);
	  modelPoints.push_back(keypoints2[m.queryIdx].pt);
	}
	Mat mask;
	Mat H = findHomography( imgPoints, modelPoints, RANSAC, 4, mask );
	for (int k = 0; k < mask.rows; k++){
	  //std::cout << (int)mask.at<uchar>(0, k) << std::endl;
	  if ( (int)mask.at<uchar>(0, k) != 0 ) // Inliers
	    good_matches.push_back(matches[k]);
	}
	if ( maximo < good_matches.size() ) maximo = good_matches.size();
	maximos.push_back(good_matches.size());
	//std::cout << i << "  " << j << "  " << (100*good_matches.size())/matches.size() << std::endl;
	//std::cout << i << "  " << j << "  " << good_matches.size() << std::endl;
      }
      else{
	maximos.push_back(0);
        //std::cout << i << "  " << j << "  " << 0 << std::endl;
      }
      
      good_matches.clear();
    }
  }
  
  //std::cout << maximo << std::endl;
  for (unsigned int k = 0; k < maximos.size(); k++)
    std::cout << pointx[k] << " " << pointy[k] << " " << (float)maximos[k]/maximo << std::endl;
  
  // drawing the results
  /*namedWindow("matches", 1);
  Mat img_matches;
  drawMatches(image, keypoints1, temp, keypoints2, matches, img_matches);
  imshow("matches", img_matches);
  waitKey(0);*/

  return 0;
}
