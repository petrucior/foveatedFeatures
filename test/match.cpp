#include "../foveatedHessianDetector.h"
#include <stdio.h>
#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"

using namespace std;
using namespace cv;
using namespace cv::xfeatures2d;

static void help()
{
    printf("\nThis program demonstrates using features2d detector, descriptor extractor and simple matcher\n"
            "Using the SURF desriptor:\n"
            "\n"
            "Usage:\n matcher_simple <image1> <image2>\n");
}

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        help();
        return -1;
    }

    Mat img1 = imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
    Mat img2 = imread(argv[2], CV_LOAD_IMAGE_GRAYSCALE);
    if(img1.empty() || img2.empty())
    {
        printf("Can't read one of the images\n");
        return -1;
    }

    // detecting keypoints
    Ptr<SURF> detector = SURF::create(400);
    vector<KeyPoint> keypoints1, keypoints2;
    //SurfFeatureDetector detector(400);
    //vector<KeyPoint> keypoints1, keypoints2;
    

    FoveatedHessianDetectorParams params(img1.cols, img1.rows, "fovea1.yml");
    foveatedHessianDetector(img1, Mat(), keypoints1, params);
    foveatedHessianDetector(img2, Mat(), keypoints2, params);
    //detector.detect(img1, keypoints1);
    //detector.detect(img2, keypoints2);
    
    // computing descriptors
    //Ptr<SURF> extractor = SURF::create();
    //SurfDescriptorExtractor extractor;
    Mat descriptors1, descriptors2;
    detector->detectAndCompute( img1, noArray(), keypoints1, descriptors1, true );
    // Para fazer o match nao eh necessario fovear a segunda imagem e por esta razao
    // o ultimo argumento ta dizendo que nao utilize os keypoints, ou seja, ele precise
    // detectar e computar.
    detector->detectAndCompute( img2, noArray(), keypoints2, descriptors2, false );
    //extractor->compute(img1, keypoints1, descriptors1);
    //extractor->compute(img2, keypoints2, descriptors2);
    //extractor.compute(img1, keypoints1, descriptors1);
    //extractor.compute(img2, keypoints2, descriptors2);

    // matching descriptors
    Ptr<BFMatcher> matcher = BFMatcher::create(NORM_L2);
    //Ptr<DescriptorMatcher> matcher = DescriptorMatcher::create(DescriptorMatcher::BRUTEFORCE);
    vector< DMatch > matches;
    matcher->match( descriptors1, descriptors2, matches);
    //BFMatcher matcher(NORM_L2);
    //vector<DMatch> matches;
    //matcher.match(descriptors1, descriptors2, matches);

    // drawing the results
    namedWindow("matches", 1);
    Mat img_matches;
    //drawKeypoints(img1, keypoints1, img_matches, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    //drawKeypoints(img2, keypoints2, img_matches, Scalar::all(-1), DrawMatchesFlags::DRAW_RICH_KEYPOINTS);
    drawMatches(img1, keypoints1, img2, keypoints2, matches, img_matches);
    imshow("matches", img_matches);
    waitKey(0);

    return 0;
}
