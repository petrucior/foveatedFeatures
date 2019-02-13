/**
 * \file multiFoveation.h
 *
 * \brief This file contains the prototype and implementation of structure multifoveated.
 *
 * \author 
 * Petrucio Ricardo Tavares de Medeiros \n
 * Universidade Federal do Rio Grande do Norte \n
 * Departamento de Computacao e Automacao Industrial \n
 * petrucior at gmail (dot) com
 *
 * \version 0.1
 * \date February 2016
 *
 * \copyright
 * Copyright (C) 2016, Petrúcio Ricardo <petrucior@gmail.com>
 * If you use this software for academic purposes, consider citing the related
 * paper: Rafael Beserra Gomes, Bruno Motta de Carvalho, Luiz Marcos Garcia 
 * Gonçalves, Visual attention guided features selection with foveated images,
 * Neurocomputing, Volume 120, 23 November 2013, Pages 34-44, ISSN 0925-2312,
 * http://dx.doi.org/10.1016/j.neucom.2012.10.033.
 *
 * This file is part of foveatedFeatures software.
 * This program is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later 
 * version. This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details. You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MULTIFOVEATION_H
#define MULTIFOVEATION_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <utility>
#include <vector>
#include "foveatedHessianDetector.h"
#include "opencv2/core/core.hpp"
//#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/xfeatures2d/nonfree.hpp"

/**
 * \struct MultiFoveation 
 *
 * \brief Struct for multiples foveae.
 */
struct MultiFoveation{

  //
  // Variables
  //
  std::vector<FoveatedHessianDetectorParams> params;

  //
  // Methods
  //

  /**
   * \fn MultiFoveation( )
   *
   * \brief Constructor multi foveae.
   */
  MultiFoveation();

  /**
   * \fn MultiFoveation(int foveae, Mat image, std::vector<String> ymlFile)
   *
   * \brief Constructor multi foveae.
   *
   * \param foveae - Number of foveae
   *        image - Image to be foveated
   *        ymlFile - Vector with yaml names.
   */
  MultiFoveation(int foveas, Mat image, std::vector<String> ymlFile);
  
  /**
   * \fn void updateParams(int fovea)
   *
   * \brief Update the params of fovea.
   *
   * \param fovea - The number of fovea
   */
  void updateParams(int fovea);
  
  /**
   * \fn void addFovea(Mat image, String ymlFile)
   *
   * \brief Add new struture of fovea.
   *
   * \param image - Image to be foveated 
   *        ymlFile - Configuration file of fovea
   */
  void addFovea(Mat image, String ymlFile);

  /**
   * \fn bool removeFovea(int fovea)
   *
   * \brief Remove the struture of fovea with id fovea.
   *
   * \param fovea - The number of fovea
   *
   * \return True if remove fovea and false otherwise.
   */
  bool removeFovea(int fovea);
  
  /**
   * \fn Point intersection(float k, int m, Size R, int fovea1, int fovea2);
   *
   * \brief Function for calculate the intersection between foveae.
   *
   * \param k - Level of fovea
   *        m - Number levels of fovea
   *        R - Size of image
   *        fovea1 - Fovea before processed
   *        fovea2 - Fovea to be processed
   *
   * \return Point of intersection between foveae.
   */
  Point intersection(float k, int m, Size R, int fovea1, int fovea2);
  
  /**
   * \fn Point directionPoint(float k, Point pointIntersection, int fovea1, int fovea2);
   *
   * \brief Function for calculate the vector.
   *
   * \param k - Level of fovea
   *        pointIntersection - The point of intersection to find with intersection function
   *        fovea1 - Fovea before processed
   *        fovea2 - Fovea to be processed
   *
   * \return Vector of direction of intersection.
   */
  Point directionPoint(float k, Point pointIntersection, int fovea1, int fovea2);
  
  /**
   * \fn void retirar(std::vector<Point>& region, std::vector<Point> sizeLevel, std::vector<Point> vectorPointIntersection, std::vector<Point> vectorDirectionIntersection);
   *
   * \brief Function for calculate the region to be processed.
   *
   * \param region - Pointer of region to be processed.
   *        sizeLevel - The size in y axis of level to be processed
   *        vectorPointIntersection - The vector of point of intersection to find with intersection function
   *        vectorDirectionIntersection - The vector point of intersection direction
   *
   * \return Vector of regions to be processed.
   */
  void retirar(std::vector<Point>& region, std::vector<Point> sizeLevel, std::vector<Point> vectorPointIntersection, std::vector<Point> vectorDirectionIntersection);
  
    
  // Publics Functions
  /**
   * \fn void extractKeypoints(Mat image, std::vector<KeyPoint>& _keypoint)
   *
   * \brief Function to extract keypoints
   *
   * \param image - Image processed
   *        _keypoint - address for keypoint pointer
   *
   * \return The keypoints extract of image
   */  
  void extractKeypoints(Mat image, std::vector<KeyPoint>& _keypoint);
  
  /**
   * \fn void drawLevels(Mat& image, std::vector<Scalar> colors)
   *
   * \brief Function to paint gride of levels 
   *
   * \param image - Address for image pointer
   *        colors - Vector of colors of grid
   *
   * \return An image with grid of levels
   */  
  void drawLevels(Mat& image, std::vector<Scalar> colors);
  
};

#endif

/**
 * \fn MultiFoveation( )
 *
 * \brief Constructor multi foveae.
 */
MultiFoveation::MultiFoveation(){
  params.clear();
}

/**
 * \fn MultiFoveation(int foveae, Mat image, std::vector<String> ymlFile)
 *
 * \brief Constructor multi foveae.
 *
 * \param foveae - Number of foveae
 *        image - Image to be foveated
 *        ymlFile - Vector with yaml names.
 */
MultiFoveation::MultiFoveation(int foveas, Mat image, std::vector<String> ymlFile){
  std::vector<int> delta;
  std::vector<int> size;
  std::vector<Point> pointsIntersection;
  std::vector<Point> pointsDirection;
  std::vector<Point> region;
  params.clear();
  for (int i = 0; i < foveas; i++){
    FoveatedHessianDetectorParams p(image.cols, image.rows, ymlFile[i]);
    params.push_back(p);
    params[i].foveaModel.init();
    int m = params[i].foveaModel.m;
    if ( i != 0 ){
      // Loop to processing levels
      for (int k = 0; k < m+1; k++){
	pointsIntersection.clear();
	pointsDirection.clear();
	// Loop to processing foveae processed
	for (int j = 0; j < i; j++){
	  Point point = Point(0, 0);
	  Point direction = Point(2, 2);
	  if ( k != 0 ){
	    point = intersection(k, m, Size(params[i].foveaModel.ux, params[i].foveaModel.uy), j, i);
	    direction = directionPoint(k, point, j, i);
	  }
	  pointsIntersection.push_back(point);
	  pointsDirection.push_back(direction);
	}
        // bubblesort
	for (unsigned int b1 = 0; b1 < pointsIntersection.size(); b1++){
	  for (unsigned int b2 = b1 + 1; b2 < pointsIntersection.size(); b2++){
	    if ( pointsIntersection[b2].x < pointsIntersection[b1].x ){
	      // swap
	      std::swap(pointsIntersection[b1], pointsIntersection[b2]);
	      std::swap(pointsDirection[b1], pointsDirection[b2]);
	    }
	  }
	}
	region.clear();
	std::vector<Point> s;
	s.push_back(Point(params[i].foveaModel.getDeltax(k), params[i].foveaModel.getDeltay(k)));
	s.push_back(Point(params[i].foveaModel.getDeltax(k) + params[i].foveaModel.getSizex(k), params[i].foveaModel.getDeltay(k) + params[i].foveaModel.getSizey(k)));
	retirar(region, s, pointsIntersection, pointsDirection);
	
	for (unsigned int v = 0; v < region.size(); v+=2){
	  delta.push_back(region[v].x); delta.push_back(region[v].y);
	  size.push_back(region[v+1].x); size.push_back(region[v+1].y);
	}
	params[i].foveaModel.setMultiFoveation(k, delta, size, i);
	delta.clear();
	size.clear();
      }
    }
  }
}

/**
 * \fn void updateParams(int fovea)
 *
 * \brief Update the params of fovea.
 *
 * \param fovea - The number of fovea
 */
void 
MultiFoveation::updateParams(int fovea){  
  std::vector<int> delta;
  std::vector<int> size;
  std::vector<Point> pointsIntersection;
  std::vector<Point> pointsDirection;
  std::vector<Point> region;
  for (unsigned int i = 1; i < params.size(); i++){
    int m = params[i].foveaModel.m;
    // Loop to processing levels
    for (int k = 0; k < m+1; k++){
      pointsIntersection.clear();
      pointsDirection.clear();
      // Loop to processing foveae processed
      for (unsigned int j = 0; j < i; j++){
	Point point = Point(0, 0);
	Point direction = Point(2, 2);
	if ( k != 0 ){
	  point = intersection(k, m, Size(params[i].foveaModel.ux, params[i].foveaModel.uy), j, i);
	  direction = directionPoint(k, point, j, i);
	}
	pointsIntersection.push_back(point);
	pointsDirection.push_back(direction);
      }
      // bubblesort
      for (unsigned int b1 = 0; b1 < pointsIntersection.size(); b1++){
	for (unsigned int b2 = b1 + 1; b2 < pointsIntersection.size(); b2++){
	  if ( pointsIntersection[b2].x < pointsIntersection[b1].x ){
	    // swap
	    std::swap(pointsIntersection[b1], pointsIntersection[b2]);
	    std::swap(pointsDirection[b1], pointsDirection[b2]);
	  }
	}
      }
      region.clear();
      std::vector<Point> s;
      s.push_back(Point(params[i].foveaModel.getDeltax(k), params[i].foveaModel.getDeltay(k)));
      s.push_back(Point(params[i].foveaModel.getDeltax(k) + params[i].foveaModel.getSizex(k), params[i].foveaModel.getDeltay(k) + params[i].foveaModel.getSizey(k)));
      retirar(region, s, pointsIntersection, pointsDirection);
	
      for (unsigned int v = 0; v < region.size(); v+=2){
	delta.push_back(region[v].x); delta.push_back(region[v].y);
	size.push_back(region[v+1].x); size.push_back(region[v+1].y);
      }
      params[i].foveaModel.setMultiFoveation(k, delta, size, i);
      delta.clear();
      size.clear();
    }
  }
}

/**
 * \fn void addFovea(Mat image, String ymlFile)
 *
 * \brief Add new struture of fovea.
 *
 * \param image - Image to be foveated 
 *        ymlFile - Configuration file of fovea
 */
void
MultiFoveation::addFovea(Mat image, String ymlFile){
  FoveatedHessianDetectorParams p(image.cols, image.rows, ymlFile);
  params.push_back(p);
  params[(int)params.size() - 1].foveaModel.init();
  updateParams((int)params.size() - 1);
}

/**
 * \fn bool removeFovea(int fovea)
 *
 * \brief Remove the struture of fovea with id fovea.
 *
 * \param fovea - The number of fovea
 *
 * \return True if remove fovea and false otherwise.
 */
bool 
MultiFoveation::removeFovea(int fovea){
  if ( ( (int)params.size() == 0 ) ||
       ( (int)params.size() - 1 < fovea ) ||
       ( fovea < 0 ) ) return false;
  else{
    params.erase(params.begin()+fovea);
    return true;
  }
}

/**
 * \fn Point intersection(float k, int m, Size R, int fovea1, int fovea2);
 *
 * \brief Function for calculate the intersection between foveae.
 *
 * \param k - Level of fovea
 *        m - Number levels of fovea
 *        R - Size of image
 *        fovea1 - Fovea before processed
 *        fovea2 - Fovea to be processed
 *
 * \return Point of intersection between foveae.
 */
Point 
MultiFoveation::intersection(float k, int m, Size R, int fovea1, int fovea2){
  FoveatedHessianDetectorParams f1, f2;
  f1 = params[fovea1];
  f2 = params[fovea2];
  Point p = Point(-1, -1);
  //
  // Implementação considera m1 == m2 ( necessita encontrar a equação para não depender desta condição )
  //
  int wmax = 0, wmin = 0; // Limit of fovea in projections
  int p1, p2; // fovea in projections
  //
  // Component x of intersection point
  //
  // wmax e wmin ( conditional ternary )
  max(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) == f2.foveaModel.fx+(R.width/2) ? wmax = f2.foveaModel.wx/2 : wmax = f1.foveaModel.wx/2;
  min(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) == f2.foveaModel.fx+(R.width/2) ? wmin = f2.foveaModel.wx/2 : wmin = f1.foveaModel.wx/2;
  // p1 e p2
  p1 = max( max(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) - wmax , min(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) + wmin );
  p2 = min( max(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) - wmax , min(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) + wmin );
  
  // Intersection between limits
  if ( ( (f1.foveaModel.fx+(R.width/2)-(f1.foveaModel.getSizex(m)/2) < f2.foveaModel.fx+(R.width/2)-(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)-(f1.foveaModel.getSizex(m)/2) < f2.foveaModel.fx+(R.width/2)+(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)+(f1.foveaModel.getSizex(m)/2) > f2.foveaModel.fx+(R.width/2)-(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)+(f1.foveaModel.getSizex(m)/2) < f2.foveaModel.fx+(R.width/2)+(f2.foveaModel.getSizex(m)/2)) ) ||
       ( (f1.foveaModel.fx+(R.width/2)-(f1.foveaModel.getSizex(m)/2) > f2.foveaModel.fx+(R.width/2)-(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)-(f1.foveaModel.getSizex(m)/2) < f2.foveaModel.fx+(R.width/2)+(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)+(f1.foveaModel.getSizex(m)/2) > f2.foveaModel.fx+(R.width/2)-(f2.foveaModel.getSizex(m)/2)) &&
	 (f1.foveaModel.fx+(R.width/2)+(f1.foveaModel.getSizex(m)/2) > f2.foveaModel.fx+(R.width/2)+(f2.foveaModel.getSizex(m)/2)) ) )
    std::swap(p1, p2);
  
  
  // DEBUG
  /*std::cout << "k = " << k << ", m = " << m << std::endl;
  std::cout << "f1.x = " << f1.foveaModel.fx+(R.width/2) << ", f2.x = " << f2.foveaModel.fx+(R.width/2) << std::endl;
  std::cout << "p1 = " << p1 << ", p2 = " << p2 << std::endl;
  std::cout << "wmax = " << wmax << ", wmin = " << wmin << std::endl;*/
  // Component x axis
  if ( min(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) == f2.foveaModel.fx+(R.width/2) ){ // Fóvea 2 mais próxima da origem
    p.x = ( k * p1 )/m;
  }
  else{ // Fóvea 1 mais próxima da origem
    p.x = ( (R.width * m) - (R.width * k) + (p2 * k) )/ m;
  }
  
  //
  // Component y of intersection point
  //
  // wmax e wmin ( conditional ternary )
  max(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) == f2.foveaModel.fy+(R.height/2) ? wmax = f2.foveaModel.wy/2 : wmin = f1.foveaModel.wy/2;
  min(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) == f2.foveaModel.fy+(R.height/2) ? wmin = f2.foveaModel.wy/2 : wmax = f1.foveaModel.wy/2;
  // p1 e p2
  p1 = max( max(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) - wmax , min(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) + wmin );
  p2 = min( max(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) - wmax , min(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) + wmin );

  // Intersection between limits
  if ( ( (f1.foveaModel.fy+(R.width/2)-(f1.foveaModel.getSizey(m)/2) < f2.foveaModel.fy+(R.width/2)-(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)-(f1.foveaModel.getSizey(m)/2) < f2.foveaModel.fy+(R.width/2)+(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)+(f1.foveaModel.getSizey(m)/2) > f2.foveaModel.fy+(R.width/2)-(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)+(f1.foveaModel.getSizey(m)/2) < f2.foveaModel.fy+(R.width/2)+(f2.foveaModel.getSizey(m)/2)) ) ||
       ( (f1.foveaModel.fy+(R.width/2)-(f1.foveaModel.getSizey(m)/2) > f2.foveaModel.fy+(R.width/2)-(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)-(f1.foveaModel.getSizey(m)/2) < f2.foveaModel.fy+(R.width/2)+(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)+(f1.foveaModel.getSizey(m)/2) > f2.foveaModel.fy+(R.width/2)-(f2.foveaModel.getSizey(m)/2)) &&
	 (f1.foveaModel.fy+(R.width/2)+(f1.foveaModel.getSizey(m)/2) > f2.foveaModel.fy+(R.width/2)+(f2.foveaModel.getSizey(m)/2)) ) )
    std::swap(p1, p2);

  // DEBUG
  /*std::cout << "f1.y = " << f1.foveaModel.fy << ", f2.y = " << f2.foveaModel.fy << std::endl;
  std::cout << "p1 = " << p1 << ", p2 = " << p2 << std::endl;
  std::cout << "wmax = " << wmax << ", wmin = " << wmin << std::endl;*/
  // Component y axis
  if ( min(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) == f2.foveaModel.fy+(R.height/2) ){ // Fóvea 2 mais próxima da origem
    p.y = ( k * p1 )/m;
  }
  else{ // Fóvea 1 mais próxima da origem
    p.y = ( (R.height * m) - (R.height * k) + (p2 * k) )/m;
  }

  
  // Verify if the fovea1 is equal fovea2
  if ( f1.foveaModel.fx+(R.width/2) == f2.foveaModel.fx+(R.width/2) ) p.x = f1.foveaModel.getDeltax(k);
  if ( f1.foveaModel.fy+(R.height/2) == f2.foveaModel.fy+(R.height/2) ) p.y = f1.foveaModel.getDeltay(k);

  // Verify if the intersection is between layers
  if ( ( !(p.x < f1.foveaModel.fx+(R.width/2) + (f1.foveaModel.getSizex(k)/2)) &&
	 !(p.x > f2.foveaModel.fx+(R.width/2) - (f2.foveaModel.getSizex(k)/2)) ) ||
       ( !(p.x > f1.foveaModel.fx+(R.width/2) - (f1.foveaModel.getSizex(k)/2)) &&
	 !(p.x < f2.foveaModel.fx+(R.width/2) + (f2.foveaModel.getSizex(k)/2)) ) ||
       ( !(p.y < f1.foveaModel.fy+(R.height/2) + (f1.foveaModel.getSizey(k)/2)) &&
	 !(p.y > f2.foveaModel.fy+(R.height/2) - (f2.foveaModel.getSizey(k)/2)) ) ||
       ( !(p.y > f1.foveaModel.fy+(R.height/2) - (f1.foveaModel.getSizey(k)/2)) &&
	 !(p.y < f2.foveaModel.fy+(R.height/2) + (f2.foveaModel.getSizey(k)/2)) ) ){
    // Not exist intersection between layers
    p.x = -1;
    p.y = -1;
  }
  
  //std::cout << "fovea1 = " << fovea1 << " e fovea2 = " << fovea2 << std::endl;
  //std::cout << p.x << " | " << p.y << std::endl;
  
  return p;
}

/**
 * \fn Point directionPoint(float k, Point pointIntersection, int fovea1, int fovea2);
 *
 * \brief Function for calculate the vector.
 *
 * \param k - Level of fovea
 *        pointIntersection - The point of intersection to find with intersection function
 *        fovea1 - Fovea before processed
 *        fovea2 - Fovea to be processed
 *
 * \return Vector of direction of intersection.
 */
Point
MultiFoveation::directionPoint(float k, Point pointIntersection, int fovea1, int fovea2){
  /*
    ------------------------------
    (-1, 1)  -- (1, 1)
    (-1, -1) -- (1, -1)
    without intersection = (0, 0)
    Without shifting = (2, 2)
    ------------------------------
   */
  FoveatedHessianDetectorParams f1, f2;
  f1 = params[fovea1];
  f2 = params[fovea2];
  // Center of levels
  Point centerF1 = Point( (int)(f1.foveaModel.getDeltax(k) + (f1.foveaModel.getDeltax(k)+f1.foveaModel.getSizex(k)))/2, 
			  (int)(f1.foveaModel.getDeltay(k) + (f1.foveaModel.getDeltay(k)+f1.foveaModel.getSizey(k)))/2 );
  Point centerF2 = Point( (int)(f2.foveaModel.getDeltax(k) + (f2.foveaModel.getDeltax(k)+f2.foveaModel.getSizex(k)))/2, 
			  (int)(f2.foveaModel.getDeltay(k) + (f2.foveaModel.getDeltay(k)+f2.foveaModel.getSizey(k)))/2 );
  // Point central in each level
  Point v = Point( centerF2.x - centerF1.x, centerF2.y - centerF1.y);
  Point direcao = Point(0, 0);
  // Without intersection
  if ( (pointIntersection.x == -1) && (pointIntersection.y == -1) ){
    return direcao;
  }
  // Horizontal shifting
  if ( (v.x > 0) && (v.y == 0) ){ // Shift out of origin
    //std::cout << "sentido horizontal negativo" << std::endl;
    direcao = Point(-1, -1);
  }
  if ( (v.x < 0) && (v.y == 0) ){ // Shift in of origin
    //std::cout << "sentido horizontal positivo" << std::endl;
    direcao = Point(1, -1);
  }

  // Vertical shifting
  if ( (v.x == 0) && (v.y > 0) ){ // Shift out of origin
    //std::cout << "sentido vertical cima" << std::endl;
    direcao = Point(1, 1);
  }
  if ( (v.x == 0) && (v.y < 0) ){ // Shift in of origin
    //std::cout << "sentido vertical baixo" << std::endl;
    direcao = Point(1, -1);
  }
  
  // Diagonal shifting
  if ( (v.x > 0) && (v.y < 0) ){ // Deslocamento para o sentido sudoeste
    //std::cout << "sentido sudoeste" << std::endl;
    direcao = Point(-1, -1);
  }
  if ( (v.x > 0) && (v.y > 0) ){ // Deslocamento para o sentido noroeste
    //std::cout << "sentido noroeste" << std::endl;
    direcao = Point(-1, 1);
  }
  if ( (v.x < 0) && (v.y < 0) ){ // Deslocamento para o sentido sudeste
    //std::cout << "sentido sudeste" << std::endl;
    direcao = Point(1, -1);
  }
  if ( (v.x < 0) && (v.y > 0) ){ // Deslocamento para o sentido nordeste
    //std::cout << "sentido nordeste" << std::endl;
    direcao = Point(1, 1);
  }
  // Without shifting
  if ( (v.x == 0) && (v.y == 0) ){
    direcao = Point(2, 2);
  }
  
  return direcao;
}

/**
 * \fn void retirar(std::vector<Point>& region, std::vector<Point> sizeLevel, std::vector<Point> vectorPointIntersection, std::vector<Point> vectorDirectionIntersection);
 *
 * \brief Function for calculate the region to be processed.
 *
 * \param region - Pointer of region to be processed.
 *        sizeLevel - The size in x, y axis of level to be processed 
 *        vectorPointIntersection - The vector of point of intersection to find with intersection function
 *        vectorDirectionIntersection - The vector point of intersection direction
 *
 * \return Vector of regions to be processed.
 */
void 
MultiFoveation::retirar(std::vector<Point>& region, std::vector<Point> sizeLevel, std::vector<Point> vectorPointIntersection, std::vector<Point> vectorDirectionIntersection){
  // Loop for detect region without shifting
  for (unsigned int p = 0; p < vectorDirectionIntersection.size(); p++){
    if ( ( vectorDirectionIntersection[p].x == 2 ) && ( vectorDirectionIntersection[p].y == 2 ) ){
      region.push_back(Point(-1, -1));
      region.push_back(Point(-1, -1));
      return;
    }
  }
  
  std::vector<int> minLimit(vectorPointIntersection.size()+1, sizeLevel[0].y);
  std::vector<int> maxLimit(vectorPointIntersection.size()+1, sizeLevel[1].y); // sizeLevel[0].y + sizeLevel[1].y
  for ( unsigned int p1 = 0; p1 < vectorPointIntersection.size(); p1++ ){    
    // Sentido noroeste
    if ( ( vectorDirectionIntersection[p1].x == -1 ) && ( vectorDirectionIntersection[p1].y == -1 ) ){
      maxLimit[p1] = vectorPointIntersection[p1].y;
      for ( unsigned int p2 = 0; p2 < p1; p2++ ){
	if ( maxLimit[p1] < maxLimit[p2] ) maxLimit[p2] = maxLimit[p1];
      }
    }
    // Sentido sudoeste
    if ( ( vectorDirectionIntersection[p1].x == -1 ) && ( vectorDirectionIntersection[p1].y == 1 ) ){
      minLimit[p1] = vectorPointIntersection[p1].y;
      for ( unsigned int p2 = 0; p2 < p1; p2++ ){
	if ( minLimit[p1] > minLimit[p2] ) minLimit[p2] = minLimit[p1];
      }
    }
    // Sentido nordeste
    if ( ( vectorDirectionIntersection[p1].x == 1 ) && ( vectorDirectionIntersection[p1].y == -1 ) ){
      maxLimit[p1+1] = vectorPointIntersection[p1].y;
      for ( unsigned int p2 = p1; p2 < vectorPointIntersection.size(); p2++ ){
	if ( maxLimit[p1+1] < maxLimit[p2+1] ) maxLimit[p2+1] = maxLimit[p1+1];
      }
    }
    // Sentido sueste
    if ( ( vectorDirectionIntersection[p1].x == 1 ) && ( vectorDirectionIntersection[p1].y == 1 ) ){
      minLimit[p1+1] = vectorPointIntersection[p1].y;
      for ( unsigned int p2 = p1; p2 < vectorPointIntersection.size(); p2++ ){
	if ( minLimit[p1+1] > minLimit[p2+1] ) minLimit[p2+1] = minLimit[p1+1];
      }
    }
  }

  // DEBUG
  /*for (unsigned int i = 0; i < vectorPointIntersection.size()+1; i++){
    std::cout <<  minLimit[i] << " - " << maxLimit[i] << std::endl;
  }*/
  
  /*for (unsigned int i = 0; i < vectorPointIntersection.size(); i++){
    std::cout << "(" << vectorPointIntersection[i].x << ", " << vectorPointIntersection[i].y << ")" << std::endl;
    std::cout << "(" << vectorDirectionIntersection[i].x << ", " << vectorDirectionIntersection[i].y << ")" << std::endl;
  }*/

  // reading left to right
  // Inicio
  // Delta
  region.push_back(Point(sizeLevel[0].x, minLimit[0]));
  // Size
  region.push_back(Point(vectorPointIntersection[0].x - sizeLevel[0].x, maxLimit[0] - minLimit[0]));
  // Meio
  for ( unsigned int i = 1; i < minLimit.size() - 1; i++ ){
    // Delta
    region.push_back(Point(vectorPointIntersection[i-1].x, minLimit[i]));
    // Size
    region.push_back(Point(vectorPointIntersection[i].x - vectorPointIntersection[i-1].x, maxLimit[i] - minLimit[i])); 
  }
  // Fim
  // Delta
  region.push_back(Point(vectorPointIntersection[vectorPointIntersection.size() - 1].x, minLimit[minLimit.size() - 1]));
  // SIze
  region.push_back(Point(sizeLevel[1].x - vectorPointIntersection[vectorPointIntersection.size() - 1].x, maxLimit[maxLimit.size() - 1] - minLimit[minLimit.size() - 1]));

}


/**
 * \fn void extractKeypoints(Mat image, std::vector<KeyPoint>& _keypoint);
 *
 * \brief Function to extract keypoints
 *
 * \param image - Image processed
 *        _keypoint - address for keypoint pointer
 *
 * \return The keypoints extract of image
 */
void 
MultiFoveation::extractKeypoints(Mat image, std::vector<KeyPoint>& _keypoint){
  for (unsigned int i = 0; i < params.size(); i++){
    vector<KeyPoint> keypoints;
    foveatedHessianDetector(image, Mat(), keypoints, params[i]);
    for (unsigned int k = 0; k < keypoints.size(); k++)
      _keypoint.push_back(keypoints[k]);
  }
}

/**
 * \fn void drawLevels(Mat& image, std::vector<Scalar> colors)
 *
 * \brief Function to paint gride of levels 
 *
 * \param image - Address for image pointer
 *        colors - Vector of colors of grid
 *
 * \return An image with grid of levels
 */  
void 
MultiFoveation::drawLevels(Mat& image, std::vector<Scalar> colors){
  for (unsigned int i = 0; i < params.size(); i++)
    drawMultiFoveatedLevels(image, params[i], i, colors[i]);
}
