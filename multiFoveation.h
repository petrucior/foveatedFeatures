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
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/nonfree/nonfree.hpp"

/**
 * \struct MultiFoveation 
 *
 * \brief Struct for múltiples foveae.
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
   * \fn void updateParams(int fovea, Size R)
   *
   * \brief Update the params of fovea.
   *
   * \param fovea - The number of fovea
   *        R - Size of image
   */
  void updateParams(int fovea, Size R);
  
  /**
   * \fn FoveatedHessianDetectorParams getParams(int fovea)
   *
   * \brief Get params of fovea.
   *
   * \param fovea - The number of fovea
   *
   * \return The structure of params foveated.
   */
  FoveatedHessianDetectorParams getParams(int fovea);

  /**
   * \fn std::vector<FoveatedHessianDetectorParams> getVectorParams()
   *
   * \brief Get the vector with fovea params.
   *
   * \return The vector structure of params foveated.
   */
  std::vector<FoveatedHessianDetectorParams> getVectorParams();
  
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
   * \fn std::vector<Point> limitProcessing(float k, Point pointIntersection, int fovea1, int fovea2);
   *
   * \brief Function for calculate the limits of processing.
   *
   * \param k - Level of fovea
   *        pointIntersection - The point of intersection to find with intersection function
   *        fovea1 - Fovea before processed
   *        fovea2 - Fovea to be processed
   *
   * \return Vector of limits of processing
   */
  std::vector<Point> limitProcessing(float k, Point pointIntersection, int fovea1, int fovea2);
  
  /**
   * \fn bool verifyRegion(int position, std::vector<Point> region, Point v);
   *
   * \brief Verify if vector v is in delimited region
   *
   * \param position - Position in region vector
   *        region - Vectors os delta and size of region
   *        v - Vector analised.
   *
   * \return True if v is inside region and false otherwise
   */
  bool verifyRegion(int position, std::vector<Point> region, Point v);
  
  /**
   * \fn std::vector<Point> updateLimit(std::vector<Point> limits);
   *
   * \brief Function for calculate the limit of processing
   *
   * \param limits - Vector of limits for update
   *
   * \return Vector of limits for region to be processed
   */
  std::vector<Point> updateLimit(std::vector<Point> limits);

};

#endif

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
  std::vector<Point> limits;
  std::vector<Point> limit;
  params.clear();
  for (int i = 0; i < foveas; i++){
    FoveatedHessianDetectorParams p(image.cols, image.rows, ymlFile[i]);
    params.push_back(p);
    params[i].foveaModel.init();
    int m = params[i].foveaModel.m;
    if ( i != 0 ){
      // Loop to processing levels
      for (int k = 0; k < m+1; k++){
	limits.clear();
	limit.clear();
	// Loop to processing foveae processed
	for (int j = 0; j < i; j++){
	  Point pontos = intersection(k, m, image.size(), j, i);
	  limits = limitProcessing(k, pontos, j, i);
	}
	if (k == 0){
	  for (unsigned int pts = 0; pts < limits.size(); pts+=2){
	    limits[pts] = Point(0, 0);
	    limits[pts+1] = Point(0, 0);
	  }
	}
	limit = updateLimit(limits);
        // Update delta and size
        for (unsigned int v = 0; v < limit.size(); v+=2){
	  delta.push_back(limit[v].x); delta.push_back(limit[v].y);
	  size.push_back(limit[v+1].x); size.push_back(limit[v+1].y);
	  //std::cout << "(" << limit[v].x << ", " << limit[v].y << ")" << std::endl;
	  //std::cout << "(" << limit[v+1].x << ", " << limit[v+1].y << ")" << std::endl;
	}
	params[i].foveaModel.setMultiFoveation(k, delta, size);
	delta.clear();
	size.clear();
      }
    }
  }
  
  // DEBUG
  /*for (unsigned int i = 0; i < limits.size(); i++){
    Point p = limits[i];
    std::cout << p.x << " :::: " << p.y << std::endl;
  }*/
}


/**
 * \fn void updateParams(int fovea, Size R)
 *
 * \brief Update the params of fovea.
 *
 * \param fovea - The number of fovea
 *        R - Size of image
 */
void 
MultiFoveation::updateParams(int fovea, Size R){
  std::vector<int> delta;
  std::vector<int> size;
  std::vector<Point> limits;
  std::vector<Point> limit;
  // Loop to processing params
  for (unsigned int p = fovea; p < params.size(); p++){
    int m = params[p].foveaModel.m;
    // Loop to processing levels
    for (int k = 0; k < m+1; k++){
      limits.clear();
      limit.clear();
      // Loop to processing foveae processed
      for (unsigned int j = 0; j < p; j++){
	Point pontos = intersection(k, m, R, j, p);
	limits = limitProcessing(k, pontos, j, p);
      }
      if (k == 0){
	for (unsigned int pts = 0; pts < limits.size(); pts+=2){
	  limits[pts] = Point(0, 0);
	  limits[pts+1] = Point(0, 0);
	}
      }
      limit = updateLimit(limits);
      // Update delta and size
      for (unsigned int v = 0; v < limit.size(); v+=2){
	delta.push_back(limit[v].x); delta.push_back(limit[v].y);
	size.push_back(limit[v+1].x); size.push_back(limit[v+1].y);
	//std::cout << "(" << limit[v].x << ", " << limit[v].y << ")" << std::endl;
	//std::cout << "(" << limit[v+1].x << ", " << limit[v+1].y << ")" << std::endl;
      }
      params[p].foveaModel.setMultiFoveation(k, delta, size);
      delta.clear();
      size.clear();
    }
  }
}

/**
 * \fn FoveatedHessianDetectorParams getParams(int fovea)
 *
 * \brief Get params of fovea.
 *
 * \param fovea - The number of fovea
 *
 * \return The structure of params foveated.
 */
FoveatedHessianDetectorParams 
MultiFoveation::getParams(int fovea){
  return params[fovea];
  //return params[fovea];
}

/**
 * \fn std::vector<FoveatedHessianDetectorParams> getVectorParams()
 *
 * \brief Get the vector with fovea params.
 *
 * \return The vector structure of params foveated.
 */
std::vector<FoveatedHessianDetectorParams> 
MultiFoveation::getVectorParams(){
  return params;
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
       ( !(p.y < f1.foveaModel.fy+(R.height/2) + (f1.foveaModel.getSizey(k)/2)) &&
	 !(p.y > f2.foveaModel.fy+(R.height/2) - (f2.foveaModel.getSizey(k)/2)) ) ){
    // Not exist intersection between layers
    p.x = -1;
    p.y = -1;
  }
  
  //std::cout << "fovea1 = " << fovea1 << " e fovea2 = " << fovea2 << std::endl;
  //std::cout << p.x << " | " << p.y << std::endl;
  
  return p;
}

/**
 * \fn std::vector<Point> limitProcessing(float k, Point pointIntersection, int fovea1, int fovea2);
 *
 * \brief Function for calculate the limits of processing.
 *
 * \param k - Level of fovea
 *        pointIntersection - The point of intersection to find with intersection function
 *        fovea1 - Fovea before processed
 *        fovea2 - Fovea to be processed
 *
 * \return Vector of limits of processing
 */
std::vector<Point> 
MultiFoveation::limitProcessing(float k, Point pointIntersection, int fovea1, int fovea2){
  FoveatedHessianDetectorParams f1, f2;
  f1 = params[fovea1];
  f2 = params[fovea2];
  // Modificar v para o centro da região
  Point v = Point(f2.foveaModel.fx - f1.foveaModel.fx, f2.foveaModel.fy - f1.foveaModel.fy);
  std::vector<Point> limits;
  Point delta = Point(0, 0);
  Point size = Point(0, 0);
  // Without intersection
  if ( (pointIntersection.x == -1) && (pointIntersection.y == -1) ){
    limits.push_back(delta);
    limits.push_back(size);
    return limits;
  }
  // Horizontal shifting
  if ( (v.x > 0) && (v.y == 0) ){ // Shift out of origin
    //std::cout << "sentido horizontal positivo" << std::endl;
    delta = Point(pointIntersection.x, f2.foveaModel.getDeltay(k));
    size = Point(f2.foveaModel.getSizex(k) - (pointIntersection.x - f2.foveaModel.getDeltax(k)), f2.foveaModel.getSizey(k));
    limits.push_back(delta);
    limits.push_back(size);
  }
  if ( (v.x < 0) && (v.y == 0) ){ // Shift in of origin
    //std::cout << "sentido horizontal negativo" << std::endl;
    delta = Point(f2.foveaModel.getDeltax(k), f2.foveaModel.getDeltay(k));
    size = Point(pointIntersection.x - f2.foveaModel.getDeltax(k), f2.foveaModel.getSizey(k));
    limits.push_back(delta);
    limits.push_back(size);
  }

  // Vertical shifting
  if ( (v.x == 0) && (v.y > 0) ){ // Shift out of origin
    //std::cout << "sentido vertical baixo" << std::endl;
    delta = Point(f2.foveaModel.getDeltax(k), pointIntersection.y);
    size = Point(f2.foveaModel.getSizex(k), f2.foveaModel.getSizey(k) - (pointIntersection.y - f2.foveaModel.getDeltay(k)));
    limits.push_back(delta);
    limits.push_back(size);
  }
  if ( (v.x == 0) && (v.y < 0) ){ // Shift in of origin
    //std::cout << "sentido vertical cima" << std::endl;
    delta = Point(f2.foveaModel.getDeltax(k), f2.foveaModel.getDeltay(k));
    size = Point(f2.foveaModel.getSizex(k), pointIntersection.y - f2.foveaModel.getDeltay(k));
    limits.push_back(delta);
    limits.push_back(size);
  }
  
  // Diagonal shifting
  if ( (v.x > 0) && (v.y < 0) ){ // Deslocamento para o sentido nordeste
    //std::cout << "sentido nordeste" << std::endl;
    // Figura 1
    delta = Point(f2.foveaModel.getDeltax(k), f2.foveaModel.getDeltay(k));
    size = Point(f2.foveaModel.getSizex(k), pointIntersection.y - f2.foveaModel.getDeltay(k));
    limits.push_back(delta);
    limits.push_back(size);
    // Figura 2
    delta = Point(pointIntersection.x, pointIntersection.y);
    size = Point(f2.foveaModel.getSizex(k) - (pointIntersection.x - f2.foveaModel.getDeltax(k)), 
		 f2.foveaModel.getSizey(k) - (pointIntersection.y - f2.foveaModel.getDeltay(k)));
    limits.push_back(delta);
    limits.push_back(size);
  }
  if ( (v.x > 0) && (v.y > 0) ){ // Deslocamento para o sentido sudeste
    //std::cout << "sentido sudeste" << std::endl;
    // Figura 1
    delta = Point(pointIntersection.x, f2.foveaModel.getDeltay(k));
    size = Point(f2.foveaModel.getSizex(k) - (pointIntersection.x - f2.foveaModel.getDeltax(k)), 
		 pointIntersection.y - f2.foveaModel.getDeltay(k));
    limits.push_back(delta);
    limits.push_back(size);
    // Figura 2
    delta = Point(f2.foveaModel.getDeltax(k), pointIntersection.y);
    size = Point(f2.foveaModel.getSizex(k), f2.foveaModel.getSizey(k) - (pointIntersection.y - f2.foveaModel.getDeltay(k)));
    limits.push_back(delta);
    limits.push_back(size);
  }
  if ( (v.x < 0) && (v.y < 0) ){ // Deslocamento para o sentido norte
    //std::cout << "sentido norte" << std::endl;
    // Figura 1
    delta = Point(f2.foveaModel.getDeltax(k), f2.foveaModel.getDeltay(k));
    size = Point(f2.foveaModel.getSizex(k), pointIntersection.y - f2.foveaModel.getDeltay(k));
    limits.push_back(delta);
    limits.push_back(size);
    // Figura 2
    delta = Point(f2.foveaModel.getDeltax(k), pointIntersection.y);
    size = Point((pointIntersection.x - f2.foveaModel.getDeltax(k)), 
		 f2.foveaModel.getSizey(k) - (pointIntersection.y - f2.foveaModel.getDeltay(k)));
    limits.push_back(delta);
    limits.push_back(size);
  }
  if ( (v.x < 0) && (v.y > 0) ){ // Deslocamento para o sentido centro-oeste
    //std::cout << "sentido centro-oeste" << std::endl;
    // Figura 1
    delta = Point(f2.foveaModel.getDeltax(k), f2.foveaModel.getDeltay(k));
    size = Point(pointIntersection.x - f2.foveaModel.getDeltax(k), pointIntersection.y - f2.foveaModel.getDeltay(k));
    limits.push_back(delta);
    limits.push_back(size);
    // Figura 2
    delta = Point(f2.foveaModel.getDeltax(k), pointIntersection.y);
    size = Point(f2.foveaModel.getSizex(k), f2.foveaModel.getSizey(k) - (pointIntersection.y - f2.foveaModel.getDeltay(k)));
    limits.push_back(delta);
    limits.push_back(size);
  }

  // Without shifting
  if ( (v.x == 0) && (v.y == 0) ){
    limits.push_back(delta);
    limits.push_back(size);
  }
  
  return limits;
}
  
/**
 * \fn bool verifyRegion(int position, std::vector<Point> region, Point v);
 *
 * \brief Verify if vector v is in delimited region
 *
 * \param position - Position in region vector
 *        region - Vectors os delta and size of region
 *        v - Vector analised.
 *
 * \return True if v is inside region and false otherwise
 */
bool 
MultiFoveation::verifyRegion(int position, std::vector<Point> region, Point v){
  Point startPoint = region[position];
  Point finishPoint = region[position+1];
  // Update finish point
  finishPoint.x += startPoint.x;
  finishPoint.y += startPoint.y;
  // Vector v - startPoint
  Point v_startPoint = Point(v.x - startPoint.x, v.y - startPoint.y);
  // Vector v - finishPoint
  Point v_finishPoint = Point(v.x - finishPoint.x, v.y - finishPoint.y);
  if ( ( v_startPoint.x > 0 ) && ( v_startPoint.y > 0 ) &&
       ( v_finishPoint.x < 0 ) && ( v_finishPoint.y < 0 ) ) return true;
  else return false;
}

/**
 * \fn std::vector<Point> updateLimit(std::vector<Point> limits);
 *
 * \brief Function for calculate the limit of processing
 *
 * \param limits - Vector of limits for update
 *
 * \return Vector of limits for region to be processed
 */
std::vector<Point>
MultiFoveation::updateLimit(std::vector<Point> limits){
  std::vector<Point> region;
  std::vector<Point> regionSearch;
  bool atualized = false;
  std::vector<int> intercepted(limits.size(), -1);
  for (unsigned int i = 0; i < limits.size(); i+=2){
    if ( i == 0 ) intercepted[i] = 0;
    // Region Analised
    regionSearch.push_back(limits[i]);
    regionSearch.push_back(limits[i+1]);
    for (unsigned int j = i+2; j < limits.size(); j+=2){ 
      if ( intercepted[j] == -1 ){ // Not intercepted
	// Build vertices ( clockwise direction )
	Point v1 = limits[j];
	Point v2 = Point((limits[j]).x+(limits[j+1]).x, (limits[j]).y);
	Point v3 = Point((limits[j]).x+(limits[j+1]).x, (limits[j]).y+(limits[j+1]).y);
	Point v4 = Point((limits[j]).x, (limits[j]).y+(limits[j+1]).y);
	// Regions update
	if ( verifyRegion(i, regionSearch, v1) ){
	  intercepted[j] = i;
	  regionSearch[i+1].x -= (v1.x - regionSearch[i].x);
	  regionSearch[i+1].y -= (v1.y - regionSearch[i].y);
	  regionSearch[i] = v1;
	  atualized = true;
        }
	if ( verifyRegion(i, regionSearch, v2) ){
	  intercepted[j] = i;
	  regionSearch[i+1].x = (v2.x - regionSearch[i].x);
	  regionSearch[i+1].y -= (v2.y - regionSearch[i].y);
	  regionSearch[i].y = v2.y;
	  atualized = true;
        }
        if ( verifyRegion(i, regionSearch, v3) ){
	  intercepted[j] = i;
	  regionSearch[i+1].x = (v3.x - regionSearch[i].x);
	  regionSearch[i+1].y = (v3.y - regionSearch[i].y);
	  atualized = true;
        }
        if ( verifyRegion(i, regionSearch, v4) ){
	  intercepted[j] = i;
	  regionSearch[i+1].x -= (v4.x - regionSearch[i].x);
	  regionSearch[i+1].y = (v4.y - regionSearch[i].y);
	  regionSearch[i].x = v4.x;
	  atualized = true;
	}
      }
    }
    int control = 0;
    // Region without intersection to set in final vector
    for (unsigned int k = i+2; k < limits.size(); k+=2){
      if ( ( intercepted[k] == -1 ) && ( k < limits.size()-(control+2) ) ){
        if ( intercepted[limits.size()-(control+2)] != -1 ){
	  std::swap(limits[k], limits[limits.size()-(control+2)]);
	  std::swap(limits[k+1], limits[limits.size()-(control+1)]);
	  std::swap(intercepted[k], intercepted[limits.size()-(control+2)]);
	  std::swap(intercepted[k+1], intercepted[limits.size()-(control+1)]);
	}
	else{ 
	  control += 2;
	  k += 2;
	}
      }
    }

    if ( atualized ){
      region.push_back(regionSearch[i]);
      region.push_back(regionSearch[i+1]);
      atualized = false;
    }
    
    if ( (i == 0) && (!atualized) ) intercepted[i] = -1;
    
  }
  
  // Add region without intersection
  for (unsigned int i = 0; i < limits.size(); i+=2){
    if ( ( intercepted[i] == -1 ) || ( limits.size() == 2 ) ){
      region.push_back(regionSearch[i]);
      region.push_back(regionSearch[i+1]);
    }
  }
  
  return region;
}
