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
   * \fn void updateParams(int fovea)
   *
   * \brief Update the params of fovea.
   *
   * \param fovea - The number of fovea
   */
  void updateParams(int fovea);
  
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
   * \fn Point intersection(int k, int m, Size R, int fovea1, int fovea2);
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
  Point intersection(int k, int m, Size R, int fovea1, int fovea2);
  
  /**
   * \fn std::vector<Point> limitProcessing(Point pointIntersection, int fovea1, int fovea2);
   *
   * \brief Function for calculate the limits of processing.
   *
   * \param pointIntersection - The point of intersection to find with intersection function
   *        fovea1 - Fovea before processed
   *        fovea2 - Fovea to be processed
   *
   * \return Vector of limits of processing
   */
  std::vector<Point> limitProcessing(Point pointIntersection, int fovea1, int fovea2);
  
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
  std::vector<int> pontox;
  std::vector<int> pontoy;
  std::vector<Point> limits;
  std::vector<int> limit;
  params.clear();
  for (int i = 0; i < foveas; i++){
    FoveatedHessianDetectorParams p(image.cols, image.rows, ymlFile[i]);
    params.push_back(p);
    int m = params[i].foveaModel.m;
    if ( i != 0 ){
      for (int k = 0; k < m + 1; k++){
	Point pontos = intersection(k, m, image.size(), i-1, i);
	// Clear points
	pontox.clear();
	pontoy.clear();
	pontox.push_back(pontos.x);
	pontoy.push_back(pontos.y);
	// Clear limits
	limits.clear();
	limit.clear();
	limits = limitProcessing(pontos, i-1, i);
	for (unsigned int k = 0; k < limits.size(); k++){
	  pontos = limits[k];
	  limit.push_back(pontos.x);
	  limit.push_back(pontos.y);
	}
      }
      params[i].foveaModel.setIntersection(pontox, pontoy, limit);
    }
  }
  // DEBUG
  /*for (int i = 0; i < limits.size(); i++){
    Point p = limits[i];
    std::cout << p.x << " :::: " << p.y << std::endl;
  }*/
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
  if ( fovea != 0 ){
    // Update of fovea params
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
}

/**
 * \fn Point intersection(int k, int m, Size R, int fovea1, int fovea2);
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
MultiFoveation::intersection(int k, int m, Size R, int fovea1, int fovea2){
  FoveatedHessianDetectorParams f1, f2;
  f1 = params[fovea1];
  f2 = params[fovea2];
  Point p = Point(0, 0);
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
  
  // DEBUG
  /*std::cout << "f1.x = " << f1.foveaModel.fx << ", f2.x = " << f2.foveaModel.fx << std::endl;
  std::cout << "p1 = " << p1 << ", p2 = " << p2 << std::endl;
  std::cout << "wmax = " << wmax << ", wmin = " << wmin << std::endl;*/
  // Component x axis
  if ( min(f1.foveaModel.fx+(R.width/2), f2.foveaModel.fx+(R.width/2)) == f1.foveaModel.fx+(R.width/2) ){ // Fóvea 1 mais próxima da origem
    p.x = ( k * p1 )/m;
  }
  else{ // Fóvea 2 mais próxima da origem
    p.x = ( (R.width * m) - (R.width * k) + (p2 * k) )/m;
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
  // DEBUG
  /*std::cout << "f1.y = " << f1.foveaModel.fy << ", f2.y = " << f2.foveaModel.fy << std::endl;
  std::cout << "p1 = " << p1 << ", p2 = " << p2 << std::endl;
  std::cout << "wmax = " << wmax << ", wmin = " << wmin << std::endl;*/
  // Component y axis
  if ( min(f1.foveaModel.fy+(R.height/2), f2.foveaModel.fy+(R.height/2)) == f1.foveaModel.fy+(R.height/2) ){ // Fóvea 1 mais próxima da origem
    p.y = ( k * p1 )/m;
  }
  else{ // Fóvea 2 mais próxima da origem
    p.y = ( (R.height * m) - (R.height * k) + (p2 * k) )/m;
  }
  
  std::cout << p.x << " | " << p.y << std::endl;
  
  return p;
}

/**
 * \fn std::vector<Point> limitProcessing(Point pointIntersection, int fovea1, int fovea2);
 *
 * \brief Function for calculate the limits of processing.
 *
 * \param pointIntersection - The point of intersection to find with intersection function
 *        fovea1 - Fovea before processed
 *        fovea2 - Fovea to be processed
 *
 * \return Vector of limits of processing
 */
std::vector<Point> 
MultiFoveation::limitProcessing(Point pointIntersection, int fovea1, int fovea2){
  FoveatedHessianDetectorParams f1, f2;
  f1 = params[fovea1];
  f2 = params[fovea2];
  Point v = Point(f2.foveaModel.fx - f1.foveaModel.fx, f2.foveaModel.fy - f1.foveaModel.fy);
  std::vector<Point> limits;
  Point start = Point(0, 0);
  Point finish = Point(0, 0);
  Point flag = Point(-1, -1);
  // Horizontal shifting
  if ( (v.x > 0) && (v.y == 0) ){ // Shift out of origin
    start.x = pointIntersection.x;
  }
  if ( (v.x < 0) && (v.y == 0) ){ // Shift in of origin
    finish.x = pointIntersection.x;
  }

  // Vertical shifting
  if ( (v.x == 0) && (v.y > 0) ){ // Shift out of origin
    finish.y = pointIntersection.y;
  }
  if ( (v.x == 0) && (v.y < 0) ){ // Shift in of origin
    start.y = pointIntersection.y;
  }
  
  // Diagonal shifting
  if ( (v.x > 0) && (v.y > 0) ){ // Deslocamento para o sentido nordeste
    //std::cout << "sentido nordeste" << std::endl;
    flag.x = 0;
    flag.y = 0;
  }
  if ( (v.x > 0) && (v.y < 0) ){ // Deslocamento para o sentido sudeste
    //std::cout << "sentido sudeste" << std::endl;
    flag.x = 0;
    flag.y = 1;
  }
  if ( (v.x < 0) && (v.y > 0) ){ // Deslocamento para o sentido norte
    //std::cout << "sentido norte" << std::endl;
    flag.x = 1;
    flag.y = 0;
  }
  if ( (v.x < 0) && (v.y < 0) ){ // Deslocamento para o sentido centro-oeste
    //std::cout << "sentido centro-oeste" << std::endl;
    flag.x = 1;
    flag.y = 1;
  }

  // Without shifting
  if ( (v.x == 0) && (v.y == 0) ){
    start = Point(0, 0);
    finish = Point(0, 0);
  }
  
  limits.push_back(start);
  limits.push_back(finish);
  limits.push_back(flag);
  return limits;
}
  
