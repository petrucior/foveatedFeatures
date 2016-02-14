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
 * Copyright (C) 2014, Petrúcio Ricardo <petrucior@gmail.com>
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
  /**
   * \fn void MultiFoveation(int foveae)
   *
   * \brief Constructor multi foveae.
   *
   * \param foveae - Number of foveae
   */
  MultiFoveation(int foveas);
  
  /**
   * \fn Point intersection(int fovea1, int fovea2);
   *
   * \brief Function for calculate the intersection between foveae.
   *
   * \param fovea1 - Fovea before processed.
   *        fovea2 - Fovea to be processed.
   *
   * \return Point of intersection between foveae.
   */
  Point intersection(int fovea1, int fovea2);
  
};

#endif

/**
 * \fn void MultiFoveation(int foveae)
 *
 * \brief Constructor multi foveae.
 *
 * \param foveae - Number of foveae
 */
MultiFoveation::MultiFoveation(int foveas){
  printf("Estou iniciando o foveamento");
}

/**
 * \fn Point intersection(int fovea1, int fovea2);
 *
 * \brief Function for calculate the intersection between foveae.
 *
 * \param fovea1 - Fovea before processed.
 *        fovea2 - Fovea to be processed.
 *
 * \return Point of intersection between foveae.
 */
Point 
MultiFoveation::intersection(int fovea1, int fovea2){
  Point p = Point(1, 2);
  return p;
}

