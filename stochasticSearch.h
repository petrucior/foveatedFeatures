/**
 * \file stochasticSearch.h
 *
 * \brief This file contains the prototype and implementation of stochastic search.
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

#ifndef STOCHASTICSEARCH_H
#define STOCHASTICSEARCH_H

/**
 * \struct StochasticSearch 
 *
 * \brief Struct for stochastic search.
 */
struct StochasticSearch{

  //
  // Variables
  //
  int samples;

  /**
   * \fn StochasticSearch( int _samples )
   *
   * \brief Constructor to search for structure foveated
   *
   * \param _samples - The samples of hipotheses
   */
  StochasticSearch( int samples );
  
  /**
   * \fn Point search( Mat image, Mat template )
   *
   * \brief Method for search the template in image
   *
   * \param image - Image for search
   *        template - Template image to be search
   */
  Point search( Mat image, Mat template );
  
};

#endif


/**
 * \fn StochasticSearch( int _samples );
 *
 * \brief Constructor to search for structure foveated.
 *
 * \param _samples - The samples of hipotheses
 */
StochasticSearch::StochasticSearch( int _samples ){
  samples = _samples;
}

/**
 * \fn Point search( Mat image, Mat template )
 *
 * \brief Method for search the template in image
 *
 * \param image - Image for search
 *        template - Template image to be search
 */
Point 
StochasticSearch::search( Mat image, Mat template ){
  Point stimation;
  Point p = Point(0, 0);
  return p;
}
