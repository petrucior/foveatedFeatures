
/*
  Copyright (C) 2014, Rafael Beserra <rafaelufrn@gmail.com>
  If you use this software for academic purposes, consider citing the related paper: Rafael Beserra Gomes, Bruno Motta de Carvalho, Luiz Marcos Garcia Gonçalves, Visual attention guided features selection with foveated images, Neurocomputing, Volume 120, 23 November 2013, Pages 34-44, ISSN 0925-2312, http://dx.doi.org/10.1016/j.neucom.2012.10.033.

  This file is part of foveatedFeatures software.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef LINEAR_FOVEATION
#define LINEAR_FOVEATION

struct LinearFoveation {

  LinearFoveation() {
    fx = fy = growthfactor = 0;
    wx = wy = ux = uy = 0;
    m = 0;
    flag = false;
    id = 0;
    deltax.clear();
    deltay.clear();
    sizex.clear();
    sizey.clear();
  }
  
  inline int getDeltax(int k) {
    return (k*(ux - wx + 2*fx))/(2*m);
  }

  inline int getDeltay(int k) {
    return (k*(uy - wy + 2*fy))/(2*m);
  }

  inline int getSizex(int k) {
    return (k*wx - k*ux + m*ux)/m;
  }
  
  inline int getSizey(int k) {
    return (k*wy - k*uy + m*uy)/m;
  }

  // Ativation flag of a multifoveation system
  inline bool getFlag(){
    return flag;
  }

  inline int getId(){
    return id;
  }

  inline int getQuantityOfIntersections(int k){
    return (int)deltax[k].size();
  }

  inline int getDeltax(int indice, int jump) {
    return deltax[indice][jump];
  }

  inline int getDeltay(int indice, int jump) {
    return deltay[indice][jump];
  }

  inline int getSizex(int indice, int jump) {
    return sizex[indice][jump];
  }
  
  inline int getSizey(int indice, int jump) {
    return sizey[indice][jump];
  }

  bool positionCalculated(int linha, int coluna, int k){
    //std::cout << "Deltax: " << getDeltax(k, 0) << " Deltay: " << getDeltay(k, 0) << " --- linha: " << linha << ", coluna: " << coluna << std::endl;
    for (int i = 0; i < getQuantityOfIntersections(k); i++){
      // When the region haven't intersections
      if ( ( k != 0 ) &&
	   ( getDeltax(k, i) == 0 ) &&
	   ( getDeltay(k, i) == 0 ) &&
	   ( getSizex(k, i) == 0 ) &&
	   ( getSizey(k, i) == 0 )  ) return false;
      
      // Position delimited
      if ( ( getDeltax(k, i) < linha ) &&
	   ( getDeltay(k, i) < coluna ) &&
	   ( getDeltax(k, i) + getSizex(k, i) > linha ) &&
	   ( getDeltay(k, i) + getSizey(k, i) > coluna ) ) return false;

    }
    return true;
  }
  
  void setMultiFoveation(int indice, std::vector<int> _delta, std::vector<int> _size, int _id){
    deltax[indice] = std::vector<int>((int)(_delta.size()/2));
    deltay[indice] = std::vector<int>((int)(_delta.size()/2));
    sizex[indice] = std::vector<int>((int)(_size.size()/2));
    sizey[indice] = std::vector<int>((int)(_size.size()/2));
    int count = 0;
    for (unsigned int i = 0; i < _delta.size(); i+=2){
      deltax[indice][count] = _delta[i];
      deltay[indice][count] = _delta[i+1];
      sizex[indice][count] = _size[i];
      sizey[indice][count] = _size[i+1];
      count++;
    }
    flag = true;
    id = _id;
    /*std::cout << "Fóvea " << _id << ", indice: " << indice << std::endl;
    std::cout << "Delta" << std::endl;
    std::cout << "(" << deltax[indice][0] << ", " << deltay[indice][0] << ")" << std::endl;
    std::cout << "Size" << std::endl;
    std::cout << "(" << sizex[indice][0] << ", " << sizey[indice][0] << ")" << std::endl;*/
  }
  
  //fix the fovea position: if fovea is outsite image domain, snap it to the closest valid position independently for each coordinate
  inline void fixFovea() {
    fx = MIN((ux - wx)/2 - growthfactor, fx);
    fx = MAX((wx - ux)/2 + growthfactor, fx);
    fy = MIN((uy - wy)/2 - growthfactor, fy);
    fy = MAX((wy - uy)/2 + growthfactor, fy);
  }

  void setFovea(int imgX, int imgY) {
    fx = imgX - ux/2;
    fy = imgY - uy/2;
    fixFovea();
  }

  void check() {
    assert(wx > 0 && wx < ux);
    assert(wy > 0 && wy < uy);
    assert(ux > 0 && uy > 0);
    assert(m >= 1);
    assert(beta.size() == eta.size());
    assert(eta.size() == level.size());
    for(unsigned int i = 0; i < beta.size(); i++) {
      assert(beta[i] == 1 || beta[i] == 0);
      assert(eta[i] >= 1);
      assert(level[i] >= 0 && level[i] <= m);
    }
    assert(growthfactor >= 0);
  }

  void init(){
     // Iniciando os vetores
    deltax = std::vector<std::vector<int> >(m+1, std::vector<int>(0));
    deltay = std::vector<std::vector<int> >(m+1, std::vector<int>(0));
    sizex = std::vector<std::vector<int> >(m+1, std::vector<int>(0));
    sizey = std::vector<std::vector<int> >(m+1, std::vector<int>(0));
  }

  int wx, wy; //smallest level size
  int ux, uy; //image size
  int m; //numberOfLevels - 1
  int fx, fy; //fovea position
  int growthfactor;
  std::vector<int> beta;
  std::vector<int> eta;
  std::vector<int> level;
  bool flag;
  int id;
  std::vector<std::vector<int> > deltax, deltay, sizex, sizey;
};


#endif

