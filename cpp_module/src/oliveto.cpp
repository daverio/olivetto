#ifndef OLIVETTO_CPP
#define OLIVETTO_CPP

#include <iostream>
#include "oliveto.hpp"
#include <float.h>
#include <cmath>



void Oliveto::set_border(boost::python::numpy::ndarray const & array)
{
  //check is the array is 2d!
  if(array.get_nd() != 2){
    throw std::invalid_argument("Oliveto::set_border argument must a 2d numpy array");
  }

  //allocate the border, fill it, and compute offsets.
  borderSize = array.shape(0);;

  double * dt = reinterpret_cast<double*>(array.get_data());
  border = new Point[borderSize];
  for(size_t i = 0;i<borderSize;i++){
    border[i].x =  dt[i*2];
    border[i].y = dt[i*2+1];
  }
  Point lowestBorderPoint;
  lowestBorderPoint.x = DBL_MAX;
  lowestBorderPoint.y = DBL_MAX;
  for(size_t i = 0;i<borderSize;i++){
    if(border[i].x < lowestBorderPoint.x)lowestBorderPoint.x = border[i].x;
    if(border[i].y < lowestBorderPoint.y)lowestBorderPoint.y = border[i].y;
  }
  PointOffset = lowestBorderPoint;
  PointOffset.x = floor(PointOffset.x);
  PointOffset.y = floor(PointOffset.y);
  std::cout<<PointOffset<<std::endl;
  //Rescale the border;
  for(size_t i = 0;i<borderSize;i++){
    border[i].x -= PointOffset.x;
    border[i].y -= PointOffset.y;
  }
}

void Oliveto::add_tree_to_border(double average_dist, double spread)
{ 
  //compute the length of the border;
  double borderLength = 0;
  for(size_t i=0;i<borderSize-1;i++)
  {
    borderLength += dist(border[i+1],border[i]);
  }
  num_borderTrees = floor(borderLength/average_dist);
  std::cout<<borderLength<<" , "<< num_borderTrees << std::endl;

  double cad = borderLength / num_borderTrees;


  borderTrees = new Point[num_borderTrees];

  borderTrees[0] = border[0];
  size_t bindex = 0;
  Point ref = borderTrees[0];
  for(size_t i=1;i<num_borderTrees;i++)
  {
    double travel = cad;
    
    while(travel>0)
    {
      Point dir = border[bindex] - border[bindex+1];
      dir /= dir.norm();
      double dist_to_next_bp = dist(ref,border[bindex+1]);
      if(dist_to_next_bp > travel)
      {
        borderTrees[i] = ref + (dir*travel);
        ref = borderTrees[i];
        travel = 0;
      }
      else
      {
        bindex++;
        ref = border[bindex];
        travel -= dist_to_next_bp;
      }


    }
  }
}

boost::python::numpy::ndarray get_border_trees();
{
  
}

#endif
