#ifndef OLIVETTO_CPP
#define OLIVETTO_CPP

#include <iostream>
#include "oliveto.hpp"
#include <float.h>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "polygone_tests.hpp"


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
  //std::cout<<PointOffset<<std::endl;
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
  //std::cout<<borderLength<<" , "<< num_borderTrees << std::endl;

  double cad = (borderLength / num_borderTrees)*1.001;


  borderTrees = new Point[num_borderTrees];

  borderTrees[0] = border[0];
  size_t bindex = 0;
  Point ref = borderTrees[0];
  for(size_t i=1;i<num_borderTrees;i++)
  {
    double travel = cad;
    
    while(travel>0)
    {
      Point dir = border[bindex+1] - border[bindex];
      dir /= dir.norm();
      double dist_to_next_bp = dist(ref,border[bindex+1]);
      if(dist_to_next_bp > travel)
      {
        borderTrees[i] = ref + (dir*travel);
        ref = borderTrees[i];
        travel = 0;
        //std::cout<<i<<" : "<< borderTrees[i]<< " ; "<<dist(borderTrees[i-1],borderTrees[i])<<std::endl;
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

boost::python::numpy::ndarray Oliveto::get_border_trees()
{
  Py_intptr_t shape[2] = {(long)num_borderTrees,2};
  boost::python::numpy::ndarray res = boost::python::numpy::zeros(2, shape, boost::python::numpy::dtype::get_builtin<double>());
  double * arr = reinterpret_cast<double*>(res.get_data());
  for(size_t i=0;i<num_borderTrees;i++)
  {
    arr[i*2] = borderTrees[i].x;
    arr[i*2+1] = borderTrees[i].y;
  }  

  return res;
}

boost::python::numpy::ndarray Oliveto::get_iner_trees()
{
  Py_intptr_t shape[2] = {(long)num_trees,2};
  boost::python::numpy::ndarray res = boost::python::numpy::zeros(2, shape, boost::python::numpy::dtype::get_builtin<double>());
  double * arr = reinterpret_cast<double*>(res.get_data());
  for(size_t i=0;i<num_trees;i++)
  {
    arr[i*2] = trees[i].x;
    arr[i*2+1] = trees[i].y;
  }  

  return res;
}

void Oliveto::add_interior_tree(size_t total_tree)
{
  num_trees = total_tree - num_borderTrees;
  trees = new Point[num_trees];

  //find recangle
  Point offset,span;
  offset.x = DBL_MAX;
  offset.y = DBL_MAX;
  span.x = DBL_MIN;
  span.y = DBL_MIN;
  for(size_t i=0;i<borderSize;i++)
  {
    if(border[i].x<offset.x)offset.x = border[i].x;
    if(border[i].y<offset.y)offset.y = border[i].y;
    if(border[i].x>span.x)span.x = border[i].x;
    if(border[i].y>span.y)span.y = border[i].y;
  }
  span = span - offset;
  std::cout<<offset<<"  "<<span<<std::endl;
  unsigned long seed = 123;


  //generate the random points:
  const gsl_rng_type * T;
  gsl_rng * rx, * ry;
  gsl_rng_env_setup();

  T = gsl_rng_default; // Generator setup
  rx = gsl_rng_alloc (T);
  ry = gsl_rng_alloc (T);
  gsl_rng_set(rx, seed);
  gsl_rng_set(ry, seed+1);

  size_t n = 0;
  Point newTree;
  while(n<num_trees)
  {
    newTree.x = gsl_rng_uniform(rx)*span.x + offset.x;
    newTree.y = gsl_rng_uniform(ry)*span.y + offset.y;

    if(InsidePolygon(border,borderSize-1,newTree))
    {
      trees[n] = newTree;
      n++;
    }

  }

  
  
  gsl_rng_free(rx);
  gsl_rng_free(ry);


}

#endif
