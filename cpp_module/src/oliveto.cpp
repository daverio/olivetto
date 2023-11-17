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
  offset = PointOffset;
}

void Oliveto::add_tree(size_t total_tree, bool flagBorder, double average_dist /* =7 */, double spread /*= 1*/)
{
  if(flagBorder) {
    add_tree_to_border(average_dist,spread);
  }
  add_interior_tree(total_tree);
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

boost::python::numpy::ndarray Oliveto::get_inner_trees()
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
  trees_vel = new Point[num_trees];
  force = new Point[num_trees];

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
  //std::cout<<offset<<"  "<<span<<std::endl;
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
      trees_vel[n].x = 0;
      trees_vel[n].y = 0;
      n++;
    }

  }
  gsl_rng_free(rx);
  gsl_rng_free(ry);

}
void Oliveto::make_step(double delta_t, double coupling,double range, double viscosity, double damping)
{
  if(viscosity<0 || viscosity>1 || damping<0 || damping>1)
  {
    std::cout<<"viscosity and damaping needs to be in [0,1]"<<std::endl;
  }
  dt = delta_t;
  compute_force(coupling,range);
  update_vel(damping);
  displace_trees(viscosity);
}

Point Oliveto::compute_pair_force(const Point & t1, const Point & t2, double coupling, double range)
{

  Point force = t2 - t1;

  double dist = force.norm();
  //if(dist<0.50)dist = 0.5;

  if(dist < range)
  {
    force *= coupling / (dist*dist*dist);
  }
  else
  {
    force.x=0;
    force.y=0;
  }
  return force;
}

void Oliveto::compute_force(double coupling, double range)
{ 
  //first compute the force from the borders
  for(int i = 0; i< num_trees; i++)
  {
    force[i].x = 0;
    force[i].y = 0;
    for(int j=0;j<num_borderTrees;j++)
    {
      force[i] += compute_pair_force(borderTrees[j],trees[i],coupling,range);
    }
  }

  for(int i = 0; i< num_trees-1; i++)
  {
    for(int j = i+1; j< num_trees; j++)
    {
      Point f = compute_pair_force(trees[j],trees[i],coupling,range);
      force[i] += f;
      force[j] -= f;
    }
  }
}
void Oliveto::update_vel(double damping)
{
  for(size_t i=0; i<num_trees; i++)
  {
    trees_vel[i] *= 1.0-damping;
    trees_vel[i] += force[i] *dt;
  }
}


void Oliveto::displace_trees(double viscosity)
{
  for(size_t i=0; i<num_trees; i++)
  {
    Point disp = trees_vel[i] * dt * (1.0-viscosity);
    move_tree(trees[i],trees_vel[i],disp,border);
  }
}

void Oliveto::move_tree(Point & Tree, Point & vel, Point & displacement, Point const * const border, int border_nocol)
{

  Point col;
  Point new_pos = Tree + displacement;
  int bnocol;
  //check if colisions:

  // if not
  // Tree = Tree + displacement
  // if yes
  // move the tree to the colision point, reset displacement call move_tree again with new displacement
  // reset the direction of velocity
  // to do that we need colision detection to return the point of colision
  bool col_flag = false;
  for(int i = 0; i < borderSize - 1;i++)
  {
    if(i!=border_nocol)
    {
      col_flag = check_colision(Tree, new_pos, border[i], border[i+1],col);
      if(col_flag)
      {
        bnocol = i;
        break;  
      }
    }
  }

  if(col_flag)
  {
    double vel_norm = vel.norm();
    double dist_to_do = displacement.norm() - (col - Tree).norm();
    

    //find new direction (after colision)

    Point normal,bvect;
    bvect =  border[bnocol+1] -  border[bnocol];
    normal.x = -bvect.y;
    normal.y = bvect.x;
    normal /= normal.norm();

    Tree = col;
    vel =  vel - normal*2*(vel.x*normal.x + vel.y*normal.y);
    displacement = vel * dist_to_do / vel_norm;

    if( displacement.norm()>0) move_tree(Tree, vel, displacement, border, bnocol);
  }
  else
  {
    Tree = new_pos;
  }


}

#endif
