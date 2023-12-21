#include <iostream>
#include "oliveto.hpp"
#include <float.h>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "polygone_tests.hpp"
#include <chrono>


void Oliveto::set_borders(boost::python::numpy::ndarray const & hedge,
                          boost::python::numpy::ndarray const & limit)
{
    if(hedge.get_nd() != 2){
        throw std::invalid_argument("Oliveto::set_border argument hedge must be a 2d numpy array");
    }
    if(limit.get_nd() != 2){
        throw std::invalid_argument("Oliveto::set_border argument limit must be a 2d numpy array");
    }
    set_hedge(hedge);
    set_limit(limit);
}

void Oliveto::set_hedge(boost::python::numpy::ndarray const & array)
{
    hedgeBorderSize_ = array.shape(0);; 
    double * dt = reinterpret_cast<double*>(array.get_data());
    hedgeBorder_ = new Point[hedgeBorderSize_];
    for(size_t i = 0;i<hedgeBorderSize_;i++){
        hedgeBorder_[i].x = dt[i*2];
        hedgeBorder_[i].y = dt[i*2+1];
    }   
    offset_.x = DBL_MAX;
    offset_.y = DBL_MAX;
    for(size_t i = 0;i<hedgeBorderSize_;i++){
        if(hedgeBorder_[i].x < offset_.x)offset_.x = hedgeBorder_[i].x;
        if(hedgeBorder_[i].y < offset_.y)offset_.y = hedgeBorder_[i].y;
    }
    offset_.x = floor(offset_.x);
    offset_.y = floor(offset_.y);   
    for(size_t i = 0;i<hedgeBorderSize_;i++){
        hedgeBorder_[i].x -= offset_.x;
        hedgeBorder_[i].y -= offset_.y;
    }
};
void Oliveto::set_limit(boost::python::numpy::ndarray const & array)
{
    treeLimitSize_ = array.shape(0);; 
    double * dt = reinterpret_cast<double*>(array.get_data());
    treeLimit_ = new Point[treeLimitSize_];
    for(size_t i = 0;i<treeLimitSize_;i++){
        treeLimit_[i].x = dt[i*2] - offset_.x;
        treeLimit_[i].y = dt[i*2+1] - offset_.y;
    }   
};

void Oliveto::add_tree(size_t total_tree, bool flagHedge, double average_dist /* =7 */, double spread /*= 1*/)
{
  if(flagHedge) {
    add_tree_to_hedge(average_dist,spread);
  }
  add_interior_tree(total_tree);
}

void Oliveto::add_tree_to_hedge(double average_dist, double spread)
{ 
  //compute the length of the border;
  double borderLength = 0;
  for(size_t i=0;i<hedgeBorderSize_-1;i++)
  {
    borderLength += dist(hedgeBorder_[i+1],hedgeBorder_[i]);
  }
  numHedgeTrees_ = floor(borderLength/average_dist);
  //std::cout<<borderLength<<" , "<< num_borderTrees << std::endl;

  double cad = (borderLength / numHedgeTrees_)*1.001;


  hedgeTrees_ = new Point[numHedgeTrees_];

  hedgeTrees_[0] = hedgeBorder_[0];
  size_t bindex = 0;
  Point ref = hedgeTrees_[0];

  for(size_t i=1;i<numHedgeTrees_;i++)
  {
    double travel = cad;
    
    while(travel>0)
    {
      Point dir = hedgeBorder_[bindex+1] - hedgeBorder_[bindex];
      dir /= dir.norm();
      double dist_to_next_bp = dist(ref,hedgeBorder_[bindex+1]);

      if(dist_to_next_bp > travel)
      {
        hedgeTrees_[i] = ref + (dir*travel);
        ref = hedgeTrees_[i];
        travel = 0;
      }
      else
      {
        bindex++;
        ref = hedgeBorder_[bindex];
        travel -= dist_to_next_bp;
      }
    }
  }
}


void Oliveto::add_interior_tree(size_t numTrees)
{
  numTrees_ = numTrees;
  trees_ = new Point[numTrees_];
  treesVel_ = new Point[numTrees_];
  force_ = new Point[numTrees_];

  //find recangle
  Point offset,span;
  offset.x = DBL_MAX;
  offset.y = DBL_MAX;
  span.x = DBL_MIN;
  span.y = DBL_MIN;
  for(size_t i=0;i<treeLimitSize_;i++)
  {
    if(treeLimit_[i].x<offset.x)offset.x = treeLimit_[i].x;
    if(treeLimit_[i].y<offset.y)offset.y = treeLimit_[i].y;
    if(treeLimit_[i].x>span.x)span.x = treeLimit_[i].x;
    if(treeLimit_[i].y>span.y)span.y = treeLimit_[i].y;
  }
  span = span - offset;
  //std::cout<<offset<<"  "<<span<<std::endl;
  auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(std::chrono::system_clock::now());
  auto value = now_ms.time_since_epoch();
  unsigned long seed = now_ms.time_since_epoch().count();


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
  while(n<numTrees_)
  {
    newTree.x = gsl_rng_uniform(rx)*span.x + offset.x;
    newTree.y = gsl_rng_uniform(ry)*span.y + offset.y;

    if(InsidePolygon(treeLimit_,treeLimitSize_-1,newTree))
    {
      trees_[n] = newTree;
      treesVel_[n].x = 0;
      treesVel_[n].y = 0;
      n++;
    }

  }
  gsl_rng_free(rx);
  gsl_rng_free(ry);

}

boost::python::numpy::ndarray Oliveto::get_border_trees()
{
  Py_intptr_t shape[2] = {(long)numHedgeTrees_,2};
  boost::python::numpy::ndarray res = boost::python::numpy::zeros(2, shape, boost::python::numpy::dtype::get_builtin<double>());
  double * arr = reinterpret_cast<double*>(res.get_data());
  for(size_t i=0;i<numHedgeTrees_;i++)
  {
    arr[i*2] = hedgeTrees_[i].x;
    arr[i*2+1] = hedgeTrees_[i].y;
  }  
  return res;
}

boost::python::numpy::ndarray Oliveto::get_inner_trees()
{
  Py_intptr_t shape[2] = {(long)numTrees_,2};
  boost::python::numpy::ndarray res = boost::python::numpy::zeros(2, shape, boost::python::numpy::dtype::get_builtin<double>());
  double * arr = reinterpret_cast<double*>(res.get_data());
  for(size_t i=0;i<numTrees_;i++)
  {
    arr[i*2] = trees_[i].x;
    arr[i*2+1] = trees_[i].y;
  }  
  return res;
}

void Oliveto::make_step(double delta_t, double coupling,double range, double viscosity, double damping)
{
  if(viscosity<0 || viscosity>1 || damping<0 || damping>1)
  {
    std::cout<<"viscosity and damaping needs to be in [0,1]"<<std::endl;
  }
  dt_ = delta_t;
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
  for(int i = 0; i< numTrees_; i++)
  {
    force_[i].x = 0;
    force_[i].y = 0;
    for(int j=0;j<numHedgeTrees_;j++)
    {
      force_[i] += compute_pair_force(hedgeTrees_[j],trees_[i],coupling,range);
    }
  }

  for(int i = 0; i< numTrees_-1; i++)
  {
    for(int j = i+1; j< numTrees_; j++)
    {
      Point f = compute_pair_force(trees_[j],trees_[i],coupling,range);
      force_[i] += f;
      force_[j] -= f;
    }
  }
}
void Oliveto::update_vel(double damping)
{
  for(size_t i=0; i<numTrees_; i++)
  {
    treesVel_[i] += force_[i] *dt_;
    treesVel_[i] *= 1.0-damping;
  }
}


void Oliveto::displace_trees(double viscosity)
{
  for(size_t i=0; i<numTrees_; i++)
  {
    Point disp = treesVel_[i] * dt_ * (1.0-viscosity);
    move_tree(trees_[i],treesVel_[i],disp,treeLimit_,treeLimitSize_);
  }
}

void Oliveto::move_tree(Point & Tree, Point & vel, Point & displacement, Point const * const border, size_t const borderSize, int border_nocol)
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
      col_flag = check_colision(Tree, new_pos, treeLimit_[i], treeLimit_[i+1],col);
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
    bvect =  treeLimit_[bnocol+1] -  treeLimit_[bnocol];
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
