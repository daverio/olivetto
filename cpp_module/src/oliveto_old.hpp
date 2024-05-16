#ifndef OLIVETO_HPP
#define OLIVETO_HPP

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <iostream>
#include "point.hpp"


class Oliveto{
private:
  Point PointOffset;
  Point * border;
  size_t borderSize;
  Point * borderTrees;
  size_t num_borderTrees;
  Point * trees;
  Point * trees_vel;
  Point * force;
  size_t num_trees;
  double dt;
  Point offset;
public:

  Oliveto(){
    border = nullptr;
    borderTrees = nullptr;
    trees = nullptr;
    trees_vel = nullptr;
    force = nullptr;
    num_borderTrees = 0;
  }

  ~Oliveto()
  {
    if(border != nullptr) delete[] border;
    if(borderTrees != nullptr) delete[] borderTrees;
    if(trees != nullptr) delete[] trees;
    if(trees_vel != nullptr) delete[] trees_vel;
    if(force != nullptr) delete[] force;
  }

  void set_border(boost::python::numpy::ndarray const & array);
  void add_tree_to_border(double average_dist, double spread);
  boost::python::numpy::ndarray get_border_trees();
  boost::python::numpy::ndarray get_inner_trees();

  void add_interior_tree(size_t total_tree);

  void add_tree(size_t total_tree, bool flagBorder, double average_dis = 7, double spread = 1, bool flagBorderReal = true);


  void update_vel(double damping);
  void displace_trees(double viscosity);
  void move_tree(Point & Tree, Point & vel, Point & displacement, Point const * const  border, int border_nocol = -1);

  void make_step(double delta_t, double coupling, double range, double viscosity, double damping);
  void compute_force(double coupling, double range);
  Point compute_pair_force(const Point & t1, const Point & t2, double coupling, double range);

  int get_number_border_tree(){return num_borderTrees;}
  int get_number_inner_tree(){return num_trees;}
  boost::python::tuple get_offset(){return boost::python::make_tuple(offset.x,offset.y);}


};

#endif
