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
  size_t num_trees;
public:

  Oliveto(){
    border = nullptr;
    borderTrees = nullptr;
    trees = nullptr;
  }

  ~Oliveto()
  {
    if(border != nullptr) delete[] border;
    if(borderTrees != nullptr) delete[] borderTrees;
    if(trees != nullptr) delete[] trees;
  }

  void set_border(boost::python::numpy::ndarray const & array);
  void add_tree_to_border(double average_dist, double spread);
  boost::python::numpy::ndarray get_border_trees();
  boost::python::numpy::ndarray get_iner_trees();

  void add_interior_tree(size_t total_tree);

};

#endif
