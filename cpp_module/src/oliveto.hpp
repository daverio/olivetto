#ifndef OLIVETO_HPP
#define OLIVETO_HPP

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>


class Oliveto
{
private:

public:
  void set_border(boost::python::numpy::ndarray const & array);
  /*tupple get_border();

  add_tree_to_boarder(double average_dist, double spread);

  tuple get_border_tree();
  */

};

#endif
