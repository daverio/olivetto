#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "oliveto.hpp"

BOOST_PYTHON_MODULE(oliveto)
{
  Py_Initialize();
  boost::python::numpy::initialize();

  boost::python::class_<Oliveto>("Oliveto")
        .def("set_border", &Oliveto::set_border)
        .def("add_tree_to_border", &Oliveto::add_tree_to_border)
    ;
}
