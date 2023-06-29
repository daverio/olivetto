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
        .def("get_border_trees", &Oliveto::get_border_trees)
        .def("add_interior_tree", &Oliveto::add_interior_tree)
        .def("get_inner_trees", &Oliveto::get_inner_trees)
        .def("make_step",&Oliveto::make_step)
        .def("get_number_border_tree",&Oliveto::get_number_border_tree)
        .def("get_number_inner_tree",&Oliveto::get_number_inner_tree)
        .def("get_offset",&Oliveto::get_offset)
    ;
}
