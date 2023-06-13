#ifndef OLIVETTO_CPP
#define OLIVETTO_CPP

#include <iostream>
#include "oliveto.hpp"




void Oliveto::set_border(boost::python::numpy::ndarray const & array)
{
  std::cout<<"set_border called"<< std::endl;
  std::cout<<array.get_nd()<<","<<array.shape(1)<< std::endl;
  //check is the array is 2d!
  if(array.get_nd() != 2){
    throw std::invalid_argument("Oliveto::set_border argument must a 2d numpy array");
  }

}

#endif
