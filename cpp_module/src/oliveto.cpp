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
        treeLimit_[i].x = dt[i*2] - offset.x;
        treeLimit_[i].y = dt[i*2+1] - offset.y;
    }   
};
