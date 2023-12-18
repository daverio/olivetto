#pragma once

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <iostream>
#include "point.hpp"


class Oliveto{
private:

    Point * hedgeBorder_;
    size_t hedgeBorderSize_;
    Point * treeLimit_;
    size_t treeLimitSize_;
    Point * hedgeTrees_;
    size_t numHedgeTrees_;
    Point * trees_;
    size_t numTrees_;
    Point * treesVel_;
    Point * force;
    double dt_;
    Point offset_;

public:

    Oliveto(){
        hedgeBorder_ = nullptr;
        treeLimit_ = nullptr;
        hedgeTrees_ = nullptr;
        trees_ = nullptr;
        treesVel_ = nullptr;
        force_ = nullptr;
        hedgeBorderSize_ = 0;
        treeLimitSize_ = 0;
        numHedgeTrees_ = 0;
        numTrees_ = 0;
    };

    ~Oliveto(){
        if(hedgeBorder_ != nullptr) delete hedgeBorder_;
        if(treeLimit_ != nullptr) delete treeLimit_;
        if(hedgeTrees_ != nullptr) delete hedgeTrees_;
        if(trees_ != nullptr) delete trees_;
        if(treesVel_ != nullptr) delete treesVel_;
        if(force_ != nullptr) delete force_;
    };
    
    void set_borders(boost::python::numpy::ndarray const & hedge,
                     boost::python::numpy::ndarray const & limit);

private:
    void set_hedge(boost::python::numpy::ndarray const & array);
    void set_limit(boost::python::numpy::ndarray const & array);
};