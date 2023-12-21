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
    Point * force_;
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

    void add_tree(size_t total_tree, bool flagHedge, double average_dist =7, double spread = 1);
    void make_step(double delta_t, double coupling, double range, double viscosity, double damping);

    boost::python::numpy::ndarray get_border_trees();
    boost::python::numpy::ndarray get_inner_trees();
    boost::python::numpy::ndarray get_trees_dist();
    int get_number_border_tree(){return numHedgeTrees_;}
    int get_number_inner_tree(){return numTrees_;}
    boost::python::tuple get_offset(){return boost::python::make_tuple(offset_.x,offset_.y);}

private:
    void set_hedge(boost::python::numpy::ndarray const & array);
    void set_limit(boost::python::numpy::ndarray const & array);
    void add_tree_to_hedge(double average_dist, double spread);
    void add_interior_tree(size_t numTrees);
    
    void update_vel(double damping);
    void displace_trees(double viscosity);
    void move_tree(Point & Tree, Point & vel, Point & displacement, Point const * const  border, size_t const borderSize, int border_nocol = -1);
    void compute_force(double coupling, double range);
    Point compute_pair_force(const Point & t1, const Point & t2, double coupling, double range);
};