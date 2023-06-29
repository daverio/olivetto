import numpy as np
import oliveto
import matplotlib.pyplot as plt
import matplotlib.animation as animation






class Simulation(object):
    def __init__(self,
                 border_file,
                 num_tree, 
                 border_tree_flag,
                 dt,
                 force_coupling,
                 force_range,
                 viscosity,
                 damping,
                 border_dist,
                 border_spread):
        self.data_path = '/home/daverio/Documents/oliveto_data/'
        

        self.border_8m = np.loadtxt(self.data_path + border_file[2], delimiter=',')
        self.border_5m = np.loadtxt(self.data_path + border_file[1], delimiter=',')
        self.border = np.loadtxt(self.data_path + border_file[0], delimiter=',')
        
        self.num_tree = num_tree
        self.border_tree_flag = border_tree_flag
        self.dt = dt
        self.force_coupling = force_coupling
        self.force_range = force_range
        self.viscosity = viscosity
        self.damping = damping

        #setup the simulation:
        self.sim = oliveto.Oliveto()
        self.sim.set_border(self.border_8m)
        self.offset = self.sim.get_offset()
        self.border[:,0] -= self.offset[0]
        self.border[:,1] -= self.offset[1]
        self.border_5m[:,0] -= self.offset[0]
        self.border_5m[:,1] -= self.offset[1]
        self.border_8m[:,0] -= self.offset[0]
        self.border_8m[:,1] -= self.offset[1]
        if(self.border_tree_flag == True):
            self.sim.add_tree_to_border(border_dist,border_spread)
            self.border_trees = self.sim.get_border_trees() 

        self.num_border_tree = self.sim.get_number_border_tree()
        self.sim.add_interior_tree(self.num_tree)
        self.num_inner_tree = self.num_tree - self.num_border_tree
        self.inner_tree = self.sim.get_inner_trees()

    def update_inner_tree(self):
        self.sim.make_step(self.dt,self.force_coupling,self.force_range,self.viscosity,self.damping)
        self.inner_tree = self.sim.get_inner_trees()

    def run_no_anim(self,num_steps,dt):
        self.dt = dt
        for i in range(num_steps):
            self.sim.make_step(self.dt,self.force_coupling,self.force_range,self.viscosity,self.damping)
        self.inner_tree = self.sim.get_inner_trees()

    def plot(self):
        plt.ion()
        plt.figure()
        plt.plot(self.border[:,0],self.border[:,1],color='grey',linestyle='dashed')
        plt.scatter(self.inner_tree[:,0],self.inner_tree[:,1],color='blue')
        if(self.border_tree_flag == True):
            plt.scatter(self.border_trees[:,0],self.border_trees[:,1],color='blue')


    def run_with_anim(self,num_steps,dt):
        self.dt = dt
        self.fig_anim, self.ax_anim = plt.subplots()

        plt.plot(self.border[:,0],self.border[:,1],color='black')
        plt.plot(self.border_5m[:,0],self.border_5m[:,1],color='grey',linestyle='dashed')
        plt.plot(self.border_8m[:,0],self.border_8m[:,1],color='grey',linestyle='dashed')
        if(self.border_tree_flag == True):
            self.ax_anim.scatter(self.border_trees[:,0],self.border_trees[:,1],color='blue')

        self.scat_anim = self.ax_anim.scatter(self.inner_tree[:,0],self.inner_tree[:,1],color='blue')

        self.animation = animation.FuncAnimation(self.fig_anim,self.update_anim,frames=10000,interval = 5)
        plt.show()

    def update_anim(self,i):
        self.update_inner_tree()
        self.scat_anim.set_offsets(np.c_[self.inner_tree[:,0],self.inner_tree[:,1]])


       
        