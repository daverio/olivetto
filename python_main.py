import numpy as np
import python_utils as utils
import oliveto
import matplotlib.pyplot as plt

data_path = '/home/daverio/Documents/oliveto_data/'

def main():
    sim = oliveto.Oliveto()
    sim.set_border(np.loadtxt(data_path + 'border_8m.csv', delimiter=','))
    sim.add_tree_to_border(7,1)

    borderTrees = sim.get_border_trees()

    sim.add_interior_tree(130)

    inerTrees = sim.get_iner_trees()
    for i in range(1000):
        sim.make_step(0.2,10,20,0.05,0.1)
    inerTrees_onestep = sim.get_iner_trees()

    plt.figure()
    plt.scatter(borderTrees[:,0],borderTrees[:,1])
    #plt.scatter(inerTrees[:,0],inerTrees[:,1])
    plt.scatter(inerTrees_onestep[:,0],inerTrees_onestep[:,1])
    plt.show()


if __name__ == "__main__":
    main()
