import numpy as np
import python_utils as utils
import oliveto

data_path = '/home/daverio/Documents/oliveto_data/'

def main():
    sim = oliveto.Oliveto()
    sim.set_border(np.loadtxt(data_path + 'border_8m.csv', delimiter=','))
    sim.add_tree_to_border(7,1)
if __name__ == "__main__":
    main()
