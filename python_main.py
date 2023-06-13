import numpy as np
import python_utils as utils
import oliveto

data_path = '/home/daverio/Documents/oliveto_data/'

def main():
    sim = oliveto.Oliveto()
    sim.set_border(np.loadtxt(data_path + 'border.csv', delimiter=','))

if __name__ == "__main__":
    main()
