from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np

problem = {
    'num_vars': 3,
    'names': ['x1', 'x2', 'x3'],
    'bounds': [[-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359],
               [-3.14159265359, 3.14159265359]]
}

param_values = saltelli.sample(problem, 1000, calc_second_order=False)
'''
Y = np.zeros([param_values.shape[0]])

for i, X in enumerate(param_values):
    Y[i] = Ishigami.evaluate(X)'''
Y = Ishigami.evaluate(param_values)
    
Si = sobol.analyze(problem, Y)

print(Si['S1'])