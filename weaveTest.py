# weaveTest.py

import math
import numpy as np
import scipy.weave as weave
import scipy.weave.converters as converters

def quantum_cat(N, kappa):
    phi = np.zeros((N,N), dtype='float64')
    pi = math.pi
    support = "#include <math.h>"
    code = """
float alpha = 2.0*pi/N;
float kap_al = kappa / alpha;

for (int k=0; k<N; ++k)
	for (int l=0; l<N; ++l)
		phi(k, l) = alpha*(k*k - k*l + l*l) + kap_al*sin(alpha*l);
"""
    weave.inline(code, ['N', 'kappa', 'pi', 'phi'],
                 type_converters = converters.blitz,
                 support_code = support, libraries = ['m'])
    return (1.0 / math.sqrt(N))*np.exp(1j * phi)