import sys
from TopoIndex import Acupuncture

filename = sys.argv[1]
data = Acupuncture(filename)

data.compute_index_beta()
data.compute_index_alpha()
data.compute_index_c(4)

print(f'Bifurcation loop index: {data.beta}, angle index: {data.alpha}, curvature index: {data.c}')
