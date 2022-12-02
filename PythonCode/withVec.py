from timeit import default_timer as timer
import numpy as np
arr = np.ones((100,1))
arr2 = np.zeros((100,1))
arr4 = np.zeros((100,1))
start = timer()
arr2 = arr*2
arr4 = arr*4
end = timer()
timetot = end-start
print(timetot)
