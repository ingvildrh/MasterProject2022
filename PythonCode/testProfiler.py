import pstats
import cProfile

prof = cProfile.Profile()

prof.enable()
a = 0
for i in range(1000000):
    a=a+i
prof.disable()

prof.print_stats()
