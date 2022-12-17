from line_profiler import LineProfiler
from MPC_run import main

lp = LineProfiler()


main = lp(main)
main()
lp.print_stats()


