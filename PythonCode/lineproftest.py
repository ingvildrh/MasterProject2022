from line_profiler import LineProfiler

lp = LineProfiler()

from MPC_run import main
main = lp(main)
main()
lp.print_stats()
lp.dump_stats(filename='testlineprof.txt')

