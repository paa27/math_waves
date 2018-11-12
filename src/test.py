#Time tests

from standing import Standing

if __name__ == '__main__':
	import timeit
	from line_profiler import LineProfiler
	lp = LineProfiler()
	lp_wrapper = lp(Standing(max_order=2,harmonics=2).solve(N=100,method='homemade'))
	lp_wrapper.print_stats()