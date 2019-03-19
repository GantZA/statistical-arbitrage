# Test for Statistical Arbitrage 
<br />
Test for Statistical Arbitrage (modified Semi-Variance definition) following the methodology in Hogan et al. (2004) and Jarrow et al. (2012) coded up in Julia 1.1. Usage requires running functions, not very user friendly at the moment.
<br />
The CM test can be used by opening the test_for_cm_stat_arb.jl file and using the run block functionality to load all functions. The gen_min_t_dist_cm() function allows one to generate the empirical Min-t distribution for a given size of the observed incremental profit process with mu, sigma and lambda parameters using a Monte Carlo method with MLE for a given number of iterations. 
<br />

<br />
References <br />
<br />
HOGAN, S., JARROW, R., TEO, M.,AND WARACHKA, M. (2004). Testing market efficiency using statistical arbitrage with applications to momentum and value strategies. Journal of Financial economics,73(3), 525–565.

JARROW, R., TEO, M., TSE, Y. K.,AND WARACHKA, M. (2012). An improved test for statistical arbitrage. Journal of Financial Markets, 15 (1), 47–80.
