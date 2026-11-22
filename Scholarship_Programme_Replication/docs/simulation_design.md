# Simulation design

Each programme is one cluster and awards scholarships to the top 10 percent of student applicants. Normalized pre-award academic scores are independent Uniform[0,1], the reference cutoff is 0.9, and the realized cutoff is the midpoint between descending ranks k and k+1, where k=m/10.

The outcome is a subsequent academic performance index. In the main model, the index without a scholarship is 40+30X+epsilon, with independent Normal(0,9) noise. Receiving a scholarship adds 5 index points. Both potential outcomes share the same student-specific noise. The realized outcome is observed for all applicants.

Score/cutoff sizes: [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]. RD starts at m=100. The bandwidth remains 0.5*m^(-1/4); the kernel is triangular. Fit separate weighted lines, subtract intercepts, and use the squared-WLS-weight HC0 sandwich. No bandwidth or seed search is performed. The small experiment holds m=10 and k=1 and increases G over [20, 50, 100, 200, 500, 1000, 2000].

For ascending order positions r=m-k and r+1, the midpoint mean is (r+0.5)/(m+1). The order-statistic covariance formula supplies its exact variance. RMSE relative to c=0.9 includes the squared finite-sample bias. In particular, the expected midpoint at m=10 is 86.363636 score points, even though the population reference is 90; this is not an implementation error.

An adjacent uniform spacing is Beta(1,m), with mean 1/(m+1) and variance m/((m+1)^2*(m+2)). Thus the small-cluster outcome contrast has mean 5+30/11 and variance 30^2*10/(11^2*12)+18. Its G-cluster sample mean has this variance divided by G. Scores shown in the figures use a factor of 100; outcome points have their own scale.

The nonlinear variant and analytic slope sensitivity are numerical supplements only. The three figures have distinct purposes: cutoff mechanics; the actual large-cluster estimator; and precision versus causal interpretation in repeated small clusters.
