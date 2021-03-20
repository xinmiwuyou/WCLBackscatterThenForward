# WCLBackscatterThenForward
Simulation codes for Lyu et al., “Backscatter then Forward: A Relaying Scheme for Batteryless IoT Networks”, IEEE Wireless Communications Letters vol. 9, no. 4, pp. 562-566, April 2020.

An introduction to reproduce Fig2. (a) and Fig2. (b) of this paper.

(1) Inside the “Results for Fig. 2(a) and Fig. 2(b)” folder, we create two folders, namely Fig. 2(a) and Fig. 2(b), and each folder contains results for Fig. 2(a) and Fig. 2(b).

(2) In each folder (i.e., Fig. 2(a) and Fig. 2(b)), please run the “Main.m” file to generate the exactly same figures in our paper, i.e., Fig. 2(a) and Fig. (b).

(3) In each sub-folders, we provide a brief guide, i.e., the “readme” file. Since the “Main.m” file needs to solve six independent problems (three schemes for both $M=5$ and $M=10$) in each iteration, running it usually takes a lot of time. If you want to quickly check the codes, please set a small number for random channel generations (e.g., 10), which is defined by “Count” for all codes. If you want to obtain more accurate results, “Count” should be sufficiently large.

(4) For the proposed scheme, we prepare the codes based on Algorithm 1 as presented in our paper. For more details, please refer to the paper. For the two benchmark schemes (the benchmark with random energy beamforming and the benchmark with equal time allocation), the corresponding problems are both convex optimization problems. So, we use the CVX tool (http://cvxr.com/cvx/) to quickly obtain the results. Note that the CVX tool should be installed in MATLAB before running the codes for the benchmark schemes. Otherwise, these codes cannot work.
