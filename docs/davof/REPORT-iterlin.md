# The static test of the iterated face-fraction corrector with the ring linear fit (branch davof-iter-linfit)

Executor: Claude (a sub-session executing the plan of the designing session). Workspace /work/scratch/tm83tomy/davof/agents/iterlin.

## Pre-registration

Pre-registered 2026-11-25 (Claude, before any run). Decision norm: L1 of 1 - n.n_e on the common cell set; orders from N 64 to 128.

PCL0 (port check): on the exact state, pass 0, the L1 of 1 - n.n_e of ringLinear, of ringArea and of plicRDF agrees with the fork's gate table (ringLinear 1 pass, ringSum, plicRDF) within 1 % at each of the 6 states. If not, the difference is located before any corrector pass is run.

PCL1 (contraction): from starts A and B at N 64 and N 128 the L1 change of alpha_f falls from pass 2 to pass 6 by a factor of at least 2 per pass in the geometric mean, and does not grow over two consecutive passes. Predicted to hold.

PCL2 (one fixed point): at pass 20 the normals from starts A, B and E agree to max |n_A - n_B|, max |n_A - n_E| <= 1e-6 at N 64 and N 128. Predicted to hold if PCL1 holds.

PCL3 (order at the fixed point): the L1 of |n - n_e| at pass 20 is first order, its order from N 64 to 128 lies in [0.7, 1.4]. Reason: with the mean of the two clips, exact normals are not reproduced; in 2D a cell with a cut of length L_0 between neighbours with cuts of length L_1 and L_2 gets an area vector tilted by kappa (L_2^2 - L_1^2) / (24 L_0), which is O(kappa h) wherever the cut lengths of neighbouring cells differ.

PCL4 (the floor): row R5 at pass 0 (the fit on the fractions that exact-normal planes induce) has an L1 of |n - n_e| at least 10 times that of row R3 at pass 0 (the fit on the exact fractions), at N 64 and N 128.

PCL5 (parity): at pass 20 the ratio to plicRDF is above 1.05 at N 128 for both placements, so parity is NOT reached there; at all 6 states the ringLinear fixed point is below the ringArea fixed point. Falsified in its first part if the ratio is at most 1.05 at N 128 for either placement.

Decision rule (the goal is plicRDF's accuracy at lower cost, Tomislav 2026-11-25; the numbers 1.05 and 2 are Claude's): the corrector with the linear fit is a lead if, from start B, after at most 2 passes the ratio to plicRDF is at most 1.05 at N 64 and N 128 for both placements (4 of 4), and the time of those passes is below the time of plicRDF's reconstruction on the same state. Otherwise it is not a lead and no transport run follows. Row R5 then says whether the re-cut rule or the estimator sets the error.

## Pre-registration of the small transport test (isoAdvector, shear3D)

Pre-registered 2026-11-25 (Claude, before any transport run). Case shear3D (leia isoShear3D), hex, N 32 (1 rank) and N 64 (8 ranks), OpenFOAM Co 0.2 and 0.5, advection isoAdvector, the same initial alpha in every row. Rows: PLIC = plicRDF; AREA0 = davof, ringArea, no corrector; LIN0 = davof, ringLinear, no corrector; LIN1 = LIN0 with faceFractionCorrector mean, correctorIterations 1, correctorTolerance 0; LIN2 = the same with correctorIterations 2. Ratios are to PLIC of the same job.

PP (parallel gate, static): on the exact sphere N 64 standard, davofTestCorrector with ringLinear on 8 ranks gives the L1 of |n - n_e| and of 1 - n.n_e at passes 0, 1, 2 equal to the serial values to a relative 1e-10. If not, N 64 runs on 1 rank.

PS1 (the gap exists in leia): E_shape(LIN0)/E_shape(PLIC) > 1.05 at N 64 at both Co.

PS2 (the corrector helps): E_shape(LIN1) < E_shape(LIN0) at N 64 at both Co.

PS3 (parity not reached): E_shape(LIN1)/E_shape(PLIC) > 1.05 at N 64 at both Co. Falsified if it is at most 1.05 at either Co.

PS4 (one pass suffices): E_shape(LIN2) within 2 % of E_shape(LIN1) at N 64 at both Co.

PS5 (cost): wall(LIN0)/wall(PLIC) < 1 and wall(LIN1)/wall(PLIC) > 1 at N 64 at both Co.

Decision rule (goal: plicRDF's accuracy at significantly lower cost, Tomislav 2026-11-25; the numbers 1.05 and 0.8 are Claude's): the corrector is a lead for the isoAdvector track if at N 64 at both Co E_shape(LIN1)/E_shape(PLIC) <= 1.05 and wall(LIN1)/wall(PLIC) <= 0.8. If only the accuracy clause holds, it is an accuracy lead whose cost is the open question. If the accuracy clause fails, the face fractions are not what separates the fit from plicRDF on this flow, and the corrector is closed for this round. N 32 is reported and not part of the rule.

(Run before this text was committed: the parallel gate PP (job 57252508; serial and 8-rank values agree to a relative 1.3e-11 at most, so PP held and N 64 runs on 8 ranks) and one probe run, PLIC at N 64, Co 0.5, 8 ranks (job 57252878, ClockTime 39 s).)
