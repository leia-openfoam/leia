# Literature survey: error amplification of semi-Lagrangian level-set schemes on unstructured FV meshes

Date: 2026-11-12. Method: web search and page/abstract fetch in this session.
Marks: VERIFIED = read in the source or its abstract online in this session (URL given).
UNVERIFIED = from memory or from a secondary source (search-engine summary). No reference is invented.
Status of the file: written incrementally, topic by topic.

## Topic 1. Stability of semi-Lagrangian schemes with Lagrange interpolation

### 1.1 Besse and Mehrenberger 2008
N. Besse, M. Mehrenberger, "Convergence of classes of high-order semi-Lagrangian schemes for the
Vlasov-Poisson system", Math. Comp. 77 (261), 93-123, 2008. PII S 0025-5718(07)01912-6.
Full text read: https://www.i2m.univ-amu.fr/perso/mehrenberg.m/BeMe.pdf

- VERIFIED. The paper proves L2 convergence at rate O(dt^2 + h^(m+1) + h^(m+1)/dt), m = polynomial degree,
  on uniform periodic grids.
- VERIFIED (introduction). Lagrange interpolation on finite elements (not centred on the foot) gives
  "unstable schemes in any Lp norm for orders greater than two". Symmetric (centred) Lagrange interpolation on
  uniform grids recovers L2 stability.
- VERIFIED (Corollary 1, Remark 3). sup |theta_(2n-1)| = 1 on [-1/2, 1/2] and sup |theta_(2n)| = 1 on [-1, 1].
  theta_m is the amplification factor of degree-m centred Lagrange interpolation. The remark says that this
  "was already proven by Strang" (Strang 1962). The stability zone is exactly [-1, 1] in the even-degree case
  and [-1/2, 1/2] in the odd-degree case, at least for n <= 82 (computer proof of Hong and Steinberg, preprint 2001).
  Meaning: even degree (odd number of nodes, centred on a NODE): the foot must lie within one cell of the central
  node. Odd degree (even number of nodes, centred on an INTERVAL): the foot must lie in the central interval.
- VERIFIED (Remark 4). "The analysis ... on uniform grids does not apply if we take the same scheme on
  unstructured meshes (triangulation) by using high-order Lagrange interpolation on triangles. The determination
  of the stability region in that case remains an interesting open question." Also: when the foot falls where the
  amplification factor exceeds one, the scheme is unstable; this explains the bad behaviour of standard
  finite-element reconstruction "whether it is on uniform grids or unstructured meshes".

Relevance. Result (a) in 1D (degree 2, centred on the cell, |G| = 1 max for Co <= 1, = Lax-Wendroff) is the
n = 1 even case of Strang's theorem. Result (b) "linear fit = FTCS" is NOT the standard odd-degree case: a
standard linear SL scheme interpolates in the upwind interval and is stable. The least-squares linear fit is
centred on the cell (symmetric 3-point data, degree 1), which is a different object (see Topic 8).
Remark 4 states that the stability region on unstructured meshes was open in 2008.

### 1.2 Strang 1962
G. Strang, "Trigonometric polynomials and difference methods of maximum accuracy", J. Math. and Phys. 41,
147-154, 1962. DOI 10.1002/sapm1962411147 (DOI UNVERIFIED, from memory).
- UNVERIFIED in the primary source (the PDF was located at
  https://people.csail.mit.edu/wangliang/Papers/Multi-points_Scheme.pdf but only the search summary was read).
  The summary confirms: schemes of maximum accuracy 2N use the Lagrange coefficients, and they are stable
  "under quite reasonable restrictions".
- VERIFIED through Besse-Mehrenberger Remark 3: Strang proved sup|theta| = 1 on [-1,1] (even degree) and
  [-1/2,1/2] (odd degree).

### 1.3 Courant, Isaacson, Rees 1952
R. Courant, E. Isaacson, M. Rees, "On the solution of nonlinear hyperbolic differential equations by finite
differences", Comm. Pure Appl. Math. 5, 243-255, 1952. UNVERIFIED (from memory, not fetched).
The CIR scheme is the linear-interpolation SL scheme (first-order upwind for Co <= 1). It is monotone and
stable for every Courant number when the interpolation interval contains the foot.

### 1.4 Ferretti 2010, 2013; Ferretti and Mehrenberger 2020; Falcone and Ferretti book
- R. Ferretti, "Equivalence of semi-Lagrangian and Lagrange-Galerkin schemes under constant advection speed",
  J. Comput. Math. 28(4), 461-473, 2010. DOI 10.4208/jcm.1003-m0012.
  VERIFIED (abstract, https://www.global-sci.com/index.php/JCM/article/view/11995): "SL schemes are equivalent
  to exact LG schemes via a suitable definition of the basis functions". Hypotheses: constant advection speed,
  uniform grid, symmetric and translation-invariant cardinal basis. Gives a simpler stability proof for
  constant coefficients. UNVERIFIED detail: the proof covers centred Lagrange interpolation of odd degree up
  to 13 (stated in the abstract of the 2020 paper, see below, which is VERIFIED).
- R. Ferretti, "On the relationship between semi-Lagrangian and Lagrange-Galerkin schemes", Numer. Math. 124,
  31-56, 2013. DOI 10.1007/s00211-012-0505-5. UNVERIFIED (search summary only): area-weighted LG schemes are a
  small perturbation of exact LG schemes; this gives stability of SL schemes with variable coefficients.
- R. Ferretti, M. Mehrenberger, "Stability of semi-Lagrangian schemes of arbitrary odd degree under constant
  and variable advection speed", Math. Comp. 89 (324), 1783-1805, 2020. DOI 10.1090/mcom/3494.
  VERIFIED (Crossref abstract, https://api.crossref.org/works/10.1090%2FMCOM%2F3494): generalises the
  SL = LG equivalence "to an arbitrary odd degree, for both the case of constant- and variable-coefficient
  equations" and "the same holds for spline interpolations".
- M. Falcone, R. Ferretti, "Semi-Lagrangian Approximation Schemes for Linear and Hamilton-Jacobi Equations",
  SIAM, 2014. DOI 10.1137/1.9781611973051. UNVERIFIED (from memory). Chapter on stability collects the
  von Neumann analysis for Lagrange interpolation and the LG equivalence.

Courant limits. Odd degree with the stencil re-centred on the interval of the foot: stable for EVERY Courant
number (the stencil moves with the foot). Even degree: stable when the foot is within half a cell of the
central node after re-centring; with a FIXED stencil centred on the arrival cell: |Co| <= 1.
The LG-equivalence theory is for ODD degree only. Our quadratic (even degree) is covered by Strang 1962 only.

Relevance to (a)-(g). (a) in 1D: KNOWN (Strang 1962; Lax-Wendroff). The multi-dimensional statement for a
least-squares fit through the full one-ring is not in these sources. (b): the instability of Lagrange
interpolation that is not centred on the foot is KNOWN (Besse-Mehrenberger introduction and Remark 4).

## Topic 2. SL schemes with least-squares / RBF / WENO / element reconstruction on unstructured meshes

### 2.1 Bradley 2024 (the most relevant single source for (e) and (f))
A. M. Bradley, "Stabilized bases for high-order, interpolation semi-Lagrangian, element-based tracer
transport", J. Comput. Phys., 2024. DOI 10.1016/j.jcp.2024.113034. Accepted manuscript read in full text:
https://www.osti.gov/servlets/purl/2345299
Companion: A. M. Bradley, P. A. Bosler, O. Guba, "Islet: interpolation semi-Lagrangian element-based
transport", Geosci. Model Dev. 15, 6285-6310, 2022, DOI 10.5194/gmd-15-6285-2022
(https://gmd.copernicus.org/articles/15/6285/2022/; VERIFIED Sect. 2.1.2: "The ISL method using the natural GLL
basis is unstable for the test problem of uniform flow on a uniform grid").

- VERIFIED. Method of analysis: A = space-time evolution matrix for periodic uniform flow. Necessary condition:
  maximum eigenvalue magnitude <= 1 for all translation distances. For an element-based method "Fourier modes
  are not eigenvectors ... and thus Fourier amplification factors are not eigenvalues"; with N elements of
  degree d one solves d x d eigenvalue problems per wave number (this is a Bloch analysis).
- VERIFIED. Review of earlier work: Baptista (1987) and Oliveira and Baptista (1995) "find an ISL method using
  a compact natural interpolant is unstable on the test problem for elements having degree three and higher".
  Bradley makes "a similar, but more general, statement" (Sect. 3.8). So the natural element interpolant of
  degree 1 and 2 passes; degree >= 3 fails.
- VERIFIED. Quote of Ferretti 2010: the LG-equivalence analysis "cannot handle ... reconstructions performed by
  finite elements, which violate the assumptions of symmetry and translation invariance". Quote of
  Ferretti-Mehrenberger 2020: "a partial result of stability for the case of piecewise quadratic finite element
  interpolations ... will be the object of a forthcoming study", but "a rigorous extension to the
  multidimensional case still seems out of reach".
- VERIFIED. Nonuniform grids: quote of Ferretti 2013 Sect. 5: "it is known that using [symmetric Lagrange
  interpolation] with unevenly spaced nodes results in general in an unstable scheme". Bradley gives a
  reproducible example: periodic 1D grid {0, 0.11242, 0.44817, 0.78392, 0.88737}, classical cubic interpolation,
  shift 0.33575: maximum eigenvalue magnitude > 1 + 1e-3.
- VERIFIED. "Sparse instability": unstable SL methods often have maximum eigenvalue 1 + eps with eps <~ 1e-2,
  "with no particular dependence on time step", in isolated spikes over the shift. Bradley speculates that ALL
  high-order ISL methods have at least a sparse pattern of lambda_max > 1 on some nonuniform grid or flow.
- VERIFIED. Perturbed uniform grid heuristic (Sect. 3.5.2, Fig. 4): element boundaries perturbed by a
  fraction delta of the subelement size; empirically lambda_max - 1 is proportional to delta^4 for the stabilised
  bases; some bases fall off more slowly.

Relevance. (e) "weak amplification on perturbed meshes" is KNOWN IN A DIFFERENT FORM: for node-based
tensor-product element interpolation in 1D/2D, with the scaling (lambda_max - 1) ~ delta^4. Not for a
cell-centred least-squares fit. (f) the Q2 point-value SL scheme is the degree-2 "natural interpolant"
element-based ISL method; its stability on the uniform test problem is KNOWN (Baptista 1987, Oliveira-Baptista
1995, Bradley 2024, Ferretti-Perrone 2007), but only as "eigenvalue <= 1", by computation, not as a proof.

### 2.2 Ferretti and Perrone 2007
R. Ferretti, G. Perrone, "On the stability of semi-Lagrangian advection schemes under finite element
interpolations", in Applied and Industrial Mathematics in Italy II, World Scientific, 2007, pp. 339-350.
- Existence VERIFIED (reference [15] of Elarif-Mehrenberger-Navoret 2026, https://arxiv.org/pdf/2607.00619).
  Content UNVERIFIED (not online). Elarif et al. cite it for: "Regarding nodal SL schemes on unstructured
  triangular grids, one popular solution is to use P2 interpolation [15], but P3 or higher order interpolation
  are unstable [11]" ([11] = Crouseilles, Mehrenberger, Vecil, ESAIM Proc. 32, 211-230, 2011).
- This is the paper to obtain for (f) and for the P2/Q2 symbol on triangles and quadrilaterals.

### 2.3 Elarif, Mehrenberger, Navoret 2026
"Hermite Semi-Lagrangian schemes on triangular meshes for advection equations", arXiv:2607.00619 (1 July 2026).
VERIFIED (full text). Statements: SL extension to unstructured meshes "is more challenging due to stability and
computational issues"; three routes for nodal SL on triangles: (i) P2 interpolation, (ii) "a larger stencil
together with least square reconstruction [8, 3]", (iii) Hermite interpolation with transported gradients
(Besse-Sonnendruecker 2003). Stability is assessed only numerically. The Bell element develops instabilities
for some time steps; Argyris does not.

### 2.4 Least-squares SL on unstructured cell-centred meshes
- G. Carre, S. Del Pino, K. Pichon Gostaf, E. Labourasse, A. V. Shapeev, "Polynomial least-squares
  reconstruction for semi-Lagrangian cell-centered hydrodynamic schemes", ESAIM: Proc. 28, 100-116, 2009.
  DOI 10.1051/proc/2009041. Existence VERIFIED (search record and reference [8] of Elarif et al.). Content
  UNVERIFIED (403 on the publisher page). No stability analysis is reported in the abstract (search summary).
- M. Bergami, W. Boscheri, G. Dimarco, "A high-order conservative semi-Lagrangian solver for 3D free surface
  flows with sediment transport on Voronoi meshes", Commun. Appl. Math. Comput., 2020. Existence VERIFIED
  (reference [3] of Elarif et al.). Content UNVERIFIED.
- W. Boscheri, M. Dumbser and co-workers: high-order SL on Voronoi/unstructured meshes with CWENO/least-squares
  reconstruction. UNVERIFIED (from memory). I know of no von Neumann analysis in these papers.
- A. Iske, M. Kaeser, "Conservative semi-Lagrangian advection on adaptive unstructured meshes", Numer. Methods
  PDE 20(3), 388-411, 2004, DOI 10.1002/num.10100. VERIFIED (search record): finite volumes on a Voronoi
  diagram, adaptive. No stability analysis reported.
- J. Behrens, A. Iske, "Grid-free adaptive semi-Lagrangian advection using radial basis functions", Comput.
  Math. Appl. 43(3), 319-327, 2002, DOI 10.1016/S0898-1221(01)00289-9. VERIFIED (search record): SL + local
  RBF (thin-plate spline) interpolation on scattered nodes. No stability analysis reported.
- E. Carlini, R. Ferretti: SL level-set scheme with RBF interpolation, usable on unstructured grids.
  UNVERIFIED (search summary only; exact title not confirmed).
- E. Carlini, R. Ferretti, G. Russo, "A weighted essentially nonoscillatory, large time-step scheme for
  Hamilton-Jacobi equations", SIAM J. Sci. Comput. 27(3), 1071-1091, 2005/2006. VERIFIED (search record;
  preprint https://www.dmi.unict.it/~russo/ps/Russo/CaFeRu-HJ-WENO-SISC05.pdf): WENO interpolation in an SL
  scheme, positivity of the WENO weights for any order. Structured grids.
- N. Besse, E. Sonnendruecker, "Semi-Lagrangian schemes for the Vlasov equation on an unstructured mesh of
  phase space", J. Comput. Phys. 191(2), 341-376, 2003. VERIFIED through Besse-Mehrenberger 2008 (introduction):
  on triangles, Hermite interpolation gives "high-order, stable and slightly diffusive schemes", Lagrange
  interpolation gives "too diffusive schemes for orders smaller than two and unstable schemes in any Lp norm
  for orders greater than two".
- M. Restelli, L. Bonaventura, R. Sacco, "A semi-Lagrangian discontinuous Galerkin method for scalar advection
  by incompressible flows", J. Comput. Phys. 216, 195-215, 2006, DOI 10.1016/j.jcp.2005.11.030.
  UNVERIFIED (search summary): unconditional von Neumann stability proven in 1D only.
- M. Mehrenberger, L. Mendoza, C. Prouveur, E. Sonnendruecker, "Solving the guiding-center model on a regular
  hexagonal mesh", ESAIM: Proc. 53, 149-176, 2016, DOI 10.1051/proc/201653010. UNVERIFIED (search summary):
  nodal SL on the equilateral-triangle (hexagonal) lattice with box splines or Hermite elements. Node-based.
  No Courant limit of sqrt(3)/2 is reported.
- P. H. Lauritzen, R. D. Nair, P. A. Ullrich, CSLAM, J. Comput. Phys. 229, 1401-1424, 2010. Existence VERIFIED
  (reference list). Cell-integrated, not point-value; not directly comparable.
- O. Desjardins, V. Moureau, H. Pitsch, "An accurate conservative level set/ghost fluid method ...", J. Comput.
  Phys. 2008: the search summary attributes "semi-Lagrangian transport + least-squares curvature" to
  a Desjardins paper. UNVERIFIED; the exact paper was not confirmed.

Conclusion for topic 2. No source found does a von Neumann or Bloch analysis of a cell-centred least-squares
SL scheme. The community states that stability on unstructured meshes is hard and mostly open
(Besse-Mehrenberger Remark 4; Ferretti-Mehrenberger 2020; Bradley 2024; Elarif et al. 2026).

## Topic 8. Von Neumann / matrix analysis of least-squares reconstruction on unstructured or perturbed meshes

- F. Haider, J.-P. Croisille, B. Courbet, "Stability analysis of the cell centered finite-volume MUSCL method on
  unstructured grids", Numer. Math. 113, 555-600, 2009, DOI 10.1007/s00211-009-0242-6. Preprint:
  https://jeanpierre-croisille.perso.math.cnrs.fr/Postscript/fhaider-article-svjour-nummat-20080327.pdf
  UNVERIFIED in the full text (search summaries of the paper and of the ICCFD5 companion
  https://jeanpierre-croisille.perso.math.cnrs.fr/Postscript/iccfd5-juillet-2008.pdf): method-of-lines
  eigenvalue stability of linear advection with linear (MUSCL) reconstruction, periodic boundaries; the
  least-squares gradient minimises every unitarily invariant norm of the reconstruction matrix; larger stencils
  should increase stability; instability can occur on 3D unstructured (tetrahedral) grids with the first
  neighbourhood. This is SEMI-DISCRETE (Eulerian flux form), not semi-Lagrangian.
- Zandsalimy and Ollivier-Gooch: eigenanalysis of the semi-discrete FV operator on unstructured meshes; vertex
  perturbation to remove unstable modes. UNVERIFIED (search summary; title not confirmed).
- "On the properties of high-order least-squares finite-volume schemes", J. Comput. Phys. 2022.
  UNVERIFIED (search record only).
- WLS-ENO (Liu, Jiao), arXiv:1602.02694: von Neumann section for a weighted-least-squares FV scheme.
  UNVERIFIED (search summary).
- Bradley 2024 (see 2.1): the only source found that computes the amplification of an SL scheme on a randomly
  perturbed periodic grid. Node-based, tensor product.

Conclusion for topic 8. The Bloch analysis of the fully discrete SL update with a one-ring quadratic
least-squares fit, on Cartesian, hexagonal and perturbed periodic meshes, was NOT FOUND.

## Topic 3. Schemes that advect point values (Active Flux, CIP, multi-moment, element-based nodal SL)

### 3.1 Element-based nodal SL (the direct counterpart of result (f))
Result (f) is NOT Active Flux. It has no cell average. It is a nodal SL scheme with the natural Q2 element
interpolant (nodes at vertices, edge midpoints, face centres, cell centre). The literature names are
"compact interpolator" (Baptista), "element-based tensor-product interpolation SL" (Bradley), "SL under
finite element interpolation" (Ferretti-Perrone).
- A. M. Baptista, PhD thesis MIT 1987; A. Oliveira, A. M. Baptista, "A comparison of integration and
  interpolation Eulerian-Lagrangian methods", Int. J. Numer. Methods Fluids 21(3), 183-204, 1995.
  VERIFIED through Bradley 2024 (Sect. 2.1 of the accepted manuscript): compact natural interpolants are
  unstable for element degree three and higher. Primary source UNVERIFIED.
- Bradley 2024, Ferretti-Perrone 2007, Ferretti-Mehrenberger 2020: see Topic 2. A proof for the quadratic
  element is announced as "a forthcoming study" in 2020; I found no published proof.
- Giraldo 1998 (spectral-element Lagrange-Galerkin/SL): cited by Bradley as a Fourier-mode analysis that
  concludes instability. UNVERIFIED in the primary source; exact title not confirmed.

### 3.2 Active Flux (cell average + point values on the cell boundary)
- T. A. Eymann, P. L. Roe, "Active flux schemes", AIAA 2011-382, 2011. VERIFIED (search record of the Deep Blue
  PDF): edge values are independent unknowns; the linear-advection case is van Leer's Scheme V.
  In 1D the point values are updated by tracing the characteristic into the quadratic reconstruction.
  UNVERIFIED (from memory): stable for CFL <= 1 in 1D.
- E. Chudzik, C. Helzel, D. Kerkmann, "The Cartesian grid Active Flux method: linear stability and bound
  preserving limiting", Appl. Math. Comput. 393, 125501, 2021. DOI 10.1016/j.amc.2020.125501.
  VERIFIED (abstract, https://ideas.repec.org/a/eee/apmaco/v393y2021ics0096300320304598.html): 2D Cartesian
  linear stability analysis; for advection, "stability for CFL <= 1 requires a more accurate flux computation
  than previously assumed"; for acoustics stability for CFL <= 1/2 is confirmed.
  Lesson: with exact point-value evolution the scheme is stable to CFL 1 only if the remaining ingredient (the
  flux quadrature for the average) is accurate enough. The point-value part alone is not the limit.
- W. Barsukow, J. Kern, C. Klingenberg, L. Lechner, "Analysis of the multi-dimensional semi-discrete Active Flux
  method using the Fourier transform", arXiv:2412.03477; Commun. Appl. Math. Comput., DOI
  10.1007/s42967-025-00512-3. VERIFIED (abstract and HTML, https://arxiv.org/html/2412.03477v2): in 3D each
  cell accesses 27 degrees of freedom: the average plus point values at 8 corners, 12 edge midpoints and
  6 face centres; reconstruction space P^{2,2,2}. Semi-discrete (Runge-Kutta) version: maximal CFL about 0.2 in
  multi-D, about 0.28 for acoustics. The 3D point set is the same as in result (f), with the cell average in
  place of the centre value.
- D. Calhoun, E. Chudzik, C. Helzel, "The Cartesian grid Active Flux method with adaptive mesh refinement",
  arXiv:2205.07572 (journal version probably J. Sci. Comput. 2023, DOI 10.1007/s10915-023-02106-8, mapping
  UNVERIFIED). VERIFIED (abstract): first Active Flux on adaptively refined Cartesian grids (ForestClaw,
  patch-based, ghost cells), third order, conservative, subcycling in time; "compact stencil in space and time".
  The abstract gives no stability analysis of the coarse-fine interface.
- R. Abgrall, J. Lin, Y. Liu, "Active flux for triangular meshes for compressible flows problems",
  arXiv:2312.11271. UNVERIFIED (search summary): on triangles, linear stability is hard to analyse; quadratic
  and cubic versions were stable in tests with the CFL limit about one half and one third of the standard value.
- Porfetye, Tang, Chu, Helzel, Lukacova-Medvid'ova, arXiv:2512.16359. UNVERIFIED (search summary): exact
  evolution operator stable to CFL 0.5 (acoustics), approximate EG2 operator CFL <= 0.279.

### 3.3 CIP, CIP-CSL, multi-moment
- T. Yabe, F. Xiao, T. Utsumi, "The constrained interpolation profile method for multiphase analysis",
  J. Comput. Phys. 169, 556-593, 2001, DOI 10.1006/jcph.2000.6625. Existence VERIFIED (search record).
  CIP advects the value AND the gradient at grid points semi-Lagrangially (cubic Hermite in the upwind cell).
- "Stability and accuracy of the Cubic Interpolated Propagation scheme" (Utsumi and co-authors; authors not
  confirmed). UNVERIFIED (search summary): Fourier analysis; stable with small gain and phase errors.
- F. Xiao, S. Ii and co-workers: CIP/multi-moment FV, VSIAM3 (cell averages + point values at vertices/face
  centres), also on triangles/tetrahedra (Xie, Jin, Xiao). Existence VERIFIED (search records). No Fourier
  stability result with a Courant limit was retrieved. UNVERIFIED.
- Besse-Sonnendruecker 2003 and Elarif et al. 2026 (Topic 2): Hermite point-value SL on triangles is the
  stable option there.

Relevance. (f) first part (full Q2 set, growth factor 1): KNOWN IN A DIFFERENT FORM. (f) second part
(the 3D set without edge midpoints amplifies) and third part (centre + vertices cannot carry a full quadratic):
NOT FOUND as stability statements. The third part is a counting fact: 9 data, 10 coefficients of P2 in 3D.

## Topic 4. BCC lattices and non-Cartesian sample sets

- T. Theussl, T. Moeller, M. E. Groeller, "Optimal regular volume sampling", IEEE Visualization 2001.
  UNVERIFIED (secondary: cited in search results): the BCC lattice needs about 30 % fewer samples than the
  Cartesian lattice for a spherically band-limited signal. A 2014 EuroVis paper is reported to say that the
  exact percentage was reported wrongly and that CC, BCC and FCC are each optimal in a different frequency
  range (UNVERIFIED, title not confirmed).
- A. Entezari, R. Dyer, T. Moeller, "Linear and cubic box splines for the body centered cubic lattice",
  IEEE Visualization 2004. VERIFIED (search record of the abstract): piecewise linear and piecewise cubic box
  spline reconstruction filters on BCC.
- A. Entezari, D. Van De Ville, T. Moeller, "Practical box splines for reconstruction on the body centered
  cubic lattice", IEEE TVCG 14(2), 2008. VERIFIED (search record of the abstract): four-direction box spline
  gives a linear C0 reconstruction; C2 quintic version; BCC is reported about twice as efficient as
  tensor-product B-splines on the Cartesian lattice.
  Note: the natural BCC reconstructions are LINEAR (C0) and QUINTIC (C2). No quadratic interpolant on the
  one-cell BCC set appears. This agrees with the counting fact in (f).
- N. Chentanez, B. E. Feldman, F. Labelle, J. F. O'Brien, J. R. Shewchuk, "Liquid simulation on lattice-based
  tetrahedral meshes", SCA 2007, DOI 10.2312/SCA/SCA07/219-228. UNVERIFIED (search summary): tetrahedral mesh
  whose interior vertices lie on a BCC lattice; surface tracked by semi-Lagrangian contouring (Bargteil,
  Goktekin, O'Brien, Strain, ACM TOG 2006). Node-based, piecewise linear on tetrahedra. No stability analysis.
- N. Molino, R. Bridson, J. Teran, R. Fedkiw, "A crystalline, red green strategy for meshing highly deformable
  objects with tetrahedra", IMR 2003: BCC tetrahedral meshes from a level set. UNVERIFIED (from memory; this
  is meshing, not level-set transport).
- T. Brochu, C. Batty, R. Bridson, "Matching fluid simulation elements to surface geometry and topology",
  ACM TOG 29(4), 2010, DOI 10.1145/1778765.1778784. VERIFIED (search record): Voronoi pressure samples;
  not BCC level-set advection.

Conclusion for topic 4. No source analyses a semi-Lagrangian level set with quadratic reconstruction on BCC
data. BCC data are used with linear interpolation on tetrahedra or with box splines.

## Topic 5. Semi-Lagrangian level sets on quadtree/octree meshes

### 5.1 Min and Gibou 2007 (read in full: UCLA CAM report 06-22)
C. Min, F. Gibou, "A second order accurate level set method on non-graded adaptive Cartesian grids",
J. Comput. Phys. 225, 300-321, 2007, DOI 10.1016/j.jcp.2006.11.034.
https://ww3.math.ucla.edu/camreport/cam06-22.pdf   All points below VERIFIED.
- Storage: phi at the NODES (cell vertices) of a non-graded quadtree/octree. The finest cells lie on the
  interface (refine when the cell diagonal times the Lipschitz constant exceeds the distance to the interface).
- T-junction (hanging) node: a node with a missing neighbour in one Cartesian direction. A ghost value is
  defined by linear interpolation along the coarse edge, then the second-order error term is cancelled with
  the transverse second derivative at the node. Result: a third-order ghost value. It "only uses the node
  values of the cells adjacent to v0". After that, a T-junction node is treated "in a same fashion as a regular
  node".
- Interpolation at the departure point. Two candidates are named. (i) "Hermite quadratic interpolation on the
  parent cell": the parent of a leaf owns 3 x 3 (3 x 3 x 3) nodes; this is exactly the tensor Q2 interpolant
  of result (f). It is rejected: "this interpolation procedure is sensitive to nearby discontinuities, e.g.
  near kinks". (ii) Adopted: bilinear/trilinear interpolation in the leaf cell plus a correction
  - phi_xx x(1-x)/2 - phi_yy y(1-y)/2, with phi_xx = the second derivative of MINIMUM ABSOLUTE VALUE among the
  cell vertices. "the choice of the smallest in absolute value enhances the numerical stability of the
  interpolation". This is a nonlinear, limited quadratic. No transverse (xy) term.
- Stability claim: semi-Lagrangian methods "are unconditionally stable thus avoiding standard CFL condition".
  This is stated, not analysed. No von Neumann or matrix analysis. No statement on amplification at a
  coarse-fine transition. Tests run at dt = dx_smallest and larger.
- Reported accuracy: second order in the maximum norm near the interface; the interface is always in the
  finest cells, so the transition is away from the zero level.
- Why nodes: the paper does not argue against cell centres. The companion paper (Min, Gibou, "A second order
  accurate projection method for the incompressible Navier-Stokes equations on fully adaptive grids",
  J. Comput. Phys. 219, 2006; CAM report 05-65, https://ww3.math.ucla.edu/camreport/cam05-65.pdf, VERIFIED
  abstract) says: "We sample all the variables at the grid nodes, producing a scheme that is straightforward
  to implement." The practical reasons visible in the text: (1) node data give a continuous multilinear
  interpolant in every leaf with no neighbour search; (2) a T-junction needs at most one ghost node in 2D
  (two in 3D) and the ghost value is built from the adjacent cells only; (3) no 2:1 balance is needed.

### 5.2 Strain 1999 and 2000
J. Strain, "Tree methods for moving interfaces", J. Comput. Phys. 151, 616-648, 1999, DOI
10.1006/jcph.1999.6205; "Semi-Lagrangian methods for level set equations", J. Comput. Phys. 151, 498-533,
1999 (pages UNVERIFIED). VERIFIED (search record of the abstracts, https://math.berkeley.edu/~strain/Publications/slt.pdf):
level set + frequent redistancing + semi-Lagrangian steps "explicit yet unconditionally stable"; quadtree
meshes; CIR predictor with trapezoidal corrector. Min-Gibou cite Strain for piecewise bilinear interpolation
as "the most natural choice" on quadtrees. Stability rests on the monotone (linear) CIR interpolation.

### 5.3 Losasso, Gibou, Fedkiw 2004; Losasso, Fedkiw, Osher 2006
F. Losasso, F. Gibou, R. Fedkiw, "Simulating water and smoke with an octree data structure", ACM TOG 23(3),
457-462, 2004, DOI 10.1145/1015706.1015745. VERIFIED (search record): unrestricted octree, semi-Lagrangian
advection of velocity, density and level set. UNVERIFIED (from memory): phi is stored at the nodes and
interpolated trilinearly; pressure at cell centres.
F. Losasso, R. Fedkiw, S. Osher, "Spatially adaptive techniques for level set methods and incompressible flow",
Comput. Fluids 35, 995-1010, 2006. UNVERIFIED (from memory, not fetched).

### 5.4 Mirzadeh, Guittet, Burstedde, Gibou 2016 (read in full text)
M. Mirzadeh, A. Guittet, C. Burstedde, F. Gibou, "Parallel level-set methods on adaptive tree-based grids",
J. Comput. Phys. 322, 345-364, 2016, DOI 10.1016/j.jcp.2016.06.017.
https://ins.uni-bonn.de/media/public/publication-media/MirzadehGuittetBursteddeEtAl16.pdf   VERIFIED:
- p4est forest of octrees, node-based, NON-graded: "we can entirely skip the 2:1 balance function".
- SL step uses "the stabilized second-order interpolation" of Min-Gibou for phi and multilinear for velocity.
- Runs at CFL 1, 10 and 100; adaptive time step with CFL = 5 in the application. "Successive application ...
  especially for large values of the CFL number, eventually degrades the signed distance property"; so
  reinitialisation every few steps is part of the method.
- No stability analysis; "unconditional stability" is asserted.

### 5.5 Guittet, Theillard, Gibou 2015; Theillard and co-workers
A. Guittet, M. Theillard, F. Gibou, "A stable projection method for the incompressible Navier-Stokes equations
on arbitrary geometries and adaptive Quad/Octrees", J. Comput. Phys. 292, 215-238, 2015, DOI
10.1016/j.jcp.2015.03.024. UNVERIFIED (search summary): MAC layout on non-graded trees; the node-based
projection of Min-Gibou "cannot be guaranteed to be stable"; SL advection. The level set stays at the nodes.
"Stable nodal projection method on octree grids", arXiv:2306.09957. UNVERIFIED (search summary).

### 5.6 Olshanskii, Terekhov, Vassilevski 2013 (read in full text) -- the one negative report
M. A. Olshanskii, K. M. Terekhov, Y. V. Vassilevski, "An octree-based solver for the incompressible
Navier-Stokes equations with enhanced stability and low dissipation", Comput. Fluids 84, 231-246, 2013,
DOI 10.1016/j.compfluid.2013.04.027. https://www.math.uh.edu/~molshan/ftp/pub/OctreeMAC.pdf   VERIFIED:
- Staggered (face-centred) unknowns on GRADED (2:1) octrees.
- "a straightforward extension is prone to produce spurious oscillatory velocity modes on the fine-to-coarse
  grids interfaces. A local linear low-pass filter is shown to reduce much of the bad influence of the
  interface modes".
- SL with LINEAR interpolation: stable but over-diffusive (first order; steady solution where a vortex street
  is expected). SL with QUADRATIC interpolation of face-centred data: second order in L2, degraded in the
  maximum norm, needed half the time step to converge, and "was found unstable" for the flow around a square
  cylinder. "an extension of higher order semi-Lagrangian methods to staggered octree grids is not available
  in the literature" (2013).
- Their remedy: drop SL, use an implicit upwind finite-difference advection.
This is the closest published evidence for question (g): a NON-node-based (face-centred) quadratic SL scheme
on a 2:1 octree was observed unstable; no analysis was given.

### 5.7 Cell-centred SL level sets on quadtrees
- Z. Wang, Z. J. Wang, "Multi-phase flow computation with semi-Lagrangian level set method on adaptive
  Cartesian grids", AIAA 2005-1390. https://cfd.ku.edu/papers/aiaa-2005-1390.pdf   VERIFIED (full text):
  cell-centred FV flow solver with least-squares reconstruction; but for the level set "phi at both vertices
  and cell centers are computed", from a distance tree (exact redistancing every step) and a piecewise LINEAR
  interpolant on a triangulation of vertices + cell centres (Strain's contouring). So the SL step never
  re-interpolates an advected cell-centre field with a quadratic. No stability analysis.
- "A parallel cell-centered adaptive level set framework for efficient simulation of two-phase flows with
  subcycling and non-subcycling", J. Comput. Phys. 448, 110740, 2022. UNVERIFIED (search summary): Eulerian
  block-structured AMR, not semi-Lagrangian.
- No paper was found that reports a cell-centre-only quadratic SL level set at a 2:1 transition as stable, and
  none that analyses it.

## Topic 6. Gerris, Basilisk, and stability theory of coarse-fine interfaces

### 6.1 Gerris and Basilisk
- S. Popinet, "Gerris: a tree-based adaptive solver for the incompressible Euler equations in complex
  geometries", J. Comput. Phys. 190, 572-600, 2003, DOI 10.1016/S0021-9991(03)00298-5. Existence and the
  use of the Bell-Colella-Glaz second-order upwind (Godunov) advection VERIFIED (search records). Graded (2:1)
  quadtree: stated by Min-Gibou 2007 ("Graded Cartesian grids are also used ... see for example the work of
  Popinet"). VERIFIED there.
- Basilisk source, read in this session. VERIFIED:
  https://basilisk.fr/src/grid/multigrid-common.h : default for a scalar is `s.prolongation = refine_bilinear`
  and `s.restriction = restriction_average`. `bilinear()` in 1D is (3*coarse + coarse_neighbour)/4; the 2D and
  3D weights are over 16 and 64. This is linear interpolation of cell-CENTRE values. `refine_linear` (coarse
  value + limited gradient, conservative) is available and not the default.
  https://basilisk.fr/src/grid/tree-common.h : halo (ghost) cells at a resolution boundary are filled level by
  level by the prolongation operator; coarse cells under fine ones get the restriction; face values on the
  coarse side are the average of the fine face values; 2:1 grading is enforced on refinement and coarsening.
- So the Gerris/Basilisk answer to a 2:1 transition is: every cell sees a locally UNIFORM Cartesian stencil at
  its own level. Missing fine neighbours are ghost cells filled by (bi)linear prolongation. Missing coarse
  neighbours are restricted averages. The advection scheme itself (BCG, dimensionally unsplit upwind) is not
  changed. Courant limit: the explicit one, on the finest cell.
- S. Popinet, "An accurate adaptive solver for surface-tension-driven interfacial flows", J. Comput. Phys. 228,
  5838-5866, 2009; "A quadtree-adaptive multigrid solver for the Serre-Green-Naghdi equations", J. Comput.
  Phys. 302, 336-358, 2015. UNVERIFIED (from memory, not fetched). I found no von Neumann or matrix stability
  analysis of the coarse-fine interface in the Gerris/Basilisk papers or pages that I read.

### 6.2 Theory of interfaces
- L. N. Trefethen, "Stability of finite-difference models containing two boundaries or interfaces",
  Math. Comp. 45 (172), 279-300, 1985. https://people.maths.ox.ac.uk/~trefethen/publication/PDF/1985_25.pdf
  VERIFIED (abstract): stability is governed by the propagation and reflection of parasitic waves (group
  velocity). "(1) Dissipativity does not ensure stability when three or more formulas are concatenated at a
  boundary or internal interface. (2) Algebraic 'GKS instabilities' can be converted by a second boundary to
  exponential instabilities only when an infinite numerical reflection coefficient is present. (3)
  'GKS-stability' ... can be established ... by showing that the numerical reflection coefficient matrices
  have norm less than one." Introduction: a boundary that can emit a wave with group velocity into the
  interior in the absence of incident waves is unstable.
- M. J. Berger, "Stability of interfaces with mesh refinement", Math. Comp. 45 (172), 301-318, 1985,
  DOI 10.1090/S0025-5718-1985-0804925-4. UNVERIFIED (secondary summaries of the abstract; publisher page 403):
  refinement in space and time; Lax-Wendroff is stable with all the interface conditions considered; leapfrog
  with interpolation interface conditions is stable when the grids overlap.
  This matters for us: the scheme of result (a) is Lax-Wendroff in 1D.
- G. Browning, H.-O. Kreiss, J. Oliger, "Mesh refinement", Math. Comp. 27, 29-39, 1973; M. Ciment, "Stable
  difference schemes with uneven mesh spacings", Math. Comp. 25, 219-227, 1971 (pages of Ciment UNVERIFIED).
  Existence VERIFIED (citation records). Content UNVERIFIED.
- B. Gustafsson, H.-O. Kreiss, A. Sundstroem, Math. Comp. 26, 649-686, 1972 (GKS theory). UNVERIFIED
  (from memory).
- D. F. Martin, P. Colella, D. Graves, J. Comput. Phys. 227, 2008: UNVERIFIED (search summary) report of weak
  instabilities at coarse-fine interfaces with Crank-Nicolson.
- Bradley 2024 Sect. 2.2.1 (VERIFIED, Topic 2): SL with cubic symmetric Lagrange interpolation on a NONUNIFORM
  1D grid has eigenvalues above 1 + 1e-3 for some shifts. A 2:1 transition is a nonuniform grid. This is the
  strongest verified warning for (g).

Conclusion for topic 6. Theory exists for 1D finite-difference interfaces (GKS, Trefethen, Berger). No source
analyses a semi-Lagrangian or least-squares scheme at a 2:1 quadtree/octree transition.

## Topic 7. Inflow boundaries

- Mirzadeh et al. 2016 and Min-Gibou 2007 (VERIFIED, read): they do not treat a physical inflow boundary for
  the level set. Their tests keep the interface inside the domain.
- Elarif, Mehrenberger, Navoret 2026 (VERIFIED, read): the advection problem is "supplemented ... with Dirichlet
  conditions at the inflow boundary". That is the standard position: an inflow boundary needs a datum.
- Semi-Lagrangian practice when the departure point leaves the domain. UNVERIFIED (search summaries):
  (i) restrict the time step so that it does not happen (arXiv:2606.21191, ESBGK model); (ii) stop the
  trajectory at the boundary and take the boundary datum at the crossing time (truncation; arXiv:2301.04923
  treats the part of the trajectory outside the domain with an explicit rule); (iii) for Hamilton-Jacobi SL
  schemes, "interpolate the solution outside the domain or ... set it to a given value" (arXiv:1310.6121).
- Level-set specific: a wildfire level-set paper (arXiv:0712.3965; UNVERIFIED, search summary) extrapolates
  one layer beyond the boundary but does not let the extrapolation decrease below the boundary value, because
  "a decrease ... at boundary node would be amplified by the extrapolation".
- GKS results on extrapolation. UNVERIFIED (search summaries of Goldberg and Tadmor's papers,
  https://www.math.umd.edu/~tadmor/pub/boundary-conditions/): extrapolation at an OUTFLOW boundary keeps a
  dissipative two-level scheme stable; inflow closures are more restricted. Trefethen 1985 (VERIFIED
  introduction): a boundary closure that can radiate a wave into the interior without an incident wave is
  unstable.
- S. Tan, C.-W. Shu, inverse Lax-Wendroff procedure for inflow boundaries; F. Vilar, C.-W. Shu, ESAIM M2AN 49,
  2015 (GKS analysis of extrapolation at outflow + inverse Lax-Wendroff at inflow). UNVERIFIED (search summary).

Conclusion for topic 7. That extrapolation without a datum at an inflow boundary is unstable is KNOWN
qualitatively (well-posedness needs an inflow datum; GKS/Trefethen wave argument). The formula
max|G| = 1 + 2 Co for the one-ring least-squares fit was NOT FOUND.

## Verdicts on our results (a)-(g)

(a) Anchored quadratic least-squares fit through the one-ring cell centres, max|G| = 1 for Co <= 1 on a
    Cartesian lattice; Lax-Wendroff in 1D; weights independent of h.
    KNOWN IN A WEAKER/DIFFERENT FORM.
    - 1D: KNOWN. Strang 1962 through Besse-Mehrenberger 2008, Corollary 1 and Remark 3 (VERIFIED): centred
      even-degree Lagrange interpolation has sup|theta| = 1 exactly on [-1, 1]. Degree 2 is Lax-Wendroff.
    - Multi-D with the TENSOR-PRODUCT Lagrange interpolant on a uniform grid: known (product of 1D symbols;
      Bradley 2024 reports that Ferretti and co-workers prove the uniform-grid test problem "in any dimension"
      for symmetric odd-degree interpolants; VERIFIED as a statement of Bradley).
    - The least-squares P2 fit through 8 (2D) / 26 (3D) cell centres is not the tensor interpolant. Its symbol
      was NOT FOUND in any source. One search hit mentions a 2D Lax-Wendroff-type scheme built from "a least
      square piecewise polynomial approximation" on 5-, 9- and 13-point stencils (cited inside
      arXiv:2210.05352 / 2405.19844, Coulombel and co-author, on the 2D Lax-Wendroff scheme; UNVERIFIED, the
      original reference was not identified). Check this one.

(b) Linear fit = FTCS (unstable); cubic singular on the one-ring.
    KNOWN IN A DIFFERENT FORM. FTCS instability is textbook. That Lagrange-type reconstruction which is not
    centred on the foot is unstable is stated by Besse-Mehrenberger 2008 (VERIFIED). The specific statements
    for the least-squares fit were NOT FOUND; they are elementary.

(c) Boundary cell without datum: max|G| = 1 + 2 Co.
    KNOWN IN A WEAKER/DIFFERENT FORM. Qualitative: an inflow boundary needs a datum; a closure that radiates
    into the interior is unstable (Trefethen 1985, VERIFIED; GKS; Goldberg-Tadmor, UNVERIFIED). Extrapolation
    in SL methods is known to lose the stability that interpolation gives (lecture notes found by search state
    that stability holds "provided that interpolation is used ... instead of extrapolation"; UNVERIFIED
    source: twister.caps.ou.edu CFD notes). The formula 1 + 2 Co: NOT FOUND.

(d) Honeycomb lattice: stable up to Co = sqrt(3)/2.
    NOT FOUND. The only SL work on the hexagonal lattice that I found (Mehrenberger, Mendoza, Prouveur,
    Sonnendruecker 2016) is node-based on equilateral triangles with box splines or Hermite elements, and no
    Courant limit of this kind is reported (UNVERIFIED, search summary).

(e) Weak amplification on perturbed meshes.
    KNOWN IN A DIFFERENT FORM. Bradley 2024 (VERIFIED): node-based tensor-product SL on randomly perturbed
    periodic element grids has lambda_max - 1 > 0, small, "sparse" in the shift, scaling like delta^4 for his
    bases. Ferretti 2013 Sect. 5 as quoted by Bradley: symmetric Lagrange interpolation on unevenly spaced
    nodes is "in general ... unstable". Bradley speculates that every high-order interpolation SL method has
    this property. For a cell-centred least-squares fit on perturbed polygonal meshes: NOT FOUND.
    Related, semi-discrete only: Haider-Croisille-Courbet 2009 (least-squares MUSCL, eigenvalue stability).

(f) Q2 point-value SL scheme (vertices, face centres, edge midpoints, cell centre): growth factor 1; the 3D
    subset without edge midpoints amplifies; centre + vertices cannot carry a full quadratic in 3D.
    - Full Q2 set stable on the uniform test problem: KNOWN IN A DIFFERENT FORM. Baptista 1987 and
      Oliveira-Baptista 1995 (compact interpolants unstable only from degree three), Ferretti-Perrone 2007
      (P2 interpolation; content UNVERIFIED), Bradley 2024 (VERIFIED). These are computed eigenvalue bounds.
      A proof for the quadratic element was announced in Ferretti-Mehrenberger 2020 and was not found.
      Min-Gibou 2007 (VERIFIED) write this exact interpolant ("Hermite quadratic interpolation" on the 3 x 3
      (x 3) nodes of the parent cell) and reject it for the level set because it is "sensitive to nearby
      discontinuities, e.g. near kinks". That is the practical objection to (f), not linear instability.
      The same 27-point set with a cell average in place of the centre value is the 3D Cartesian Active Flux
      (Barsukow et al., VERIFIED).
    - Subset without edge midpoints amplifies: NOT FOUND.
    - Centre + vertices (BCC data) cannot carry a quadratic: NOT FOUND as a stability statement; consistent
      with the BCC reconstruction literature, which uses linear or quintic box splines (Entezari et al.).

(g) 2:1 transitions and hybrid own-cell-point + neighbour-centre schemes.
    NOT FOUND (no analysis). Evidence that exists:
    - Negative: Olshanskii-Terekhov-Vassilevski 2013 (VERIFIED): quadratic SL interpolation of non-nodal
      (face-centred) data on a graded octree "was found unstable"; spurious modes sit on the fine-to-coarse
      interfaces. Bradley 2024 (VERIFIED): symmetric high-order interpolation on a nonuniform grid has
      eigenvalues above one for some shifts.
    - Positive, for comparison: Berger 1985 (UNVERIFIED summary): Lax-Wendroff in 1D stays stable at a
      refinement interface for all interface conditions that she considers.
    - Hybrid schemes (own-cell point values + neighbour cell centres): NOT FOUND. The closest are Active Flux
      (average + boundary points; AMR version by Calhoun-Chudzik-Helzel, no interface stability analysis in
      the abstract) and Wang-Wang 2005 (phi at vertices AND cell centres, but linear on a triangulation and
      re-distanced exactly every step).

## The 10 sources to read in full, ranked

1. Bradley 2024, J. Comput. Phys., DOI 10.1016/j.jcp.2024.113034 (open: https://www.osti.gov/servlets/purl/2345299).
   The method to copy: maximum eigenvalue of the space-time operator over all shifts; the reduction to small
   eigenproblems per wave number; the perturbed-uniform-grid metric; "sparse instability". Sections 2.2, 3.3,
   3.5, 3.7, 3.8, Appendix A.
2. Ferretti and Perrone 2007, in Applied and Industrial Mathematics in Italy II, World Scientific, 339-350.
   Stability of SL under finite-element (P2) interpolation. Not online; get it from the library.
3. Besse and Mehrenberger 2008, Math. Comp. 77, 93-123 (https://www.i2m.univ-amu.fr/perso/mehrenberg.m/BeMe.pdf).
   Section 5.1 and Appendix A: closed form of the symbol theta_m, the exact stability zones, Remarks 3 and 4.
   Together with Strang 1962.
4. Min and Gibou 2007, J. Comput. Phys. 225, 300-321 (https://ww3.math.ucla.edu/camreport/cam06-22.pdf).
   Sections 4, 5, 7: T-junction ghost values, the limited quadratic interpolation.
5. Olshanskii, Terekhov, Vassilevski 2013, Comput. Fluids 84, 231-246
   (https://www.math.uh.edu/~molshan/ftp/pub/OctreeMAC.pdf). The interface modes on graded octrees and the
   unstable quadratic SL scheme.
6. Ferretti and Mehrenberger 2020, Math. Comp. 89, 1783-1805, with Ferretti 2010 (J. Comput. Math. 28) and
   Ferretti 2013 (Numer. Math. 124, especially Sect. 5 on unevenly spaced nodes).
7. Trefethen 1985, Math. Comp. 45, 279-300, and Berger 1985, Math. Comp. 45, 301-318. The wave/reflection
   criterion for interfaces; Lax-Wendroff at a refinement interface.
8. Chudzik, Helzel, Kerkmann 2021, Appl. Math. Comput. 393, 125501, and Calhoun, Chudzik, Helzel
   (arXiv:2205.07572). Linear stability of the Cartesian Active Flux and its AMR version (ghost-cell
   treatment of point values at coarse-fine patch boundaries).
9. Haider, Croisille, Courbet 2009, Numer. Math. 113, 555-600. Eigenvalue stability of least-squares
   reconstruction on unstructured grids (semi-discrete).
10. Elarif, Mehrenberger, Navoret 2026, arXiv:2607.00619, and its references Carre et al. 2009 (least-squares SL
    on cell-centred meshes) and Besse-Sonnendruecker 2003. The current state on triangles; Hermite data as the
    stable route. Also Oliveira and Baptista 1995 for the original compact-interpolator result.

## Lessons for (g)

What the octree level-set community does (VERIFIED in Min-Gibou 2007 and Mirzadeh et al. 2016):
1. It stores phi at the cell VERTICES. Every leaf then owns a complete multilinear interpolant. No neighbour
   search and no 2:1 balance are needed.
2. At a hanging (T-junction) node it does not create an unknown with a special update. It builds a ghost value
   from the nodes of the adjacent cells: linear interpolation along the coarse edge, corrected to third order
   with the transverse second difference. All difference formulas then look regular.
3. At the departure point it uses multilinear interpolation in the leaf plus a quadratic correction with the
   second derivative of minimum modulus over the cell vertices. It does not use the full Q2 (parent-cell)
   interpolant, because that interpolant oscillates near kinks of the distance function. The stencil never
   extends beyond the leaf and the second differences at its vertices.
4. It keeps the finest cells on the interface. The 2:1 (or larger) transition lies in the far field, where
   phi is re-distanced every few steps. Reinitialisation is part of the method and removes far-field errors.
5. It asserts unconditional stability. It gives no von Neumann or matrix analysis, also not at transitions.

What the cell-centred tree codes do (VERIFIED in the Basilisk source):
6. Gerris/Basilisk enforce 2:1 grading, fill ghost cells of the fine level by bilinear prolongation and the
   coarse cells under fine ones by averaging. Each cell uses a uniform stencil at its own level. The advection
   scheme is Eulerian (BCG) with the explicit Courant limit of the finest cell. No stability analysis of the
   interface was found.

Are cell-centre-only schemes at 2:1 transitions reported stable?
7. For a QUADRATIC semi-Lagrangian reconstruction: no positive report was found. The one report that I found is
   negative (Olshanskii et al. 2013, face-centred data, graded octree, quadratic interpolation unstable, modes
   on the fine-to-coarse interface). For LINEAR interpolation the scheme is stable and first order
   (CIR, monotone weights).
8. Theory says that a nonuniform node set breaks the symmetry on which the |G| <= 1 result rests
   (Ferretti 2013 as quoted by Bradley; Bradley's 5-point example with lambda_max > 1 + 1e-3). Expect
   lambda_max - 1 > 0 at a 2:1 transition for a least-squares quadratic fit, small and dependent on the shift.
   Whether it is "sparse" (Bradley's term) or present for all Courant numbers must be computed.
9. Concrete steps that follow from the literature:
   a. Do the Bloch analysis on a periodic patch that contains a 2:1 transition, as Bradley does for his
      perturbed grid: build the space-time operator for the patch, scan the shift over (0, 1] times the fine
      cell, report lambda_max - 1 against the Courant number of the fine and of the coarse cell.
   b. Test the Basilisk option first: ghost values by prolongation, so that every cell fits on a uniform
      one-ring at its own level. Then the fit weights are the Cartesian ones and only the ghost values carry
      the transition. Bilinear ghost values are second order; Min-Gibou show how to make a ghost value third
      order with a transverse second difference.
   c. Test Trefethen's criterion: no mode with group velocity away from the transition may exist without an
      incident wave; the reflection coefficients at the transition must be below one.
   d. For the hybrid scheme, the Active Flux AMR paper (Calhoun, Chudzik, Helzel) is the nearest template for
      point values on a coarse-fine boundary; read its ghost-cell interpolation of the point values.
   e. Keep the transition away from the interface and state that in the method; that is what every published
      octree level-set method does, and it is why none of them needed the analysis.

## Limits of this survey
- Publisher pages of AMS, ESAIM and BIRS returned HTTP 403. Those abstracts were read through Crossref, author
  PDFs or citing papers, and are marked accordingly.
- Not obtained: Ferretti-Perrone 2007, Strang 1962 (full text), Oliveira-Baptista 1995, Carre et al. 2009,
  Popinet 2003/2009 (full text), Losasso et al. 2004/2006 (full text), Eymann's thesis.
- "NOT FOUND" means not found in this session's searches. It is not a proof that no publication exists.
