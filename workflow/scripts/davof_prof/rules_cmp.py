# Profile round 3 (2026-11-05): the categories of hpan.py for the comparison of the two advection loops
# function by function, applied alike to the swept-face advection (fork 407f969, rules_swept3.py's
# line ranges) and to isoAdvector + plicRDF (the same binary's libVoF): the reconstruction split into
# the plane positioning (surfaceIteratorPLIC / cutCellPLIC, shared code) and the normal (plicRDF's RDF
# and its iterations; davof's ring normal), the flux into its sub-parts (ours: the trace, the pre-test,
# the clips; theirs: cutFaceAdvect's plane sweep per face), the point velocity and the stencil
# machinery named alike for both. A sample belongs to the NEAREST ancestor context matching a rule.
S = r"isoAdvectionSwept\.C$"
P = r"davofSweptPieces\.H$"
RULES = [
    ("flux ours: pieces, the clip by the cell planes and the PLIC plane (sweptPieceVolume)", "line", P, (337, 485)),
    ("flux ours: pieces, the cap closure and the volume (sweptPieceVolume)", "line", P, (486, 587)),
    ("flux ours: pieces, sweptPieceVolume other", "fn", r"davof::sweptPieceVolume"),
    ("flux ours: the face update, the pre-image clip (sweptPieceArea)", "fn", r"davof::sweptPieceArea"),
    ("flux ours: the face update, the donor plane (donorWet)", "fn", r"isoAdvection::donorWet"),
    ("flux ours: the ring search (vertex ring, bounding boxes)", "line", S, (749, 842)),
    ("flux ours: the cell planes, state, wisp rule", "line", S, (843, 961)),
    ("flux ours: the pre-test", "line", S, (962, 1027)),
    ("flux ours: the shortcut and the check (off)", "line", S, (1028, 1122)),
    ("flux ours: the pre-image fan, areas", "line", S, (657, 748)),
    ("flux ours: sweptFaceCellsFlux other", "fn", r"isoAdvection::sweptFaceCellsFlux"),
    ("flux ours: the one-donor clip (sweptFaceFlux, fallback)", "fn", r"isoAdvection::sweptFaceFlux\b"),
    ("flux ours: the trace, the cell search (ring, octree)", "line", S, (297, 331)),
    ("the velocity at a point (interpolationCellPoint::interpolate, cellPointWeight, findTetrahedron)", "fn", r"interpolationCellPoint<[^>]*>::interpolate\(|cellPointWeight|polyMeshTetDecomposition"),
    ("flux ours: the trace, the velocity lookup's own lines", "line", S, (332, 345)),
    ("flux ours: the trace, the midpoint trace and the pre-image", "line", S, (346, 419)),
    ("flux ours: the trace, the triangle list (sweptTriangles other)", "fn", r"isoAdvection::sweptTriangles"),
    ("the point velocity (volPointInterpolation of U, the whole mesh)", "fn", r"interpolationCellPoint<[^>]*>::interpolationCellPoint\(|volPointInterpolation::"),
    ("flux ours: syncProcPatches (the all-to-all of the face fractions)", "fn", r"isoAdvection::syncProcPatches"),
    ("flux ours: curvedFlux (the band, the wisp surfaces, the exchanges, the face loop)", "fn", r"isoAdvection::curvedFlux"),
    ("flux ours: the remote exchange (sweptExchange, pack/unpack)", "fn", r"isoAdvection::sweptExchange|isoAdvection::packSweptCell|davof::unpackRemoteCell|davof::packCellSurface"),
    ("the stencil machinery (zoneDistribute, zoneCPCStencil, mapDistribute, cellToCellStencil)", "fn", r"zoneDistribute|zoneCPCStencil|mapDistribute|cellToCellStencil"),
    ("surface build ours (buildCellSurfaces: fit, anchor, volume rule)", "fn", r"isoAdvection::buildCellSurfaces|isoAdvection::anchorTensorCentroids|davof::kromerBothe|davof::kb[A-Z]|davof::solveOffset|cellSurface::"),
    ("reconstruction: the plane positioning (surfaceIteratorPLIC, cutCellPLIC; ours: ccs on the same cut)", "fn", r"surfaceIteratorPLIC::|cutCellPLIC::|cutFacePLIC::|cutCellIso|cutFaceIso|cutCell::|cutFace::|davofCcs|ccsPositioning"),
    ("reconstruction: the normal, plicRDF (the RDF, its interpolation, the gradient, the residual, the initial normals)", "fn", r"reconstructedDistanceFunction::|plicRDF::centreAndNormalBC|plicRDF::gradSurf|leastSquareGrad|plicRDF::calcResidual|plicRDF::setInitNormals|plicRDF::interpolateNormal|multiDimPolyFitter"),
    ("reconstruction: the normal, davof (the ring normal model)", "fn", r"davofNormalModels::"),
    ("reconstruction: the rest (plicRDF::reconstruct / davof::reconstruct own lines, field algebra)", "fn", r"reconstruction::davof::|reconstruction::plicRDF::|reconstructionSchemes"),
    ("bounding (limitFluxes, boundCellContent, second-layer sweep)", "fn", r"isoAdvection::limitFluxes|isoAdvection::boundCellContent|isoAdvection::applyBruteForceBounding|isoAdvection::sweepSecondRing|isoAdvection::boundFlux|isoAdvection::checkIfOnProcPatch"),
    ("flux theirs: the plane sweep per face (cutFaceAdvect::timeIntegratedFaceFlux, timeIntegratedArea)", "fn", r"cutFaceAdvect::"),
    ("flux theirs: timeIntegratedFlux other (the downwind faces, the face loop)", "fn", r"isoAdvection::timeIntegratedFlux|isoAdvection::timeIntegratedFaceFlux|isoAdvection::setDownwindFaces|isoAdvection::netFlux|isoAdvection::publishWetFaces"),
    ("isoAdvector shared: defaultAreaFractions, initAreaFractions", "fn", r"isoAdvection::defaultAreaFractions|isoAdvection::initAreaFractions|isoAdvection::setDefaultAreaFraction"),
    ("isoAdvector shared: advect other (alpha update, counters, I/O of the step)", "fn", r"isoAdvection::advect|isoAdvection::reportCurved|isoAdvection::sweptConsistencyProbe|isoAdvection::classifyPreImage"),
    ("MPI library", "mod", r"^lib(mpi|open-pal|open-rte|pmix|uc[xpts]|ucm|hwloc|event|psm|fabric|ibverbs|rdmacm|mlx|hcoll)"),
    ("Pstream wrappers", "fn", r"Foam::(UPstream|UIPstream|UOPstream|Pstream|IPstream|OPstream|PstreamBuffers)::"),
    ("shared: setFlow (the function object)", "fn", r"functionObjects::setFlow"),
    ("shared: the Courant number, other function objects", "fn", r"functionObjects::|CourantNo"),
    ("shared: I/O and start-up (mesh, fields read/write)", "fn", r"Foam::polyMesh::|Foam::fvMesh::|Foam::primitiveMesh::|Foam::regIOobject::|Foam::IOobject::|Foam::Time::|Foam::OFstream|Foam::IFstream|Foam::ISstream|Foam::OSstream|Foam::fileOperation|Foam::Istream|Foam::Ostream|Foam::ITstream|Foam::dictionary::|Foam::entry::|Foam::primitiveEntry|Foam::Field<[^>]*>::writeEntry|Foam::GeometricField<[^>]*>::writeData|::readField|::writeData\("),
    ("shared: the solver's own loop (main)", "fn", r"^main$"),
]
ALLOC = r"^(malloc|free|cfree|_int_malloc|_int_free|realloc|calloc|__libc_malloc|__libc_free|_mid_memalign|unlink_chunk|malloc_consolidate|sysmalloc|operator new|operator delete|__GI___libc_malloc|tcache_)"
MPI_SITE_EXCLUDE = ("MPI library", "Pstream wrappers")
