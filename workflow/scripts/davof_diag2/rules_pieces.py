# The attribution of the own samples (hpan.py) for the pieces path of leia, branch davof-ccs-pieces at d50f9d96 (library code of
# 112c284a): a sample belongs to the category of the NEAREST ancestor context matching a rule (a function node by name, a source line by
# file and line range, or a load module); libc's malloc/free go to the caller (the allocation share is listed separately).
# Lines refer to src/davof/geometry/davofSweptPieces.H and src/davof/advection/davofSweptAdvectionFlux.C of that commit.
H = r"davofSweptPieces\.H$"
F = r"davofSweptAdvectionFlux\.C$"
RULES = [
    # sweptPieceExactVolume (signed tets clipped by the cell planes and the PLIC plane), by source range
    ("ev: clipOne (distances, crossing, copy into out)", "line", H, (437, 476)),
    ("ev: areaOf", "line", H, (477, 485)),
    ("ev: coverage clips of the triangle", "line", H, (498, 515)),
    ("ev: tet set-up (poly array, set3)", "line", H, (516, 549)),
    ("ev: clipAll polygon loop (poly out, q = out)", "line", H, (550, 565)),
    ("ev: cap frame and centroid", "line", H, (566, 584)),
    ("ev: cap atan2", "line", H, (585, 590)),
    ("ev: cap insertion sort", "line", H, (591, 601)),
    ("ev: cap copy and orientation", "line", H, (602, 622)),
    ("ev: plane loops, dump branch", "line", H, (623, 648)),
    ("ev: volume sums (areaOf, v3)", "line", H, (649, 665)),
    ("ev: other lines of sweptPieceExactVolume", "fn", r"davof::sweptPieceExactVolume"),
    ("area clips (sweptPieceArea)", "fn", r"davof::sweptPieceArea|davof::clipWith|davof::planeDistances|davof::clipByHalfSpace"),
    ("closure path / one-plane clips (sweptPieceVolume, sweptPieceOnePlane)", "fn", r"davof::sweptPieceVolume|davof::sweptPieceOnePlane"),
    # sweptFaceCellsFlux around them
    ("flux: pre-image fan, areas, ring search, donor fast path", "line", F, (668, 835)),
    ("flux: ring state, planes, wisp rule, pre-test", "line", F, (836, 1165)),
    ("flux: piece loop other lines", "line", F, (1166, 1459)),
    ("flux: assembly (scale, donor wet, alpha_f)", "line", F, (1460, 1708)),
    ("flux: other (sweptFaceCellsFlux)", "fn", r"davofSweptAdvection::sweptFaceCellsFlux"),
    # the trace
    ("trace: velocity interpolation", "fn", r"velocityAt|interpolationCellPoint|cellPointWeight|polyMeshTetDecomposition|volPointInterpolation"),
    ("trace: sweptTriangles", "fn", r"davofSweptAdvection::sweptTriangles"),
    # curvedFlux and its parts
    ("curvedFlux: sweptExchange, remote cells, zoneDistribute", "fn", r"davofSweptAdvection::sweptExchange|davofSweptAdvection::packSweptCell|davof::unpackRemoteCell|zoneDistribute|zoneCPCStencil|mapDistribute"),
    ("curvedFlux: syncProcPatches", "fn", r"davofSweptAdvection::syncProcPatches"),
    ("curvedFlux: surface build (buildCellSurfaces, anchor, volume rule)", "fn", r"davofSweptAdvection::buildCellSurfaces|davofSweptAdvection::anchorTensorCentroids|davof::kromerBothe|davof::kb[A-Z]|davof::solveOffset|cellSurface::"),
    ("curvedFlux: faces loop and other", "fn", r"davofSweptAdvection::curvedFlux|davofSweptAdvection::onePlane|davofSweptAdvection::sweptOnePlaneFlux"),
    # the reconstruction
    ("reconstruction: davof (area normal, ring normal, CCS positioning)", "fn", r"reconstruction::davof::|davofNormalModel|ringAreaNormal|davofCcs|cutCellPLIC|cutFacePLIC|cutCell::|cutFace::|davofCutCellPLIC"),
    ("reconstruction: plicRDF / isoAlpha (OpenFOAM)", "fn", r"reconstruction::plicRDF::|reconstruction::isoAlpha::|reconstructionSchemes|cutCellIso|cutFaceIso|multiDimPolyFitter|leastSquareGrad|reconstructedDistanceFunction|zoneDistribute::setUpCommforZone|isoAlpha"),
    # bounding and the advection step
    ("limiter (limitFluxes, boundFlux, bounding)", "fn", r"davofSweptAdvection::limitFluxes|davofSweptAdvection::boundFlux|davofSweptAdvection::applyBruteForceBounding|davofSweptAdvection::boundedFaceFraction|isoAdvection::limitFluxes|isoAdvection::boundFlux|isoAdvection::applyBruteForceBounding|isoAdvection::setDownwindFaces"),
    ("isoAdvection: timeIntegratedFlux, face fractions", "fn", r"isoAdvection::timeIntegratedFlux|isoAdvection::timeIntegratedFaceFlux|isoAdvection::calcIsoFace|isoAdvection::netFlux|cutFaceAdvect"),
    ("advect other (alpha update, face fractions, reports)", "fn", r"davofSweptAdvection::advect|davofSweptAdvection::reportCurved|davofSweptAdvection::defaultAreaFractions|davofSweptAdvection::publishWetFaces|davofSweptAdvection::reinitFaceFractions|isoAdvection::advect"),
    ("MPI library", "mod", r"^lib(mpi|open-pal|open-rte|pmix|uc[xpts]|ucm|hwloc|event|psm|fabric|ibverbs|rdmacm|mlx|hcoll)"),
    ("Pstream wrappers", "fn", r"Foam::(UPstream|UIPstream|UOPstream|Pstream|IPstream|OPstream|PstreamBuffers)::"),
    ("function objects (error measures, Courant number, flow setup)", "fn", r"functionObjects::|CourantNo|davofInit::"),
    ("I/O and start-up (mesh, fields)", "fn", r"Foam::polyMesh::|Foam::fvMesh::|Foam::primitiveMesh::|Foam::regIOobject::|Foam::IOobject::|Foam::Time::|Foam::OFstream|Foam::IFstream|Foam::ISstream|Foam::OSstream|Foam::fileOperation|Foam::Istream|Foam::Ostream|Foam::ITstream|Foam::dictionary::|Foam::entry::|Foam::primitiveEntry|writeData|readField"),
    ("main: the solver's own loop", "fn", r"^main$"),
]
MPI_SITE_EXCLUDE = ("MPI library", "Pstream wrappers")
ALLOC = r"^(malloc|free|cfree|_int_malloc|_int_free|realloc|calloc|__libc_malloc|__libc_free|_mid_memalign|unlink_chunk|malloc_consolidate|sysmalloc|operator new|operator delete|__GI___libc_malloc|tcache_)"
