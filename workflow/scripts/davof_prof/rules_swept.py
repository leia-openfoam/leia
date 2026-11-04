# The attribution of the own samples (hpan.py): a sample belongs to the category
# of the NEAREST ancestor context matching a rule (a function node by name, a
# source line by file and line range, or a load module). Inlined helpers
# (clipByHalfSpace, triNegativePolygon, planeOf) are therefore attributed to
# the call site's range; libc's malloc/free to the caller's category (the
# allocation share is listed separately, ALLOC). Lines refer to fork 720cfc6.
S = r"isoAdvectionSwept\.C$"
C = r"isoAdvectionCurved\.C$"
P = r"davofSweptPieces\.H$"
RULES = [
    # the flux of a face through the ring cells' pieces (sweptFaceCellsFlux)
    ("pieces: clip by cell planes + PLIC (sweptPieceVolume)", "line", P, (296, 394)),
    ("pieces: cap closure + volume (sweptPieceVolume)", "line", P, (395, 496)),
    ("pieces: sweptPieceVolume other", "fn", r"davof::sweptPieceVolume"),
    ("face update: pre-image clip (sweptPieceArea)", "fn", r"davof::sweptPieceArea"),
    ("face update: donor plane (donorWet)", "line", S, (947, 974)),
    ("swept cells: ring search (vertex ring, bounding boxes)", "line", S, (694, 731)),
    ("swept cells: cell planes, state, wisp rule", "line", S, (774, 846)),
    ("swept cells: pre-image fan, areas", "line", S, (671, 692)),
    ("swept cells: other (sweptFaceCellsFlux)", "fn", r"isoAdvection::sweptFaceCellsFlux"),
    ("one-donor clip (sweptFaceFlux, fallback)", "fn", r"isoAdvection::sweptFaceFlux\b"),
    # the swept triangles and the backward trace
    ("trace: cell search (ring, octree)", "line", S, (297, 328)),
    ("trace: velocity interpolate (cellPointWeight, findTetrahedron)", "fn", r"interpolationCellPoint<[^>]*>::interpolate\(|cellPointWeight|polyMeshTetDecomposition"),
    ("trace: velocity interpolate", "line", S, (329, 344)),
    ("trace: midpoint trace + pre-image", "line", S, (345, 392)),
    ("trace: triangle list + other (sweptTriangles)", "fn", r"isoAdvection::sweptTriangles"),
    # the driver curvedFlux
    ("curvedFlux: point velocity (volPointInterpolation of U, with its patch exchange)", "fn", r"interpolationCellPoint<[^>]*>::interpolationCellPoint\(|volPointInterpolation::"),
    ("curvedFlux: point velocity (volPointInterpolation)", "line", C, (1176, 1204)),
    ("curvedFlux: syncProcPatches (all-to-all of the face fractions)", "fn", r"isoAdvection::syncProcPatches"),
    ("curvedFlux: processor-patch surface exchanges", "line", C, (1205, 1253)),
    ("curvedFlux: wisp surfaces (neighbour plane, anchoring)", "line", C, (1254, 1383)),
    ("curvedFlux: candidate band (interface cells + neighbours + vertex rings)", "line", C, (1384, 1453)),
    ("curvedFlux: second exchange (candidate flags)", "line", C, (1454, 1553)),
    ("curvedFlux: second ring / surfaceOf", "line", C, (1625, 1705)),
    ("curvedFlux: face loop + other", "fn", r"isoAdvection::curvedFlux"),
    # the remote exchange of the ring cells
    ("remote exchange (sweptExchange, zoneDistribute, pack/unpack)", "fn", r"isoAdvection::sweptExchange|isoAdvection::packSweptCell|davof::unpackRemoteCell|davof::packCellSurface|zoneDistribute|zoneCPCStencil|mapDistribute"),
    # the surfaces and the reconstruction
    ("surface build (buildCellSurfaces: fit, anchor, volume rule)", "fn", r"isoAdvection::buildCellSurfaces|isoAdvection::anchorTensorCentroids|davof::kromerBothe|davof::kb[A-Z]|davof::solveOffset|cellSurface::"),
    ("reconstruction (davof::reconstruct, ring normal, plane positioning)", "fn", r"reconstruction::davof::|davofNormalModels::|reconstruction::plicRDF::|reconstructionSchemes|cutCellPLIC|cutFacePLIC|cutCellIso|cutFaceIso|cutCell::|cutFace::|multiDimPolyFitter|leastSquareGrad|zoneDistribute::setUpCommforZone|reconstructedDistanceFunction"),
    # the bounding and isoAdvector's own parts
    ("bounding (limitFluxes, boundCellContent, second-layer sweep)", "fn", r"isoAdvection::limitFluxes|isoAdvection::boundCellContent|isoAdvection::applyBruteForceBounding|isoAdvection::sweepSecondRing|isoAdvection::boundFlux|isoAdvection::checkIfOnProcPatch|isoAdvection::syncProcPatches"),
    ("isoAdvector: timeIntegratedFlux (serial PLIC flux)", "fn", r"isoAdvection::timeIntegratedFlux|isoAdvection::timeIntegratedFaceFlux|isoAdvection::setDownwindFaces|isoAdvection::netFlux|isoAdvection::publishWetFaces|isoAdvection::defaultAreaFractions|isoAdvection::initAreaFractions"),
    ("isoAdvector: advect other (alpha update, counters, I/O of the step)", "fn", r"isoAdvection::advect|isoAdvection::reportCurved|isoAdvection::sweptConsistencyProbe|isoAdvection::classifyPreImage"),
    ("MPI library", "mod", r"^lib(mpi|open-pal|open-rte|pmix|uc[xpts]|ucm|hwloc|event|psm|fabric|ibverbs|rdmacm|mlx|hcoll)"),
    ("Pstream wrappers", "fn", r"Foam::(UPstream|UIPstream|UOPstream|Pstream|IPstream|OPstream|PstreamBuffers)::"),
    ("flow setup (setFlow function object)", "fn", r"functionObjects::setFlow"),
    ("Courant number, other function objects", "fn", r"functionObjects::|CourantNo"),
    ("I/O and start-up (mesh, fields read/write)", "fn", r"Foam::polyMesh::|Foam::fvMesh::|Foam::primitiveMesh::|Foam::regIOobject::|Foam::IOobject::|Foam::Time::|Foam::OFstream|Foam::IFstream|Foam::ISstream|Foam::OSstream|Foam::fileOperation|Foam::Istream|Foam::Ostream|Foam::ITstream|Foam::dictionary::|Foam::entry::|Foam::primitiveEntry|Foam::Field<[^>]*>::writeEntry|Foam::GeometricField<[^>]*>::writeData|::readField|::writeData\("),
    ("main: the solver's own loop", "fn", r"^main$"),
]
# the categories whose MPI time is broken down by call site (the MPI and Pstream rules are dropped for that view)
MPI_SITE_EXCLUDE = ("MPI library", "Pstream wrappers")
ALLOC = r"^(malloc|free|cfree|_int_malloc|_int_free|realloc|calloc|__libc_malloc|__libc_free|_mid_memalign|unlink_chunk|malloc_consolidate|sysmalloc|operator new|operator delete|__GI___libc_malloc|tcache_)"
