# The categories of hpan.py for the fork at 407f969 (2026-11-04, the reference setting): the
# attribution of rules_swept.py / rules_swept2.py with the line ranges of the current sources
# (isoAdvectionSwept.C after the trace cache 557fab7 and the check 09a6c07; davofSweptPieces.H
# unchanged since 1e279d4). A sample belongs to the NEAREST ancestor context matching a rule.
S = r"isoAdvectionSwept\.C$"
C = r"isoAdvectionCurved\.C$"
P = r"davofSweptPieces\.H$"
RULES = [
    ("pieces: clip by cell planes + PLIC (sweptPieceVolume)", "line", P, (337, 485)),
    ("pieces: cap closure + volume (sweptPieceVolume)", "line", P, (486, 587)),
    ("pieces: sweptPieceVolume other", "fn", r"davof::sweptPieceVolume"),
    ("face update: pre-image clip (sweptPieceArea)", "fn", r"davof::sweptPieceArea"),
    ("face update: donor plane (donorWet)", "fn", r"isoAdvection::donorWet"),
    ("swept cells: ring search (vertex ring, bounding boxes)", "line", S, (749, 842)),
    ("swept cells: cell planes, state, wisp rule", "line", S, (843, 961)),
    ("swept cells: the pre-test", "line", S, (962, 1027)),
    ("swept cells: the state shortcut", "line", S, (1028, 1069)),
    ("swept cells: the pre-test check", "line", S, (1070, 1122)),
    ("swept cells: pre-image fan, areas", "line", S, (657, 748)),
    ("swept cells: other (sweptFaceCellsFlux)", "fn", r"isoAdvection::sweptFaceCellsFlux"),
    ("one-donor clip (sweptFaceFlux, fallback)", "fn", r"isoAdvection::sweptFaceFlux\b"),
    ("trace: cell search (ring, octree)", "line", S, (297, 331)),
    ("trace: velocity interpolate (cellPointWeight, findTetrahedron)", "fn", r"interpolationCellPoint<[^>]*>::interpolate\(|cellPointWeight|polyMeshTetDecomposition"),
    ("trace: velocity interpolate", "line", S, (332, 345)),
    ("trace: midpoint trace + pre-image", "line", S, (346, 419)),
    ("trace: triangle list + other (sweptTriangles)", "fn", r"isoAdvection::sweptTriangles"),
    ("curvedFlux: point velocity (volPointInterpolation of U, with its patch exchange)", "fn", r"interpolationCellPoint<[^>]*>::interpolationCellPoint\(|volPointInterpolation::"),
    ("curvedFlux: syncProcPatches (all-to-all of the face fractions)", "fn", r"isoAdvection::syncProcPatches"),
    ("curvedFlux: wisp surfaces, candidate band, exchanges, face loop", "fn", r"isoAdvection::curvedFlux"),
    ("remote exchange (sweptExchange, zoneDistribute, pack/unpack)", "fn", r"isoAdvection::sweptExchange|isoAdvection::packSweptCell|davof::unpackRemoteCell|davof::packCellSurface|zoneDistribute|zoneCPCStencil|mapDistribute"),
    ("surface build (buildCellSurfaces: fit, anchor, volume rule)", "fn", r"isoAdvection::buildCellSurfaces|isoAdvection::anchorTensorCentroids|davof::kromerBothe|davof::kb[A-Z]|davof::solveOffset|cellSurface::"),
    ("reconstruction (davof::reconstruct, ring normal, plane positioning)", "fn", r"reconstruction::davof::|davofNormalModels::|reconstruction::plicRDF::|reconstructionSchemes|cutCellPLIC|cutFacePLIC|cutCellIso|cutFaceIso|cutCell::|cutFace::|multiDimPolyFitter|leastSquareGrad|zoneDistribute::setUpCommforZone|reconstructedDistanceFunction|davofCcs|ccsPositioning"),
    ("bounding (limitFluxes, boundCellContent, second-layer sweep)", "fn", r"isoAdvection::limitFluxes|isoAdvection::boundCellContent|isoAdvection::applyBruteForceBounding|isoAdvection::sweepSecondRing|isoAdvection::boundFlux|isoAdvection::checkIfOnProcPatch"),
    ("isoAdvector: timeIntegratedFlux (serial PLIC flux)", "fn", r"isoAdvection::timeIntegratedFlux|isoAdvection::timeIntegratedFaceFlux|isoAdvection::setDownwindFaces|isoAdvection::netFlux|isoAdvection::publishWetFaces|isoAdvection::defaultAreaFractions|isoAdvection::initAreaFractions"),
    ("isoAdvector: advect other (alpha update, counters, I/O of the step)", "fn", r"isoAdvection::advect|isoAdvection::reportCurved|isoAdvection::sweptConsistencyProbe|isoAdvection::classifyPreImage"),
    ("MPI library", "mod", r"^lib(mpi|open-pal|open-rte|pmix|uc[xpts]|ucm|hwloc|event|psm|fabric|ibverbs|rdmacm|mlx|hcoll)"),
    ("Pstream wrappers", "fn", r"Foam::(UPstream|UIPstream|UOPstream|Pstream|IPstream|OPstream|PstreamBuffers)::"),
    ("flow setup (setFlow function object)", "fn", r"functionObjects::setFlow"),
    ("Courant number, other function objects", "fn", r"functionObjects::|CourantNo"),
    ("I/O and start-up (mesh, fields read/write)", "fn", r"Foam::polyMesh::|Foam::fvMesh::|Foam::primitiveMesh::|Foam::regIOobject::|Foam::IOobject::|Foam::Time::|Foam::OFstream|Foam::IFstream|Foam::ISstream|Foam::OSstream|Foam::fileOperation|Foam::Istream|Foam::Ostream|Foam::ITstream|Foam::dictionary::|Foam::entry::|Foam::primitiveEntry|Foam::Field<[^>]*>::writeEntry|Foam::GeometricField<[^>]*>::writeData|::readField|::writeData\("),
    ("main: the solver's own loop", "fn", r"^main$"),
]
ALLOC = r"^(malloc|free|cfree|_int_malloc|_int_free|realloc|calloc|__libc_malloc|__libc_free|_mid_memalign|unlink_chunk|malloc_consolidate|sysmalloc|operator new|operator delete|__GI___libc_malloc|tcache_)"
MPI_SITE_EXCLUDE = ("MPI library", "Pstream wrappers")
