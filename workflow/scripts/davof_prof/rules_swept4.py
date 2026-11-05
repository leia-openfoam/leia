# The categories of hpan.py for the fork at a11ebc1 (2026-11-05, profile round 3b: the class
# statistics, the donor fast path, the face update from the volume clip, the pre-test's coplanar
# side): rules_swept3.py's attribution with the line ranges of the current sources. A sample belongs
# to the NEAREST ancestor context matching a rule.
S = r"isoAdvectionSwept\.C$"
P = r"davofSweptPieces\.H$"
RULES = [
    ("pieces: clip by cell planes + PLIC (sweptPieceVolume)", "line", P, (323, 505)),
    ("pieces: cap closure + volume (sweptPieceVolume)", "line", P, (506, 607)),
    ("pieces: sweptPieceVolume other", "fn", r"davof::sweptPieceVolume"),
    ("face update: pre-image clip (sweptPieceArea)", "fn", r"davof::sweptPieceArea"),
    ("face update: donor plane (donorWet)", "fn", r"isoAdvection::donorWet"),
    ("swept cells: the donor fast path test", "line", S, (753, 795)),
    ("swept cells: ring search (vertex ring, bounding boxes)", "line", S, (796, 885)),
    ("swept cells: cell planes, state, wisp rule + the pre-test (first pass)", "line", S, (886, 1048)),
    ("swept cells: the state shortcut", "line", S, (1049, 1109)),
    ("swept cells: the class statistics", "line", S, (1266, 1280)),
    ("swept cells: pre-image fan, areas", "line", S, (659, 752)),
    ("swept cells: other (sweptFaceCellsFlux)", "fn", r"isoAdvection::sweptFaceCellsFlux"),
    ("one-donor clip (sweptFaceFlux, fallback)", "fn", r"isoAdvection::sweptFaceFlux\b"),
    ("trace: cell search (ring, octree)", "line", S, (297, 331)),
    ("trace: velocity interpolate (cellPointWeight, findTetrahedron)", "fn", r"interpolationCellPoint<[^>]*>::interpolate\(|cellPointWeight|polyMeshTetDecomposition"),
    ("trace: velocity interpolate", "line", S, (332, 345)),
    ("trace: midpoint trace + pre-image", "line", S, (346, 430)),
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
