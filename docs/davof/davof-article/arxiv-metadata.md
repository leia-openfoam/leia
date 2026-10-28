# arXiv metadata for the DAVOF preprint (2026-10-28; Tomislav submits)

**Title**
Dual Area-Volume-of-Fluid (DAVOF): co-advecting volume and face-area fractions for a closed-form, cell-local interface normal and position on unstructured meshes

**Author**
Tomislav Marić, Mathematical Modeling and Analysis, Department of Mathematics, Technical University of Darmstadt

**Primary category** physics.flu-dyn. **Cross-lists** math.NA, physics.comp-ph.

**Comments** 45 pages, 9 figures, 18 tables. Code: https://github.com/tmaric/TwoPhaseFlow (branch davof) and https://github.com/leia-openfoam/leia (branch development).

**License** CC BY 4.0 (choose the arXiv non-exclusive licence instead if the journal version may go to a publisher that objects to CC BY).

**Abstract for the metadata field** (arXiv limit 1920 characters; the PDF keeps the full abstract of ms.tex). Draft, shortened from ms.tex without new claims:

Geometric Volume-of-Fluid (VOF) on unstructured polyhedral meshes reconstructs the interface orientation from the volume fractions of a stencil of neighbours; on irregular cells and in the presence of wisps it degrades, and the curvature built on it drives the parasitic currents of capillary flows. We transport, consistently to the volume fraction of every cell, the wetted fraction of every face area. From the area and volume fractions, closed-form equations recover the interface normal and position of a cell, for the first time in VOF: the normal from one sum over the cell's own faces, the divergence theorem on the liquid sub-volume, and the plane's position from a second, without neighbours or iteration. Both identities are exact for a plane on hexahedral, perturbed and polyhedral cells. On a sphere and an ellipsoid the normal has the order of the face fractions, the position is second order for every state and mesh, and a quadric through one ring of recovered normals gives a converging curvature where the reference models show none. Transported with the plane clipping, the ring normal, the same identity on the union of the cell's ring, stabilises the state and carries a sphere over 3.3 diameters bounded and conservative, first order in the face fractions as the clipping law predicts, at plicRDF's position and shape errors. A preliminary transport implementation with one paraboloid per cell, its frame from the normal, its curvature tensor from the ring and its position from the volume fraction, carries the sphere bounded and conservative over the resolution ladder with face fractions of order 2.8 per period with the exact normal and 2.1 with the ring quadric's normal over N = 32 to 128. The loop has the order of the normal it is fed; a transported normal of second order is the research goal. This document is the basis on which the further development of the method is documented.

**Data and code (decision of 2026-10-27)** No figshare deposit for version 1: the tables stay in the preprint and as CSV files in the leia repository; both repositories are cited at the release tag `davof-v1.0` (leia: branch development with the manuscript folder; TwoPhaseFlow fork: branch davof at the commit of the tables, with the branches davof-fluxint and davof-tiltdiag named for the appendix results).

**Resolved (2026-10-27, late)**
- leia organisation name: `leia-openfoam` is correct (https://github.com/leia-openfoam/leia answers 200, `leia-methods/leia` 404; leia's own remote is leia-openfoam). The manuscript, the bibliography and the tag URLs use it.
- Puckett 1991: title, proceedings (Fourth International Symposium on Computational Fluid Dynamics, Davis, CA, 1991) and pages 933–938 agree with the citation used in the literature (e.g. the reference list of arXiv:1703.02202); the editor field (H. Dwyer) is the symposium's organiser and is kept. No DOI exists for these proceedings.

- Review comment 16 (heading "Plane: exactness"): withdrawn by Tomislav on 2026-10-27; the heading stays.

**Nothing is open. The archive of tag `davof-v1.0` is the submission.**
2. If the manuscript changes after the tag, move the leia tag to the final commit (`git tag -f davof-v1.0 && git push -f origin davof-v1.0`) or tag `davof-v1.0.1`; the fork tag stays.
3. Run `bash make_arxiv_package.sh`; upload `_build/davof-preprint-arxiv.tar.gz`; paste the metadata above.


**2026-10-28 (proposal session)** The transport section was restructured (short results, eight new tables, the ladder figure, the implementation record as Appendix B), the limitations are stated as research goals, and the Elsevier declaration of generative AI was added after the acknowledgments; all of it is committed on development (last commit 3195655e). The tag davof-v1.0 (24e13fa8) predates these commits: before the upload, move it to the final commit as item 2 above says, or tag davof-v1.0.1 and change the two bibliography entries. The package in _build/arxiv was built from the final manuscript.
