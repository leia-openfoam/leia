# plan-library-split-and-build-policy.md

Per-clone binaries, version stamps in every library, and the split of
`libleiaLevelSet` into a core and seven method libraries. Working document for a
coding agent -- **v0.1**, 2026-09-22, written from the question series of that
day. Style: ASD-STE100. Nothing in this document has been executed.

Scope: how leia is built and installed on the laptop and on Lichtenberg, and how
the library is cut. Out of scope: any change to a numerical method, any study
result, the pending SDPLS science items (the R arm's noise amplification, the
suspended parallel ladders), and the two pending briefs (`interfaceDefects`,
`sdplsCombined`).

---

## 0. Decisions taken (user, 2026-09-22)

| # | question | decision |
|---|---|---|
| 1 | goal of this plan | finish the cluster hygiene: build policy, study safety, and the library split |
| 2 | why two sessions collided | not between methods: every method is runtime-selected in one library; the collision was between two code VERSIONS installed into one shared folder |
| 3 | cut lines | core; `sdplsSource`; `semiLagrangian`; `velocityExtension`; `redistancer`; `surfaceTensionForce` (with the `fvOptions` capillary term); `volumeCorrection`; all separate and combinable at run time |
| 4 | the advection hub | its own thin library `libleiaAdvection`, linking `sdplsSource`, `velocityExtension` and `semiLagrangian`, so the Eulerian solver keeps pairing a velocity extension with a source term |
| 5 | build policy | binaries inside the clone: `WM_PROJECT_USER_DIR` is the clone root, binaries land in `<clone>/platforms/`; one committed environment file sourced after OpenFOAM's `etc/bashrc` in every launch path |
| 6 | gates for the split | the full standing set (section 3.3), compared at tolerance 0 |
| 7 | loading of runtime-selected models | each solver links every library it can select from |
| 8 | layout | in place: one `Make/` per part under `src/leiaLevelSet/`; the core keeps `src/leiaLevelSet/Make` |
| 9 | order | environment file, then version stamps, then the split; one gated commit each |
| 10 | cluster work | this session, both clones, authorized for this plan; nothing while a job of either session runs |
| 11 | study safety | `snakemake --touch` every finished study once, and add a guard that refuses a launch that would re-run a finished case |
| 12 | version stamp | `git describe --always` plus `-dirty` per library; printed in the solver banner; written to `leia.version` in the case; copied into `case_params.json` next to `gitCommit` |
| 13 | documents | CLAUDE.md and AGENTS.md build line; CLUSTER.md and SLURM.md; workflow/README.md and STATUS.md; the `docs/IMPROVEMENTS.md` ownership table (committed with the split) |
| 14 | old install folders and `$HOME/.leia_env` | rename with a dated suffix; delete after one week of runs from the new folders |
| 15 | deliverable | this document, uncommitted, for review |

---

## 1. Ground truth -- measured 2026-09-22, do not re-litigate

**Runtime selection.** `libleiaLevelSet` (88 translation units) holds eleven
runtime-selection families: `semiLagrangian` (8 tables), `sdplsSource` (8),
`advection`, `narrowBand`, `phaseIndicator`, `profile`, `redistancer`,
`surfaceTensionForce`, `velocityExtension`, `velocityModel`, `volumeCorrection`
(2 each). One solver binary carries every method; the case dictionaries select.

**Compile-time dependencies between the parts** (header includes; the only ones):

| part | files | lines | includes headers from |
|---|---|---|---|
| `semiLagrangian` | 44 | 9980 | none |
| `sdplsSource` | 32 | 5388 | none (finds `NarrowBand` and `phi` by name at run time) |
| `surfaceTensionForce` | 27 | 4177 | `phaseIndicator` (`levelSetPlaneReconstruction.H`) |
| `velocityExtension` | 16 | 3195 | `semiLagrangian` (`closestPoint.H` -> `slReconstruction.H`), `velocityModel` (`fluxCorrection.H`) |
| `redistancer` | 10 | 2015 | `phaseIndicator` (`levelSetPlaneReconstruction.H`) |
| `phaseIndicator` | 12 | 1677 | none |
| `volumeCorrection` | 4 | 1177 | `phaseIndicator` |
| `narrowBand` | 12 | 1172 | none |
| `velocityModel` | 6 | 1023 | none |
| `advection` | 6 | 630 | `velocityExtension`, `sdplsSource`, `velocityModel`, `semiLagrangian` |
| `fvOptions` | 2 | 563 | none (`semiImplicitCapillaryForce`) |
| `profile` | 6 | 482 | none |
| `schemes` | 2 | 403 | none (`levelSetBlended`) |

**Consumers.** 17 binaries link `-lleiaLevelSet` and include
`src/leiaLevelSet/lnInclude`: 5 solvers, 10 test applications, 2 utilities. The
other leia libraries (`levelSetImplicitSurfaces`, `finiteVolume`,
`functionObjects`) do not link it. By header use: `leiaLevelSetFoam` uses the
advection hub, `redistancer`, `volumeCorrection`; `leiaLevelSetTwoPhaseFoam` uses
`sdplsSource`, `surfaceTensionForce`, `redistancer`, `volumeCorrection`;
`leiaRedistancedLevelSetFoam` uses `redistancer`, `sdplsSource`;
`leiaSemiLagrangeLevelSetFoam` uses `slAdvection`;
`leiaSemiLagrangianLevelSetTwoPhaseFoam` uses `slAdvection`,
`velocityExtension`, `sdplsSource`, `surfaceTensionForce`, `redistancer`,
`volumeCorrection`.

**Launch paths that source OpenFOAM.** `Allwmake` assumes a sourced shell.
`run-studies.sbatch` sources `etc/bashrc` (line 30) and, since b16fc41, a
clone-local `.leia_env` if present (line 31, to be replaced). The committed
`profiles/slurm/config.yaml` line 132 passes an `env_preamble` that sources
`etc/bashrc` in every job; the study configs carry the same preamble for the
local profile. The Snakefile builds every job's shell in ONE helper, `sh()`
(`workflow/Snakefile` line 119): `set +eu` + PREAMBLE + `set -e` + the commands.
That helper is the single place where the environment file reaches every job,
local or SLURM, whatever profile or config supplied the preamble.

**The incident.** Both sessions install into `$HOME/OpenFOAM/tm83tomy-v2512` by
default. On 2026-09-09 the curvature clone compiled there; the SDPLS clone, at
b3aa65e, then held binaries about 200 commits ahead of its source. Nothing in
any log showed it: the `gitCommit` column is stamped by `materialize.py` from
the CLONE's git when a case is built, and the binaries carry no stamp.
Measured on 2026-09-22 (job 54823248): a job submitted from a pinned driver
shell resolves `leiaLevelSetFoam` to the default folder, because the job's
preamble re-sources `etc/bashrc`. `.gitignore` already ignores `platforms/`
(lines 42, 46) and `lnInclude` (line 34). `git describe --always --dirty` on
the laptop gives `shared-method-config-2026-09-01-132-gd613e27-dirty`.

**Study safety.** After the 2026-09-22 mtime sweep on `/work/scratch`, a dry
run with the Makefile's flags (`--rerun-triggers mtime`) would re-run 6 of 42
cases of `sdplsConv2Dvortex` and 2 of 14 of `sdplsExpSource2Dvortex`; the
`solve` rule deletes the case CSV before it runs. Without that flag every case
would re-run.

---

## 2. Work packages, in order

Each work package is one commit on `development`, gated before the push, then
pulled into both cluster clones and rebuilt there (section 4). Every commit is
staged by explicit path. The thread-C entries in `git status` are not touched.

### WP1 -- the environment file: binaries inside the clone

**Files.** New `etc/leia-env.sh` (new directory `etc/`). Edits: `Allwmake`,
`Allwclean`, `workflow/Snakefile` (the `sh()` helper), `run-studies.sbatch`,
`run-sink-decay.sbatch` if it launches a leia binary, and every other script
that a grep for `etc/bashrc` finds launching a leia binary. NOT edited:
`profiles/slurm/config.yaml`, the 34 study configs (their preambles stay; the
Snakefile hook follows them).

**Content of `etc/leia-env.sh`** (sourced, never executed):

```bash
# Source AFTER OpenFOAM's etc/bashrc. Installs and finds THIS clone's binaries in
# <clone>/platforms/<WM_OPTIONS>/{bin,lib}. Idempotent.
if [ -z "$WM_OPTIONS" ]; then
    echo "leia-env.sh: WM_OPTIONS is empty; source OpenFOAM etc/bashrc first" >&2
    return 1 2>/dev/null || exit 1
fi
LEIA_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export LEIA_ROOT
export WM_PROJECT_USER_DIR="$LEIA_ROOT"
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
# drop every OpenFOAM user dir that etc/bashrc or an earlier source put first
_strip() { echo "$1" | tr ':' '\n' | grep -v '/OpenFOAM/[^/]*-v2512/platforms/' | grep -v "^$LEIA_ROOT/platforms/" | paste -sd:; }
export PATH="$FOAM_USER_APPBIN:$(_strip "$PATH")"
export LD_LIBRARY_PATH="$FOAM_USER_LIBBIN:$(_strip "$LD_LIBRARY_PATH")"
unset -f _strip
```

**Hooks.**

1. `Allwmake`, first lines: `. "$(cd "$(dirname "$0")" && pwd)/etc/leia-env.sh" || exit 1`.
   Same in `Allwclean`.
2. `workflow/Snakefile`, `sh()`: after PREAMBLE insert
   `. "{REPO}/etc/leia-env.sh"` (REPO is the Python variable at line 16).
3. `run-studies.sbatch`: replace line 31 (the `.leia_env` hook of b16fc41) by
   `. "$PWD/etc/leia-env.sh" || exit 1`.
4. The user's interactive shells: documented in the CLAUDE.md build line (WP6):
   `source $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && . ./etc/leia-env.sh`.

**Gate WP1** (laptop, then cluster):

1. `source etc/bashrc; . ./etc/leia-env.sh; ./Allwmake`. Check: binaries appear
   under `<clone>/platforms/linux64GccDPInt32Opt/{bin,lib}`; `which
   leiaLevelSetFoam` prints the clone path; `ls platforms/*/bin | wc -l` >= 17.
2. Bit identity of the relocation: re-run `studies/sdpls1Dstretch/1Dstretch_00001`
   (R arm, serial, ran on this laptop) as in the 2026-09-22 gate; compare both
   CSVs with `compare_metrics_csv.py --tol 0 --skip ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`.
   Expected: PASS (same code, new location).
3. Through the workflow (exercises the `sh()` hook): `snakemake
   --workflow-profile profiles/local --configfile config/sdpls1Dstretch.yaml
   --config studies_dir=$PWD/tmp/gate-wp1 --nolock`. Check: every case log's
   `Exec` line names `<clone>/platforms`; the 12 CSVs equal the published ones
   at tolerance 0; `git status --short docs/` unchanged (the report rule
   regenerates identical tables).
4. Cluster, after the pull into clone A (section 4): `sbatch
   --export=ALL,STUDY=sdpls1Dstretch run-studies.sbatch` (that study does not
   exist there, so nothing is overwritten). Check: `grep -m1 '^Exec'
   studies/sdpls1Dstretch/1Dstretch_00000/log.leiaLevelSetFoam` names
   `/work/scratch/tm83tomy/leia/platforms/`; the 12 CSVs equal the laptop's at
   tolerance 1e-10 or better (cross-machine, so not byte-exact by rule).

**Rollback WP1.** `git revert <hash>`; rebuild; the old folders still exist.

### WP2 -- version stamps

**Files.** New `etc/leia-stamp.sh`; new `src/leiaLevelSet/leiaVersionRegistry.{H,C}`
(core); generated, git-ignored `src/leiaLevelSet/leiaStamp_leiaLevelSet.C` (one
file now, one per library after WP3; `.gitignore` gets `leiaStamp_*.C`); two
lines in each of the 5 solvers after `createMesh.H`; `workflow/Snakefile`
(`solve` rule) and `workflow/scripts/aggregate.py` (one column).

**Mechanism.**

- `etc/leia-stamp.sh <library-dir> <library-name>` computes
  `git -C "$LEIA_ROOT" describe --always` and appends `-dirty` when
  `git status --porcelain -- <the same code paths materialize.py uses>` is not
  empty (so a docs edit does not mark the code dirty, exactly like `gitCommit`).
  It writes `leiaStamp_<name>.C` only when the content changed, so wmake
  recompiles exactly when the stamp changes. `Allwmake` calls it before each
  `wmake`.
- The generated file: one static object
  `Foam::leia::versionRegistrar reg_("libleia<Name>", "<describe>")`. The
  registry (core) collects (library, stamp) pairs at load time.
- Each solver, after `createMesh.H`: `leia::reportVersions(Info)` (one banner
  line per loaded leia library) and `leia::writeVersions(runTime)` (file
  `<case>/leia.version`, one line per library). Test applications: the same
  banner call, optional.
- Workflow: the `solve` rule appends the lines of `leia.version` into
  `case_params.json` under `libStamps`; `aggregate.py` writes them as the
  column `libStamps` in `<study>_database.csv` and `_errors.csv`.

**Gate WP2.** Banner shows `leia library libleiaLevelSet : <describe>`;
`leia.version` exists after a run; `case_params.json` and the database carry
`libStamps`; a pulled-but-not-rebuilt clone shows a stamp that differs from
`gitCommit` (test it once by hand: touch a source, do not rebuild, run, read).
Bit identity: the WP1 gate cases at tolerance 0 (the stamp code adds no
arithmetic).

**Rollback WP2.** `git revert <hash>`; the banner lines vanish; nothing else.

### WP3 -- the split

**Libraries** (in place; each part gets `Make/files` and `Make/options`; the
core keeps `src/leiaLevelSet/Make`):

| library | parts | links (leia) |
|---|---|---|
| `libleiaCore` | `profile`, `narrowBand`, `phaseIndicator`, `velocityModel`, `schemes`, the root file, `leiaVersionRegistry` | `levelSetImplicitSurfaces` |
| `libleiaSdplsSource` | `sdplsSource` | `Core` |
| `libleiaSemiLagrangian` | `semiLagrangian` | `Core`, `levelSetImplicitSurfaces` |
| `libleiaVelocityExtension` | `velocityExtension` | `Core`, `SemiLagrangian` |
| `libleiaRedistancer` | `redistancer` | `Core` |
| `libleiaVolumeCorrection` | `volumeCorrection` | `Core` |
| `libleiaSurfaceTension` | `surfaceTensionForce`, `fvOptions` | `Core` |
| `libleiaAdvection` | `advection` | `Core`, `SdplsSource`, `SemiLagrangian`, `VelocityExtension` |

`libleiaLevelSet` is retired: `Allwmake` deletes a stale
`$FOAM_USER_LIBBIN/libleiaLevelSet.so` when it finds one, so no binary can load
the monolith next to the new libraries.

**Link matrix** (every binary also links `-lleiaCore`; wmake links executables
with `--no-as-needed`, so a linked library is loaded and its models register):

| binary | links |
|---|---|
| `leiaLevelSetFoam` | `Advection` (and, through it, `SdplsSource`, `SemiLagrangian`, `VelocityExtension`), `Redistancer`, `VolumeCorrection` |
| `leiaLevelSetTwoPhaseFoam` | `SdplsSource`, `SurfaceTension`, `Redistancer`, `VolumeCorrection` |
| `leiaRedistancedLevelSetFoam` | `Redistancer`, `SdplsSource` |
| `leiaSemiLagrangeLevelSetFoam` | `SemiLagrangian` |
| `leiaSemiLagrangianLevelSetTwoPhaseFoam` | `SemiLagrangian`, `VelocityExtension`, `SdplsSource`, `SurfaceTension`, `Redistancer`, `VolumeCorrection` |
| test applications | per their header use (section 1); `leiaTestSdplsSource` links `SdplsSource` |
| `leiaSetFields`, `leiaPerturbMesh` | `Core` |

**Compile-time discipline.** wmake's `lnInclude` for `src/leiaLevelSet/` links
every header below it, so a solver's include line
`-I../../../src/leiaLevelSet/lnInclude` stays valid and hides no boundary. The
boundaries are enforced in two other ways: (1) every library's `Make/options`
links with `-Wl,--no-undefined`, so a symbol used from a library that is not in
its `LIB_LIBS` fails the library link instead of surfacing at run time; (2) the
dependency script of 2026-09-22 (header includes per part) runs as a gate and
must print exactly the table of section 1 and nothing more.

**Build order in `Allwmake`.** `levelSetImplicitSurfaces`, `Core`,
`SdplsSource`, `SemiLagrangian`, `VelocityExtension`, `Redistancer`,
`VolumeCorrection`, `SurfaceTension`, `Advection`, `finiteVolume`,
`functionObjects`, then `wmake all applications`. `Allwclean` mirrors it.
Every library gets its own stamp file (WP2) and appears in the banner.

**Gate WP3 -- the full standing set, at tolerance 0.** Baselines are produced
with the WP2 binaries BEFORE the split, on the same machine, and preserved
(`studies/<study>` copied to `studies/<study>.pre-split-<date>`); every gate
compares case CSVs pairwise with `compare_metrics_csv.py --tol 0 --skip
ELAPSED_CPU_TIME,ELAPSED_CLOCK_TIME`, never aggregates alone.

| rung | what | where |
|---|---|---|
| unit | `leiaTestSdplsSource` (89 assertions) and `cases/sdplsSourceUnit/Allrun.sh` | laptop |
| exact 1D | `config/sdpls1Dstretch.yaml`, 12 serial cases | laptop |
| advection regression, hex 2D | `config/advConv2Dvortex.yaml`, `config/advConv2Dtranslation.yaml`, three resolutions | laptop |
| advection regression, hex 3D | `config/advConv3DshearHex.yaml` | laptop or cluster by size |
| advection regression, poly 3D | `config/advConv3DshearPoly.yaml` on the shared mesh (`workflow/scripts/advect_bound_arm.sh`) | cluster |
| coupled, SL two-phase | `config/stationaryDroplet3DbitIdentity.yaml` | cluster |
| coupled, Eulerian two-phase | the N = 32 arm of `config/sdplsPsiBudgetDroplet2D.yaml` | laptop |
| seam | `config/seamConsistency3Dserial.yaml` and `config/seamConsistency3Dpar4.yaml` | cluster |
| four ranks, every solver | one np = 4 case per solver (`sdplsExpSource2Dvortex` N = 32 for `leiaLevelSetFoam`; the bit-identity droplet for each two-phase solver; a 2D vortex arm for `leiaSemiLagrangeLevelSetFoam`) | laptop |
| loading | `ldd` of all 17 binaries shows the expected leia libraries and no `libleiaLevelSet` | laptop and cluster |

PASS = every pair identical. Any difference stops the push; a refactor that
changes one number is not a refactor.

**Rollback WP3.** `git revert <hash>`; `./Allwclean && ./Allwmake`; the
pre-split baselines stay preserved.

### WP4 -- study safety: touch once, then guard

1. On each cluster clone, for every study directory with a matching config:
   `snakemake --workflow-profile profiles/local --configfile config/<study>.yaml
   --touch --nolock` (marks outputs up to date, runs nothing). Check: the dry
   run with the Makefile's flags then lists no `solve` job. Record the list of
   touched studies in STATUS.
2. Guard `workflow/scripts/guard_finished_cases.py <config.yaml> --profile <dir>`:
   runs `snakemake -n` with the Makefile's flags, collects the `solve` jobs'
   `idx` wildcards, and for each checks `studies/<study>/<case>/<solver>.csv`
   with `foam_log_state.sh`. If any would-be solve targets a COMPLETED case, it
   prints them and exits 2. Override only by `LEIA_ALLOW_RERUN=1`, whose use is
   documented as "preserve the study directory first".
3. `Makefile`: `studies-one` and the study loops call the guard before
   `$(SNAKE)`. This is a shared file; the change is additive and both clones
   have no local edit of it.

### WP5 -- docs (with each work package, same commit)

- CLAUDE.md and AGENTS.md: the build line becomes `source
  $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc && . ./etc/leia-env.sh && ./Allwmake`,
  plus one sentence: binaries live in `<clone>/platforms/`; two clones never
  share binaries. `diff CLAUDE.md AGENTS.md` prints nothing.
- CLUSTER.md and SLURM.md: the `$HOME/.leia_env` overlay recipe and the
  private profile edit are retired in favour of `etc/leia-env.sh`; the measured
  incidents (2026-09-01, 2026-09-09, 2026-09-22) stay as history.
- workflow/README.md: the library table, the link matrix, and the guard.
- STATUS.md: a section 10 with the migration state of both clones, the gate
  verdicts, and the retired folders with their deletion date.
- docs/IMPROVEMENTS.md: the ownership table names libraries, not directories;
  committed with the split (it was untracked).

### WP6 -- retirement of the old folders

After the first real run from each clone shows the new stamps and the clone
path in its `Exec` line: rename `$HOME/OpenFOAM/tm83tomy-v2512`,
`curvature-v2512`, `sdpls-v2512` and `$HOME/.leia_env` to
`<name>.retired-<date>`. Delete them one week later, after both sessions have
run a study. Record both dates in STATUS section 10.

---

## 3. Cluster migration (this session, both clones, authorized)

Preconditions before every pull or build: `squeue -u tm83tomy -h -o '%.10i %.30j %Z'`
shows no job with either clone as WorkDir; `sacct` since the previous step shows
none either. Every remote step is a script piped through `ssh host 'bash -s' <
step.sh` that writes to `/work/scratch/tm83tomy/leia-sync-<date>/log/<step>.log`;
a second call reads the log. Job ids go to `.my_jobs`.

Per work package, per clone:

1. Clone B first, once, before the WP1 pull: `git stash push -m "private
   profile overlay, superseded by etc/leia-env.sh" -- profiles/slurm/config.yaml`.
   The stash is recoverable. `$HOME/.leia_env` is renamed at WP6, not now.
2. `git pull --ff-only origin development`; verify HEAD.
3. Clone A only, at WP1: delete the clone-local `.leia_env` and its line in
   `.git/info/exclude` (superseded).
4. `module purge; module load gcc/11.5.0-z7mc openmpi/4.1.8-6xzv; source
   $HOME/OpenFOAM/OpenFOAM-v2512/etc/bashrc; . ./etc/leia-env.sh; ./Allwmake`.
   Check: binaries under `<clone>/platforms/`; `strings` finds the stamp.
5. The cluster rungs of the work package's gate (section 2).
6. At WP4: the `--touch` pass over every finished study of that clone.

Rollback of a clone: `git reset --keep <previous hash>` (never `--hard`), then
rebuild in the same shell setup, then `snakemake -n` before any launch.

---

## 4. Do not

- Do not change `profiles/slurm/config.yaml`: the Snakefile hook makes it
  unnecessary, and the file is the one clone B edits locally.
- Do not `git add -A`, `-a`, `.`; do not stash, checkout or reset anything on
  the laptop; do not commit the 14 thread-C modified files.
- Do not build into `$HOME/OpenFOAM/*-v2512` again; do not source
  `$HOME/.leia_env` in any shell of this plan.
- Do not run `make decks` or `docs/build-decks.sh` (they rewrite thread-C data
  mirrors); build one deck with `export_html.py` when needed.
- Do not touch `src/leiaLevelSet/semiLagrangian/**`,
  `surfaceTensionForce/**`, `velocityExtension/**` beyond moving files between
  `Make/files` lists: no code inside a method changes in this plan.
- Do not delete a retired folder before its week has passed.
- Do not launch or re-run any study on the cluster except the gate rungs.

---

## 5. Time budget

| work package | wall time |
|---|---|
| WP1 environment file, laptop gate, both clones | half a day |
| WP2 stamps, gate, both clones | half a day |
| WP3 split, baselines, full gate set, both clones | one and a half days |
| WP4 touch and guard | half a day |
| WP5 docs | with each package |
| WP6 retirement | ten minutes, twice |

---

## Appendix A -- commands used to measure section 1

```bash
# dependency map (parts, files, lines, cross-part includes)
python3 - <<'PY'
import os,re,collections
os.chdir("src/leiaLevelSet"); hdr={}
for dp,dn,fn in os.walk("."):
    if "/Make" in dp or "lnInclude" in dp: continue
    for f in fn:
        if f.endswith(".H"): hdr[f]=dp.split("/")[1] if dp.count("/")>=1 else "(root)"
for part in sorted({v for v in hdr.values()}):
    out=collections.defaultdict(set)
    for dp,dn,fn in os.walk(part if part!="(root)" else "."):
        if "/Make" in dp or "lnInclude" in dp: continue
        for f in fn:
            if not f.endswith((".C",".H")): continue
            for m in re.finditer(r'#include\s+"([^"]+)"', open(os.path.join(dp,f),errors="ignore").read()):
                d=hdr.get(os.path.basename(m.group(1)))
                if d and d!=part: out[d].add(f)
    print(part, dict((k,sorted(v)) for k,v in out.items()) or "-")
PY
# consumers
grep -l 'lleiaLevelSet' applications/*/*/Make/options | wc -l
# runtime-selection families
grep -rlE 'declareRunTimeSelectionTable|defineRunTimeSelectionTable' src/leiaLevelSet | cut -d/ -f3 | sort | uniq -c
```
