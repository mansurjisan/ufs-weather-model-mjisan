# One-way NWM River Forcing into ADCIRC through the UFS Coastal Coupler

Status: implementation in progress on branch `feature/nwm_river_forcing`.
Nothing here has been compiled or run yet -- see "Open items".

This mirrors the already-validated NWM -> SCHISM route documented in
`NWM_SCHISM_COUPLING_DOCUMENTATION.md`. Read that document first: the CDEPS
`dnwm` data component, the `fd_ufs.yaml` field, the driver registration of the
`NWM` component and the direct-connector run sequence are all shared and are
NOT changed by the ADCIRC work.

---

## 1. Objective

Drive ADCIRC's normal-flux open-boundary node strings (IBTYPE 22/32) from
National Water Model discharge delivered through the NUOPC/CDEPS coupling stack,
replacing the legacy `fort.20` file as the data source. One-way only: ADCIRC
exports nothing back to NWM.

Route (identical to the validated SCHISM route):

```
DATM (CDEPS)  --MED (CMEPS)-->  ADCIRC (OCN)
NWM (CDEPS dnwm) --------------> ADCIRC (OCN)     direct connector,
                                                 :remapMethod=redist
```

`advance_to_next_time=.true.` on the dnwm side, paired with the NWM-last run
sequence (see section 5, runSeq finding): the export consumed at the start of
window k is stamped at that same time, and the discharge is held constant over
each coupling window (zero-order hold; QNIN1=QNIN2 in the cap). Note this
differs from the SCHISM direct-connector setup, which uses `.false.` because
the SCHISM cap's CheckImportRiver specialization exempts the river field from
the strict NUOPC timestamp check.

---

## 2. Architecture decisions

These are locked design decisions (D1-D8 of the implementation spec).

**D1 -- dnwm is unchanged.** `CDEPS-interface/CDEPS/dnwm/nwm_comp_nuopc.F90` is
reused as-is. It exports `river_volume_flux` (m3 s-1), `datamode=copyall`, one
value per mesh entity. The `redist` connector has `mapfcopy` semantics: no area
weighting, and the sequence index of an entity in the dnwm mesh/stream is the
receiver's global entity id. Everything below follows from that.

**D2 -- ID contract for ADCIRC.** One dnwm mesh entity per ADCIRC GLOBAL NODE, in
natural global node order `1..NP` (`NP` from `fort.14`). The stream carries
per-node volumetric shares `Q_i` in m3/s: a full-size sparse array, zero
everywhere except at river-boundary nodes. This is the direct analogue of the
SCHISM contract, which carried 45,167 elements with 6 nonzero values.

**D3 -- River entry point and the width conversion.** The coupled discharge lands
on IBTYPE 22/32 normal-flux boundary node strings, the same arrays legacy ADCIRC
fills from `fort.20`: `QNIN1`/`QNIN2` (declared `src/global.F:324`, read
`src/cstart.F:220-248`, hot-started `src/hstart.F:924-980`, re-read per step in
`src/timestep.F:877-894`). ADCIRC's `QNIN` is a flux PER UNIT WIDTH, and it is
indexed by the flux-boundary list index `J = 1..NVEL`, not by node number
(`NBV(J)` gives the node, `LBCODEI(J)` the boundary code -- `src/boundaries.F`).

The reach-to-node SPLIT of NWM discharge is computed OFFLINE (Phase 0, static
fractions). At runtime the cap only does the unit conversion:

```
QNIN_i = Q_i / w_i
w_i    = effective width of node i [m]
       = half the sum of the lengths of the boundary edges adjacent to node i
         along its own string
```

`w_i` is computed ONCE at init from the model's own mesh geometry, using
spherical/CPP distances consistent with how ADCIRC treats the mesh coordinates.
Conservation holds by construction for any width formula:
`sum_i QNIN_i * w_i = sum_i Q_i`.

**D4 -- Zero-order hold.** Each coupling window the cap sets
`QNIN1 = QNIN2 = current value` and brackets the window with
`QTIME1`/`QTIME2`, so ADCIRC's internal time interpolation
(`QN2 = RampExtFlux*(QNIN1 + QTRATIO*(QNIN2-QNIN1))`) becomes a no-op. All
`fort.20` file I/O is skipped while the coupled river is active.

**D5 -- Compile flag.** The same name as the SCHISM work: `USE_NUOPC_RIVER`,
plumbed to the ADCIRC targets when the superproject is configured with
`-DUSE_NUOPC_RIVER=ON` and `ADCIRC=ON`.

**D6 -- Runtime activation is connection-driven.** Flag ON but
`river_volume_flux` NOT connected must behave EXACTLY like legacy ADCIRC:
`fort.20` path intact, plain `coastal_ike_shinnecock_atm2adc` unaffected. The
gate is `NUOPC_IsConnected`, as in the SCHISM cap. A constant-discharge stub
mode is available for cheap runtime checks, opt-in via the OCN attributes
`river_stub` (explicit `true`) and `river_stub_q` (total m3 s-1) -- the same
pair the SCHISM cap uses. With `nwm_coupling=true`, no provider and no stub the
cap aborts rather than silently injecting zero (see open item 8).

**D7 -- One-way, 2DDI focus.** No temperature/salinity handling (unlike SCHISM's
`ath3`, which carries tracers).

**D8 -- Input hygiene.** `dnwm` already aborts the run on NaN/Inf and clamps
negatives and fill values to zero upstream, so the cap does not need to re-guard.
A cheap defensive `max(0, .)` at the `QNIN` fill is acceptable.

---

## 3. Hazards

**H1 -- Parallel decomposition / ghost nodes.** `adcprep` splits `fort.20` per
subdomain (`src/boundaries.F`, around line 500). The ESMF import field covers
OWNED nodes only, but a river-boundary node can also exist as a GHOST in a
neighbouring subdomain that still expects its `QNIN` entry filled. The SCHISM
work hit exactly this (owner-only read, then replication). ADCIRC already has the
right tool: the met-field import in `adc_cap.F90` fills owned nodes from the
import pointer via `mdataOut%owned_to_present_nodes` and then calls
`UPDATER(WVNX1, WVNY1, PRN1, 3)` (`src/messenger.F:663`) to fill ghosts. The
same pattern applies here -- scatter into an `NP`-sized node array, `UPDATER` it,
then map node -> `J` through `NBV`/`LBCODEI` to fill `QNIN`. A full-`NP`
allreduce is not needed.

**H2 -- Advertise/realize lifecycle.** An unconnected `river_volume_flux` must
not abort or be silently dropped in a way that breaks the rest of the state. The
SCHISM cap had a bug here (`SCHISM_RemoveUnconnectedFields` removed the field).
`adc_cap.F90:ADCIRC_RealizeFields` already only realizes connected fields and
records `%connected`, which is the hook D6 needs.

**H3 -- Boundary formulation is unchanged.** `FluxSettlingIT` ramping,
`NFLUXF`/`NFFR` bookkeeping and the IBTYPE 22/32 formulation keep working; only
the DATA SOURCE for `QNIN` changes.

**H4 -- Sequence-index base.** Verified: the ADCIRC ESMF mesh node ids ARE the
global `fort.14` node ids. `adc_mod.F90:EXTRACT_MSG_TABLE_FORT18` reads the
local-to-global node map out of `fort.18` into `the_data%NdIds`, and
`create_parallel_esmf_mesh_from_meshdata` passes those (abs-valued) as
`nodeIDs`. So the D2 contract "sequence index i == global node i" is sound.

**H5 -- The node width `w_i` is a SUBDOMAIN quantity and must be healed across
decomposition cuts.** This was found as a real defect in the first np=3 coupled
demo and is fixed; the detail is in section 3.1 below because it is the single
most subtle part of the cap.

### 3.1 H5 in full: the subdomain-local width defect, the fix, and the self-checks

#### The defect

`ADCIRC_RiverInit` builds the effective node width from ADCIRC's own boundary
segment array:

```fortran
prevlen = 0
do j = 1, NVEL
   w = 0.75d0*(prevlen + BNDLEN2O3(j))     ! = 0.5*(L_{j-1} + L_j)
   prevlen = BNDLEN2O3(j)
   ...
   river_width(j) = w
```

`mesh.F:2320` sets `BNDLEN2O3(j) = (2/3)|P_{j+1} - P_j|` and, because
`NBVV(k,NVELL(k)+1) = NBVV(k,NVELL(k))`, the LAST entry of every boundary string
gets zero. That terminator is what makes the recursion restart correctly at each
string, and it is the same convention the GWCE relies on.

`BNDLEN2O3` is built from the SUBDOMAIN mesh, though. `adcprep` splits a global
river string across subdomains and extends each piece by the neighbouring
subdomain's node, so a node that is INTERIOR to the global string but sits at
the end of a subdomain's truncated copy of it gets the zero terminator on one
side and therefore HALF its true width -- and `ADCIRC_ImportRiver`, which sets
`QNIN(j) = Q_i / (w_j * n_j)`, gives it exactly 2x the correct per-unit-width
flux.

Those entries are always GHOSTS, which is why the defect looks harmless and is
not. The GWCE specified-normal-flow load assembly (`src/gwce.F:1845-1852`) is
segment based:

```fortran
BndLenO6NC = NCI*NCJ*BndLen2O3(J-1)/4.d0            ! = L/6
GWCE_LV(NBDI) = GWCE_LV(NBDI) + BndLenO6NC*(2*QForceI + QForceJ)
GWCE_LV(NBDJ) = GWCE_LV(NBDJ) + BndLenO6NC*(2*QForceJ + QForceI)
```

so each node's load depends on `QN2` at its NEIGHBOUR. A doubled ghost `QNIN`
therefore corrupts the OWNED node on the other side of the cut.

Measured on the `adcirc_rivers` demo mesh at np=3 (7-node string
`[1532,1452,1372,1295,1219,1143,1142]` at 50 m3/s, 5-node string
`[517,466,419,375,373]` at 200 m3/s, total 250 m3/s): global nodes 1219 (ghost
on PE0000) and 1143 (ghost on PE0002) each got half their width, the effective
injection came out 252.777778 m3/s instead of 250 (+1.11% overall, +5.56% on the
split string), and the water level differed from standalone ADCIRC by
max 1.84e-4 m, 1.25% of the 1.47e-2 m signal. Invisible at np=1 and at any np
that keeps every river string inside one subdomain.

Both of the cap's original diagnostics are structurally blind to it:

- `injected river discharge [m3 s-1]` sums `qflux*w` over owned entries, and
  `qflux` is `Q/w`, so `w` cancels identically. It printed 2.50000000E+02
  throughout the buggy run.
- `global river-boundary width [m]` sums `river_width` over OWNED entries only,
  and every OWNED entry does have its correct full width -- the half widths live
  exclusively on ghosts. It printed the correct 5.12338715E+03 m.
- the "no usable boundary width" warning only fires for `w <= river_wmin`; a
  half width is not small.

#### The fix

`adcprep` extends every local boundary string by the neighbouring subdomain's
node, so the OWNER of a river node always sees both of its boundary segments and
always holds the global-string width. Only ghost copies at the ends of a
truncated local piece are short. `ADCIRC_RiverInit` therefore does, once at init:

1. build the subdomain-local per-entry width exactly as before;
2. project it onto a node-indexed vector of length `np` (LOCAL node count);
3. `call UPDATER(wnode_own, river_scratch, river_scratch, 1)` -- one
   nearest-neighbour ghost exchange, the same ADCIRC messenger idiom the met
   import and `ADCIRC_ImportRiver` already use -- which overwrites every ghost
   slot with the owner's value;
4. write the owner value back into `river_width(j)` for every local entry.

Properties:

- Correct at any `np`, wherever the cuts land, and for a string cut any number of
  times: each cut node is owned by exactly one PET, and that PET has the whole
  stencil.
- One-time, init-only. Nothing is added to the per-window path, which still has
  exactly the one `UPDATER` and the one `ESMF_VMAllReduce` it had before.
- Memory `O(np_local)`; no global-node-sized array and no full-`NP` allreduce.
  (A global reconstruction -- deposit each string edge into a global-node-indexed
  vector and `ESMF_VMAllReduce(MAX)` -- also works and needs no assumption about
  `adcprep`, but costs `2*NP_global` reals per rank at init, which is 138 MB on a
  10M-node mesh. The self-checks below verify the cheap route instead of trusting
  it.)
- Collective-safe: every ADCIRC PET calls `UPDATER` unconditionally, including
  PETs with zero river nodes (`NVEL == 0` just leaves the vector at zero). The
  two `ESMF_VMAllReduce` calls that follow are likewise reached by every PET
  before any abort decision is taken, which is the existing discipline in this
  routine.
- Nodes that carry MORE THAN ONE flux-boundary entry keep their per-entry widths
  (a node-indexed heal cannot represent two different strings). This is the
  already-warned duplicate case and the warning now says so.

#### The self-checks (new, init-only)

Three lines were added, all cheap and all reduced across every PET:

```
(adc_cap:ADCIRC_RiverInit) local string 1: river entries = 3, global nodes 1219..1142, width sum [m] =   1.20985352E+03, healed = 1
(adc_cap:ADCIRC_RiverInit) river-boundary width [m]: owned view =   5.12338715E+03, all local views (owned+ghost) =   6.09126996E+03
(adc_cap:ADCIRC_RiverInit) H5 decomposition width heal: entries healed = 2, max relative deviation =   5.00000005E-01
(adc_cap:ADCIRC_RiverInit) GWCE conservation identity: expected =   1.20000000E+01, assembled =   1.20000000E+01, rel err =   0.0000E+00
```

1. **Per-string totals.** One line per LOCAL boundary string that carries river
   entries: entry count, first/last GLOBAL node id, width sum, and how many of
   its entries had to be healed. Local strings are delimited by the `BNDLEN2O3`
   zero terminator, i.e. by the same rule the GWCE uses, so a decomposition cut
   is exactly where a local string ends.
2. **Owned view vs all local views.** The owned-only width total is blind to a
   truncated ghost by construction, so the total over EVERY local entry (owned
   plus ghost) is reported next to it, together with the number of entries healed
   and the largest owner-vs-local relative deviation seen before healing.
3. **GWCE conservation identity.** The specified-normal-flow assembly is re-run
   at init with a unit probe, `q_j = 1/(w_j*n_j)`, which puts exactly 1 m3/s on
   every river NODE. Summed over owned nodes the assembly must return the number
   of owned river nodes. This is the ONLY check that can see a wrong width on a
   ghost entry, because the assembly at an owned node consumes its neighbour's
   `QN`. It is an exact identity in real arithmetic; a relative error above
   1e-10 (`river_wtol`) raises a loud WARNING. `NodeCode` is deliberately left
   out -- every node is wet at init and this is a geometric identity. Note the
   probe's relative error is NOT the discharge error: it is measured with a
   uniform 1 m3/s per node, so it reports sensitivity, not the actual
   over-injection. Any nonzero value means the run does not conserve.

Two further WARNINGs cover the cases the heal cannot repair: an entry WIDER than
the copy on its owning subdomain (the owner's string is truncated -- see the
`adcprep` hanging-node drop in section 5), and an entry whose owning subdomain
has no flux-boundary entry for that node at all.

Rebuilt with the heal disabled (one line, container only) the same np=3 demo
reproduces the pre-fix run bit for bit and the checks read

```
river-boundary width [m]: owned view =   5.12338715E+03, all local views (owned+ghost) =   5.60732855E+03
H5 decomposition width heal: entries healed = 2, max relative deviation =   5.00000005E-01
GWCE conservation identity: expected =   1.20000000E+01, assembled =   1.23333333E+01, rel err =   2.7778E-02
WARNING: river discharge will NOT be conserved: the assembled GWCE boundary load is off by   2.7778E-02 ...
```

while `injected river discharge [m3 s-1] = 2.50000000E+02` stays reassuringly
wrong. That is the class of bug these three lines exist to catch.

#### Verification

With the fix, the np=3 coupled run reproduces standalone ADCIRC driven by a
`fort.20` built from GLOBAL `BNDLEN2O3` widths to **max |dz| = 2.1e-14 m**
(RMS 7.1e-15 m, 1.4e-12 of the signal) over 12 records x 6509 nodes, and the
np=1 coupled control is BIT-IDENTICAL to the np=3 run (max |dz| = 0.0). Full
numbers: `logs/adcirc_nwm/fix_width/`.

---

## 4. Files changed / added

Split across three workstreams. Only section 4.2 is complete and in the tree as
described; 4.1 and 4.3 are stated as INTENT (owned by the other workstreams) and
should be re-read from the code once they land.

### 4.1 ADCIRC cap + core guards (intent -- `ADCIRC-interface/**`, root `CMakeLists.txt`)

| Path | Intended change |
| --- | --- |
| `ADCIRC-interface/ADCIRC/thirdparty/nuopc/adc_cap.F90` | Add `river_volume_flux` to `fldsToAdc` via `fld_list_add` (next to the met fields near lines 534-541); advertise/realize node-located, wrapped in `#ifdef USE_NUOPC_RIVER`. |
| `ADCIRC-interface/ADCIRC/thirdparty/nuopc/adc_cap.F90` or `adc_mod.F90` | Import routine: once per coupling window, if connected, read the field, use cached `w_i`, fill `QNIN1`/`QNIN2` + `QTIME1`/`QTIME2` per D3/D4, replicate to ghosts per H1. If not connected, do nothing. |
| `ADCIRC-interface/ADCIRC/src/cstart.F` | Skip the `fort.20` open + the two initial `QNIN1`/`QNIN2` read loops when the coupled river is active. FIXED-FORM Fortran. |
| `ADCIRC-interface/ADCIRC/src/timestep.F` (and `src/hstart.F` for hot starts) | Skip the aperiodic `READ(20,...)` rotate branch when the coupled river is active. FIXED-FORM Fortran. |
| `ADCIRC-interface/CMakeLists.txt` | `if(USE_NUOPC_RIVER) target_compile_definitions(adcirc PRIVATE USE_NUOPC_RIVER)` on the single `adcirc` target, which carries both the fixed-form core guards and the free-form cap (D5), plus a warning when `USE_NUOPC_RIVER` is given without `COUPLED=ON`. Mirrors `SCHISM-interface/CMakeLists.txt`. |
| root `CMakeLists.txt` | Reports `USE_NUOPC_RIVER` and warns when it is ON with neither `ADCIRC` nor `SCHISM` enabled. Deliberately NOT declared as a `set(... CACHE BOOL ...)` option, so that an unset flag stays undefined everywhere and does not collide with SCHISM's `define_opt()` declaration; `-DUSE_NUOPC_RIVER=ON` on the cmake command line creates the untyped cache entry both interface files test. |

Init-time logging should mirror dnwm's style: number of river-boundary nodes
found, total width per string, and per-window total injected discharge
(`sum Q_i`).

### 4.2 Superproject config + tests (this workstream)

| Path | Change |
| --- | --- |
| `tests/tests/coastal_ike_shinnecock_atm2adc_nwm_cdeps` | NEW. `coastal_ike_shinnecock_atm2adc` plus the dnwm block; merges the ADCIRC side of `coastal_ike_shinnecock_atm2adc` with the NWM side of `coastal_ike_shinnecock_atm2sch_nwm_cdeps`. |
| `tests/rt_coastal.conf` | Added a COMMENTED `COMPILE` line (`atm2adc_nwm`, the `atm2adc` flags plus `-DUSE_NUOPC_RIVER=ON`) and a COMMENTED `RUN` line for `coastal_ike_shinnecock_atm2adc_nwm_cdeps`, next to the existing ADCIRC block, with a note about the staging gap. |
| `tests/parm/ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` | REUSED, not forked. Every directive was already `@[]`-templated and OCN-model agnostic; only the header comment changed, to document both the SCHISM (per-element) and ADCIRC (per-node) id contracts and to note that `river_stub_q`/`nwm_coupling` are read only by caps that query them. |
| `tests/parm/dnwm_in.IN`, `tests/parm/dnwm.streams.coastal.IN` | UNCHANGED. Already generic: `@[MESH_NWM]`, `@[NWM_NX_GLB]`, `@[nwm_stream_files]`, `@[nwm_stream_variables]`. |
| `tests/parm/fd_ufs.yaml` | `river_volume_flux` entry already present with the right `standard_name` and `canonical_units: m3 s-1`; only its free-text `description` was widened from "per SCHISM source element" to cover the ADCIRC per-node case. No functional change. |

Variables the ADCIRC NWM test adds relative to `coastal_ike_shinnecock_atm2adc`:

```
NWM_CDEPS=true
DNWM_STREAM_CONFIGURE=dnwm.streams.coastal.IN
MESH_NWM="adcirc_node_ESMFmesh.nc"      # per-NODE mesh (SCHISM route: per-element)
NWM_NX_GLB=<NP>                         # ADCIRC global NODE count (SCHISM route: nelem)
nwm_stream_files="\"INPUT/nwm_discharge.nc\""
nwm_stream_variables="\"river_discharge river_volume_flux\""
river_stub_q=100.0                      # unused on a connected run
nwm_model=dnwm                          # MUST be exactly "dnwm"
nwm_petlist_bounds="0 10"
nwm_omp_num_threads=1
UFS_CONFIGURE=ufs.configure.coastal_datm_ocn_nwm_cdeps.IN   # was ..._datm_ocn.IN
```

and one variable it CHANGES relative to that sibling:

```
meshloc=node                            # sibling uses element; see open item 3
```

Everything else (dates 2008-09-04 +36 h, `HURR=Ike`, `CPL_CONF=atm2adc`,
`ATM_compute_tasks=11`, `OCN_tasks=11`, `LIST_FILES`, `FV3_RUN`) is
byte-identical to `coastal_ike_shinnecock_atm2adc`.

Note on PET bounds: `rt_utils.sh:compute_petbounds_and_tasks_esmf_threading`
unsets and recomputes `atm_/ocn_/med_petlist_bounds` from the `*_tasks`
variables, so exporting them in a test file has no effect. It has no `NWM`
branch, so `nwm_petlist_bounds` and `nwm_omp_num_threads` MUST be exported by
the test. With `ATM_compute_tasks=11` / `OCN_tasks=11` the harness produces
ATM `0 10`, OCN `11 21`, MED `0 10`, `TASKS=22`; dnwm is placed on `0 10` so it
adds no tasks.

### 4.3 Phase 0 offline tool (intent -- `logs/adcirc_nwm/`)

| Path | Intended content |
| --- | --- |
| `logs/adcirc_nwm/gen_adcirc_nwm_phase0.py` | Parse `fort.14` (nodes, elements, open/land boundary tables), find the IBTYPE 22/32 strings, compute per-node widths `w_i` (haversine for lon/lat, planar for cartesian, `--coords {sphere,cart}`), split each reach's `Q` across its string proportional to a conveyance weight `w_i * max(h_i, hmin)` (`--hmin` default 0.1 m) normalized so `sum_i Q_i == Q` exactly. |
| `logs/adcirc_nwm/gen_adcirc_nwm_phase0.py` outputs | `adcirc_node_ESMFmesh.nc` (ESMF mesh file, ONE ELEMENT PER GLOBAL NODE in natural node order -- only the ordering matters for `redist`); `nwm_discharge.nc` with `river_discharge(time, nnode)` full-size sparse; and a mapping sidecar JSON (string -> nodes -> fractions -> widths) for auditability. |
| `logs/adcirc_nwm/README.md` | Usage plus the D2 contract statement. |

Modelled on the SCHISM Phase-0 generator `logs/gen_nwm_phase0.py`.

---

## 5. How to build and run

### Build

```
# from the superproject root
cmake -DAPP=CSTLA -DADCIRC_CONFIG=PADCIRC -DCOUPLED=ON -DMPI=ON -DUSE_NUOPC_RIVER=ON ...
```

Same as the existing `atm2adc` compile line in `tests/rt_coastal.conf` with
`-DUSE_NUOPC_RIVER=ON` appended. With the flag OFF the build and the results
must be byte-identical to today's `atm2adc`.

**`-DMPI=ON` is mandatory and is easy to lose.** `ADCIRC-interface/CMakeLists.txt`
gates `add_mpi(adcirc)` -- which is what defines `CMPI` -- on the bare CMake
variable `MPI`, and `find_package(MPI REQUIRED)` in the root `CMakeLists.txt`
does NOT set it. The RT harness always supplies it
(`tests/compile.sh:86: CMAKE_FLAGS+=" -DMPI=ON"`), so a hand-written cmake line
is the only place it goes missing. Without it ADCIRC is compiled SERIALLY inside
a parallel executable: every OCN rank reports `myProc == 0`, takes the non-`CMPI`
branch at `src/sizes.F:372` (`LOCALDIR = ROOTDIR = '.'`), ignores the `PE0000/`
subdomain directories, reads the GLOBAL `fort.14`, and all ranks then race to
create `./fort.63.nc`:

```
ERROR: check_err: Permission denied
```

The failure is silent at build time. Confirm `-DCMPI` appears in
`ADCIRC-interface/CMakeFiles/adcirc.dir/flags.make` before trusting a run.

### Phase 0 (offline, once per mesh)

```
python3 logs/adcirc_nwm/gen_adcirc_nwm_phase0.py \
    --fort14 /path/to/fort.14 \
    --spec   /path/to/reach_spec.json \
    --coords sphere \
    --outdir phase0_adcirc_nwm_out
```

Produces `adcirc_node_ESMFmesh.nc`, `nwm_discharge.nc` and the mapping sidecar.
Take `NP` (the `nnode` dimension) from the sidecar and put it in the test as
`NWM_NX_GLB`.

### Run

RT auto-staging of the dnwm inputs is NOT wired yet, so a first connected run is
manual, exactly as for the SCHISM `_nwm_cdeps` test:

1. Let the harness build the run directory for
   `coastal_ike_shinnecock_atm2adc_nwm_cdeps` (or copy a working `atm2adc` run
   directory).
2. Copy `phase0_adcirc_nwm_out/adcirc_node_ESMFmesh.nc` and
   `phase0_adcirc_nwm_out/nwm_discharge.nc` into `INPUT/`.
3. Generate `dnwm_in` and `dnwm.streams` by hand from
   `tests/parm/dnwm_in.IN` and `tests/parm/dnwm.streams.coastal.IN`
   (substituting `MESH_NWM`, `NWM_NX_GLB`, `dtlimit`, `SYEAR`,
   `nwm_stream_files`, `nwm_stream_variables`).
4. Make sure `ufs.configure` is the `..._nwm_cdeps` variant, i.e. that it has
   the `NWM` component, `NWM_model: dnwm`, and the
   `NWM -> OCN :remapMethod=redist` line in `runSeq`.
5. Make sure `fort.14` really has an IBTYPE 22/32 normal-flux string and that
   `fort.15` (`NFLUXF`, `NFFR`) agrees with it.
6. Uncomment the two lines in `tests/rt_coastal.conf` only once the above is
   automated and a baseline exists.

### Runtime gotchas found in the first end-to-end runs

**The shipped `runSeq` order is wrong for the DIRECT (mediator-less) connector,
and the ADCIRC cap is stricter than the SCHISM cap here.**
`tests/parm/ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` orders the slot as
`NWM` / `NWM -> OCN :remapMethod=redist` / `OCN`. On the direct connector that
aborts on the very first window:

```
WARNING PET2  OCN: Field 'river_volume_flux' in the importState is not at the expected time.
ERROR   PET2  OCN: .../NUOPC_ModelBase.F90:2473 Invalid argument
        - NUOPC INCOMPATIBILITY DETECTED: Import Fields not at current time
```

NUOPC stamps a model's export state at the END of its advance and there is no
mediator on this route to re-stamp it, so `dnwm`'s own
`NUOPC_SetTimestamp(exportState, clock)` (currTime) is superseded and OCN sees a
field 1800 s in its future. Config-only workaround -- advance NWM at the END of
the slot:

```
runSeq::
@1800
  NWM -> OCN :remapMethod=redist
  OCN
  NWM
@
::
```

Window 1 then consumes the export `dnwm` stamps in `InitializeRealize` (start
time) and window k>1 consumes the export stamped at the end of window k-1, which
is the start of window k. Pair it with `advance_to_next_time = .true.` in
`dnwm_in` so the VALUE `dnwm` interpolates is the one valid at the timestamp it
carries (zero lag). This contradicts the note in
`CDEPS-interface/CDEPS/dnwm/nwm_comp_nuopc.F90` that `advance_to_next_time=.false.`
plus NWM-first is "the validated direct NWM->OCN route": either the SCHISM cap
tolerates the future timestamp where the ADCIRC cap (a plain `NUOPC_Model` with
the default `checkImport` -- its `label_CheckImport` specialization is commented
out in `adc_cap.F90`) does not, or the SCHISM validation used the mediator
route. The shipped template should be fixed once that is settled, rather than
leaving every run to rediscover it.

**`adcprep` can silently DROP a river node at an unlucky `np`.** METIS
partitions can cut a river string so that one subdomain gets a single node of
it, which `adcprep` then removes:

```
ERROR: The land boundary number 1 is only one node long in subdomain 3.
INFO: Eliminating hanging boundary node from this subdomain.
```

On the `adcirc_rivers` demo mesh this happens at np=4 and np=6 for the 7-node
string; np=2, 3 and 5 are clean. The consequence is an owned river node with no
boundary entry at all (`w = 0`), so its share of the discharge is lost and the
GWCE never assembles that segment. No width fix can repair it -- the entry does
not exist. The cap now detects it: `ADCIRC_ImportRiver`'s "node(s) carried
discharge but have no usable boundary width" warning fires, and the H5 checks add
"river entry/entries are WIDER than the copy on the owning subdomain" plus the
GWCE conservation identity, which will show a shortfall. Grep the `adcprep`
prepall log for `hanging boundary node` before trusting a decomposition.

**`fort.20` must still exist for `adcprep` even though the run never reads it.**
`prep/prep.F::PREP20` opens unit 20 whenever `NFLUXF=1` and drops into an
interactive prompt (then dies on EOF) if the file is missing. Supply a dummy
zero `fort.20`; `run.log` will confirm
`NORMAL FLOW INFORMATION SUPPLIED BY THE NUOPC COUPLER (UNIT 20 NOT USED)`.
Also note `adcprep --prepall` creates the GLOBAL `fort.63.nc` / `fort.64.nc` /
`maxele.63.nc` / `maxvel.63.nc` skeletons; ADCIRC refuses to create them itself,
so they must be regenerated if they are deleted. And CDEPS wants
`dnwm.streams.xml`, not `dnwm.streams` (`nwm_comp_nuopc.F90` appends `.xml`
unless `DISABLE_FOX`); a missing file gives an `ESMF_Config` open error followed
by a segfault in `dshr_strdata_mod.F90`.

### What to check in the logs

- The cap logs `river_volume_flux ... is connected.` (not "not connected").
- Init log: number of IBTYPE 22/32 nodes found and the total width per string.
- Per window: total injected discharge `sum Q_i` matches what
  `nwm_discharge.nc` holds for that time, and `sum_i QNIN_i * w_i` equals it.
  Treat this line as NECESSARY BUT NOT SUFFICIENT: `w` cancels out of it, so it
  cannot see a wrong width (see H5 / section 3.1).
- `GWCE conservation identity: ... rel err = 0.0000E+00`. This is the check that
  actually proves the discharge will be conserved. Anything above 1e-10 also
  raises a WARNING.
- `H5 decomposition width heal: entries healed = N`. `N > 0` is normal and
  expected whenever a river string is split across subdomains; it means the heal
  did its job. `max relative deviation` is the pre-heal error (0.5 = a halved
  ghost width).
- No `WARNING: river entry/entries are WIDER than the copy on the owning
  subdomain` and no `whose OWNING subdomain has no flux-boundary entry` -- both
  mean `adcprep` mangled the boundary string for that `np`.
- `dnwm.log` shows no clamping/NaN messages.

---

## 6. Open items

1. ~~**No compile verification.**~~ CLOSED. The cap and the fixed-form `fort.20`
   guards have been built (`-DAPP=CSTLA -DADCIRC_CONFIG=PADCIRC -DCOUPLED=ON
   -DMPI=ON -DUSE_NUOPC_RIVER=ON`, gfortran 13, spack-stack) and run end to end
   against a real CDEPS `dnwm` component on the `adcirc_rivers` test mesh at
   np=3 and np=1, and validated against standalone `padcirc` driven through the
   native `fort.20` path to 2.1e-14 m. See `logs/adcirc_nwm/demo_rivers/` and
   `logs/adcirc_nwm/fix_width/`. What is still NOT verified: a hot start
   (`src/hstart.F` guard), IBTYPE 32, more than one coupling component (no ATM,
   no MED in these runs), and any mesh other than the demo one.
2. **`USE_NUOPC_RIVER` is an untyped command-line-only flag.** By design (see
   the note in the root `CMakeLists.txt`) it is not a cache option, so it does
   not appear in `cmake -L` output or in `ufs-weather-model` build docs; the only
   signal is the `USE_NUOPC_RIVER .. ON/OFF (not set)` line the root file prints.
   Verify that line in the build log before trusting a coupled-river run.
3. **`meshloc` is global per component in the ADCIRC cap, so the NWM test runs
   with `meshloc=node`.** `ADCIRC_RealizeFields` creates every field with the
   single `meshloc` OCN attribute (`ESMF_FieldCreate(..., meshloc=meshloc)`),
   and every other coastal test exports `meshloc=element`. The D2 contract is
   per-NODE, and the cap enforces it: a CONNECTED `river_volume_flux` with
   `meshloc /= node` is rejected with `ESMF_RC_NOT_VALID` in
   `ADCIRC_RealizeFields`, and `ADCIRC_ImportRiver` asserts
   `size(ptr) == mdataIn%NumOwnedNd`. `coastal_ike_shinnecock_atm2adc_nwm_cdeps`
   therefore exports `meshloc=node` (also the cap's default when the attribute
   is absent, and consistent with the met-field import, which already indexes
   owned NODES via `owned_to_present_nodes` / `NumOwnedNd`). Consequence: the
   ATM -> OCN half of the NWM test is NOT bit-comparable to the plain `atm2adc`
   test. Note `eliminate_ghosts` is called only on the element path, so the node
   path is a genuinely different (though cap-default) code path and needs its
   own first run before it is trusted.
4. **RT auto-staging gap.** `NWM_CDEPS`, `DNWM_STREAM_CONFIGURE`, `MESH_NWM`,
   `NWM_NX_GLB`, `nwm_stream_files`, `nwm_stream_variables` are consumed only by
   the `tests/parm` templates; nothing in `default_vars.sh` / `run_test.sh`
   generates `dnwm_in` / `dnwm.streams` or stages the Phase-0 NetCDF files (the
   `datm` equivalents are `run_test.sh:353-354`). Closing this needs an
   `export_dnwm_cdeps` in `default_vars.sh`, two `atparse` calls in
   `run_test.sh` gated on `NWM_CDEPS`, and a `fv3_conf` staging snippet. It is
   deliberately out of scope here and is the same gap the SCHISM `_nwm_cdeps`
   test has.
5. **No real river-capable `fort.14`.** The Shinnecock test mesh has no
   IBTYPE 22/32 string (the SCHISM/Duck counterpart had zero river boundaries
   too and needed a synthetic carved channel). A real pairing needs either a
   mesh with a genuine river boundary or a synthetic string added to `fort.14`
   plus matching `fort.15` bookkeeping.
6. **`NWM_NX_GLB` in the new test is a placeholder** (`3070`, the classic
   Shinnecock Inlet `NP`; `NE=5780`). It must be re-verified against the
   `fort.14` actually staged for the case before any run.
7. **No baseline.** The commented `RUN` line has an empty baseline field on
   purpose.
8. **Stub mode is not reachable from the RT templates.** Both caps gate the
   constant-discharge stub on an OCN attribute `river_stub` (`river_stub_q`
   alone is only the magnitude), but neither
   `ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` nor
   `ufs.configure.coastal_datm_ocn_nwm.IN` emits a `river_stub` line, so the
   stub can only be exercised by hand-editing the generated `ufs.configure`.
   Adding `river_stub = @[river_stub]` to the shared template also needs a
   default in `tests/default_vars.sh`, otherwise the existing SCHISM
   `_nwm_cdeps` test leaves the token unresolved. Deliberately not done here
   because it changes a template two SCHISM tests already depend on. Note this
   is a fail-loud gap, not a silent one: with `nwm_coupling=true`, no provider
   and no stub, `ADCIRC_RiverInit` aborts with `ESMF_RC_NOT_VALID`.

---

## 7. References

- `NWM_SCHISM_COUPLING_DOCUMENTATION.md` -- the completed SCHISM route.
- `NWM_RIVER_INCREMENT1_VALIDATION.md` -- SCHISM stub/connected validation
  procedure, including the manual staging recipe this document's section 5
  follows.
- `CDEPS-interface/CDEPS/dnwm/nwm_comp_nuopc.F90` -- the data component.
- `SCHISM-interface/SCHISM-ESMF/src/schism/schism_nuopc_cap.F90` -- the proven
  cap-side ingest idiom (`USE_NUOPC_RIVER`, `SCHISM_ImportRiver`,
  `NUOPC_IsConnected` fallback).
- `ADCIRC-interface/ADCIRC/thirdparty/nuopc/adc_cap.F90`,
  `adc_mod.F90` -- the ADCIRC cap and its mesh construction.
- `logs/gen_nwm_phase0.py` -- the SCHISM Phase-0 generator the ADCIRC one is
  modelled on.
- `logs/adcirc_nwm/demo_rivers/NOTES.md` -- the first end-to-end
  `dnwm -> ADCIRC` run (np=3, 250 m3/s, 6 h) and every deviation it needed.
- `logs/adcirc_nwm/demo_rivers/compare_standalone/COMPARISON.md` -- the coupled
  vs standalone-`padcirc` comparison that isolated the H5 defect.
- `logs/adcirc_nwm/fix_width/FIXNOTES.md` -- the H5 fix, its self-checks and the
  post-fix verification numbers.
