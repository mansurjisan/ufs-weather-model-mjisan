# NWM → SCHISM one-way river forcing — Increment 1 (de-risking core)

Validates the **SCHISM receiving side** of online NWM river forcing, independent
of CDEPS: a `USE_NUOPC_RIVER` build injects a uniform constant discharge into
SCHISM source elements via the NUOPC cap, exercising
`init → cap fills ath3 → schism_step → vsource → continuity → output`.

Branch: `feature/nwm_river_forcing` (superproject + the `SCHISM-interface/SCHISM`
and `SCHISM-interface/SCHISM-ESMF` submodules).

## What changed (all `#ifdef USE_NUOPC_RIVER`; OFF ⇒ byte-identical to baseline)

| File | Change |
|---|---|
| `SCHISM-interface/SCHISM/src/Hydro/schism_init.F90` | `#error` guard vs `USE_BMI`; reuse sparse `source_sink.in` read, set `th_dt3=3600` |
| `SCHISM-interface/SCHISM/src/Hydro/misc_subs.F90` (`other_hot_init`) | skip `.th`/`.nc` opens; cold-start init `ath3=0`, `th_time3=[time,time+th_dt3]` |
| `SCHISM-interface/SCHISM/src/Hydro/schism_step.F90` | fold `USE_NUOPC_RIVER` into the `USE_BMI` skip-read/rotate branch |
| `SCHISM-interface/SCHISM-ESMF/src/schism/schism_nuopc_cap.F90` | module var `river_stub_q` + attribute read; `SCHISM_ImportRiver` fills replicated `ath3` (both time levels, zero-order hold); called in `ModelAdvance` before the step loop |
| `SCHISM-interface/SCHISM/src/CMakeLists.txt` | `define_opt(USE_NUOPC_RIVER ... OFF)` → `-D` for SCHISM core (init/step/misc_subs) |
| `SCHISM-interface/CMakeLists.txt` | `target_compile_definitions(schism PRIVATE USE_NUOPC_RIVER)` → `-D` for the cap |
| `tests/parm/ufs.configure.coastal_datm_ocn_nwm.IN` | base config + OCN attr `river_stub_q` |
| `tests/tests/coastal_ike_shinnecock_atm2sch_nwm` | test clone (uses `_nwm` config, exports `river_stub_q`) |
| `tests/parm/coastal_shinnecock_source_sink.in.template` | sparse `source_sink.in` template |
| `tests/rt_coastal.conf` | `atm2sch_nwm` COMPILE/RUN entries (commented until data+baseline exist) |

**Design:** new runtime `if_source` value was ruled out (`iabs(if_source)>1` aborts),
so a compile flag is used, mutually exclusive with `USE_BMI`. `ath3` is the full,
rank-replicated source array, so the cap fills it identically on every rank; the
existing apply loop (`iegl(ieg_source(i))%rank==myrank`) handles per-rank
ownership. msource stays at the ambient default (S=0, T=−9999), so only `vsource`
is driven. The river **import field** (advertise/realize/read) is **deferred to
Phase 5** — an unconnected field would be removed by `SCHISM_RemoveUnconnectedFields`,
so Increment 1 fills `ath3` directly from the constant `river_stub_q`.

## Build (on your spack-stack HPC)

Add `-DUSE_NUOPC_RIVER=ON` to the SCHISM coastal build, e.g.:

```
-DAPP=CSTLS -DUSE_ATMOS=ON -DNO_PARMETIS=OFF -DOLDIO=ON -DUSE_NUOPC_RIVER=ON
```

(or uncomment the `atm2sch_nwm` COMPILE line in `tests/rt_coastal.conf`). With the
flag OFF the executable is unchanged from the stock `atm2sch` build — verify the
existing `coastal_ike_shinnecock_atm2sch` regression still passes.

## Stage input data (NOT in this repo — required to run)

The Shinnecock `param.nml`/`hgrid`/forcing are staged from input data:

1. In the staged SCHISM **`param.nml`** set **`if_source = 1`** (stock ships `0`)
   and keep `dramp_ss > 0` (smooth source ramp-up).
2. Place a **`source_sink.in`** in the SCHISM run dir. Start from
   `tests/parm/coastal_shinnecock_source_sink.in.template`, but **replace the
   placeholder element ids (1,2,3) with real wet land-boundary / river-mouth
   global element ids** for the Shinnecock mesh (interior/wet elements give a
   clean signal; a real pairing comes from
   `SCHISM/src/Utility/Pre-Processing/NWM/NWM_coupling`). Format: `nsources`,
   then one global elem id per line, a blank line, then `nsinks` (0 here).
3. No `vsource.th`/`msource.th`/`source.nc` are needed.

Set the discharge via `river_stub_q` (m³/s) in the test (default 100.0).

## Run & inspect (no baseline yet — inspect, don't byte-compare)

Run `coastal_ike_shinnecock_atm2sch_nwm` (or the base case with the staged
`param.nml`/`source_sink.in` and the `_nwm` ufs.configure). Confirm:

1. **Init**: SCHISM logs read `nsources` from the sparse `source_sink.in` (== your
   count), not `ne_global`; no allocation error.
2. **Cap**: ESMF log shows `... river_stub_q [m^3/s] = 100.000` (or your value).
3. **Cold-start**: completes the first step with **no** `STEP: wrong sign vsource`
   or `STEP: rat out in vsource.th` abort.
4. **Effect**: in `outputs/schout_*` (elevation/velocity), a localized freshwater
   inflow at the source elements; magnitude scales ~linearly with `river_stub_q`
   (try 0 vs 100 vs 500). `river_stub_q=0` ≈ the no-river run.
5. **Parallel**: with `OCN_tasks=3`, the source elements are split across SCHISM
   PETs — confirm the inflow appears regardless of which rank owns each element.

## Status of later phases
- **Phase 3–5 — DONE (code, compile-verified in the container):** `dnwm` CDEPS
  component + driver/CMake registration; the SCHISM cap advertises/realizes
  `river_volume_flux` and reads it per source when connected (replacing the
  constant); `river_volume_flux` in `fd_ufs.yaml`; full config
  `ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` (+ `NWM -> OCN` redist connector);
  `dnwm.streams.coastal.IN` / `dnwm_in.IN` templates.
- **Phase 0 — DONE (synthetic, generated 2026-05-29):** `logs/gen_nwm_phase0.py`
  parsed the Duck/"shinnecock" `hgrid.ll` → `phase0_nwm_out/`:
  `schism_elem_ESMFmesh.nc` (45167 elems, **natural order ⇒ ESMF seq index ==
  SCHISM global element id**, the redist contract — verified: file element 0 ==
  hgrid elem 1, centers match), `source_sink.in` (6 wet-shoreline source elems:
  3505/12140/12656/14892/25500/29582), `nwm_discharge.nc`
  (`river_discharge(time,nelem)` = 100 m³/s at those 6, 0 elsewhere, hourly).
  Real NWM streamflow would replace the synthetic discharge (via
  `NWM_coupling/coupling_nwm.f90` + `NWM_shp_ll.nc` + CHRTOUT). Mesh format mirrors
  the production ATM ESMFmesh; `ESMF_RegridWeightGen` fails identically on BOTH my
  mesh and the known-good ATM mesh in this container (an RWG/env issue, not a mesh
  issue) — definitive read test is the coupled run.
- **Phase 6** (optional): CMEPS `comprof` mediator route for an upstream PR.

## Phase 5 connected-run prerequisites & known limitations (from sanity review)
The connected (`_nwm_cdeps`) path compiles and is decomposition-corrected, but has
NOT been runtime-exercised. Before a real connected run:
- **CRITICAL — element-id contract (silent-failure risk):** the `NWM -> OCN`
  connector uses `remapMethod=redist`, which routes by global sequence index. The
  dnwm `MESH_NWM` element ESMFmesh MUST have element ids identical to SCHISM's
  global element ids **and** the same element count (`ne_global`); otherwise
  discharge lands on the WRONG elements with no error. Phase 0 must generate the
  mesh from the SCHISM hgrid accordingly. (Differing PET decomposition is fine;
  mismatched id *values* are not.) There is no runtime assertion for this yet —
  treat it as a release gate; consider an element-id round-trip check.
- **Connection is required, not optional:** if the NWM component is present but the
  `river_volume_flux` field fails to connect, the cap **silently falls back to the
  constant `river_stub_q`** (a broken coupling can look "working"). Verify in the
  ESMF logs that the field is connected (not stubbed) for `_nwm_cdeps` runs.
- **`nwm_model` must be exactly `dnwm`** (the driver matches `trim(model)=="dnwm"`).
- **No connected RT test wired yet:** the `atm2sch_nwm` entries in `rt_coastal.conf`
  are commented; the `_nwm_cdeps` config + `nwm_*`/`MESH_NWM` test vars need a test.
- **meshloc / if_source:** now guarded — the cap aborts if `meshloc/=element`, and
  `schism_init` aborts on `if_source=-1`, under `USE_NUOPC_RIVER`.
- **Decomposition:** the connected path replicates `ath3` via `ESMF_VMAllReduce`
  (owner reads its element, then all-reduce) so owner+ghost source copies apply the
  same value — matching the standard bcast path. SH_MEM_COMM not yet exercised.

## Fixed in sanity review (2026-05-29)
- **Temperature bug (HIGH):** source elements were getting T=0 °C (msource overwrote
  the `-9999` ambient sentinel). Fixed: `other_hot_init` inits `ath3(:,1,:,3)=-9999`
  (T ambient); salinity stays 0 (fresh). Affected even the stub run.
- **Connected-path OOB / ghost replication (HIGH):** per-rank field read indexed
  past `farrayPtr1(1:ne)` for ghost source elements and didn't replicate. Fixed with
  the owner-only (`%id<=ne`) read + `ESMF_VMAllReduce`.

## Connected run — manual staging recipe (Phase 0 artifacts)
Until the RT auto-staging is wired, run the connected path
(`coastal_ike_shinnecock_atm2sch_nwm_cdeps`) by hand from a copy of the working
Duck run dir, with the `USE_NUOPC_RIVER`+dnwm `ufs_model`:
1. `cp phase0_nwm_out/schism_elem_ESMFmesh.nc phase0_nwm_out/nwm_discharge.nc <rundir>/INPUT/`
2. `cp phase0_nwm_out/source_sink.in <rundir>/` and set **`if_source = 1`** in `<rundir>/param.nml`.
3. Use `ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` (substitute `@[...]`) — `EARTH_component_list: ATM OCN NWM MED`, `NWM_model=dnwm`, the `NWM -> OCN :remapMethod=redist` connector.
4. Provide `dnwm_in` (from `dnwm_in.IN`: `model_meshfile="INPUT/schism_elem_ESMFmesh.nc"`, `nx_global=45167`, `datamode="copyall"`) and `dnwm.streams` (from `dnwm.streams.coastal.IN`: `stream_mesh_file="INPUT/schism_elem_ESMFmesh.nc"`, `stream_data_files="INPUT/nwm_discharge.nc"`, `stream_data_variables="river_discharge river_volume_flux"`, `mapalgo=redist`).
5. Run on the same PE layout (ATM/OCN/MED/NWM all `0 2`).
**Expect & verify:** ESMF log shows `river_volume_flux` **connected** (not the stub
fallback); freshwater inflow appears at the 6 source elements
(3505/12140/12656/14892/25500/29582 — centers near lon −75.7…−75.82); magnitude
matches the constant-stub run at `river_stub_q=100` (cross-check). A `0`-discharge
file ⇒ no inflow.
