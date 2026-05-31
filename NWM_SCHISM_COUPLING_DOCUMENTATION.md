# One-way NWM River Forcing into SCHISM through the UFS Coastal Coupler

**Online NUOPC / CDEPS / CMEPS coupling of National Water Model (NWM) streamflow
into SCHISM as one-way river/freshwater forcing.**

- Repository (forks): `mansurjisan/ufs-weather-model-mjisan`, branch `feature/nwm_river_forcing`
- Submodule forks (same branch): `mansurjisan/schism`, `mansurjisan/schism-esmf`, `mansurjisan/CDEPS`, `mansurjisan/CMEPS`
- Tracking issue: `mansurjisan/My-Workplan` #257
- Status: the **direct-connector** route is implemented, compile-verified, and runtime-validated in a local Docker container (connected run completes; river produces a physical, localized water-level + current response). The mediator route is scaffolded but not yet runnable (see section 9.4). Full DATM+SCHISM+dnwm case on a real mesh is a Hercules task.

---

## 1. Objective and design

Deliver NWM streamflow to SCHISM **inside the UFS coupled framework** (not via offline
`vsource.th`/`msource.th` preprocessing). One-way only: NWM forces SCHISM; there is no
backwater feedback to NWM.

**Design decisions (locked):**

1. **A new CDEPS data component `dnwm`** (cloned from `drof`) reads NWM discharge from a
   NetCDF stream and exports a single field `river_volume_flux` (m3 s-1, volumetric
   per-element discharge), `datamode=copyall`.
2. **SCHISM ingests it via a new compile flag `USE_NUOPC_RIVER`.** The SCHISM NUOPC cap
   fills SCHISM's `ath3` source/sink buffer each `ModelAdvance` (the same injection point
   the NextGen/BMI work uses), so no `.th`/`.nc` source files are read. A new runtime
   `if_source` value was ruled out because `schism_init` aborts on `iabs(if_source)>1`.
3. **Reach-to-element pairing is offline (Phase 0).** The runtime exchange is a **redist**
   (element-id -> element-id) of per-element discharge, avoiding any gridded
   m3/s-vs-kg/m2/s area-flux mismatch.
4. **Two delivery routes:**
   - **Direct connector** `NWM -> OCN :remapMethod=redist` (primary; the only
     runtime-validated route).
   - **CMEPS mediator route** (NWM as the ROF role, `ROF -> MED -> OCN`): the field
     advertise/map/merge is scaffolded in `esmFldsExchange_coastal_mod.F90` and compiles,
     but it does **not** yet run end-to-end (coastal-mode rof fractions and a
     `med_phases_post_rof` init are still needed to pass the realize phase — see section
     9.4). It is the recommended target for an upstream-clean PR, not a working route today.

**Data flow (direct connector):**

```
NWM discharge file (nwm_discharge.nc)
   -> dnwm (CDEPS) reads stream, exports river_volume_flux [m3/s] on element mesh
   -> NUOPC connector  NWM -> OCN  (remapMethod=redist, element-id -> element-id)
   -> SCHISM cap SCHISM_ImportRiver: fills ath3(:,1,2,1) [vsource, m3/s]
   -> schism_step: vsource applied at source elements -> continuity
   -> water-level / current response
```

**The element-id contract (critical):** `redist` routes by global sequence index, and
CMEPS `redist` **silently ignores unmatched indices**, so the dnwm mesh element ids MUST
equal SCHISM global element ids and have the same count (`ne_global`); a mismatch
misplaces discharge with **no run-time error**. Phase 0 builds the element ESMFmesh in
natural element order to guarantee this. Enforce it with the release-gate script
`logs/verify_id_contract.py` (asserts mesh count == `ne`, connectivity in natural order so
seqindex == global id, valid sources, `nsinks==0`) before any coupled run. A stronger
in-framework round-trip (dnwm exports its element ids; the cap compares them to `ielg`
after import) remains future work (section 9.5).

---

## 2. Files modified / added

### 2.1 SCHISM core (`SCHISM-interface/SCHISM`, submodule fork `mansurjisan/schism`)

| File | Change |
|---|---|
| `src/CMakeLists.txt` | `define_opt(USE_NUOPC_RIVER ... OFF)` |
| `src/Hydro/schism_init.F90` | `#error` guard vs `USE_BMI`; reuse sparse `source_sink.in` read; abort on `if_source=-1`; set `th_dt3=3600` |
| `src/Hydro/misc_subs.F90` (`other_hot_init`) | skip `.th`/`.nc` opens; cold-start init `ath3=0`, msource T=-9999 ambient, `th_time3=[time,time+th_dt3]` |
| `src/Hydro/schism_step.F90` | fold `USE_NUOPC_RIVER` into the `USE_BMI` skip-read/rotate branch |

### 2.2 SCHISM NUOPC cap (`SCHISM-interface/SCHISM-ESMF`, fork `mansurjisan/schism-esmf`)

| File | Change |
|---|---|
| `src/schism/schism_nuopc_cap.F90` | module vars `river_stub_q`, `keep_river_field`, `river_stub_mode`; read `river_stub_q` + `nwm_coupling` + `river_stub` attributes; advertise/realize `river_volume_flux` on the element mesh; `SCHISM_ImportRiver` (fills `ath3` from the connected field; aborts if `nsinks>0`; uses the constant stub only in explicit `river_stub` mode, else aborts); `CheckImportRiver` (re-checks all imports at currTime, exempts only `river_volume_flux`); exempt `river_volume_flux` from `SCHISM_RemoveUnconnectedFields` |
| `src/schism/schism_esmf_util.F90` | `SCHISM_MeshCreateElement`: read/loop the `elementIds` and connectivity export fields to `ne` (owned), not `nea` — the mesh is owned-only, so `nea` (ghosts) overran the `ne`-sized arrays |

### 2.3 CDEPS `dnwm` data component (`CDEPS-interface/CDEPS`, fork `mansurjisan/CDEPS`)

| File | Change |
|---|---|
| `dnwm/nwm_comp_nuopc.F90` | NEW (+544 lines): module `cdeps_dnwm_comp`, cloned from drof; exports `river_volume_flux`; registers it in the NUOPC field dictionary; guards `cpl_scalars` SetScalar; stamps export timestamp |
| `dnwm/CMakeLists.txt` | NEW (+33): standalone build of dnwm |
| `CMakeLists.txt` | add `dnwm` to the component `foreach` (standalone build only) |

### 2.4 CMEPS mediator (`CMEPS-interface/CMEPS`, fork `mansurjisan/CMEPS`)

| File | Change (+44 lines) |
|---|---|
| `mediator/esmFldsExchange_coastal_mod.F90` | import `comprof`; `rof_present` flag; read `ROF_model`; advertise + map(`mapfcopy`) + merge(`copy`) `river_volume_flux` rof->ocn; register the field in the dictionary on MED PETs |

### 2.5 Superproject (`ufs-weather-model`, fork `mansurjisan/ufs-weather-model-mjisan`)

| File | Change |
|---|---|
| `.gitmodules` | repoint SCHISM/SCHISM-ESMF/CDEPS/CMEPS to the `mansurjisan` forks on `feature/nwm_river_forcing` |
| `CMakeLists.txt` | `FRONT_CDEPS_DNWM=cdeps_dnwm_comp` in the `if(CDEPS)` block |
| `driver/UFSDriver.F90` | `FRONT_CDEPS_DNWM` use + `NUOPC_DriverAddComp` for `model=="dnwm"` (mirrors docn) |
| `SCHISM-interface/CMakeLists.txt` | `target_compile_definitions(schism PRIVATE USE_NUOPC_RIVER)` for the cap |
| `CDEPS-interface/CMakeLists.txt` | `add_library(dnwm OBJECT ...)` + `$<TARGET_OBJECTS:dnwm>` in the `cdeps` aggregate (the UFS build path) |
| `CDEPS-interface/cdeps_files.cmake` | `cdeps_dnwm_files` list |
| `tests/parm/fd_ufs.yaml` | `river_volume_flux` entry (m3 s-1) |
| `tests/parm/ufs.configure.coastal_datm_ocn_nwm.IN` | NEW: stub config (river_stub_q, no NWM comp) |
| `tests/parm/ufs.configure.coastal_datm_ocn_nwm_cdeps.IN` | NEW: connected config (NWM comp + `NWM -> OCN` connector + `nwm_coupling=true`) |
| `tests/parm/ufs.configure.coastal_datm_ocn_nwm_mediator.IN` | NEW: mediator-route config (ROF role) |
| `tests/parm/dnwm_in.IN`, `dnwm.streams.coastal.IN` | NEW: dnwm namelist + stream templates |
| `tests/tests/coastal_ike_shinnecock_atm2sch_nwm[_cdeps]` | NEW: RT test cases |
| `phase0_nwm_out/` | Phase-0 artifacts (element ESMFmesh, source_sink.in, nwm_discharge.nc) |

---

## 3. Key code (actual committed diffs)

### 3.1 New compile flag (SCHISM `src/CMakeLists.txt`)

```cmake
define_opt(USE_NUOPC_RIVER "Enable one-way NWM river forcing via the NUOPC coupler (SCHISM source/sink filled by the cap)" OFF)
```

### 3.2 Cold-start init of the source buffer (`misc_subs.F90`, `other_hot_init`)

```fortran
#elif defined(USE_NUOPC_RIVER)
    ! NUOPC river coupling: the SCHISM cap fills ath3 every ModelAdvance (which
    ! runs before schism_step), so there are no .th/.nc files to open here.
    ! Cold-start-initialize ath3 and the window [time, time+th_dt3] so the first
    ! step is valid (rat in [0,1]) and the ath3>=0 sign check passes at t0.
    if(if_source/=0) then
      ath3=0 !cold start; cap overwrites the new vsource level (>=0) before first step
      ath3(:,1,:,3)=-9999.   !msource temperature -> 'use ambient' sentinel (salinity stays 0 = fresh)
      th_time3(1,:)=time
      th_time3(2,:)=time+th_dt3(:)
    endif !if_source
```

### 3.3 Skip-read / rotate branch (`schism_step.F90`)

```fortran
#if defined(USE_BMI) || defined(USE_NUOPC_RIVER)
    !USE_BMI (NextGen) and USE_NUOPC_RIVER (UFS coupler) both fill ath3 at the
    !new time level externally before this point, so bypass the .th/.nc reads and
    !only rotate the time levels.
    if(nsources>0) then
      if(time>th_time3(2,1)) then
        ath3(:,1,1,1)=ath3(:,1,2,1)
        ...
```

### 3.4 The cap: filling `ath3` from the connected field (`schism_nuopc_cap.F90`, `SCHISM_ImportRiver`)

```fortran
! nsinks>0 is unsupported (the cap fills only vsource) -> abort early.
if (nsinks > 0) call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg="nsinks>0 unsupported by USE_NUOPC_RIVER", ...)

! Is 'river_volume_flux' present AND connected? (query presence first so a
! removed/unconnected field is routed to the explicit-mode decision below, not a crash.)
connected = .false.
call ESMF_StateGet(importState, itemName="river_volume_flux", itemType=itemtype_river, rc=localrc)
if (localrc == ESMF_SUCCESS .and. itemtype_river == ESMF_STATEITEM_FIELD) then
  call ESMF_StateGet(importState, itemName="river_volume_flux", field=field, rc=localrc)
  if (localrc == ESMF_SUCCESS) connected = NUOPC_IsConnected(field, rc=localrc)
end if

if (connected) then
  ! Connected to dnwm. The field is on the element mesh (farrayPtr1 indexed by
  ! local owned element ie). SCHISM applies a source at owner AND ghost copies and
  ! needs ath3 rank-REPLICATED, so each rank reads only the sources it OWNS
  ! (%id<=ne) then VMAllReduce(SUM) replicates the full source list to all ranks.
  call ESMF_FieldGet(field, farrayPtr=farrayPtr1, rc=localrc)
  allocate(q_local(nsources), q_global(nsources)); q_local = 0.0_ESMF_KIND_R8
  do i = 1, nsources
    if (iegl(ieg_source(i))%rank == myrank .and. iegl(ieg_source(i))%id <= ne) then
      ie = iegl(ieg_source(i))%id
      q_local(i) = max(0.0_ESMF_KIND_R8, farrayPtr1(ie))   ! owner reads its element [m^3/s]
    end if
  end do
  call ESMF_VMAllReduce(vm, q_local, q_global, nsources, ESMF_REDUCE_SUM, rc=localrc)
  ath3(1:nsources,1,2,1) = real(q_global, 4)        ! new vsource, replicated on all ranks
  ath3(1:nsources,1,1,1) = ath3(1:nsources,1,2,1)   ! zero-order hold (old=new)
else if (keep_river_field) then
  ! nwm_coupling=true: a real NWM provider was expected but is absent/unconnected.
  ! ABORT rather than silently fabricate flow (do not mask a mis-wired connector).
  call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg="nwm_coupling=true but river field unconnected", ...)
else if (river_stub_mode) then
  ! Explicit stub mode (river_stub=true, no NWM component): inject the constant river_stub_q.
  ath3(1:nsources,1,2,1) = real(max(0.0_ESMF_KIND_R8, river_stub_q), 4)
  ath3(1:nsources,1,1,1) = ath3(1:nsources,1,2,1)
else
  ! Neither a provider nor an explicit stub selected -> ABORT (previously silently stubbed).
  call ESMF_LogSetError(ESMF_RC_NOT_VALID, msg="no NWM provider and river_stub not set", ...)
end if
```

### 3.5 dnwm export-field registration (`CDEPS/dnwm/nwm_comp_nuopc.F90`)

```fortran
! river_volume_flux is not a CF/standard NUOPC field; register it (m3 s-1
! volumetric discharge) so it can be advertised, realized and connected.
if (.not. NUOPC_FieldDictionaryHasEntry(trim(fldname_river), rc=rc)) then
   call NUOPC_FieldDictionaryAddEntry(trim(fldname_river), "m3 s-1", rc=rc)
end if
```

### 3.6 Driver registration of dnwm (`driver/UFSDriver.F90`)

```fortran
#ifdef FRONT_CDEPS_DNWM
      use FRONT_CDEPS_DNWM, only: DNWM_SS  => SetServices
#endif
...
#ifdef FRONT_CDEPS_DNWM
          if (trim(model) == "dnwm") then
            call NUOPC_DriverAddComp(driver, trim(prefix), DNWM_SS, &
              petList=petList, comp=comp, rc=rc)
            found_comp = .true.
          end if
#endif
```

### 3.7 fd_ufs.yaml field definition

```yaml
    - standard_name: river_volume_flux
      canonical_units: m3 s-1
      description: NWM (dnwm) data-river export - volumetric river discharge per SCHISM source element (one-way NWM -> SCHISM forcing)
```

### 3.8 CMEPS mediator route (`mediator/esmFldsExchange_coastal_mod.F90`)

```fortran
! advertise (rof -> ocn); register the custom field on the mediator PETs
if (coastal_attr%rof_present .and. coastal_attr%ocn_present) then
  if (.not. NUOPC_FieldDictionaryHasEntry('river_volume_flux', rc=rc)) then
     call NUOPC_FieldDictionaryAddEntry('river_volume_flux', canonicalUnits='m3 s-1', rc=rc)
  end if
  call addfld_from(comprof, 'river_volume_flux')
  call addfld_to(compocn, 'river_volume_flux')
end if

! map + merge (rof -> ocn): mapfcopy = redistribution (no area weighting), correct
! for per-element volumetric m3/s under the element-id contract
if (coastal_attr%rof_present .and. coastal_attr%ocn_present) then
  fldname = 'river_volume_flux'
  if (fldchk(is_local%wrap%FBExp(compocn),trim(fldname),rc=rc) .and. &
      fldchk(is_local%wrap%FBImp(comprof,comprof),trim(fldname),rc=rc)) then
     call addmap_from(comprof, trim(fldname), compocn, mapfcopy, coastal_attr%mapnorm, 'unset')
     call addmrg_to(compocn, trim(fldname), mrg_from=comprof, mrg_fld=trim(fldname), mrg_type='copy')
  end if
end if
```

---

## 4. Build

Add `-DUSE_NUOPC_RIVER=ON` to the coastal SCHISM build flags:

```
-DAPP=CSTLS -DUSE_ATMOS=ON -DNO_PARMETIS=OFF -DOLDIO=ON -DUSE_NUOPC_RIVER=ON
```

With `-DUSE_NUOPC_RIVER` OFF the **SCHISM side** is unchanged: every SCHISM core/cap code
path added for river forcing is `#ifdef USE_NUOPC_RIVER`-gated, so the existing coastal
`atm2sch` regression behaves exactly as before. The build is **not** byte-identical to a
pre-NWM tree, however: the `dnwm` component and its driver registration are gated by
`if(CDEPS)` (CMakeLists.txt: `FRONT_CDEPS_DNWM=cdeps_dnwm_comp`), **not** by
`USE_NUOPC_RIVER`, so they compile into any CDEPS-enabled coastal build. dnwm is simply
never instantiated unless a run sequence adds an `NWM`/`ROF` component, so the runtime
behavior of existing configurations is unaffected.

**Pull (Hercules):**

```bash
git clone -b feature/nwm_river_forcing https://github.com/mansurjisan/ufs-weather-model-mjisan.git
cd ufs-weather-model-mjisan
git submodule update --init --recursive     # NOT --depth 1 (CMEPS commit unreachable shallow)
```

---

## 5. Runtime configurations

Three `ufs.configure` variants are provided:

- **Stub** (`...coastal_datm_ocn_nwm.IN`): `if_source=1` + `source_sink.in` + OCN attribute
  `river_stub_q`, NO dnwm component. Validates the SCHISM receiving side alone.
- **Connected / direct connector** (`...coastal_datm_ocn_nwm_cdeps.IN`):
  `EARTH_component_list: ATM OCN NWM MED`, `NWM_model=dnwm`, `OCN` attribute
  `nwm_coupling=true`, runseq `NWM -> OCN :remapMethod=redist`.
- **Mediator route** (`...coastal_datm_ocn_nwm_mediator.IN`): NWM as the ROF role,
  runseq `ROF -> MED :remapMethod=redist`, `MED med_phases_post_rof`.

dnwm needs `dnwm_in` (namelist: `model_meshfile`, `nx_global`, `datamode=copyall`) and
**both** `dnwm.streams` and `dnwm.streams.xml` staged (dnwm appends `.xml`). The discharge
file's time epoch must cover the model run date.

---

## 6. Experiments performed

All validation was done in a local Docker container (`ufsnwm:built`, ESMF 8.8.0,
OpenMPI 4.1.6) because the full DATM+era5 atmosphere OOMs a 12 GB box. The NWM coupling
is independent of the atmosphere, so a memory-fitting **OCN + dnwm** configuration
(SCHISM tides + river, no ATM, built `USE_ATMOS=OFF` to allow `nws=0`) was used for the
connected-path science. Idealized **carved-channel** runs demonstrate channelized flow.

### 6.1 Compile + stub
- Full `ufs_model` builds with `-DUSE_NUOPC_RIVER=ON` (SCHISM core + cap + dnwm + driver).
- STUB run: `if_source=1` + 6 sources + constant `river_stub_q=100`. SCHISM ran a full
  hour, exit 0, `river_stub_q` logged on all OCN PETs — validates the receiving side.

### 6.2 Connected direct-connector run (the key coupling validation)
Resolved a cascade of bring-up blockers (see Section 7), then:
```
SCHISM_ImportRiver: CONNECTED to NWM provider; nsources=6, vsource(1) [m^3/s] = 100.000   (all OCN PETs)
dnwm_comp_run :ES: river_volume_flux  0.0  100.0  400.0                                    (reads the FILE, not the stub)
EXIT=0 ; PROGRAM ufs HAS ENDED ; fatal.error empty ; schout written per rank
```
Exercises the full chain: dnwm reads discharge -> connector redistributes by element id ->
cap fills `ath3` -> vsource -> SCHISM timesteps -> output. **CONNECTED, not the stub.**

### 6.3 With-river vs without-river water level (idealized carved channel)
A river channel was carved into the test mesh (DEPTH-ONLY edit of `hgrid.gr3`;
connectivity untouched, so decomposition / ESMFmesh / element-id contract stay valid). The
channel runs from a (formerly dry-land) head down to the shelf; the head element is the
single river source. Matched runs (identical except discharge):

| Run | discharge | result |
|---|---|---|
| 1 h, 300 m3/s point source | 300 | head WL rises monotonically; difference = 0 at t=0 -> discharge is the cause; effect local |
| 1 h, channel, 1000 m3/s | 1000 | head +46 cm, mid/mouth ~0 (channel conveys flow) |
| 3 h, channel | 1000 | head clean +8.5 cm (after fully-wet-head fix); mid/mouth ~0 over the full 3 h |

**Conclusion:** the river produces a localized backwater at the constricted head; the
channel conveys the discharge seaward without ponding, so no stage signal accumulates
downstream. An off-channel control showed ~0, confirming the effect is localized, not a
domain-wide artifact. Discharge conservation check: head conveyed flow (9.9 m2/s) x channel
width (~1.26 km) ~= 993 m3/s, recovering the injected 1000 m3/s.

### 6.4 Tide-current interaction (9 h run, ~most of the M2 cycle)
At the channel head: the tide alone drives ~0 current (dead-end head), while with river a
steady 0.9-1.5 m/s seaward current appears and is **tidally modulated** (strongest near
low water / ebb, weakest near high water). With-river WL rides ~5 cm above the tidal curve
through the whole cycle. At mid-channel, the conveyed flow is the tidal oscillation lifted
by a steady river offset (tide + river superposition).

---

## 7. Debugging journey (blockers found and fixed)

The connected path required clearing a cascade of NUOPC/CDEPS bring-up issues. Recorded
here because they are the non-obvious parts and will recur for similar data components.

1. **Field dropped at realize.** The SCHISM cap is legacy IPDv00; its single realize pass
   ran `SCHISM_RemoveUnconnectedFields` before the connector connected `river_volume_flux`.
   Fix: `nwm_coupling=true` OCN attribute exempts the field from removal.
2. **`SCHISM_ImportRiver` aborted** "no ESMF_Field found" when the field was absent. Fix:
   presence (`itemType`) + `NUOPC_IsConnected` check with a graceful stub fallback.
3. **dnwm could not open `dnwm.streams.xml`** (it appends `.xml`). Fix: stage the `.xml`.
4. **`shr_stream_findbounds` crash** — discharge time epoch (2008) didn't cover the model
   date (2012). Fix: repoint the discharge `time:units` to the model start.
5. **dnwm `dshr_state_SetScalar` "Bad Object"** — the direct connector doesn't consume
   `cpl_scalars`, so it stays unrealized. Fix: guard SetScalar with `NUOPC_IsConnected`.
6. **Run-phase timestamp incompatibility** — NUOPC default `CheckImport` requires the
   import field at currTime; the direct connector has no mediator to broker time. Fix: the
   cap specializes `label_CheckImport` to re-check every connected import at currTime while
   **exempting only** `river_volume_flux` (one-way zero-order hold), and dnwm stamps its
   export with `NUOPC_SetTimestamp`. (Originally a blanket no-op; narrowed per review #1.)
7. **Field-dictionary error (mediator route)** — `river_volume_flux is not a StandardName`
   on the MED PETs (which don't host the dnwm/SCHISM caps that register it). Fix: the
   mediator adds the entry (guarded) before advertising.
8. **Source-cell drying/wetting instability** — a source element with a dry node went
   unstable under sustained discharge (WL plunged to -3 m and oscillated). Fix: a river
   source element must be FULLY wet (all nodes deep). Carve wider / pick a fully-wet head.

**Two corrections made during the work (recorded for honesty):**
- An early "carved channel" was a near no-op: the carve only *deepened* nodes, but the
  chosen axis was already deep water, so it changed almost nothing. Fix: OVERWRITE the
  depth along a corridor that crosses real land, and always diff vs the original to confirm
  the carve took.
- An early "transfer-realize segfault" attributed to the connector was actually DATM
  crashing on the era5 mesh read (a container-memory artifact), not the NWM code.

---

## 8. Scope and limitations (important)

- **This is a synthetic per-element forcing prototype, not native NWM ingestion.** What is
  validated is the COUPLING MECHANISM (file -> dnwm -> redist -> `ath3` -> SCHISM
  continuity) and SCHISM's physical response. The hydrology pipeline a production system
  needs is NOT built: NWM reach -> element aggregation and placement, sink handling, and a
  full-mesh time-varying stream all remain. Phase 0 here is a synthetic stand-in; the real
  pairing is `NWM_coupling/coupling_nwm.f90`.
- **Volume + freshwater only, no river momentum.** The exchange injects a volume source
  (`vsource`, m3/s) with fresh/ambient tracers (S=0, T=-9999 ambient); it does NOT impart
  inflow momentum. Plume/jet and momentum-driven dynamics are therefore NOT represented or
  validated.
- **Sinks are not supported.** The cap fills only `vsource`; `source_sink.in` must have
  `nsinks=0`. The cap now ABORTS when `nsinks>0` rather than silently dropping the
  withdrawals.
- **One-way only, by construction.** dnwm is a CDEPS DATA model (reads file, exports); it
  has no inlet for SCHISM state. Two-way (backwater feedback) is NOT built and would need
  the live NWM model as a NUOPC component or a bidirectional cap.
- **The validation source locations are SYNTHETIC.** The Duck/"shinnecock" test mesh has
  no river boundaries and ships no NWM reach data (`NWM_shp_ll.nc` is a git-LFS stub), so
  a hydrologically-correct source set cannot be built for it.
- **Real NWM forcing is not yet exercised.** A constant/idealized discharge was used. A
  *ramping* hydrograph is needed to confirm coupling-interval timing (see the dnwm
  currTime fix, section 9.5) — a constant discharge cannot reveal a one-interval offset.
- **The full atmosphere case (DATM+era5, USE_ATMOS=ON, nws=4) was not run locally** — only
  the memory limit prevented it; the river path is independent of `USE_ATMOS`.

---

## 9. Remaining work

1. Run the full DATM + SCHISM + dnwm connected case on Hercules (ample memory).
2. Real NWM discharge via `nos_utils.forcing.nwm.NWMProcessor` on an operational mesh
   (SECOFS / STOFS-3D) that ships river boundaries + a RiverConfig; derive the dnwm
   per-element `discharge.nc` from its `vsource.th` + `source_sink.in` element order.
3. RT auto-staging: `export_dnwm_cdeps` in `tests/default_vars.sh` + `fv3_conf` (must emit
   `dnwm.streams.xml`); uncomment the `atm2sch_nwm` entries in `tests/rt_coastal.conf`.
4. Finish the CMEPS mediator route (coastal-mode rof fractions + `med_phases_post_rof_init`
   so it passes the realize phase) for an upstream-clean PR toward `oceanmodeling#161`.
5. **Timing + contract validation.** (a) Validate the dnwm `currTime` timing fix with a
   *ramping* discharge run on Hercules — a constant discharge cannot expose a one-interval
   offset. (b) The offline element-id gate `logs/verify_id_contract.py` now exists; add the
   stronger in-framework round-trip (dnwm exports its element ids; the cap asserts them
   against `ielg` after the first import) so the contract is enforced at run time, not only
   pre-run.

---

## 10. Review fixes applied (2026-05-31)

A code review (tracked in memory `nwm-schism-review-fixes`) flagged nine issues; all were
verified against the tree and addressed. The SCHISM-side changes (#1, #2, #4, #9) were
compile-verified in the Docker container (`make schism` -> `libschism.a`, exit 0); the dnwm
change (#3) compiles standalone against the CDEPS modules. The timing fix (#3) still needs a
*ramping*-discharge run on Hercules to confirm the offset is gone at run time.

| # | Issue | Fix |
|---|---|---|
| 1 | `CheckImportRiver` was a blanket no-op, disabling the timestamp staleness check for ALL imports (ATM/wave/ice) | Re-implemented to re-check every connected import at `currTime`, **exempting only** `river_volume_flux` (cap) |
| 2 | An absent/unconnected river field silently fell back to `river_stub_q` and looked successful | Stub is now opt-in via the `river_stub` attribute; with `nwm_coupling=true` (or no mode set) an unconnected field ABORTS (cap) |
| 3 | dnwm read the stream one step ahead (`nextTime`) on the mediator-less direct route -> hydrograph applied one coupling interval early | dnwm now reads at `currTime`, consistent with its `currTime` export stamp (dnwm); validate with a ramp run |
| 4 | `nsinks` was read but never applied; sinks silently dropped | Cap ABORTS when `nsinks>0` (sinks unsupported by this prototype) |
| 5 | Phase-0 generator deliberately selected partly-dry source elements (`mindep<=0`) -> drying/wetting instability | Generator now requires every source node deeper than `WETMIN` (fully wet); regenerated sources verified all-wet |
| 6 | The element-id contract was a doc "release gate" with no assert | Added `logs/verify_id_contract.py` (count / natural-order / sources / `nsinks==0`), smoke-tested pass + fail; in-framework round-trip noted as follow-up |
| 7 | Branch not reproducible: dnwm/config staging files were untracked | Tracked the staging files (see git history of the superproject fork) |
| 8 | Mediator template used the wrong token `@[coupling_interval_slow]` | Corrected to `@@[coupling_interval_slow_sec]` (matches all other coastal templates) |
| 9 | `SCHISM_MeshCreateElement` read/looped element-id and connectivity export fields to `nea` from `ne`-sized arrays (ghost overread) | Read/loop to `ne` (the mesh is owned-only); consistent with the cap consuming only `%id<=ne` (esmf_util) |

Doc corrections from the same review: removed the "byte-identical with flag OFF" claim
(dnwm/driver are gated by `if(CDEPS)`, not `USE_NUOPC_RIVER`); made the mediator-route
status consistent (scaffolded, not runnable); sharpened the scope (section 8:
synthetic-prototype framing, volume+freshwater only / no momentum, sinks unsupported).

---

## 11. References

- `oceanmodeling/ufs-weather-model#161` — "NWM data ingestion through CDEPS" (the placeholder this implements)
- `schism-dev/schism-esmf#25` — "NWM BMI for SCHISM" (the NextGen/BMI ath3-injection route this mirrors)
- Companion docs in-repo: `NWM_to_SCHISM_Integration_Report.md`, `NWM_RIVER_INCREMENT1_VALIDATION.md`
