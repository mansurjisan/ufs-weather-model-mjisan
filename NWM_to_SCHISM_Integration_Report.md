# Integrating NWM Streamflow as One-Way River Forcing into SCHISM (UFS Coastal)

**Scope:** one-way NWM → SCHISM freshwater/river forcing in `oceanmodeling/ufs-weather-model` (UFS Coastal). Offline preprocessing vs. online in-framework NUOPC coupling, with a recommendation and an implementation path for *this* repo.

**Provenance / confidence:** This report was synthesized from a salvaged deep-research run (workflow `wf_cd8662c5-935`: 103 agents, ~1.78M tokens) that completed its search/fetch/extract phases but died before its Verify→Synthesize phases finished. Most findings carry *prose* verification reasoning; a subset carry formally-emitted structured verdicts. The two most decision-critical facts (GitHub issues #161 and #25) were **re-verified live via the GitHub API on 2026-05-29**. Confidence is labeled per claim. Data-quality caveats the verifiers caught (dead URLs, paraphrased "quotes," figure overreaches) are listed in the Appendix — read it before quoting specifics externally.

---

## TL;DR — Recommendation

**Use the OFFLINE path: generate SCHISM-native source/sink files from NWM with `pyschism`, and turn them on with `if_source` in `param.nml`.** It is the conservation-correct, operationally-proven, near-term-feasible mechanism, and it is already tooled in this very repo. **[Confidence: HIGH]**

Do **not** route NWM through CDEPS/CMEPS for one-way forcing today: the online plumbing does not exist (only an unstarted placeholder issue), and ESMF/NUOPC regridding is the wrong primitive for reach-based point discharge. Reserve an online component for when you need **two-way** coupling (tidal/surge backwater into the routing) — which is out of scope for one-way forcing.

---

## Part 1 — What NWM actually gives you

NWM output is **not** a gridded runoff field; it is a **1-D, reach-based point dataset**:

| Property | Value | Source / confidence |
|---|---|---|
| Channel output file | `CHRTOUT` / `channel_rt`, **Point Type** | NOAA NODD `nodd-data-docs/nwm/README.md` — *"Channel Output in the CHRTOUT files (Point Type)" / "Streamflow values at points associated with flow lines"* **[HIGH]** |
| Variable & unit | `streamflow`, **m³ s⁻¹** (volumetric discharge) | NODD README — *"streamflow \| River Flow (m³ s⁻¹)"* **[HIGH]** |
| Network size | ~**2.7 million** river reaches (CONUS) | water.noaa.gov/about/nwm — *"streamflow for over 2.7 million river reaches"* **[HIGH]** |
| Reach key / geometry | `feature_id` (NHDPlus COMID); reach lat/lon in separate `Route_Link.nc` | NODD docs + searcher **[HIGH]** |
| Cadence | hourly (`CHRTOUT` published "Every Hour") | NODD README **[HIGH]** |

**Why this matters:** `streamflow` in **m³/s** is exactly what SCHISM's `vsource` wants, and is *fundamentally different* from the **kg m⁻² s⁻¹ area-flux** that CDEPS DROF and ESMF bilinear/conservative regrid assume. This unit/representation mismatch is the crux of the whole offline-vs-online decision (Part 4). **[HIGH]**

### Access (for driving a forecast or hindcast)
- **Real-time feed:** `s3://noaa-nwm-pds` (us-east-1), **free, anonymous** (`aws s3 ls --no-sign-request s3://noaa-nwm-pds/`). Live anonymous read confirmed during research. **[HIGH]**
  - File naming: `nwm.tCCz.<config>.channel_rt[_memN].fFFF.<domain>.nc`. Also mirrored on GCP and NOMADS (48 h rolling); an SNS new-data topic + Kerchunk/Zarr companion support low-latency automation.
  - **Caveat:** the AWS registry *describes* a "rolling four-week archive," but a live bucket listing showed `IsTruncated=false` with substantially longer retention — so "data ages out in 4 weeks, must pull every cycle" is **overstated**. The real-time constraint is per-cycle *latency*, not aging-out. **[MEDIUM — verifiers split]**
- **Retrospective / hindcast:** the NWM CONUS Retrospective Dataset (Feb 1979–present) is the streamflow source the official SCHISM docs use with `if_source=1`; a separate retrospective archive bucket exists for reproducible regression testing. **[HIGH]**

### Forecast configurations (cadence / lead time for one-way forcing)
| Config | Cycle | Horizon | Notes |
|---|---|---|---|
| analysis_assim | hourly | 3–28 h lookback | nowcast |
| short_range | hourly | 18 h (CONUS), 48 h (HI/PR) | HRRR/RAP-forced |
| medium_range | 4×/day | ~10 days | GFS-forced ensemble |
| long_range | every 6 h | 30 days | 16-member ensemble |

*(Searcher-verified against AWS registry + water.noaa.gov; structured output submitted. **[HIGH]**)*

---

## Part 2 — The two mechanisms

### Mechanism A — OFFLINE (pyschism → SCHISM native source/sink) ✅ recommended

SCHISM ingests rivers through its own **source/sink** machinery, gated by `if_source` in `param.nml`:

- `if_source = 0` — off (current state of all coastal tests, e.g. `coastal_ike_shinnecock_atm2sch`). **[HIGH]**
- `if_source = 1` — **ASCII**, requires four files: `source_sink.in` (element-ID lists: sources block, then sinks block), `vsource.th` (volume inflow, m³/s), `vsink.th` (volume outflow), `msource.th` (tracer/mass = T, S, …). **[HIGH]**
- `if_source = -1` — **NetCDF** `source.nc` (single-file equivalent). **[HIGH, searcher-confirmed]**
- `msource` ordering is `T, S, <tracers>`; the sentinel **`-9999.` injects ambient concentration**; pyschism defaults rivers to **salinity = 0.0 (freshwater), temperature = -9999.0 (ambient)**. **[HIGH, source-code-confirmed]**

**The tooling (`pyschism.forcing.source_sink.nwm`):**
- `NWMElementPairings(hgrid)` — geometrically intersects NWM reaches (LineStrings keyed by `feature_id`) with the SCHISM mesh hull / land boundary using an **STR-tree** spatial index, then assigns each reach to the **nearest element centroid via KDTree**. **Flow direction decides role:** reach flowing *into* the domain → **source**; flowing *out* → **sink**. Pairings are computed once per mesh and cached as `sources.json` / `sinks.json`. **[HIGH, source-code-confirmed; the exact inline-comment "quote" the run cited could not be byte-verified — treat wording as paraphrase]**
- `NationalWaterModel(pairings=...).write(outdir, hgrid, startdate, rnday)` — reads `streamflow`/`feature_id` from `CHRTOUT` and writes the four files for a given window. Driven by `startdate` + `rnday`, cached per date → **re-run per forecast cycle**. **[HIGH]**
- `AWSDataInventory` (and subclasses `AWSForecastInventory` / `AWSHindcastInventory` / `GOOGLEHindcastInventory`) — pulls the NWM data. **[HIGH]**

**Mesh prerequisite — RiverMapper / RiverMeshTools** (Ye et al. 2023, *Env. Modelling & Software* 166:105731, DOI 10.1016/j.envsoft.2023.105731): a parallel Python toolchain (`pyDEM` extracts thalwegs → 1-D network; `RiverMapper` places river "arcs" to guide mesh generation) that **resolves river channels in the unstructured grid at continental scale** and explicitly **accepts NWM flowlines** as input. This is the *mesh-side* enabler (where rivers physically live in the grid) — distinct from, but a prerequisite for accurate, source/sink *forcing*. Used to build the STOFS-3D-Atlantic East/Gulf domain. **[HIGH — DOI confirmed via Crossref. NOTE: correct title is "A parallel Python-based tool for meshing watershed rivers at continental scale"; the research prompt's title was a paraphrase.]**

**In-repo assets you already have:**
- `SCHISM-interface/SCHISM/src/Utility/Pre-Processing/NWM/` (incl. `NWM_coupling`, `Sflux2Source` subdirs) and `gen_sourcesink_nwm.py` / `gen_sourcesink.py` (both contain `streamflow = nc['streamflow'][:]` / `feature_id = nc['feature_id'][:]`). **[HIGH, in-tree]**
- `RiverMapper` lives at `SCHISM-interface/SCHISM/src/Utility/Grid_Scripts/Compound_flooding/RiverMapper/`. **[HIGH, in-tree]**

### Mechanism B — ONLINE (NUOPC data-river + CMEPS routing into SCHISM) ❌ not viable today

To force SCHISM from NWM *through the UFS coupling layer*, you would need all of:
1. **Add a runoff/river import field** to the SCHISM NUOPC cap (`SCHISM-interface/SCHISM-ESMF/src/schism/schism_nuopc_cap.F90`, which today advertises ~39 atmosphere/wave/ice fields and **no** runoff field) and wire it to SCHISM's source/sink arrays. **[HIGH, in-tree grounding]**
2. **Add rof→ocn routing** to the coastal mediator (`CMEPS-interface/CMEPS/mediator/esmFldsExchange_coastal_mod.F90`), which today handles only ATM/OCN/WAV/ICE. The runoff machinery exists **only** in the CESM path: `esmFldsExchange_cesm_mod.F90` uses `addmap_from(comprof,'Forr_rofl',compocn,mapconsd,...)` (a **conservative area remap**) and `addmrg_to(compocn,'Foxx_rofl',...,mrg_type='sum')`. Stock CMEPS has **no coastal driver** at all. **[HIGH, source-code-confirmed]**
3. **Add a data-river component** (CDEPS DROF or bespoke) to a `ufs.configure.coastal_*` run sequence — none exists in any coastal config. **[HIGH]**
4. **Solve the representation mismatch** — turn reach-based m³/s point discharge into something the ESMF regridder can move conservatively onto SCHISM elements (Part 4).

---

## Part 3 — Honest pros/cons comparison

| Criterion | OFFLINE (pyschism source/sink) | ONLINE (NUOPC data-river + mediator) |
|---|---|---|
| **Exists today?** | ✅ Yes — tooled in-repo, operational in STOFS-3D | ❌ No — only placeholder issue #161, zero code |
| **Physical fidelity** | ✅ Volume injected at correct land-boundary elements (mesh-resolved) | ⚠️ ESMF regrid smears a sharp point source; needs nearest-neighbor (non-conservative) hacks |
| **Mass conservation** | ✅ Inherent (volume in = volume specified) | ⚠️ Conservative regrid is area-based; doesn't fit m³/s points |
| **Unit match to NWM** | ✅ m³/s → `vsource` directly | ❌ Mediator expects kg/m²/s area flux |
| **Real-time feasibility** | ✅ Per-cycle file gen is scriptable (STOFS does it daily) | — (would be in-line once built) |
| **Dev cost / maintenance** | ✅ Low — config + a preprocessing script | ❌ High — cap field + mediator route + new component + mapping |
| **UFS regression testing** | ✅ Fits existing DATM+SCHISM RT pattern | ❌ No ROF/DROF RT exists to extend |
| **Alignment w/ where NOAA is heading** | ✅ Active path (ufs-coastal-app #8, STOFS-3D-Alaska #23) | ⚠️ The *online* NWM↔SCHISM effort is going **BMI/NextGen** (#25), not UFS-NUOPC |
| **Two-way (backwater) capable** | ❌ One-way only | ✅ The reason you'd eventually build online |
| **Latency coupling to atmosphere** | ⚠️ Offline, not in-step with the coupled clock | ✅ In-step |

---

## Part 4 — The technical crux: point discharge vs. area flux

This is *why* offline wins for one-way, beyond mere convenience:

- **ESMF `LocStream` (disconnected points) can only be a regrid *destination*, never a source.** NWM reaches are naturally a point cloud, so they **cannot** be the source of a conservative regrid onto the SCHISM mesh. **[HIGH — ESMF Reference Manual]**
- **Conservative regridding is cell/area-based** — it preserves ∫(value × area) and needs cell corners/areas. A volumetric point discharge (m³/s) has no cell geometry, so it doesn't fit the conservative model. **[HIGH]**
- **Nearest-source-to-destination** is the *only* method that maps points one-to-one to nearest neighbors (structurally like reach→element) — but it is **explicitly non-conservative**. **[HIGH]**
- **The CESM/MOSART precedent confirms the paradigm gap:** MOSART routes on a regular **0.5° grid** and the coupler "automatically re-routes river runoff to the nearest ocean grid cell" as a gridded flux — the inverse of SCHISM's native point-source ingestion of NWM m³/s. CDEPS DROF inherits this gridded-runoff assumption. **[HIGH — CLM5/CTSM Tech Note]**

**Net:** routing NWM through ESMF/CMEPS would force you to fabricate a runoff mesh from reach points and convert m³/s → per-area flux, then spread a sharp point source — hydrodynamically wrong for estuary-driving rivers. The offline source/sink path sidesteps all of it and is inherently volume-conserving. **[HIGH]**

---

## Part 5 — Recommendation

**Adopt the offline pyschism source/sink path for one-way NWM → SCHISM forcing.** Rationale: it's the only mechanism that (a) exists today, (b) is physically/conservation correct for point discharge, (c) matches NWM's m³/s reach data, (d) is operationally proven (STOFS-3D-Atlantic), and (e) aligns with the active near-term direction in the UFS Coastal and SCHISM communities. **[HIGH]**

Build an online component **only** when you need two-way backwater feedback — and even then, the community signal (#25) is that NWM↔SCHISM coupling is being built on **BMI/NextGen**, a different track from UFS CMEPS/NUOPC. **[HIGH]**

---

## Part 6 — Concrete implementation path in THIS repo (offline)

Target a coastal test case (e.g. the Shinnecock domain) and add NWM rivers:

1. **Confirm the mesh resolves your rivers.** If channels aren't represented, run `RiverMapper` (`.../Grid_Scripts/Compound_flooding/RiverMapper/`) with NWM flowlines + a DEM to generate river arcs before (re)meshing. For Shinnecock's small domain this may be minimal.
2. **Compute reach→element pairings (once per mesh):**
   ```python
   from pyschism.mesh import Hgrid
   from pyschism.forcing.source_sink.nwm import NationalWaterModel, NWMElementPairings
   hgrid = Hgrid.open('hgrid.gr3', crs='epsg:4326')
   pairings = NWMElementPairings(hgrid)
   pairings.save_json(sources='sources.json', sinks='sinks.json')   # cache
   ```
3. **Generate the forcing files for your window:**
   ```python
   nwm = NationalWaterModel(pairings=pairings, cache=cache)
   nwm.write(output_directory, hgrid, startdate, rnday, overwrite=True)
   # -> source_sink.in, vsource.th, vsink.th, msource.th
   ```
   (Or use the in-repo `Utility/Pre-Processing/NWM/gen_sourcesink*.py` if you prefer the bundled scripts over the pip package.)
4. **Turn it on in `param.nml`:** set **`if_source = 1`** (ASCII) — or `-1` if you generate `source.nc`. Verify `msource.th` uses `S=0.0` and `T=-9999.` (ambient) unless you have river temperatures.
5. **Place the four files** in the run directory alongside the existing inputs; confirm SCHISM reads them at init (`schism_init.F90`) and applies them per step (`schism_step.F90`).
6. **Wire a regression test:** clone the `coastal_ike_shinnecock_atm2sch` pattern into a new `..._riv` case that flips `if_source` and ships the source/sink files (or generates them from a small cached NWM slice). This keeps it reproducible in the UFS RT framework — the one piece STOFS does operationally that the UFS test suite doesn't yet cover.
7. **For real-time:** script steps 3–4 per cycle, pulling from `s3://noaa-nwm-pds` with `--no-sign-request`, keyed to the appropriate config (short_range 18 h, etc.). Pairings (step 2) are reused.

**Effort:** config + preprocessing wiring + one RT case. No Fortran/cap/mediator changes. **[HIGH]**

---

## Part 7 — What an online path would require (future / two-way)

If/when two-way coupling is justified, the in-framework work is: add a source/sink import field to the SCHISM cap; add a rof→ocn route to the coastal mediator (you can't just switch to CESM mode without breaking coastal ATM/WAV/ICE exchanges); add a reach-aware data-river component (a CDEPS "NWM data mode" per #161, or a bespoke point-source cap); and solve point→element mapping with `mapnstod` (nearest-source-to-destination) rather than conservative regrid. Note the community is pursuing this as a **SCHISM BMI under NextGen** (#25: *"USE_NWM_BMI goes with if_source/=0"* — injecting into SCHISM's internal vsource/msource arrays, bypassing NetCDF files), not as a UFS NUOPC data-river. A real operational hazard flagged there: **cold-start** (no source/boundary data at t0). **[HIGH]**

---

## Part 8 — Operational & community landscape

- **STOFS-3D-Atlantic** (NOAA's operational 3D coastal system, SCHISM core; operational since **Jan 2023**, daily 12 UTC, 24 h nowcast + 96 h forecast; ~2.9M nodes / 5.65M elements; ~8 m river resolution) already does **one-way NWM forcing via native source/sink**. This is the production proof of the offline recommendation. **[HIGH; the specific "6000+ rivers" figure is widely cited but was NOT on the AWS registry page — treat as MEDIUM]**
- **NWM v3.0** added a coastal-routing **linkage to SCHISM** for flood guidance. True and operational — **but inside NWM's NextGen framework, not the UFS/oceanmodeling NUOPC framework** this question targets. **[HIGH — water.noaa.gov verbatim]**
- **`oceanmodeling/ufs-weather-model` #161** — "Placeholder issue: NWM data ingestion through CDEPS," OPEN, opened **2025-04-10** by janahaddad, **0 comments**, points only at offline `gen_sourcesink.py`. The UFS online path is **un-started**. **[HIGH — live-verified 2026-05-29]**
- **`schism-dev/schism-esmf` #25** — "Development of a NWM BMI for SCHISM," OPEN, opened **2023-04-21** by josephzhang8, **11 comments** (with NOAA-OWP's Jason Ducker). The real online NWM↔SCHISM effort is **BMI/NextGen**. **[HIGH — live-verified 2026-05-29]**
- **`ufs-coastal-app` #8** (maintain pyschism incl. NWM source/sink) and **STOFS-3D-Alaska-Dev #23** (rivers → `source.nc`, GLOFAS fallback outside CONUS) confirm offline reach→element generation is the active near-term path. **[MEDIUM — single-source, not re-verified live]**
- **Who to talk to:** uturuncoglu (CDEPS/NUOPC), josephzhang8 (SCHISM), saeed-moghimi-noaa (NOS coastal), janahaddad (#161), jduckerOWP / Jason Ducker (OWP BMI), platipodium / Carsten Lemmen (Hereon).
- **Two-way reference:** Zhang et al. 2024, *Hydrology* 11(9):145 (DOI 10.3390/hydrology11090145) — a two-way NWM-SCHISM research coupling (Hurricane Matthew 2016); bounds where online coupling is justified (it adds tidal/surge backwater that one-way omits). **[MEDIUM — bibliographic only; full text paywalled]**

---

## Appendix — Source quality & data caveats (read before quoting externally)

The verification phase did not formally complete, and several agents initially emitted *fabricated* quotes that were later self-corrected. Known issues:

- **Dead/wrong URLs:** `schism-dev.github.io/.../tutorials/compound-flooding.html` is a real **404** — the canonical page is `.../getting-started/pre-processing-with-pyschism/nwm.html` ("Source and Sink"). The Frontiers STOFS-3D DOI `10.3389/fmars.2023.1199463` did not resolve. The bucket `noaa-nws-stofs3d-pds` does not exist; the real one is `noaa-nos-stofs3d-pds`.
- **Paraphrase risk:** Some "verbatim quotes" were summarizer paraphrases (the pyschism inline source/sink comment; several Ye et al. abstract lines). The *facts* are corroborated across mirrors/code; the exact wording is high-confidence-not-byte-verified.
- **Title correction:** Ye et al. 2023 EMS 166:105731 = "A parallel Python-based tool for meshing watershed rivers at continental scale" (confirmed via Crossref), not the prompt's paraphrased title.
- **Overreaches downgraded:** "noaa-nwm-pds holds *only* 4 weeks" (live listing showed longer retention) and "6000+ rivers" in STOFS (not on the cited registry page) → MEDIUM.
- **Verification status:** facts labeled HIGH are either live-verified, code-grounded in this tree, or had ≥1 formally-submitted structured verdict plus corroboration. MEDIUM = single-source or prose-only reasoning without independent live confirmation. No claim here is contradicted by the evidence; the open question is wording precision on a few quotes, not substance.

*Full salvaged research digest: `logs/nwm_research_digest.md` (25 distilled claims + 21 source extractions + 5 search syntheses).*
