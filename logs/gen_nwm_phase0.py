#!/usr/bin/env python3
"""
Phase 0 generator for one-way NWM -> SCHISM coupling (Duck/"shinnecock" test mesh).

Produces, from the SCHISM hgrid.ll:
  1. <out>/schism_elem_ESMFmesh.nc  -- element ESMF mesh in NATURAL element order,
     so the ESMF file-order sequence index == SCHISM global element id. This is
     the redist contract artifact: dnwm exports on this mesh and the NWM->OCN
     redist connector then routes element-id -> element-id into SCHISM.
  2. <out>/source_sink.in           -- sparse source list (wet shoreline elements),
     SCHISM format (nsources, elem ids, blank, nsinks).
  3. <out>/nwm_discharge.nc         -- per-element discharge [m3/s] (time, nelem)
     on the element mesh, 0 except at source elements (constant test value), with
     an hourly time axis over the run window -- the dnwm stream forcing.

This is a SYNTHETIC pairing (no live NWM streamflow available locally): real
production pairing would run NWM_coupling/coupling_nwm.f90 against NWM_shp_ll.nc +
NWM CHRTOUT. The element mesh + ordering are exact; the source set + discharge are
representative so the connected path can be run and verified end-to-end.
"""
import sys, numpy as np, netCDF4 as nc

HGRID   = sys.argv[1] if len(sys.argv) > 1 else \
    "/mnt/d/SCHISM_MESH_GEN/RT_DUCK_DATM_SCH/coastal_ike_shinnecock_atm2sch_intel/hgrid.ll"
OUTDIR  = sys.argv[2] if len(sys.argv) > 2 else "/mnt/d/ufs-weather-model/phase0_nwm_out"
NSRC    = 6          # number of source elements to select
QTEST   = 100.0      # constant test discharge per source element [m3/s]
START   = "2008-08-23 00:00:00"   # run start; hourly axis spanning the 24h window+
NT      = 49         # hourly steps (covers a 48h window, > the 24h run)

import os; os.makedirs(OUTDIR, exist_ok=True)

# ---- parse hgrid.ll ---------------------------------------------------------
with open(HGRID) as f:
    f.readline()                                  # comment
    ne, np_ = (int(x) for x in f.readline().split())
    xlon = np.empty(np_); ylat = np.empty(np_); dp = np.empty(np_)
    for i in range(np_):
        t = f.readline().split()
        xlon[i] = float(t[1]); ylat[i] = float(t[2]); dp[i] = float(t[3])
    i34  = np.empty(ne, dtype=np.int32)
    conn = np.full((ne, 4), -1, dtype=np.int32)   # 1-based node ids, -1 pad
    for e in range(ne):
        t = f.readline().split()
        n = int(t[1]); i34[e] = n
        conn[e, :n] = [int(t[2+k]) for k in range(n)]   # keep 1-based
print(f"parsed hgrid: ne={ne} np={np_}")
assert ne == 45167 and np_ == 23018, "unexpected mesh size"

# ---- element centers + depths (mean of member nodes) ------------------------
cx = np.empty(ne); cy = np.empty(ne); cdep = np.empty(ne); mindep = np.empty(ne)
for e in range(ne):
    nd = conn[e, :i34[e]] - 1                      # 0-based node idx
    cx[e] = xlon[nd].mean(); cy[e] = ylat[nd].mean()
    cdep[e] = dp[nd].mean(); mindep[e] = dp[nd].min()

# ---- select source elements: wet (mean depth>0) shoreline (>=1 node depth<=0),
#      spread across the domain by longitude -------------------------------
wet_shore = np.where((cdep > 0.0) & (mindep <= 0.0))[0]
if wet_shore.size < NSRC:                          # fallback: wettest elements
    wet_shore = np.argsort(cdep)[::-1][:max(NSRC, 50)]
order = wet_shore[np.argsort(cx[wet_shore])]
pick  = order[np.linspace(0, order.size - 1, NSRC).round().astype(int)]
src_eid = np.sort(np.unique(pick)) + 1             # -> 1-based global element ids
print(f"selected {src_eid.size} source elements (1-based global ids): {list(src_eid)}")
for eid in src_eid:
    e = eid - 1
    print(f"   elem {eid:6d}  center=({cx[e]:.4f},{cy[e]:.4f})  meandep={cdep[e]:.2f} m")

# ---- 1. element ESMF mesh (natural order => seqindex == global elem id) -----
mpath = f"{OUTDIR}/schism_elem_ESMFmesh.nc"
m = nc.Dataset(mpath, "w", format="NETCDF4")
m.createDimension("nodeCount", np_)
m.createDimension("elementCount", ne)
m.createDimension("maxNodePElement", 4)
m.createDimension("coordDim", 2)
v = m.createVariable("nodeCoords", "f8", ("nodeCount", "coordDim")); v.units = "degrees"
v[:, 0] = xlon; v[:, 1] = ylat
ec = m.createVariable("elementConn", "i4", ("elementCount", "maxNodePElement"), fill_value=-1)
ec.long_name = "Node indices that define the element connectivity"; ec[:, :] = conn
nec = m.createVariable("numElementConn", "i4", ("elementCount",))
nec.long_name = "Number of nodes per element"; nec[:] = i34
cc = m.createVariable("centerCoords", "f8", ("elementCount", "coordDim")); cc.units = "degrees"
cc[:, 0] = cx; cc[:, 1] = cy
em = m.createVariable("elementMask", "i4", ("elementCount",)); em.units = "none"; em[:] = 1
m.gridType = "unstructured mesh"; m.version = "0.9"
m.inputFile = HGRID; m.note = "element order == SCHISM global element id (redist contract)"
m.close()
print(f"wrote {mpath}  (elementCount={ne}, order=global-id)")

# ---- 2. source_sink.in (sparse; nsources, ids, blank, nsinks=0) -------------
spath = f"{OUTDIR}/source_sink.in"
with open(spath, "w") as f:
    f.write(f"{src_eid.size}\n")
    for eid in src_eid: f.write(f"{eid}\n")
    f.write("\n0\n")
print(f"wrote {spath}")

# ---- 3. per-element discharge NetCDF (dnwm stream forcing) -------------------
dpath = f"{OUTDIR}/nwm_discharge.nc"
d = nc.Dataset(dpath, "w", format="NETCDF4")
d.createDimension("time", None)
d.createDimension("nelem", ne)
tv = d.createVariable("time", "f8", ("time",))
tv.units = f"hours since {START}"; tv.calendar = "standard"
tv[:] = np.arange(NT, dtype="f8")
q = d.createVariable("river_discharge", "f8", ("time", "nelem"), zlib=True)
q.units = "m3 s-1"; q.long_name = "volumetric river discharge per SCHISM element"
arr = np.zeros((NT, ne))
arr[:, src_eid - 1] = QTEST                        # constant test discharge at sources
q[:, :] = arr
d.history = f"synthetic Phase-0 NWM discharge; {QTEST} m3/s at {src_eid.size} source elements"
d.close()
print(f"wrote {dpath}  (time={NT}, nelem={ne}, Q={QTEST} m3/s at sources)")
print("PHASE0_DONE")
