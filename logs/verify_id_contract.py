#!/usr/bin/env python3
"""
Release gate for the one-way NWM -> SCHISM redist element-id contract.

The NWM->OCN connector uses remapMethod=redist, which routes element-by-element by
DISTGRID SEQUENCE INDEX, and CMEPS redist SILENTLY IGNORES unmatched indices
(med_map_mod). So if the dnwm element ESMFmesh is not in the SAME order as SCHISM's
global element numbering, discharge lands on the wrong elements (or nowhere) with NO
error at run time. The coupling doc calls an id-contract check a "release gate" but
shipped no assert; this is that assert. Run it BEFORE any coupled run:

  1. mesh elementCount == SCHISM ne (from hgrid)                       [count]
  2. mesh element i connectivity == hgrid element (i+1) connectivity   [order: seqindex == global id]
  3. source_sink.in: nsources>=1, unique ids in [1,ne], nsinks == 0    [sources valid; sinks unsupported]
  4. discharge nelem == ne, nonzero ONLY at the source elements        [forcing placed at the sources]

Exit 0 and print CONTRACT_OK if all pass; exit 1 with a FAIL line otherwise.

Usage:
  verify_id_contract.py [ESMFmesh.nc] [hgrid.ll] [source_sink.in] [nwm_discharge.nc]
(defaults match gen_nwm_phase0.py output)
"""
import sys, numpy as np, netCDF4 as nc

def fail(msg):
    print(f"FAIL: {msg}")
    sys.exit(1)

MESH  = sys.argv[1] if len(sys.argv) > 1 else "/mnt/d/ufs-weather-model/phase0_nwm_out/schism_elem_ESMFmesh.nc"
HGRID = sys.argv[2] if len(sys.argv) > 2 else "/mnt/d/SCHISM_MESH_GEN/RT_DUCK_DATM_SCH/coastal_ike_shinnecock_atm2sch_intel/hgrid.ll"
SRC   = sys.argv[3] if len(sys.argv) > 3 else "/mnt/d/ufs-weather-model/phase0_nwm_out/source_sink.in"
DIS   = sys.argv[4] if len(sys.argv) > 4 else "/mnt/d/ufs-weather-model/phase0_nwm_out/nwm_discharge.nc"

# ---- parse hgrid (the truth: ne and element connectivity in natural/global order)
with open(HGRID) as f:
    f.readline()
    ne, npn = (int(x) for x in f.readline().split())
    for _ in range(npn):
        f.readline()
    i34  = np.empty(ne, dtype=np.int32)
    conn = np.full((ne, 4), -1, dtype=np.int32)        # 1-based node ids, -1 pad
    for e in range(ne):
        t = f.readline().split()
        n = int(t[1]); i34[e] = n
        conn[e, :n] = [int(t[2 + k]) for k in range(n)]
print(f"hgrid: ne={ne}")

# ---- 1 + 2. mesh count and natural order (seqindex == global element id) --------
m = nc.Dataset(MESH)
mne = m.dimensions["elementCount"].size
if mne != ne:
    fail(f"mesh elementCount {mne} != hgrid ne {ne} (count mismatch -> redist would drop the tail)")
ec = m.variables["elementConn"][:]
mconn = ec.filled(-1) if np.ma.isMaskedArray(ec) else np.asarray(ec)
nec = m.variables["numElementConn"][:]
mnec = nec.filled(0) if np.ma.isMaskedArray(nec) else np.asarray(nec)
m.close()
if not np.array_equal(mnec.astype(np.int32), i34):
    fail("numElementConn != hgrid node-count-per-element (element order broken)")
badrows = np.where(np.any(mconn.astype(np.int32) != conn, axis=1))[0]
if badrows.size:
    fail(f"{badrows.size} mesh elements differ from hgrid natural order "
         f"(first at file index {badrows[0]}, global id {badrows[0] + 1}) -> seqindex != global id")
print(f"PASS [count]  mesh elementCount == ne == {ne}")
print(f"PASS [order]  all {ne} mesh elements match hgrid natural order (seqindex == global element id)")

# ---- 3. source_sink.in (sparse: nsources, ids, blank, nsinks) -------------------
lines = open(SRC).read().split("\n")
ns = int(lines[0])
src = [int(lines[1 + i]) for i in range(ns)]
rest = [l for l in lines[1 + ns:] if l.strip() != ""]   # skip the blank separator
nsink = int(rest[0]) if rest else 0
if ns < 1:
    fail("nsources < 1")
if any(s < 1 or s > ne for s in src):
    fail(f"a source id is outside [1,{ne}]")
if len(set(src)) != ns:
    fail("duplicate source ids in source_sink.in")
if nsink != 0:
    fail(f"nsinks={nsink} but the one-way NWM prototype does not support sinks (must be 0)")
print(f"PASS [sources] nsources={ns}, ids unique in [1,{ne}], nsinks=0")

# ---- 4. discharge placed at exactly the source elements -------------------------
d = nc.Dataset(DIS)
dne = d.dimensions["nelem"].size
q = d.variables["river_discharge"][:]
d.close()
if dne != ne:
    fail(f"discharge nelem {dne} != ne {ne}")
nz = set((np.where(np.any(q != 0.0, axis=0))[0] + 1).tolist())   # 1-based elems with any nonzero Q
extra = sorted(nz - set(src))
missing = sorted(set(src) - nz)
if extra:
    fail(f"discharge nonzero at {len(extra)} non-source elements (first: {extra[:10]})")
if missing:
    fail(f"source elements receiving zero discharge: {missing}")
print(f"PASS [forcing] discharge nonzero at exactly the {ns} source elements")

print("CONTRACT_OK: dnwm ESMFmesh order == SCHISM global element ids; sources and discharge consistent")
