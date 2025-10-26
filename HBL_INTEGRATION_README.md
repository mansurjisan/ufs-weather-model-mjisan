# HBL Integration Branch

**Branch**: `feature/hbl-integration`
**Remote**: https://github.com/mansurjisan/ufs-weather-model-mjisan
**Base**: `feature/coastal_app`
**Date**: October 26, 2025

## Purpose

This branch integrates the URI Hurricane Boundary Layer (HBL) model into UFS-Coastal through CMEPS mediator.

## HBL NUOPC Cap

- **Repository**: https://github.com/mansurjisan/HBL-NUOPC-Interface (private)
- **Location**: `/work2/noaa/nos-surge/mjisan/ufs-weather-model-test/HBL-interface`
- **Status**: Complete and validated with Intel compilers on Hercules HPC
- **Library**: `libhbl_nuopc.a` (1.3MB, Intel-compiled)

## Integration Plan

See `/work2/noaa/nos-surge/mjisan/ufs-weather-model-test/HBL-interface/UFS_COASTAL_INTEGRATION_PLAN.md` for complete details.

## Files to be Modified

1. `CMEPS-interface/CMEPS/mediator/med_internalstate_mod.F90` - Add `comphbl` component
2. `CMEPS-interface/CMEPS/mediator/esmFldsExchange_coastal_mod.F90` - Add HBL field exchanges
3. `CMakeLists.txt` - Link HBL NUOPC library
4. `CMEPS-interface/CMEPS/mediator/med.F90` - Register HBL component

## Workflow

All changes on this branch will:
- Be pushed to `mansurjisan/ufs-weather-model-mjisan`
- NOT affect the upstream `oceanmodeling/ufs-weather-model`
- Be available for pull requests when ready

## Quick Commands

```bash
# Check current branch
git status

# Push changes to personal repo
git push personal feature/hbl-integration

# View integration plan
cat ../HBL-interface/UFS_COASTAL_INTEGRATION_PLAN.md
```

## Next Steps

1. Modify CMEPS mediator files
2. Update build system
3. Test compilation
4. Run coupled test
5. Create pull request (if desired)
