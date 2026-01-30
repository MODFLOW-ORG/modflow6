# MODFLOW 6.7.0 Particle Tracking (PRT) Model Migration Guide

This document describes how to migrate PRT models from MODFLOW 6.6.x to 6.7.0, as well as significant changes in behavior.

This document is not a comprehensive list of changes. See the release notes for a complete list.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->

- [Local Z release coordinates](#local-z-release-coordinates)
- [Mass budget](#mass-budget)
- [Period block release settings](#period-block-release-settings)
- [Event timing and reporting](#event-timing-and-reporting)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Local Z release coordinates

The Particle Release Point (PRP) package supports configuring release point vertical coordinates as relative to the cell's vertical thickness via the `LOCAL_Z` option. In MF6.6.0, Newton models using `LOCAL_Z` could produce invalid negative particle tracking times when the head in the release-point cell was above the cell top or below the cell bottom. This was caused by an incorrect conversion of local release point vertical coordinates to model coordinates.

To correct this, the release point is now interpolated vertically between the *effective top* of the cell and the cell bottom:

- If the cell is **convertible**, the effective top is the head in the cell, constrained to be no higher than the geometric top and no lower than the cell bottom.
- If the cell is **confined**, the effective top is the geometric top of the cell.

To distinguish convertible from confined cells at particle release time, the user must now specify a binary grid file entry (`GWFGRID`) in the PRT FMI package when release coordinates are configured with `LOCAL_Z` and the PRT model is running in a separate simulation from the GWF model. A `GWFGRID` entry is not required if the GWF and PRT models run in the same simulation.

For instance, a sample FMI package input file with `GWFGRID`:

```
BEGIN OPTIONS
END OPTIONS

BEGIN PACKAGEDATA
  GWFBUDGET  FILEIN  gwf.bud
  GWFHEAD    FILEIN  gwf.hds
  GWFGRID    FILEIN  gwf.dis.grb
END PACKAGEDATA
```

## Mass budget

The PRT model mass budget has been rebalanced with a new term to accumulate the mass of terminated particles. Previously, particle mass remained in storage regardless of particle status, leaving the particle mass budget unbalanced at the end of the simulation. Cell-by-cell mass storage is also now reported in the binary budget file, where previously only cumulative mass storage was reported in the list file.

**Note**: Scripts or post-processing tools that parse PRT budget output may need to be updated to account for the new terminated-mass term.

## Period block release settings

The Particle Release Point (PRP) package is responsible for configuring when particles should be released. Release settings can be configured by period block, as with other MODFLOW 6 stress packages.

In MF6.6.0, the PRP package could fail to correctly load period block release settings, which could prevent particles from being released as configured, or cause excess particles to be released because empty period blocks intended to stop fill-forward at a certain stress period were not properly recognized. This has been corrected in MF6.7.0.

**Note**: If a model relies on period block release settings, results may differ from previous versions. Verify that particle releases occur as intended.

## Event timing and reporting

PRT previously reported events occurring exactly on a boundary between time steps as occurring during the *preceding* step. While this properly credited the event to the flow conditions responsible for the particle's trajectory, it was inappropriately applied to events such as particle releases, causing the stress period and time step to be reported incorrectly for some events.

In MF6.7.0, PRT reports each event as occurring in the stress period and time step in which the event is triggered. This means simultaneous releases can be credited to different time steps. For example, a particle release scheduled using an explicit time that falls on the boundary between two time steps occurs at the first opportunity, i.e., at the end of the earlier time step. On the other hand, if a period block is used to schedule a particle release at the beginning of the later time step, the release occurs in that time step, not the previous one.

**Note**: When comparing tracking output between versions, stress period and time step values for boundary events may differ.
