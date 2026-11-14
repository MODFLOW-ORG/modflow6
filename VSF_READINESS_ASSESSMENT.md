# VSF Capability Readiness Assessment

Status of the Variably Saturated Flow (VSF / Richards) capability on `vsf-dev`:
VSF (Richards flow), VSET (ET), VSSB (surface boundary), VSSP (seepage), and the
Brooks-Corey / Haverkamp / Van Genuchten soil models.

## Bottom line

The core is solid and the full VSF-family integration suite passes — Richards
flow, three soil models, kr-averaging, Newton *and* Picard, transient storage,
water-table interaction, VSSB (incl. mover and evaporation), VSSP, and VSET —
plus guard tests for every unsupported model/package combination. Remaining work
is mainly documentation (technical reference, release note) and a few deferred
numerical/feature items.

## Solid (verified)

- **Tests green.** The VSF-family suite passes, including negative/validation
  cases that confirm the guards fire cleanly (`autotest/test_gwf_vsf_*.py`,
  `test_gwf_vss{b,p}.py`, `test_gwf_vset.py`).
- **Misuse guards** with clear errors: partial-domain Richards (whole-domain
  `IVARSAT` required), convertible NPF cells under VSF, VSET-requires-VSF,
  VSSP-requires-VSF.
- **Unsupported-combination guards**, each with a test: transport via GWT/GWE
  exchanges, NPF XT3D, GWF-GWF exchanges (hence LGR and parallel), and BUY/VSC
  (`autotest/test_gwf_vsf_{transport,xt3d,gwfgwf,density}.py`).
- **Dev-only options gated** behind `IDEVELOPMODE` (`devsoilmodel`,
  `dev_soil_param`, `model_kr`); released binaries reject them.
- **Registered** in the GWF model, and mf6io input sections exist for all four
  packages (`doc/mf6io/gwf/{vsf,vset,vssb,vssp}.tex`).

## Open gaps

### Documentation
- No technical/mathematical reference for the Richards formulation or soil
  models in `doc/SuppTechInfo` (currently lives in a separate development repo).
  A standalone VSSP Picard/Newton note exists (`vssp-seepage-note.pdf`) but is
  not wired into the repo technical reference.
- No release-note item for VSF itself (only
  `doc/ReleaseNotes/items/vssb-package.toml`).
- Dev-only options still appear in the generated input tables (confusing even
  though they are runtime-gated).

### Functional / numerical
- Partial-domain Richards not supported — `IVARSAT` must be set everywhere
  (guarded).
- VSSP does not support the mover — hard error. (VSSB's mover works.)
- VSF storage Newton Jacobian (`uft_fn` in `VsfStorage.f90`) is a stub, so the
  Jacobian is approximate for transient VSF; hard problems may converge slowly.
- The VSSP seepage-face switch is a hard cutoff (not smoothed), so its Jacobian
  is discontinuous at the switch; the relative-permeability Newton term is in
  place.

## Deferred to release

- [ ] Technical reference (Richards/VSF formulation + soil models) and a VSF
      release-note item.
- [ ] VSF storage `uft_fn` Newton term; smooth the VSSP seepage-face switch for
      a continuous Jacobian.
- [ ] VSSP mover support.
- [ ] Partial-domain Richards support.
- [ ] Storage-budget entry question (`gwf-sto.f90`); exfiltration handling
      (`gwf.f90`).

## For teaching now

- Constrain exercises to single-model, GWF-only setups.
- State up front: whole-domain `IVARSAT`, non-convertible NPF cells, no VSSP
  mover.
