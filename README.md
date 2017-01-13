# genray_scan
## Instructions
First add the `$genray_scan/idl` directory to your IDL path. Likely something like

```
!PATH=expand_path("+~/code/genray_scan/idl/:<IDL_DEFAULT>")
```
in whatever file the `IDL_STARTUP` environment variable points to.

```
mkdir rundir
cd rundir
cp -r $genray_scan/template3 .
idl
IDL>genray_scan, /run
```

## Notes
Runs seem to fail for three reasons ...

1. eps_xe being too large (or small?).
2. having a particular magnetic configuration (try tweaking a current +/- 1)
3. due to being run from IDL spawn (re-run at command line and re-run genray_scan, /run) ... this may be due to the '> genray.log' bit. Maybe not. Trying a separate `rungenray.sh` script instead.
4. prmt4 = 1.0d-6 seems to fix the unfixable?

## Failures and how to fix

```
Xe is close to xe0. BEFORE transmit_coef_ox_xyz:
 xe, xe0-eps_xe, cn2-cnpar2   0.5969E-02   0.9990E+00   0.5010E-02
 AFTER: transm_ox, transm_ox_old, Nper^2=   7.0566019774174318E-039   0.0000000000000000        5.0100313371724825E-003
 transm_ox<transm_ox_old, initiate OX jump
   find_rho_X_xyz ENTER: gradne=0
```

