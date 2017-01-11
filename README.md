# genray_scan
## Instructions
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
3. due to being run from IDL spawn (re-run at command line and re-run genray_scan, /run)
