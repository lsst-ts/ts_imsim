# ts_imsim
Package to create Rubin AOS simulations with imSim

## Needed Packages

- [ts_wep](https://github.com/lsst-ts/ts_wep)
- [ts_ofc](https://github.com/lsst-ts/ts_ofc)
- [imSim](https://github.com/LSSTDESC/imSim)
- [astroplan](https://github.com/astropy/astroplan)
- [black](https://github.com/psf/black) (optional)
- [documenteer](https://github.com/lsst-sqre/documenteer) (optional)
- [plantuml](http://plantuml.com) (optional)
- [sphinxcontrib-plantuml](https://pypi.org/project/sphinxcontrib-plantuml/) (optional)

## Use of Module

1. Install and Setup the WEP and OFC packages first.

2. Setup the imSim environment by `eups`:

```bash
cd $ts_imsim_directory
setup -k -r .
scons
```
