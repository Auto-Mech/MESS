# Instructions

These instructions show how to build and upload the `mess-static` conda package,
containing the static executables in this directory.
The process relies in the Pixi package manager, which can be installed as follows.
```
curl -fsSL https://pixi.sh/install.sh | sh
```

### Build

Build the conda package:
```
pixi build
```

### Upload

Upload the conda package:
```
export ANACONDA_API_KEY=<API token>
pixi run upload
```
