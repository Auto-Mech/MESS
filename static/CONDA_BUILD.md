# Instructions for building & uploading the conda package

The following steps use the `rattler-build` and `anaconda-client` packages to
create a conda package for installing the static MESS executables. See
instructions below for how to install them. Once the package is built and
uploaded, it can be installed in an environment as follows.
```
conda install -c auto-mech mess-static
```

### Build the conda package

1. Update the version number in `static/recipe.yaml` with today's date: `YYYY.MM.DD`.
2. Run `rattler-build build -r static/` and look for `/path/to/file.conda` in the output. Copy that path for the steps below.

### Upload the conda package

1. Navigate to `https://anaconda.org/` to create an API token for yourself. Save it in a secure location where you can access it.
2. `export ANACONDA_OWNER=<username>`
3. `export ANACONDA_API_KEY=<API token>`
4. `anaconda upload -u auto-mech --label main /path/to/file.conda`

## Testing

### Testing the build

To test that the build works, you can create a test environment with it installed as follows.
```
conda create -n test-mess /path/to/file.conda
```
You can then test that the conda package works as follows.
```
conda activate test-mess
cd examples/
mess mess.inp
messpf messpf.inp
```
On subsequent tests, you will need to do remove the environment and clear the
cache before re-running the above commands (note that the cache must be cleared *after* removing the test environment).
```
conda env remove -n test-mess
conda clean --all   # Remove the cached version of test-mess
```

### Testing the upload

Follow the same steps as shown above for testing the build, but create the
environment using the uploaded version as follows. Don't forget to remove the
`test-mess` environment and run `conda clean --all` if it was already created.
```
conda create -n test-mess -c auto-mech mess-static
```

## Dependencies

The following steps use the `rattler-build` and `anaconda-client` packages.

These can be installed through pixi as follows:
```
curl -fsSL https://pixi.sh/install.sh | bash   # installs pixi
pixi global install rattler-build anaconda-client
```

Or, if you prefer, they can also be installed through conda:
```
conda install -c conda-forge rattler-build anaconda-client
```