# mess-static

The `mess-static` conda package holds fully static MESS executables (`mess`, `messpf`, `messabs`).
They have no runtime dependencies and run on any x86-64 Linux.
They are built from the repository root with the same `CMakeLists.txt` as the dynamic `mess` package, using `-DMESS_STATIC=ON` (see `recipe.yaml`).

This directory is a package of the root Pixi workspace. From the repository root:

```
pixi run -e static test                                # build and test the static executables
pixi publish --path static --target-dir dist           # build dist/mess-static-*.conda
```

The executables are then in `.pixi/envs/static/bin/`.
The GitHub Actions workflow uploads them as a build artifact.
Each release publishes the conda package to the `automech` channel on prefix.dev and attaches the executables to the GitHub release (see `RELEASING.md`).
