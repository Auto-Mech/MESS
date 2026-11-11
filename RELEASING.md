# Releasing MESS

These notes are for the maintainers who make releases.

`mess` and `mess-static` share one version, `YYYY.M.N`, where `N` counts the releases within the month, starting from 0.
To release, from an up-to-date `master` of this repository:
1. List the changes to the source since the last release:
   ```
   git log --no-merges $(git describe --tags --abbrev=0)..HEAD -- src examples
   ```
2. Describe them under `[Unreleased]` in `CHANGELOG.md`, and commit.
3. Run `pixi run release YYYY.M.N` (for example, `pixi run release 2026.11.0`).
   This updates the version (in `pixi.toml` and `static/recipe.yaml`) and the changelog, then commits, tags and pushes.

The tag starts the release workflow (`.github/workflows/release.yml`), which:
- stops unless the tagged commit is on `master`;
- builds and tests both packages;
- publishes them to the [`automech`](https://prefix.dev/channels/automech) channel on prefix.dev;
- makes a GitHub release with the packages, a tarball of the static executables, and the changelog entry as its notes.

Publishing uses Trusted Publishing: the channel lists this repository and `release.yml` as a trusted publisher, so there are no API keys or secrets.
If the publish job fails, fix the cause and re-run the job: files that were already uploaded are skipped.

`pixi run version` shows the current version.
