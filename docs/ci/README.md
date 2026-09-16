# CI workflow definitions

These two files belong in `.github/workflows/`. They are parked here because
the token used to create them could not write to that directory (adding a
GitHub Actions workflow needs the `workflow` OAuth scope). Move them with:

```bash
git mv docs/ci/pipeline_ci.yml .github/workflows/pipeline_ci.yml
git mv docs/ci/container.yml   .github/workflows/container.yml
git rm docs/ci/README.md
```

- **`pipeline_ci.yml`** lints the pipeline and runs it in stub mode, both the
  core workflow and the identification branch. No container or real data
  needed, so it is fast enough for every push.

- **`container.yml`** builds `containers/quandenser/Dockerfile` when it or the
  sources change, and publishes to GHCR on a `rel-*` tag. It also builds weekly
  on a schedule, which is not decoration: the image depends on ProteoWizard via
  a floating `.lastSuccessful` URL, so it can stop building with no commit to
  this repository at all. That is precisely how the existing release workflow
  decayed unnoticed. A weekly build means CI finds out first.

The existing `build_and_release.yml` is a separate matter and needs repair on
its own terms: it uses `actions/upload-artifact@v1` and `download-artifact@v3`,
which GitHub has shut down; it targets `centos:centos7` and `fedora:35`, whose
package mirrors return 404; the Maven URLs in the Windows and macOS builders
are dead; and it hard-codes `x86_64` while `macos-latest` is arm64.
