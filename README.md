[![Build Status](https://travis-ci.org/icgc-argo/dna-seq-processing-tools.svg?branch=master)](https://travis-ci.org/icgc-argo/dna-seq-processing-tools)
# DNA sequencing reads processing

This repository keeps a collect of data processing tools for DNA-Seq analysis. All tools are defined using Nextflow workflow language.

Every tool is self-sufficient, can be independently developed, tested, released and used. This clean isolation allows maximum flexibility, maintainability and portability.

These tools are building blocks to create multi-step data analysis workflows as needed, like the
workflows here: https://github.com/icgc-argo/dna-seq-processing-wfs

## Development
As tools are meant to be independent from each other, arguably a better choice is to
develop each tool using its own source control repository and container image. In
reality, it's undesirable to have to manage too many repositories, so we ended up
with using one repository for many tools. Despite sharing the same repository, in
tools development, we still want to follow good practices to ensure as much as possible
tools are independent.

Besides common software development practices such as feature branches, PRs and code reviews,
we'd like to follow additional guidelines.
1. To start development, start a branch from the current master.
2. If the planned work affects only one tool, name the branch with the tool's name as prefix,
   followed by '.release_version' (which is the same as the next new release tag, eg, `score-download.0.1.4`).
   During development, no code/file not related to this tool should be changed.
3. In the tool's Nextflow process file update docker image and tag in 'container' to include the same
   release version, for example: `quay.io/icgc-argo/score-download:score-download.0.1.4`.
4. During the tool development, when it's ready local commits should be pushed to the git server to
   trigger automated Travis CI tests. Note that first trigger of the test may fail due to missing docker
   image (which is to be or being built by Quay.io, this problem maybe solved if we use Travis
   or Github to built docker image).
5. Repeat #3 until planned features are done and all tests pass. You should also run local testing as
   well by invoking `pytest -v`, some tests are only executed locally due to security issues (at some
   point we should come up with plan to run those tests on Travis). If needed,
   you can also test the updated version of the tool in workflows that need it, just update the tool URL
   to point to tool's development branch.
6. Create PR for review when ready.
7. Review and update as needed until all tests (Travis + local tests) pass and PR approved.
8. Merge feature branch to master, and delete the branch. Create a new release from the current master branch,
   use the same tag as the tool development branch that was just merged and deleted, eg `score-download.0.1.4`.
10. Note that if the work involves changing the sharded base docker (ie, /docker/Dockerfile) with new or
   updated installed dependencies all tools that use the base docker should be updated individually going
   through the above process.
