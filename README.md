[![Build Status](https://travis-ci.org/icgc-argo/dna-seq-processing-tools.svg?branch=master)](https://travis-ci.org/icgc-argo/dna-seq-processing-tools)
# DNA sequencing reads processing

This repository keeps a collect of data processing tools for DNA-Seq analysis. All tools are defined using Common Workflow Language (CWL).
Eevery tool is self-sufficient, can be independently developed, tested, released and used. This clean isolation allows maximum flexibility, maintainability and portability.

These tools are building blocks to create multi-step data analysis workflows as needed, like the
workflows here: https://github.com/icgc-argo/dna-seq-processing-wfs

## Development
As tools are meant to be independent from each other, arguably a better choice is to
develop each tool using its own source control repository and container image. In
reality, it's undesirable to have to manage too many repositories, so we ended up
with using one repository for many tools. Despite sharing the same repository, in
tools development, we still want to follow good practices to ensure as much as possible
tools are independent.

Besides common software development practices such as feature branches and PRs,
we'd like to follow additional guidelines.
1. To start development, start a branch from the current master.
2. If the planned work affects only one tool, name the branch with the tool's name as prefix,
   followed by '.something'. During development, no code/file not related to this tool should
   not be changed.
3. During the tool development, when it's ready commits should be pushed to the server to
   trigger automated Travis CI tests. Note that first trigger of the test may fail due to missing docker
   image (which is to be or being built by Quay.io, this problem maybe solved if we use Travis
   or Github to built docker image).
4. Repeat #3 until planned features are done and all tests pass. You can run local testing as
   well by invoking `pytest -v`.
5. Create PR for review when ready.
6. Review and update as needed until all tests pass and PR approved.
7. Merge feature branch to master.
8. When ready for a new release, create a new release branch named as intented release version,
   and in the tool CWL file update docker image and tag in dockerPull to include the same
   release version, for example: 'quay.io/icgc-argo/dna-seq-processing-tools:0.1.1'. Then push
   the branch to GitHub to trigger automated build and testing.
9. Note that if the work involves changing docker file with new or updated installed dependencies
   all tools should update to use the new docker image tag in the release branch, all tests require
   to be passed.
10. With tests pass in #8, merge to release branch to master and create a new release from master.
