# wdl-hub

A centralized repository for organizing, storing, and sharing WDL tasks and workflows.

> [!NOTE]
> This repository is a work in progress. I use it to centralize my WDL tasks and workflows so that I can easily reuse them in my projects.
> However, I'm open to suggestions on how to improve this repository. If you have any ideas, please open an issue or a pull request.


## Project Structure

The project is organized into the following directories:

* `tasks`: Contains WDL tasks.
    * `tests`: Contains tests for WDL tasks.
* `pipelines`: Contains WDL workflows.
    * `<workflow-name>`: Contains a WDL workflow and its associated files.
        * `tests`: Contains tests for the WDL workflow.

## Releases

We use Pumbaa to build our pipeline releases. It packs all dependencies of a WDL workflow into a single zip file. This zip file can be used to run the workflow on Cromwell.

To build a release, run the following command:

```bash
# Install pumbaa
curl https://raw.githubusercontent.com/lmtani/pumbaa/main/assets/install.sh | bash

# Build a release
pumbaa build --wdl pipelines/benchmark_genome_assemblies/BenchmarkGenomeAssemblies.wdl --out releases/
```
