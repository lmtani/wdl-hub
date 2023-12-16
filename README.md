# WDL-hub

A centralized repository for organizing, storing, and sharing WDL tasks and workflows. Currently it supports only WDL 1.0.

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

## Tests

Here we just test that the WDL tasks and workflows can be run without errors. We don't test the correctness of the results.

## How to use it

1. Download the release .tar.gz file from the releases page.
1. Extract in the directory you will develop your WDL workflow.
1. Import the tasks and workflows you want to use in your WDL workflow.

### Importing tasks

To import a task, add the following line to the top of your WDL file:

```wdl
import "wdl-hub/tasks/<module-name>.wdl" as a_name
```

Then you can use the task in your WDL workflow:

```wdl
workflow my_workflow {
    call a_name.task_name { input: ... }
}
```

### Building a workflow

We use [Pumbaa]() to build our workflows. It packs all dependencies of a WDL workflow into a single zip file. This zip file can be used to run the workflow on Cromwell.

To build a workflow, run the following command:

```bash

# Install pumbaa
curl https://raw.githubusercontent.com/lmtani/pumbaa/main/assets/install.sh | bash

# Build a workflow
pumbaa build --wdl MyWorkflow.wdl --out releases/
```
