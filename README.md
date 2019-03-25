
# QC Slow Pipeline

Uses Snakemake to Batch-Run Maxquant on raw files and adds QC-related metrics to a database.

The Workflow, in general, looks like this:

![](img/wokflow.jpg)

## Dependencies:

  - Windows
  - Python 3.6
  - MaxQuant
  
## Installation

Download the Pipeline with git, create a virtual environment (conda version shown here), 
and install all python dependencies into that environment:

```bash
git clone https://github.com/KeepCalmAndRushB/QCSlowPipeline
cd QCSlowPipeline
conda create --name QCSlow python=3.6
conda activate QCSlow
pip install -r requirements.txt
```

## Usage

### Do a Dry Run to confirm the pipeline is correct

```bash
snakemake -n
```


### Vizualize the Workflow

Note: This example uses [GraphViz](https://graphviz.gitlab.io/_pages/Download/Download_windows.html) to convert the workflow graph to PDF.
After installing it, add the path with the "dot.exe" file (C:\Program Files (x86)\Graphviz2.38\bin on my PC) to the PATH environment variable
so the **dot** command works in your terminal.

```
snakemake -n --dag | dot -Tpdf > img/workflow.pdf
```

Details for this command can be found at the [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/executable.html#visualization)

### Run the workflow on your files

```bash
snakemake
```

 
 