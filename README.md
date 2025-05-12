# CAT

Cancer Associated TCR

## Run it using local conda environment


```
conda create --name CAT python=3.10.8
conda activate CAT
conda install -c conda-forge pandas numpy scanpy anndata seaborn matplotlib
conda install -c conda-forge python-igraph leidenalg
```
## To run on the Fred Hutch Cluster, first
```
grabnode
```
and select hardwares (some larger files it is better to use 32 CPUs and 512GB memory). Then enter the working directory:
```
cd /fh/working/sun_w/CAT
```
Load the modules:
```
ml purge 
ml SciPy-bundle/2023.02-gfbf-2022b
ml scenicplus/1.0.0-foss-2022b
ml JupyterLab/4.0.3-GCCcore-12.2.0
jupyter lab --ip=$(hostname) --port=$(fhfreeport) --no-browser
```
Copy the links, open in browser and work from there.

If any additional library is needed (for example scanpy), one can use this command to install into the user directory:
```
pip3 install --user scanpy
```
or search Fred Hutch python modules in terminal:
```
module avail xxx
```


