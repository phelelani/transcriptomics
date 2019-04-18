## Differential Expression and Pathway analysis R script

### Clone Repository:
```
git clone https://github.com/phelelani/transcriptomics.git
```

### Download/build singularity image
```
## Download
singularity pull --name "R.simg" shub://phelelani/transcriptomics:r
mv R.simg containers/

## Build
sudo singularity build containers/R.simg containers/Singularity.R
```

### Execute:
```
singularity exec --cleanenv containers/R.simg R -e 'source("scleroderma_expression.R", echo=TRUE)'
```
