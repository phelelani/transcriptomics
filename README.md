## Differential expression and pathway analysis R script

### Clone repository:
```
git clone https://github.com/phelelani/transcriptomics.git
```

### Download or build singularity image
```
## To download:
singularity pull --name "R.simg" shub://phelelani/transcriptomics:r
mv R.simg containers/

## To build:
sudo singularity build containers/R.simg containers/Singularity.R
```

### Execute:
```
singularity exec --cleanenv containers/R.simg R -e 'source("scleroderma_expression.R", echo=TRUE)'
```
