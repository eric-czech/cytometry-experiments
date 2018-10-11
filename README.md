# cytometry-experiments
Shiny app for RGLab tools


### Docker

```bash
# Build container
cd $REPO/docker
docker build -t cytoui .

```

```
# Run this and then run the shinyCyto server using shinyCyto/R/start.R (at localhost:8787)
export DATA_DIR=$HOME/data/research/hammer
export ANALYSIS_REPO_DIR=$HOME/repos/cytometry-experiments
export BOOK_REPO_DIR=$HOME/repos/hammer/learning-about-t-cells-through-data
export UI_REPO_DIR=$HOME/repos/hammer/rglab/forks/shinyCyto
export PASSWORD=cytoui
docker run --rm -p 8787:8787 -p 3838:3838 \
-e USERID=$UID -e USER=$USER -e PASSWORD=$PASSWORD \
-v $UI_REPO_DIR:/home/$USER/shinyCyto \
-v $ANALYSIS_REPO_DIR:/home/$USER/analysis \
-v $BOOK_REPO_DIR:/home/$USER/book \
-v $DATA_DIR:/lab/data \
cytoui
```