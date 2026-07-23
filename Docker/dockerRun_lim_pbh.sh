#!/bin/bash

## Initialise docker container.
docker run -d \
  -it \
  -e PASSWORD=hello1 \
  -e USER=rstudio \
  --memory=248g \
  --cpus=28 \
  -p 8888:8787 \
  --name LIM_PBH-docker \
  --mount type=bind,source=/home/<username>/LIM_PBH,target=/home/rstudio/LIM_PBH \
  ceresbarros/lim_pbh

