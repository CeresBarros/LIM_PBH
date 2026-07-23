#!/bin/bash

## Initialise docker container.
docker run -d \
  -it \
  -e PASSWORD=hello1 \
  -e USER=rstudio \
  --memory=248g \
  --cpus=28 \
  -p 8080:8787 \
  --name LIM_PBH-docker \
  --mount type=bind,source=/home/rstudio/LIM_PBH,target=/home/rstudio/LIM_PBH \
  ceresbarros/lim_pbh

