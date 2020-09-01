#!/usr/bin/env bash
docker run -d \
  -it \
  --name devtest \
  --mount type=bind,source="$(pwd)"/io_files,target=/mounted \
  dnpsoup:latest
