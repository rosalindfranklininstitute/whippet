# whippet

A quick script for simulating stuff for Mark and Nikolai.

# Installation

## Pip

This requires the same dependencies as parakeet.

```sh
  python -m pip install git+https://github.com/rosalindfranklininstitute/parakeet.git@master
```

## Apptainer (singularity)

```sh

  apptainer build whippet.sif docker://ghcr.io/rosalindfranklininstitute/whippet:main

  function container {
    apptainer exec --nv --bind=$(pwd):/mnt --pwd=/mnt whippet.sif $@
  }

  container whippet --new_config
  container whipper -c config.yaml
```
