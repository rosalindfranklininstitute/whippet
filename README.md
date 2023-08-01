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
  function container {
    apptainer exec --nv --bind=$(pwd):/mnt --pwd=/mnt \
      docker://ghcr.io/rosalindfranklininstitute/whippet:main $@
  }

  container whippet --new_config
  container whipper -c config.yaml
```
