# Launch scripts

This directory contains Slurm launch scripts running the solver and validate its output for various supercomputers.

Il faut en premier se connecter à un noeud de Ruche depuis la frontale. Voici un exemple pour pouvoir utiliser un gpu.

```bash
salloc \
      --reservation prouveurc_178 \
      --partition gpua100 \
      --gres gpu:1 \
      -t 30
```

Exemple d'output :
```bash
salloc: Granted job allocation 11192643
salloc: Waiting for resource configuration
salloc: Nodes ruche-gpu17 are ready for job
```

Il faut ensuite se connecter à la machine qui a un GPU : 

```bash
shh ruche-gpu17
```

Puis il est possible de soumettre un job :

```bash
srun \
      --reservation prouveurc_178 \
      --partition gpua100 \
      --gres gpu:1 \
      --pty path/to/exe
```


