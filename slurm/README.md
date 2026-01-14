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

Il faut ensuite relancer la préparation du noeud :

```bash
module purge \
&& module load "gcc/13.2.0/gcc-4.8.5" \
&& module load "cuda/12.8.0/gcc-13.2.0" \
&& module load "cmake/3.28.3/gcc-11.2.0" \
&& module load "python-nospack/3.12.4/gcc-11.2.0"
```

Utiliser la commande mini-run pour compiler et débuter l'exécution :

```bash
mini-run
```

Puis il est possible de soumettre un job :

```bash
srun \
      --reservation prouveurc_178 \
      --partition gpua100 \
      --gres gpu:1 \
      --pty path/to/exe
```

Un petit exemple rapide à compiler et exécuter :

```bash
mini-run -g gpu-a100 -s thermal
```