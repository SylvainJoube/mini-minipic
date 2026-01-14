# Hackathon Kokkos, 12-16 janvier 2026 à la MdlS

Nous sommes l'équipe 5.

[Version en ligne, pas à jour](https://codimd.math.cnrs.fr/1Dj1bIr9R4SyfbW7IHdriQ?both#)

## Connexion Ruche

Host : ruche.mesocentre.universite-paris-saclay.fr

Lien vers le projet du hackathon : https://github.com/CExA-project/mini-minipic


```bash
ssh cexa-hk24@ruche.mesocentre.universite-paris-saclay.fr
cd $WORKDIR
```

Aller dans le fichier https://github.com/CExA-project/mini-minipic/blob/main/doc/compilation.md pour apprendre comment compiler.



## A faire une seule fois en local

```bash
# Pour ne pas avoir à taper mon mot de passe à chaque fois
ssh-copy-id -i /home/sylvain/.ssh/id_ed25519.pub cexa-hk24@ruche.mesocentre.universite-paris-saclay.fr

mkdir /mnt/hkokkos_2026

sshfs cexa-hk24@ruche.mesocentre.universite-paris-saclay.fr:/gpfs/workdir/cexa-hk24 /mnt/hkokkos_2026

# Ajouté à mon .ssh
  # Hackathon Kokkos 2026 CEXA
  Host hkokkos2026
    User cexa-hk24
    Hostname ruche.mesocentre.universite-paris-saclay.fr
    # Pour ne pas avoir à taper mon mot de passe à chaque fois
    # ssh-copy-id -i /home/sylvain/.ssh/id_ed25519.pub cexa-hk24@ruche.mesocentre.universite-paris-saclay.fr
```


## A faire une seule fois, lorsque connecté à Ruche

```bash
# Une fois connecté au serveur Ruche
cd $WORKDIR \
&& git clone --recurse-submodules https://github.com/CExA-project/mini-minipic.git \
&& cd mini-minipic && git submodule update --init \
&& pip3 install --user .
```


## A chaque connexion à Ruche

```bash
# A chaque lancement du projet
module purge \
&& module load "gcc/13.2.0/gcc-4.8.5" \
&& module load "cuda/12.8.0/gcc-13.2.0" \
&& module load "cmake/3.28.3/gcc-11.2.0" \
&& module load "python-nospack/3.12.4/gcc-11.2.0"

# Test de compilation qui a fonctionné :
##  Configuration
cmake -B build -DCMAKE_CXX_COMPILER=g++ -DKokkos_ENABLE_OPENMP=ON 
##  Build du projet
cmake --build build --parallel 10
```

## Fichiers Python permettant d'automatiser le cmake

???


## Lancement de slurm

Pour lancer une compilation et un job, il faut demander un noeud depuis la frontale, et ensuite s'y connecter via ssh.

[Lien vers les scripts slurm](slurm/README.md)


Pour le lancement slurm, suivre ces instructions : https://github.com/CExA-project/mini-minipic/blob/main/slurm/launch_ruche_a100.slurm

Mais remplacer "thermal_heavy" de la commande `mini-run -g gpu-a100 -s thermal_heavy` par `thermal`

