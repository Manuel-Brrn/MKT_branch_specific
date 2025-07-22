#!/bin/bash
# You can (must) adapt the following variables :
GeCKO_path="/home/barrientosm/projects/GE2POP/2024_TRANS_CWR/2024_MANUEL_BARRIENTOS/03_scripts/mapping/GeCKO"
partition=agap_normal
config_file=$GeCKO_path/CONFIG/config_ReadMapping_speltoides_TrMo.fasta.yml
cluster_config=$GeCKO_path/CONFIG/RM_CLUSTER_PROFILE_SLURM_MAPPING_SPELTOIDES_TrMo/

echo "Lancement de GeCKO..."

# Run GeCKO using sbatch
sbatch --partition=$partition --wrap="${GeCKO_path}/runGeCKO.sh --workflow ReadMapping --config-file $config_file --cluster-profile $cluster_config --jobs 20 --latency-wait 60"

if [[ $? -eq 0 ]]; then
    echo "GeCKO lancé avec succès !"
else
    echo "Erreur lors du lancement de GeCKO !" >&2
    exit 1
fi

echo "Fin du script."
