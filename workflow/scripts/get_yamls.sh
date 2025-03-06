
#!/usr/bin/bash

declare -a envs=("sp_utils" "sp_R" "sp_homology" "sp_tree" "sp_python")

for env in "${envs[@]}"; do
    echo $env
    conda env export -n $env --from-history | grep -E -v "^name:|^prefix:" > workflow/envs/${env}.yaml
done