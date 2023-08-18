#!/bin/bash








for ranked in 0 1 2 3 4
do
    file_output="results/plddt/plddt_ranked_${ranked}.csv"
    echo ID,plddt_avg,plddt,aa_pos > ${file_output}
    for FILE in data/structures/ranked_${ranked}/*
    do      
         PLDDT_VEC=$(cat ${FILE} | cut -c 62- | cut -f 1 -d " " | tr '\n' ';')
         PLDDT_MEAN=$(cat ${FILE} | cut -c 62- | cut -f 1 -d " " | awk ' {total+=$1} END {print total/NR}'| tr ',' '.')
         AA_POS=$(cat ${FILE}| cut -c 24- | rev | cut -c 55- | rev | sed  -r s'/^[ ]+//g'| tr '\n' ';')  
         ID=$(echo ${FILE} | rev | cut -f 1 -d "/"| rev)
         echo ${ID}
         echo "${ID},${PLDDT_MEAN},${PLDDT_VEC},${AA_POS}" >> ${file_output}
    done
done