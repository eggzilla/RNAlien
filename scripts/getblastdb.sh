#!/bin/bash
declare -a ntfiles=("nt.00.tar.gz" "nt.01.tar.gz" "nt.02.tar.gz" "nt.03.tar.gz" 
                "nt.04.tar.gz" "nt.05.tar.gz" "nt.06.tar.gz" "nt.07.tar.gz"
                "nt.08.tar.gz" "nt.09.tar.gz" "nt.10.tar.gz" "nt.11.tar.gz"
                "nt.12.tar.gz" "nt.13.tar.gz" "nt.14.tar.gz" "nt.15.tar.gz"
                "nt.16.tar.gz" "nt.17.tar.gz" "nt.18.tar.gz" "nt.19.tar.gz"
                "nt.20.tar.gz" "nt.21.tar.gz" "nt.22.tar.gz" "nt.23.tar.gz"
                "nt.24.tar.gz" "nt.25.tar.gz" "nt.26.tar.gz" "nt.27.tar.gz"
                "nt.28.tar.gz" "nt.29.tar.gz" "nt.30.tar.gz" "nt.31.tar.gz"
                "nt.32.tar.gz" "nt.33.tar.gz" "nt.34.tar.gz" "nt.35.tar.gz"
                "nt.36.tar.gz" "nt.37.tar.gz" "nt.38.tar.gz" "nt.39.tar.gz"
                "nt.40.tar.gz" "nt.41.tar.gz" "nt.42.tar.gz" "nt.43.tar.gz"
                "nt.44.tar.gz" "nt.45.tar.gz" "nt.46.tar.gz" "nt.47.tar.gz"
                "nt.48.tar.gz" "nt.49.tar.gz" "nt.50.tar.gz" "nt.51.tar.gz"
                "nt.52.tar.gz" "nt.53.tar.gz" "nt.54.tar.gz" "nt.55.tar.gz")
for f in "${ntfiles[@]}"
do
   echo "$f"
   wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/v5/$f
   tar zxvpf $f
done

