#!/usr/bin/bash

#SBATCH --job-name=db_test
#SBATCH --output=db_test.out
#SBATCH --error=db_test.err
#SBATCH --time=159:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

#seq FIRST STEP LAST
chrom_num=($(seq 1 1 22))

chrom_num+=("X" "Y")

for i in "${chrom_num[@]}"
do
    echo $i
    echo "#!/usr/bin/bash" > /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "#SBATCH --job-name=file${i}" >> /jobs/file${i}.sh
    echo "#SBATCH --output=/jobs/file${i}.out" >> /jobs/file${i}.sh
    echo "#SBATCH --error=/jobs/file${i}.err" >> /jobs/file${i}.sh
    echo "#SBATCH --time=159:59:59" >> /jobs/file${i}.sh
    echo "#SBATCH --mem=96gb" >> /jobs/file${i}.sh
    echo "#SBATCH --nodes=1" >> /jobs/file${i}.sh
    echo "#SBATCH --open-mode=append" >> /jobs/file${i}.sh
    echo "#SBATCH --export=NONE" >> /jobs/file${i}.sh
    echo "#SBATCH --get-user-env=L" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "# Load ml" >> /jobs/file${i}.sh
    echo "ml Anaconda3/5.3.0" >> /jobs/file${i}.sh
    echo "source activate stage" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "echo 'HALLO ${i}'" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
    echo "" >> /jobs/file${i}.sh
done