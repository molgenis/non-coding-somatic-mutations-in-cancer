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

for chrom in "${chrom_num[@]}"
do
    echo ${chrom}
    # echo $i
    # echo "#!/usr/bin/bash" > file${i}.sh
    # echo "" >> file${i}.sh
    # echo "#SBATCH --job-name=file${i}" >> file${i}.sh
    # echo "#SBATCH --output=file${i}.out" >> file${i}.sh
    # echo "#SBATCH --error=file${i}.err" >> file${i}.sh
    # echo "#SBATCH --time=159:59:59" >> file${i}.sh
    # echo "#SBATCH --mem=96gb" >> file${i}.sh
    # echo "#SBATCH --nodes=1" >> file${i}.sh
    # echo "#SBATCH --open-mode=append" >> file${i}.sh
    # echo "#SBATCH --export=NONE" >> file${i}.sh
    # echo "#SBATCH --get-user-env=L" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "# Load ml" >> file${i}.sh
    # echo "ml Anaconda3/5.3.0" >> file${i}.sh
    # echo "source activate stage" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "echo 'HALLO ${i}'" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
    # echo "" >> file${i}.sh
done