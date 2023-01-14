#! /bin/bash
#SBATCH -p q_student
#SBATCH --cpu-freq=High 
#SBATCH --time=3:00 

echo "1 Node"
srun --cpu-freq=High --time=3:00 -p q_student --nodes 1 --ntasks-per-node 16 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 1 --ntasks-per-node 32 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 1 --ntasks-per-node 16 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 10, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 1 --ntasks-per-node 32 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 10, blocksize 4, generate textfile

echo "20 Nodes"
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 1 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 16 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 32 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfilesni
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 1 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 16 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 20 --ntasks-per-node 32 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile

echo "36 Nodes"
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 1 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate no textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 16 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate no textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 32 ./Ex1 -c 100 -p 2 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 1 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 16 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile
srun --cpu-freq=High --time=3:00 -p q_student --nodes 36 --ntasks-per-node 32 ./Ex1 -c 100 -p 10 -g 1 # vector length 1e7, powers of 2, blocksize 4, generate textfile