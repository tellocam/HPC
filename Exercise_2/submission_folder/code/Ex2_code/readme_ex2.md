# Run stuff on hydra
## Step 1: Compile Executable from C file on hydra
    $module load mpi/openmpiS
    $mpicc -o Ex2 Ex2.c -lm -O3

## Step 2: Run bash script on hydra
    $chmod +x Ex2_script.sh
    $./Ex2_script.sh

# Run stuff locally
## Step 1: Compile Executable
    $mpicc -o Ex2 Ex2.c -lm -O3

## Step 2: Run Executable
    $mpirun -np 1 ./Ex2 -c 1000 -p 10 -b 4 -h 1 -g 1

# Commandline args explained (not relevant on hydra, its prepared in `.sh` file)


## -c 1000
Maximal count lenght or buffer length of 1000

## -p 10
Every entry of count is a power of 10

## -b 4
Blocksize for pipelining is set to 4

## -h 1
only relevant if textgeneration is enabled, to concatenate the hydranodes number to the .txt files

## -g 1
If set to 1, a .txt file with relevant commandline arguments is generated. If set to 0, no .txt generation.

