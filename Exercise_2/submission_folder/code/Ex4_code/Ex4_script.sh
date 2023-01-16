mc=10000 # max. Count
txtgen=1 # if set to 1 -> generates .txt, if set to 0 -> does NOT generate .txt

cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"

for Nd in 4 20 36 # Hydra Nodes
do

    for TPN in 4 16 32 # Tasks Per Node
    do

        #Run Binaries with srun
        echo "Exercise 4 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
        srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex4 -c $mc -p 2 -b 100 -h $Nd -g $txtgen
        srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex4 -c $mc -p 10 -b 100 -h $Nd -g $txtgen

    done
    
done

echo "Exericse 4 Hydra Run Finished! :-)"
