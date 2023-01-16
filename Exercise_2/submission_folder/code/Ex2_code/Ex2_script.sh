mc=10000 # maximal count
txtgen=1 # if set to 1 -> generates .txt, if set to 0 -> does NOT generate .txt

cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"

for Nd in 4 20 36 # number of Nodes
do

    for TPN in 4 16 32 # number of Tasks Per Node
    do

        for bs in 100 1000 10000 100000 1000000  # blockSize
        do

            #Run Binaries with srun
            echo "Exercise 2 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex2 -c $mc -p 2 -b $bs -h $Nd -g $txtgen
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex2 -c $mc -p 10 -b $bs -h $Nd -g $txtgen
        
        done

    done
    
done

echo "Exericse 2 Hydra Run Finished! :-)"
