mc=10000 #maximal count
txtgen=1 # set to 1 -> generates txt file, set to 0 -> does NOT generate txt file

cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"

for Nd in 1 20 36 # number of hydra nodes
do

    for TPN in 1 16 32 # number of tasks per node
    do

        if [ "$TPN" -ne "1" ] || [ "$Nd" -ne "1" ];
        then

            #Run Binaries with srun
            echo "Exercise 1 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex1 -c $mc -p 2 -h $Nd -g $txtgen
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex1 -c $mc -p 10 -h $Nd -g $txtgen
            
        fi

    done
    
done

echo "Exericse 1 Hydra Run Finished! :-)"
