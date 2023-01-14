cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"
bs=4 # blockSize
mc=1000 # max. Count

for Nd in 1 20 36 # Hydra Nodes
do

    for TPN in 1 16 32 # Tasks Per Node
    do

        if [ "$TPN" -ne "1" ] || [ "$Nd" -ne "1" ];
        then
            #Run Binaries with srun
            echo "Exercise 4 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex4 -c $mc -p 2 -b $bs -h $Nd -g 1
            srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex4 -c $mc -p 10 -b $bs -h $Nd -g 1
        fi

    done
    
done

echo "Exericse 4 Hydra Run Finished! :-)"
