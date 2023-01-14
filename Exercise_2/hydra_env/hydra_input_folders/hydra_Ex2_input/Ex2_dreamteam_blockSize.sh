cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"
mc=1000000

for Nd in 1 20 36 # Nodes
do

    for TPN in 1 16 32 # Tasks Per Node
    do

        for bs in 10 100 1000 10000 100000 1000000 #blockSize
        do

            if [ "$TPN" -ne "1" ] || [ "$Nd" -ne "1" ];
            then
                #Run Binaries with srun
                echo "Exercise 2 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
                srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex2 -c $mc -p 2 -h $Nodes -b $bs -g 1
                srun $cs $ts $ss --nodes=$Nd --ntasks-per-node=$TPN ./Ex2 -c $mc -p 10 -h $Nodes -b $bs -g 1
            fi

        done

    done
    
done

echo "Exericse 2 Hydra Run Finished! :-)"
