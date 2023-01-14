cs="--cpu-freq=High"
ts="--time=5:00"
ss="-p q_student"
mc=10000

for Nd in 1 2 4 # Nodes
do

    for TPN in 4 6 8 # Tasks Per Node
    do

        for bs in 100 1000 10000 100000 1000000 10000000   # blockSize
        do

            if [ "$TPN" -ne "1" ] || [ "$Nd" -ne "1" ];
            then
                #Run Binaries with srun
                echo "Exercise 2 with $Nd Node(s) and $TPN Task(s) per Node with both powers of 2 and 10"
                mpirun -np $TPN ./Ex2 -c $mc -p 2 -b $bs -h $Nd -g 1
                mpirun -np $TPN ./Ex2 -c $mc -p 10 -b $bs -h $Nd -g 1
            fi

        done

    done
    
done

echo "Exericse 2 Hydra Run Finished! :-)"
