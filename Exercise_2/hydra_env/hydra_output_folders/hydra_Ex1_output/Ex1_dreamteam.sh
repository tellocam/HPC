cpustr="--cpu-freq=High"
timestr="--time=5:00"
studstr="-p q_student"
mcount=10000000

for Nodes in 1 20 32
do

    for tasksPerNode in 1 16 32
    do

        if [ "$tasksPerNode" -ne "1" ] || [ "$Nodes" -ne "1" ];
        then
            #Run Binaries with srun
            echo "Exercise 1 with $Nodes Node(s) and $tasksPerNode Task(s) per Node with both powers of 2 and 10"
            srun $cpustr $timestr $studstr --nodes=$Nodes --ntasks-per-node=$tasksPerNode ./Ex1 -c $mcount -p 2 -h $Nodes -g 1
            srun $cpustr $timestr $studstr --nodes=$Nodes --ntasks-per-node=$tasksPerNode ./Ex1 -c $mcount -p 10 -h $Nodes -g 1
        fi

    done
    
done

echo "Exericse 1 Hydra Run Finished! :-)"
