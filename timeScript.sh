for run in {1..50}
do
	./stressTest b c < meta/interpolatedCollisionBones >> timeLog
	echo ${run}
done