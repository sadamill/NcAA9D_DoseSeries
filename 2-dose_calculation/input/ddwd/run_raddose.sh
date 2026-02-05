MAX_JOBS=4
count=0

for i in {1..38}
do
  (
    cd ../../output/ddwd/wedge${i}
    java -jar /Applications/RADDOSE-3D-master/raddose3d.jar -i ../../../input/ddwd/wedge${i}.txt
  ) &
  
  count=$((count + 1))

  if (( count % MAX_JOBS == 0 )); then
    wait  # Waits for background jobs to finish before continuing
  fi
done

wait  # Wait for any final jobs to finish
