counter=1
for x in $(seq -16 0.5 -5); do
  for y in $(seq -6 0.5 9); do
    for z in $(seq -20 0.5 -2); do
      printf "HETATM%5d  O   HOH S   1    %8.3f%8.3f%8.3f  1.00 15.00           O\n" \
        $counter $x $y $z >> QuantificationGrid.pdb
      counter=$((counter+1))
    done
  done
done
