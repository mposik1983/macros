#bin/bash

#pTarray=(0.0 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 10.0 15.0 20.0)
#pTarray=(0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 10.0)
pTarray=(0.3 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0)
#pTarray=(7.0 8.0)
etamin=0.0
etamax=0.5
count=1
echo "$etamin < eta < $etamax"
for j in ${pTarray[@]}
  do
  echo "count: $count"
  if [ $count == 12 ]; then  
    exit
  fi
  j1=${pTarray[$count]}
  echo " $j < pT < ${j1}"
  #Output Root file name identifies the parameters of the generator
  root -b -l -q './RunAnalysis.C('$j', '$j1','$etamin', '$etamax')'
  count=${count}+1
done
