#!/bin/csh

foreach exp (2000)
  echo "Processing experiment "${exp}
  grads -bpc "run climate.gs atmsmyl${exp} ${exp}" >& /tmp/tuc
  grads -bpc "run climatezm.gs atmsmyl${exp} ${exp}" >& /tmp/tuc
  grads -bpc "run climateglob.gs atmsmyl${exp} 10" 
end
