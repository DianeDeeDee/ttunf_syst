getLimits='grep -A 5 " The computed upper limit is" */lim*'

for mass in Z400 Z500 Z750 Z1000 Z1250 Z1500 Z1750 Z2000 Z2250 Z2500 Z3000  KKg500 KKg600 KKg700 KKg800 KKg900 KKg1000 KKg1150 KKg1300 KKg1600 KKg1800 KKg2000 KKg2250 KKg2500 KKg2750 KKg3000
#for mass in Z400 Z500 Z750 KKg400 KKg500 KKg600 KKg700 KKg800 KKg900
  do
  line1=(`grep 'expected limit (median)' *${mass}/limit*.log | tail -1`)
  line2=(`grep 'expected limit (-1 sig)' *${mass}/limit*.log | tail -1`)
  line3=(`grep 'expected limit (+1 sig)' *${mass}/limit*.log | tail -1`)
  line4=(`grep 'expected limit (-2 sig)' *${mass}/limit*.log | tail -1`)
  line5=(`grep 'expected limit (+2 sig)' *${mass}/limit*.log | tail -1`)
  
  # echo $line1 $line2 $line3 $line4

  #echo $mass GeV:  x-section = $line1[6] 16% = $line2[7]   84% = $line3[7]  2.5% = $line4[7]  97.5% = $line5[7]
  echo $mass $line1[6] $line2[7]  $line3[7]  $line4[7]   $line5[7]

done
