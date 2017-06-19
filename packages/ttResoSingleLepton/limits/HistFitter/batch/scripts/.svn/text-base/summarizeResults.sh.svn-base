ls -d */ | while read dir
  do

  echo "\n\n"
  echo $dir
  grep "the signal scale factor is"  ${dir}/lim* 
  grep "fixed scan  in interval" ${dir}/lim*
  grep "Second derivative zero" ${dir}/lim* | uniq -c
  grep "MIGRAD TERMINATED WITHOUT CONVERGENCE" ${dir}/lim* | uniq -c
  #grep "Could not find histogram" ${dir}/lim* | uniq -c
  #grep "nominal sample QCDmultijet is empty for channel" ${dir}/lim* | uniq -c

  grep -A 6 "scan point(s) for hypo test inversion"  ${dir}/lim*

done

#alias getLimits='grep -A 5 "scan point(s) for hypo test inversion"  */lim*'
#alias getRemoved='grep "scan point(s) for hypo test inversion"  */lim*'
