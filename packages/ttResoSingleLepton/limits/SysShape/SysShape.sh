mkdir sysShape

for chan in e mu
do
  for cond in resolved_cat1 resolved_cat2 resolved_cat3 boosted_cat1 boosted_cat2 boosted_cat3
  do
#     for sys in JER0 JER1 JER2 JES0 JES1 JVF BoostedJES0 BoostedJES1 BoostedJES2 WRWptjmin10 WRWiqopt3 WHFC0 WHFC3 WHFC4 NormW Btag BtagC BtagL IFSR MCGen PartonShower PDF EWS topmass norm_tt norm_QCDe norm_QCDmu luminosity EleSF MuSF
     for sys in all
     do 
       root -b -q SysShape.C'("'$chan'","'$cond'","'$sys'")'
     done 
  done
done