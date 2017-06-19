

echo  '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">'
echo  '<html>' 
echo  '<head>'
echo  '   <title></title>'
echo  '</head>'


echo '<body>'
echo '<h1></h1>'
echo '<table>'

 echo '<tr>'
  echo '<td> '${point}' </td>'
  echo '</tr>'
  echo '<tr>'
  echo '<td><img border="0" src="c_corrMatrix.png" height="300"/> </td>'
  echo '<td><img border="0" src="fitresult_Boosted_MassSpectraFit_Excl_btagcat1btagcat2btagcat3_Nominal_default_.png" height="300" /> </td>'
  echo '</tr>'

  echo '<tr>'
  echo '<td> '"Before Fit"' </td>'
  echo '<td> '"After Fit"' </td>'
  echo '</tr>'

for cat in cat1ElectronsResolved cat1MuonsResolved cat2ElectronsResolved cat2MuonsResolved cat3ElectronsResolved cat3MuonsResolved cat1ElectronsBoosted cat1MuonsBoosted cat2ElectronsBoosted cat2MuonsBoosted cat3ElectronsBoosted cat3MuonsBoosted 
  do
  echo '<tr>'
  echo '<td> '${cat}' </td>'
  echo '</tr>'
  echo '<tr>'  
  echo '<td><img border="0" src="can_SRexcbtag_'${cat}'_cuts_beforeFit.png" height="300"/> </td>'
  echo '<td><img border="0" src="can_SRexcbtag_'${cat}'_cuts_afterFit.png" height="300"/> </td>'
  echo '</tr>'

done
