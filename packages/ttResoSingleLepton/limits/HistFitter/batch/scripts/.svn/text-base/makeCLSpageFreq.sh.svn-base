Zsamples=( "Z400" "Z500" "Z750" "Z1000" "Z1250" "Z1500" "Z1750" "Z2000" "Z2250" "Z2500" "Z3000" )
KKgSamples=( "KKg400" "KKg500" "KKg600" "KKg700" "KKg800" "KKg900" "KKg1000" "KKg1150" "KKg1300" "KKg1600" "KKg1800" "KKg2000" "KKg2250" "KKg2500" "KKg2750" "KKg3000" )
KKgWidthSamples1=( "KKg1000_width10pc" "KKg1000_width15pc" "KKg1000_width20pc" "KKg1000_width25pc" "KKg1000_width30pc" "KKg1000_width35pc" "KKg1000_width40pc" )
KKgWidthSamples2=( "KKg2000_width10pc" "KKg2000_width15pc" "KKg2000_width20pc" "KKg2000_width25pc" "KKg2000_width30pc" "KKg2000_width35pc" "KKg2000_width40pc" )
KKgWidthSamples3=( "KKg3000_width10pc" "KKg3000_width15pc" "KKg3000_width20pc" "KKg3000_width25pc" "KKg3000_width30pc" "KKg3000_width35pc" "KKg3000_width40pc" )
RSGSamples=("RSG400" "RSG500" "RSG600" "RSG700" "RSG800" "RSG900" "RSG1000" "RSG1200" "RSG1400" "RSG1600" "RSG1800" "RSG2000" "RSG2500")

sigSamplesCLS=(${Zsamples} ${KKgSamples} ${KKgWidthSamples1} ${KKgWidthSamples2} ${KKgWidthSamples3} ${RSGSamples})


echo  '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">'
echo  '<html>' 
echo  '<head>'
echo  '   <title></title>'
echo  '</head>'


echo '<body>'
echo '<h1></h1>'
echo '<table>'

for point in "${sigSamplesCLS[@]}"
  do
  echo '<tr>'
  echo '<td> '${point}' </td>'
  echo '</tr>'
  echo '<tr>'  
  echo '<td><img border="0" src="upperlimit_cls_poi_'${point}'_noClean_Freq_CLs_grid_ts3.root.png" height="300"/> </td>'
  echo '<td><img border="0" src="upperlimit_cls_poi_'${point}'_Freq_CLs_grid_ts3.root.png" height="300" /> </td>'
  echo '</tr>'

done
