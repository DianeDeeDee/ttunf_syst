// TTbarEWSREweighter.h
//
// provide a simple routine for weighting SM ttbar with EW sudakov corrections according to
// http://arxiv.org/abs/1201.3926
//
// To use with MC@NLO find top and anti-top with status 155 and input their invariant mass (in GeV)
//
//  James Ferrando <james.ferrando@glasgow.ac.uk> 2/7/2012
//   v1.0
//

float GetEWSWeight(float mtt /*GeV*/, int ecm=8 /*TeV*/){
  if (ecm == 8) {
    float nmtt = mtt/1000.0;
    float fitt = 1.0068099807259914 - 0.0802243110713301*nmtt + 0.022269990954561628*nmtt*nmtt - 0.004426804328171152*nmtt*nmtt*nmtt + 0.0003639444976699816*nmtt*nmtt*nmtt*nmtt;
    return fitt;
  }
  else if (ecm == 7){
    float bins[10]={0,350.0,500.0,750.0,1000.0,1500.0,2000.0,2500.0,3000.0,3500.0};
    float values[9]={1.0,0.98,0.97,0.95,0.94,0.92,0.90,0.88,0.87};
    for (int i=0; i<9;i++) {
      if (mtt < bins[i+1]) {
	return values[i];
      }
    }
    return 1.0;
  } 
  return 1.0; // because we need to return something
}
