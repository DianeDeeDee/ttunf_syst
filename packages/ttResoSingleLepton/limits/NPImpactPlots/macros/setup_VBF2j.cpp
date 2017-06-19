if (do2j) {

    //Fill channels list for 2j SRs
    vector<Channel> channels2j;//  names, ne,nm,nj
    if (!mergeSFSR) {
      if (doee) channels2j.push_back(Channel("ee","2j",2,0,2));
      if (domm) channels2j.push_back(Channel("mm","2j",0,2,2));
    } else if (doee && domm)
      channels2j.push_back(Channel("SF","2j",2,0,0));
    if (doem) channels2j.push_back(Channel("em","2j",1,1,2));
    if (splitem && dome) channels2j.push_back(Channel("me","2j",1,1,2));
    if (!channels2j.size())    {
      cout << "ERROR::No channels2j selected!" << endl;
      exit(1);
    }
    
    //Fill channels list for 2j CRs
    vector<Channel> crChannels2j;
    if (combineCRs)    {
      if (doem && dome) crChannels2j.push_back(Channel("OF","2j",1,1,2));
      if (!combineSFCRs)      {
        if (doee) crChannels2j.push_back(Channel("ee","2j",2,0,2));
        if (domm) crChannels2j.push_back(Channel("mm","2j",0,2,2));
      }
    }
    if (combineSFCRs)    {
      if (doee && domm) crChannels2j.push_back(Channel("SF","2j",2,0,2));
      if (!combineCRs)      {
        if (doem)            crChannels2j.push_back(Channel("em","2j",1,1,2));
        if (splitem && dome) crChannels2j.push_back(Channel("me","2j",1,1,2));
      }
    }





    //2j SRs definitions
    Region sr("signalLike");
    if (splitmjj || (doVBF2j && doMVAVBF2j && splitBDT))  sr.name = "signalLike1";
    sr.channels = channels2j;

    Region sr2("signalLike2");
    sr2.channels = channels2j;
    
    Region sr3("signalLike3");
    sr3.channels = channels2j;


    // LOWPT SIGNALREGION // FIXME update the region if needed
    Region sr_lo("lowPt");
    sr_lo.channels = channels2j;





    // // WW CONTROLREGION
    // Region cr("mainControl");
    // cr.channels = crChannels2j;
    // cr.isCR = true;
    // if (ZMode == 2) cr.isSFCR = false;
  

    // TOP CONTROLREGION
    Region tb("topbox");
    tb.channels = crChannels2j;
    tb.isCR = true;
    tb.isSFCR = true;

    Region tb2("topbox2");
    if ( (splitmjjtopCR && !doMVAVBF2j) || (doVBF2j && doMVAVBF2j && splitBDT) ){
      tb.name = "topbox1";
      tb2.channels = crChannels2j;
      tb2.isCR = true;
      tb2.isSFCR = true;
      tb2.name = "topbox2";
    }

    // Z CONTROL REGION FOR VBF
    Region CRE("zcrlowmet");
    CRE.channels = crChannels2j;
    CRE.isCR = true;
    CRE.isSFCR = true;

    Region CRE2("zcrlowmet2");
    if(doVBF2j && doMVAVBF2j && splitBDT){
      CRE.name = "zcrlowmet1";
      CRE2.channels = crChannels2j;
      CRE2.isCR = true;
      CRE2.isSFCR = true;
      CRE2.name = "zcrlowmet2";
    }








    //BDT SRs
    vector<Region> bdt_SRs;
    vector<double> bounds;
    if(binTopNF && doMVAVBF2j && doVBF2j){
      bounds = readBoundaryFile("config"+string(do2012?"_2012/":"_2011/")+"bdt_bins_2j.txt");
      nBinsTopNF = bounds.size()-1;
      for(int ii=0; ii < nBinsTopNF; ii++){
        stringstream srName;
	srName << "signalLike" << ii;
        Region sr4(srName.str());
        sr4.channels = channels2j;
	bdt_SRs.push_back(sr4);
      }
    }

    vector<Region> bdt_CRs;
    if(doVBF2j && doMVAVBF2j && binTopNF){
      for(int ii=0; ii < nBinsTopNF; ii++){
        stringstream tbName;
	tbName << "topbox" << ii;
        Region tb4(tbName.str());
        tb4.channels = crChannels2j;
        tb4.isCR = true;
        tb4.isSFCR = true;
	bdt_CRs.push_back(tb4);
      }
    }

    vector<Region> bdt_ZCRs;
    if(doVBF2j && doMVAVBF2j && binTopNF){
      for(int ii=0; ii < nBinsTopNF; ii++){
        stringstream CREName;
	CREName << "zcrlowmet" << ii;
        Region CRE4(CREName.str());
        CRE4.channels = crChannels2j;
        CRE4.isCR = true;
        CRE4.isSFCR = true;
	bdt_ZCRs.push_back(CRE4);
      }
    }


    // Z CONTROL REGION FOR SPIN
    Region zb("zbox");
    zb.channels = crChannels2j;
    

    // PACMAN STUFF STARTS HERE

    // HIGH MLL CONTROLREGION WITH RECOIL FAIL FOR PACMAN
    Region E("EWWCR");
    E.channels = crChannels2j;
    E.isCR = true;
    E.isSFCR = true;

    // HIGH MLL CONTROLREGION WITH RECOIL PASS FOR PACMAN
    Region E_frec("EfrecWWCR");
    E_frec.channels = crChannels2j;
    E_frec.isCR = true;
    E_frec.isSFCR = true;

    // Z WINDOW WITH RECOIL FAIL FOR PACMAN
    Region C("CZpeak");
    C.channels = crChannels2j;
    C.isCR = true;
    C.isSFCR = true;

    // Z WINDOW WITH RECOIL FAIL FOR PACMAN
    Region C_frec("CfrecZpeak");
    C_frec.channels = crChannels2j;
    C_frec.isCR = true;
    C_frec.isSFCR = true;

    // SIGNALREGION WITH RECOIL FAIL FOR PACMAN
    Region A("ASR");
    if (!mergeSFSR)
    {
      A.channels = channels2j;
      A.isCR = false;
      A.isSFCR = false;
    }
    else
    {
      A.channels = crChannels2j;
      A.isCR = true;
      A.isSFCR = true;
    }

    // SIGNALREGION WITH RECOIL PASS FOR PACMAN
    Region A_frec("AfrecSR");
    if (!mergeSFSR) 
    {
        A_frec.channels = channels2j;
        A_frec.isCR = false;
        A_frec.isSFCR = false;
    }
    else 
    {
      A_frec.channels = crChannels2j;
      A_frec.isCR = true;
      A_frec.isSFCR = true;
    }


    // DECIDE WHICH REGIONS TO CONSIDER
    // SR
    if (dosignalregion)
    {
      if (useHighPt)
      {
	if (doVBF2j && doMVAVBF2j && binTopNF)
	{
          for(int ii=0; ii < nBinsTopNF; ii++){ myRegions.push_back(bdt_SRs[ii]); }
	}
	else if (splitmjj || (doVBF2j && doMVAVBF2j && splitBDT && !binTopNF)) {
          myRegions.push_back(sr);
          myRegions.push_back(sr2);
          if (doSpin) myRegions.push_back(sr3);
        }
        else 
	  myRegions.push_back(sr);
      }
      if (useLowPt)   myRegions.push_back(sr_lo);
      if ( doPacman2j && (doee || domm)) myRegions.push_back(A);
      if ( doPacman2j && (doee || domm)) myRegions.push_back(A_frec);
    }

    // WW CR
    if (doWWCR2j)
      myRegions.push_back(cr);
    

    // TOP
    if(!doSpin && !(doVBF2j && doMVAVBF2j && binTopNF)) myRegions.push_back(tb);
    if(splitmjjtopCR && !doMVAVBF2j){
      myRegions.push_back(tb2);
    }
    else if(doMVAVBF2j && doVBF2j && splitBDT && !binTopNF){
      myRegions.push_back(tb2);
    }
    else if(doMVAVBF2j && doVBF2j && binTopNF){
      for(int ii=0; ii < nBinsTopNF; ii++){ myRegions.push_back(bdt_CRs[ii]);}
    }
   
    // Z CR
    if (doSpin && doMVAVBF2j && doZCR_MVA){
      myRegions.push_back(zb);
    }
    if (doMVAVBF2j && doVBF2j && doABCD2j && (doee || domm)){
      if (!binTopNF){
        myRegions.push_back(CRE);
        if (splitBDT) myRegions.push_back(CRE2);
      }
      else {
        for(int ii=0; ii < nBinsTopNF; ii++){ myRegions.push_back(bdt_ZCRs[ii]);}
      }
    }

    // PACMAN CR
    if ( doPacman2j && (doee || domm))
    {
      myRegions.push_back(C);
      myRegions.push_back(C_frec);
      if (!(overrideCuts && massMode == 3) && !useHighMass2) myRegions.push_back(E);
      if (!(overrideCuts && massMode == 3) && !useHighMass2) myRegions.push_back(E_frec);
    }

}



