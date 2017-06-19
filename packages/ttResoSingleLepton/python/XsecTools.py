import numpy as np

#
# XsecTools.py : tools for calculating theoretical uncertainties on cross sections
#
# Copyright James Ferrando - University of Oxford 2010 <james@ferrando.co.uk>
#                - University of Glasgow 2013
#
# v1.0 - Calculate PDF uncertainties on cross sections : cross checked against http://arxiv.org/abs/1007.5489
# v1.1 - Add in PDF uncertaitnies on asymmetries as for : http://arxiv.org/abs/1010.2130
# v1.2 - Add in PDF uncertainties on cross section ratios and correlations : used in http://arxiv.org/abs/1109.5141
# v1.3 - add in treatment of scale variations  : used in  http://arxiv.org/abs/1206.5731
# v1.4 - add in symmetric Hessian PDF error treatment
# v1.5 - add some print functions for easier output 
# v1.6 - add a fix from M. Lisovyi for symmetric error case when both sets give small downwards uncertainty
# v2.0 -  updates for PDF4LHC presciption 8/7 teV ratios for ATLAS H->bb team
# v2.1 - add Metadata member to xsecs 
# v2.2 - add treatment for ABM errors
#
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


def MakeGraphs(xsecs,name):
    from ROOT import TGraphAsymmErrors, Double
    import string
    
    stat=TGraphAsymmErrors(len(xsecs))
    total=TGraphAsymmErrors(len(xsecs))
    uppertotal=TGraphAsymmErrors(len(xsecs))
    lowertotal=TGraphAsymmErrors(len(xsecs))

    stat.SetName(name+"_stat")
    total.SetName(name+"_total")
    uppertotal.SetName(name+"_uppertotal")
    lowertotal.SetName(name+"_lowertotal")

    i=1
    for xsec in xsecs:
        ax=Double(string.atof(xsec.description))
        ay=Double(xsec.Mean())
        adummy=Double(0.0)
        upperr= Double(xsec.StatError())
        downerr=Double(xsec.StatError())
#        tupperr= Double((xsec.StatError()**2+xsec.UpperError()**2)**0.5)
        tupperr= Double((xsec.UpperError()**2)**0.5)
#        tdownerr=Double((xsec.StatError()**2+xsec.LowerError()**2)**0.5)
        tdownerr=Double((xsec.LowerError()**2)**0.5)
        auy=Double(xsec.Mean()+tupperr)
        ady=Double(xsec.Mean()-tdownerr)
#        print auy, ady
        stat.SetPointError(i,adummy,adummy,upperr,downerr)
        total.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        uppertotal.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        lowertotal.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        stat.SetPoint(i,ax,ay)
        total.SetPoint(i,ax,ay)
        uppertotal.SetPoint(i,ax,auy)
        lowertotal.SetPoint(i,ax,ady)
        i=i+1


    return [stat,total,uppertotal,lowertotal]


def MakeGraphsScale(xsecs,name):
    from ROOT import TGraphAsymmErrors, Double
    import string
    
    scale=TGraphAsymmErrors(len(xsecs))
    total=TGraphAsymmErrors(len(xsecs))
    uppertotal=TGraphAsymmErrors(len(xsecs))
    lowertotal=TGraphAsymmErrors(len(xsecs))
    uppertotal.SetName(name+"_upperscaletotal")
    lowertotal.SetName(name+"_lowerscaletotal")

    i=1
    for xsec in xsecs:
        ax=Double(string.atof(xsec.description))
        ay=Double(xsec.Mean())
        adummy=Double(0.0)
        upperr= Double(xsec.UpperScaleError())
        downerr=Double(xsec.LowerScaleError())
#        tupperr= Double((xsec.StatError()**2+xsec.UpperError()**2)**0.5)
        tupperr= Double((xsec.UpperScaleError()**2)**0.5)
#        tdownerr=Double((xsec.StatError()**2+xsec.LowerError()**2)**0.5)
        tdownerr=Double((xsec.LowerScaleError()**2)**0.5)
        auy=Double(xsec.Mean()+tupperr)
        ady=Double(xsec.Mean()-tdownerr)
        print "Scale:",xsec.Mean(), ady-xsec.Mean(), auy - xsec.Mean()
        scale.SetPointError(i,adummy,adummy,upperr,downerr)
        total.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        uppertotal.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        lowertotal.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        scale.SetPoint(i,ax,ay)
        total.SetPoint(i,ax,ay)
        uppertotal.SetPoint(i,ax,auy)
        lowertotal.SetPoint(i,ax,ady)
        i=i+1


    return [scale,total,uppertotal,lowertotal]

def MakeSymmGraphs(xsecs,name):
    from ROOT import TGraphAsymmErrors, Double
    import string
    
    stat=TGraphAsymmErrors(len(xsecs))
    total=TGraphAsymmErrors(len(xsecs))
    i=1
    for xsec in xsecs:
        ax=Double(string.atof(xsec.description))
        ay=Double(xsec.Mean())
        adummy=Double(0.0)
        upperr= Double(xsec.StatError())
        downerr=Double(xsec.StatError())
        tupperr= Double((xsec.StatError()**2+xsec.SymmError()**2)**0.5)
        tdownerr=Double((xsec.StatError()**2+xsec.SymmError()**2)**0.5)
        stat.SetPointError(i,adummy,adummy,upperr,downerr)
        total.SetPointError(i,adummy,adummy,tupperr,tdownerr)
        stat.SetPoint(i,ax,ay)
        total.SetPoint(i,ax,ay)
        i=i+1


    return [stat,total]

def PDF4LHC(mstw,ct10,nnpdf,rmup=1.0,rmdown=1.0):
    toterrorup=0.0
    toterrordo=0.0
    maxv=max(mstw.Mean()+mstw.UpperError(),ct10.Mean()+ct10.UpperError(), nnpdf.Mean()+nnpdf.UpperError() )
    minv=min(mstw.Mean()-mstw.LowerError(),ct10.Mean()-ct10.LowerError(),nnpdf.Mean()-nnpdf.LowerError())
#    print mstw.Mean()+mstw.UpperError(),ct10.Mean()+ct10.UpperError(), nnpdf.Mean()+nnpdf.UpperError()
    mean=(maxv+minv)/2.0
    pdferr=(maxv-minv)/2.0
    toterrorup=toterrorup+pdferr**2
    toterrordo=toterrordo+pdferr**2
    maxmt=0.5*(maxv+minv)*max(rmup,rmdown) - 0.5*(maxv+minv)
    minmt=0.5*(minv+maxv)-0.5*(minv + maxv)*min(rmup,rmdown)
    toterrorup=toterrorup+maxmt**2
    toterrordo=toterrordo+minmt**2

    print "Combined:"
    print (maxv+minv)*0.5, "+-",  0.5*(maxv-minv), "("+ str(100*(maxv-minv)/(maxv+minv))+"%) [PDF]",
    # currently use the mstw values for scale uncertainty (rescale to norm of central PDF4LHC)
    supper=mstw.UpperScaleError()/mstw.Mean()
    slower=mstw.LowerScaleError()/mstw.Mean()
    
    toterrorup=toterrorup+(supper*mean)**2
    toterrordo=toterrordo+(slower*mean)**2
    
    toterrorup=toterrorup**0.5
    toterrordo=toterrordo**0.5

    print  "+",str(supper*0.5*(maxv+minv)) , "("+ str(100*supper)+"%)","+",str(slower*0.5*(maxv+minv)) , "("+ str(100*slower)+"%)  [Scale]", "+", maxmt, "(",100*maxmt/(0.5*(maxv+minv)) ,"%)", "-", minmt,"(",100*minmt/(0.5*(maxv+minv)) ,"%) [m_t]"
    print mean, "+", toterrorup, "(",100.0*toterrorup/mean ,"%)" , "-", toterrordo,  "(",100.0*toterrordo/mean ,"%)" ,"[Total]"

    return 1

class CrossSection(object):
    def __init__(self):
        self.values=[]
        self.staterror=0.0
        self.description="No Description"
        self.metadata={}
    def AddValue(self,value):
        self.values.append(value)
        self.staterror=np.std(self.values)/(len(self.values)**0.5)
    def Mean(self):
        return np.mean(self.values)
    def StatError(self):
        return self.staterror
    def SetValueAndStatError(self,val,err):
        self.values=[val]
        self.staterror=err
    
class CrossSectionPDFErrors(CrossSection):
    def __init__(self):
        super(CrossSectionPDFErrors,self).__init__()
        self.ErrorTreatment="MSTW"
#       error treament can be "MSTW","NNPDF","CT10","NONE","ABM"        
        self.PDFSubsetCrossSections=[]
        self.AlphaSErrorSets=[]
        self.UpError=None
        self.SymmetricError=None
        self.DownError=None
    def Print(self):
        if self.ErrorTreatment=="NNPDF":
                    print "central:", self.Mean(), "+" , self.UpperError(), "-", self.LowerError() 
        else:
            print "central:", self.Mean(), "+" , self.UpperError(), "-", self.LowerError() 
    def AddSubsetCrossSection(self,xs):
        self.PDFSubsetCrossSections.append(xs)
    def AddAlphaSErrorSets(self,xs):
        self.AlphaSErrorSets.append(xs)
    def SetValueAndStatError(self,val,err):
        self.values=[val]
        self.staterror=err
    def UpperError(self):
        if self.UpError==None:
            self.UpError=0.0
            i=0
            if self.ErrorTreatment=="NNPDF":
                valarray=[]
                while i < len(self.PDFSubsetCrossSections):
                    valarray.append(self.PDFSubsetCrossSections[i].Mean())
                    i=i+1
                self.UpError=np.std(valarray)
                self.DownError=np.std(valarray)
# unneeded scaling of std deviation to mean error /(len(valarray)**0.5)
                if len(self.AlphaSErrorSets)>0:
                    maxas= self.AlphaSErrorSets[0].Mean()
                    minas= self.AlphaSErrorSets[0].Mean()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        maxas= max(self.AlphaSErrorSets[j].Mean(),maxas)
                        minas= min(self.AlphaSErrorSets[j].Mean(),minas)
                        j=j+1
                    self.UpError=self.UpError**2
                    self.UpError=self.UpError+((maxas-self.Mean())/1.645)**2
                    self.UpError=self.UpError**0.5
                    self.DownError=self.DownError**2
                    self.DownError=self.DownError+((minas-self.Mean())/1.645)**2
                    self.DownError=self.DownError**0.5
                return self.UpError
            if self.ErrorTreatment=="ABM":
                while i < len(self.PDFSubsetCrossSections):
                    self.UpError=self.UpError+(self.PDFSubsetCrossSections[i].Mean()-self.Mean())**2
                    i=i+1
                self.UpError=self.UpError**0.5
                return self.UpError
            while i < len(self.PDFSubsetCrossSections):
                big=max(self.PDFSubsetCrossSections[i].Mean(),self.PDFSubsetCrossSections[i+1].Mean())
 #               print big, max(0,big-self.Mean())
                self.UpError=self.UpError+max(0,big-self.Mean())**2
                i=i+2

            self.UpError=self.UpError**0.5
            if self.ErrorTreatment=="MSTW":
                if len(self.AlphaSErrorSets)>0:
                    maxas= self.AlphaSErrorSets[0].Mean()+self.AlphaSErrorSets[0].UpperError()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        maxas= max(self.AlphaSErrorSets[j].Mean()+self.AlphaSErrorSets[j].UpperError(),maxas)
                        j=j+1
                    self.UpError=max(self.UpError,maxas-self.Mean())
            if self.ErrorTreatment=="CT10":
                if len(self.AlphaSErrorSets)>0:
                    maxas= self.AlphaSErrorSets[0].Mean()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        maxas= max(self.AlphaSErrorSets[j].Mean(),maxas)
                        j=j+1
                    self.UpError=self.UpError**2
                    self.UpError=self.UpError+(maxas-self.Mean())**2
                    self.UpError=self.UpError**0.5
                self.UpError=self.UpError/1.645
            return self.UpError
        else:
            return self.UpError

    def LowerError(self):
        if self.DownError==None:
            self.DownError=0.0
            i=0
            if self.ErrorTreatment=="NNPDF":
                valarray=[]
                while i < len(self.PDFSubsetCrossSections):
                    valarray.append(self.PDFSubsetCrossSections[i].Mean())
                    i=i+1
                self.UpError=np.std(valarray)
                self.DownError=np.std(valarray)
                if len(self.AlphaSErrorSets)>0:
                    maxas= self.AlphaSErrorSets[0].Mean()
                    minas= self.AlphaSErrorSets[0].Mean()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        maxas= max(self.AlphaSErrorSets[j].Mean(),maxas)
                        minas= min(self.AlphaSErrorSets[j].Mean(),minas)
                        j=j+1
                    self.UpError=self.UpError**2
                    self.UpError=self.UpError+((maxas-self.Mean())/1.645)**2
                    self.UpError=self.UpError**0.5
                    self.DownError=self.DownError**2
                    self.DownError=self.DownError+((minas-self.Mean())/1.645)**2
                    self.DownError=self.DownError**0.5
                return self.DownError
            if self.ErrorTreatment=="ABM":
                while i < len(self.PDFSubsetCrossSections):
                    self.DownError=self.DownError+(self.PDFSubsetCrossSections[i].Mean()-self.Mean())**2
                    i=i+1
                self.DownError=self.DownError**0.5
                return self.DownError
            while i < len(self.PDFSubsetCrossSections):
                big=min(self.PDFSubsetCrossSections[i].Mean(),self.PDFSubsetCrossSections[i+1].Mean())
#                print big, min(0,big-self.Mean())
                self.DownError=self.DownError+max(0,self.Mean()-big)**2
                i=i+2
            self.DownError=self.DownError**0.5
            if self.ErrorTreatment=="MSTW":
                if len(self.AlphaSErrorSets)>0:
                    minas= self.AlphaSErrorSets[0].Mean()-self.AlphaSErrorSets[0].LowerError()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        minas= min(self.AlphaSErrorSets[j].Mean()-self.AlphaSErrorSets[j].LowerError(),minas)
                        j=j+1
                    self.DownError=max(self.DownError,self.Mean()-minas)
                return self.DownError
            if self.ErrorTreatment=="CT10":
                if len(self.AlphaSErrorSets)>0:
                    minas= self.AlphaSErrorSets[0].Mean()
                    j=1
                    while j < len(self.AlphaSErrorSets):
                        minas= min(self.AlphaSErrorSets[0].Mean(),minas)
                        j=j+1
                    self.DownError=self.DownError**2
                    self.DownError=self.DownError+(minas-self.Mean())**2
                    self.DownError=self.DownError**0.5
                self.DownError=self.DownError/1.645
            return self.DownError
        else:
            return self.DownError

    def SymmError(self):
        if self.SymmetricError==None:
            self.SymmetricError=0.0
            i=0
            while i < len(self.PDFSubsetCrossSections):
                big=max(self.PDFSubsetCrossSections[i].Mean(),self.PDFSubsetCrossSections[i+1].Mean())
                small=min(self.PDFSubsetCrossSections[i].Mean(),self.PDFSubsetCrossSections[i+1].Mean())
                #                print big, min(0,big-self.Mean())
                if (big-self.Mean())*(small-self.Mean())<0.0:
                    self.SymmetricError=self.SymmetricError+(0.5*(big-small))**2
                else:
                    self.SymmetricError=(self.SymmetricError+max((big-self.Mean())**2,(small-self.Mean())**2))
#
                i=i+2
            self.SymmetricError=self.SymmetricError**0.5
            return self.SymmetricError
        else:
            return self.SymmetricError
                                                                                                                        


    def PDFMean(self):
        i=0
        temp_array=[]
        while i < len(self.PDFSubsetCrossSections):
            temp_array.append(self.PDFSubsetCrossSections[i].Mean())
            i=i+1
        return np.mean(temp_array)

    def PDFRMS(self):
        i=0
        temp_array=[]
        while i < len(self.PDFSubsetCrossSections):
            temp_array.append(self.PDFSubsetCrossSections[i].Mean())
            i=i+1
        return np.std(temp_array)


    def AddedCrossSection(self,other):
        newxs=CrossSectionPDFErrors()
        newxs.SetValueAndStatError(self.Mean()+other.Mean(),(self.StatError()**2+other.StatError()**2)**0.5)
#        print "AddedCross:", (self.StatError()**2+other.StatError()**2)**0.5 , newxs.StatError(), newxs.staterr

        i=0
        while i < len(self.PDFSubsetCrossSections):
            merged=CrossSection()
            merged.SetValueAndStatError(self.PDFSubsetCrossSections[i].Mean()+other.PDFSubsetCrossSections[i].Mean(),(self.PDFSubsetCrossSections[i].StatError()**2+other.PDFSubsetCrossSections[i].StatError()**2)**0.5)
            newxs.AddSubsetCrossSection(merged)
            i=i+1
        return newxs
    def CrossSectionRatio(self,other):
        newxs=CrossSectionPDFErrors()
        newxs.SetValueAndStatError(self.Mean()/other.Mean(),(self.StatError()**2/(other.Mean()**2)+(other.StatError()**2)*self.Mean()**2/(other.Mean())**4)**0.5)
        i=0
        while i < len(self.PDFSubsetCrossSections):
            merged=CrossSection()
            merged.SetValueAndStatError(self.PDFSubsetCrossSections[i].Mean()/other.PDFSubsetCrossSections[i].Mean(),(self.PDFSubsetCrossSections[i].StatError()**2/(other.PDFSubsetCrossSections[i].Mean()**2)+(other.PDFSubsetCrossSections[i].StatError()**2)*self.PDFSubsetCrossSections[i].Mean()**2/(other.PDFSubsetCrossSections[i].Mean())**4)**0.5)
            newxs.AddSubsetCrossSection(merged)
            i=i+1
        i=0
        while i < len(self.AlphaSErrorSets):
            merged=self.AlphaSErrorSets[i].CrossSectionRatio(other.AlphaSErrorSets[i])
            newxs.AddAlphaSErrorSets(merged)
            i=i+1
        return newxs
    def CrossSectionAsymmetry(self,other):
        newxs=CrossSectionPDFErrors()

        a   = self.Mean();
        b   = other.Mean();
        bot = a+b;
        da=self.StatError()
        db=other.StatError()
        c2=1.0
        dc2=0.0
        totasymerror= 2*((a*a*c2*c2*db*db + c2*c2*b*b*da*da+a*a*b*b*dc2*dc2)**0.5)/(bot*bot)
        newxs.SetValueAndStatError((a-b)/(a+b),totasymerror)
        i=0
        while i < len(self.PDFSubsetCrossSections):
            merged=CrossSection()
            a   = self.PDFSubsetCrossSections[i].Mean();
            b   = other.PDFSubsetCrossSections[i].Mean();
            bot = a+b;
            da=self.PDFSubsetCrossSections[i].StatError()
            db=other.PDFSubsetCrossSections[i].StatError()
            c2=1.0
            dc2=0.0
            totasymerror= 2*((a*a*c2*c2*db*db + c2*c2*b*b*da*da+a*a*b*b*dc2*dc2)**0.5)/(bot*bot)
            merged.SetValueAndStatError((a-b)/(a+b),totasymerror)
            newxs.AddSubsetCrossSection(merged)
            i=i+1
        return newxs


    def PDFCorrelations(self,other):
        i=0
        asum=0.0
        sumerrx=0.0
        sumerry=0.0
        while i < len(self.PDFSubsetCrossSections):
            asum=asum+(self.PDFSubsetCrossSections[i].Mean()-self.PDFSubsetCrossSections[i+1].Mean())*(other.PDFSubsetCrossSections[i].Mean()-other.PDFSubsetCrossSections[i+1].Mean())
            sumerrx=sumerrx+(self.PDFSubsetCrossSections[i].Mean()-self.PDFSubsetCrossSections[i+1].Mean())**2
            sumerry=sumerry+(other.PDFSubsetCrossSections[i].Mean()-other.PDFSubsetCrossSections[i+1].Mean())**2
            i=i+2

        sumerrx=0.5*(sumerrx**0.5)
        sumerry=0.5*(sumerry**0.5)
        asum=asum/(4*(sumerrx*sumerry))
        print "Cos Phi:", asum
        print "Delta X", sumerrx
        print "Delta Y:", sumerry
        return 1

    def PDFStatCorrelations(self,other):
        i=0
        asum=0.0
        sumerrx=0.0
        sumerry=0.0

        while i < len(self.PDFSubsetCrossSections):
            asum=asum+(self.PDFSubsetCrossSections[i].Mean()-self.Mean())*(other.PDFSubsetCrossSections[i].Mean()-other.Mean())
            sumerrx=sumerrx+((self.PDFSubsetCrossSections[i].Mean()-self.Mean())**2)
            sumerry=sumerry+((other.PDFSubsetCrossSections[i].Mean()-other.Mean())**2)
            i=i+1
        sumerrx=(sumerrx**0.5)
        sumerry=(sumerry**0.5)
        asum=asum/((sumerrx*sumerry))
        print "Rho:", asum
        return 1
# read cross sections for a particular mass value from a top++-file, pdf uncertainites
    def ReadTopPPOutputPDF(self,filename,mass):
        import string
        my_file=open(filename)
        for line in my_file.readlines():
            splitline=string.split(line)
            if line[0]=="#":
                if splitline[1]=="Collider:":
                    self.description=str(string.atof(splitline[5])*1.0e-3)
                    self.metadata["ecm"]=splitline[5]
                    self.metadata["collider"]=splitline[2]
                if splitline[1]=="pdf":
                    self.metadata["pdf_pres"]=splitline[-1]
            if len(splitline)==6:
                if mass==splitline[0]:
                    if splitline[1]=="0":
                        self.SetValueAndStatError(string.atof(splitline[-1]),0.0)
                    else:
                        # add resilience to -inf
                        if splitline[-1]=="-inf" and self.metadata["pdf_pres"]=='NNPDF':
                            print "WARNING: infinity ignored"
                        else:
                            tempxs=CrossSection()
                            tempxs.SetValueAndStatError(string.atof(splitline[-1]),0.0)
                            self.AddSubsetCrossSection(tempxs)
        my_file.close()
        
class CrossSectionPDFScaleErrors(CrossSectionPDFErrors):
    def __init__(self):
        super(CrossSectionPDFScaleErrors,self).__init__()
        self.ScaleCrossSections=[]
        self.UpScaleError=None
        self.DownScaleError=None
    def AddScaleCrossSection(self,xs):
        self.ScaleCrossSections.append(xs)
    def UpperScaleError(self):
        self.UpScaleError=0.0
        i=0
        while i < len(self.ScaleCrossSections):
            self.UpScaleError=max(self.ScaleCrossSections[i].Mean()-self.Mean(),self.UpScaleError)
            i=i+1
        return self.UpScaleError
    def LowerScaleError(self):
        self.DownScaleError=0.0
        i=0
        while i < len(self.ScaleCrossSections):
            self.DownScaleError=max(-self.ScaleCrossSections[i].Mean()+self.Mean(),self.DownScaleError)
            i=i+1
        return self.DownScaleError
    def CrossSectionRatio(self,other):
        newxs=CrossSectionPDFScaleErrors()
        newxs.SetValueAndStatError(self.Mean()/other.Mean(),(self.StatError()**2/(other.Mean()**2)+(other.StatError()**2)*self.Mean()**2/(other.Mean())**4)**0.5)
        i=0
        while i < len(self.PDFSubsetCrossSections):
            merged=CrossSection()
            merged.SetValueAndStatError(self.PDFSubsetCrossSections[i].Mean()/other.PDFSubsetCrossSections[i].Mean(),(self.PDFSubsetCrossSections[i].StatError()**2/(other.PDFSubsetCrossSections[i].Mean()**2)+(other.PDFSubsetCrossSections[i].StatError()**2)*self.PDFSubsetCrossSections[i].Mean()**2/(other.PDFSubsetCrossSections[i].Mean())**4)**0.5)
            newxs.AddSubsetCrossSection(merged)
            i=i+1
        i=0
        while i < len(self.AlphaSErrorSets):
            merged=self.AlphaSErrorSets[i].CrossSectionRatio(other.AlphaSErrorSets[i])
            newxs.AddAlphaSErrorSets(merged)
            i=i+1
        i=0
        while i < len(self.ScaleCrossSections):
            merged=CrossSection()
            merged.SetValueAndStatError(self.ScaleCrossSections[i].Mean()/other.ScaleCrossSections[i].Mean(),(self.ScaleCrossSections[i].StatError()**2/(other.ScaleCrossSections[i].Mean()**2)+(other.ScaleCrossSections[i].StatError()**2)*self.ScaleCrossSections[i].Mean()**2/(other.ScaleCrossSections[i].Mean())**4)**0.5)
            newxs.AddScaleCrossSection(merged)
            i=i+1
        return newxs

    def Print(self):
        if self.ErrorTreatment=="NNPDF":
            print "central:", self.Mean(), "+" , self.UpperError(), "-", self.LowerError(), "(PDF)", "+", self.UpperScaleError(), "-", self.LowerScaleError(),  "(Scale)"
        else:
            print "central:", self.Mean(), "+" , self.UpperError(), "-", self.LowerError() , "(PDF)", "+", self.UpperScaleError(), "-", self.LowerScaleError(), "(Scale)"

    def ReadTopPPOutputScale(self,filename,mass):
        import string
        my_file=open(filename)
        for line in my_file.readlines():
            splitline=string.split(line)
            if len(splitline)==5:
                if mass==splitline[0]:
                    tempxs=CrossSection()
                    tempxs.SetValueAndStatError(string.atof(splitline[-1]),0.0)
                    self.AddScaleCrossSection(tempxs)
        my_file.close()
