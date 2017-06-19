
class TopColorParams:
    def __init__(self):
#        self.cottheta=3.18
        self.f1=1.0
        self.f2=0.0
        self.sin2thetaw=0.232
        self.cos2thetaw=1-self.sin2thetaw
        self.costhetaw=self.cos2thetaw**0.5
        self.sinthetaw=self.sin2thetaw**0.5
        self.mt=172.5
        self.cot2thetah=1.5856
# 1.58571676838731 (from supriya is good to  < 1e-5)
        self.alphaZ0=1/129.0
        self.alpha=1/129.0
        self.g=((self.alpha*4*3.141592653589793238)**0.5)/(self.sinthetaw)
#        self.mz=90.18
#        self.alpha=1/128.0
    
        

class ZpWidthCalculator:
    def __init__(self):
        self.tcpar=TopColorParams()
    def __init__(self,tcpar):
        self.tcpar=tcpar
    def Width(self,masszp):
        tc2=self.tcpar
        mzp=masszp
        mt=tc2.mt
        my_width=tc2.alpha*tc2.cot2thetah*masszp/(8*tc2.cos2thetaw)
        extra=((1-4*(mt**2)/(mzp**2))**0.5)*(2+4*mt**2/mzp**2)+4
        my_width=my_width*extra
        return 100*my_width/masszp
    def TCPars(self):
        return self.tcpar

class ZpBrCalculator:
    def __init__(self):
        self.tcpar=TopColorParams()
    def __init__(self,tcpar):
        self.tcpar=tcpar
    def Br(self,pdgid,mzp):
        partial=0.0
        mt=self.tcpar.mt
        total=((1-4*(mt**2)/(mzp**2))**0.5)*(2+4*mt**2/mzp**2)+4
        f1=self.tcpar.f1
        f2=self.tcpar.f2
        rm2=(mt**2)/(mzp**2)
        if(pdgid==6):
            prefac=(1-4*rm2)**0.5
            partial= prefac*((1+f1**2)*(1-rm2)+6*f1*rm2)
        if(pdgid==5):
            partial=1+f2**2
        if(pdgid==2):
            partial=1+f1**2
        if(pdgid==1):
            partial=1+f2**2

        ratio=partial/total
        return ratio

class ZpSigmaCalculator:
    def __init__(self):
        self.wcalc=ZpWidthCalculator()
        self.tcpar=self.wcalc.TCPars()
        self.flavours=[1,2,3,4,5]
        self.pdf="cteq6ll.LHpdf"
        self.mzp=1000.0
        self.sqrts=14000.0
        self.nsteps=10000
    def __init__(self,wcalc):
        self.wcalc=wcalc
        self.tcpar=self.wcalc.TCPars()
        self.flavours=[1,2,3,4,5]
        self.pdf="cteq6ll.LHpdf"
        self.mzp=1000.0
        self.sqrts=14000.0
        self.nsteps=10000
    def initPDF(self):
        import lhapdf
        lhapdf.initPDFSetByName(self.pdf)

    def dSdMdy(self,yb,m):
        # Gev/C to microbarns
        #hbar=1.054571726e-34
        #c=299792458
        #sifac=(hbar*c)**2
        from math import log, exp, pi
        
        import lhapdf
    # Gev/C to picobarns
        hbar=1.054571726e-34
        c=299792458
        gevj=1.602176565e-10
#        sifac=1.0e6*1.0e34*((hbar*c)**2)/(gevj**2)
#      but for numerical stability only apply 1e6 here
        sifac=1.0e6
        sum=0.0
        shat=(m**2)
#*((self.sqrts)**2)
        if m<0.0000001:
            return 0.0
        x1=exp(yb)*(m/self.sqrts)
        x2=((m**2)/(x1*(self.sqrts)**2))
#        print x1,x2, m**2, m , self.sqrts, m**2/(self.sqrts**2)
        prefac=9.0*pi*((self.tcpar.alpha)**2)*(self.tcpar.cot2thetah**2)
        prefac=prefac/(16*(self.tcpar.cos2thetaw**2))
        scale= m/2.0
        gamma=0.5*((shat**0.5)/(self.tcpar.mt))
        beta=1.0
        
        for flav in self.flavours:
            flavfac=0.0
            topfac=0.0
            if scale>self.tcpar.mt:
                topfac=1.0
                beta=((1-1/(gamma**2)))**0.5
            fstate=2*beta*(1+(beta**2)/3)+beta*(1-beta**2)
            if flav==1 or flav==5:
                flavfac=1.0
            if flav==2 or flav==6:
                flavfac=2.0
            xfx1=max(lhapdf.xfx(x1,scale,flav),0.0)
            xfx2=max(lhapdf.xfx(x2,scale, -flav),0.0)
            xfx3=max(lhapdf.xfx(x1,scale,-flav),0.0)
            xfx4=max(lhapdf.xfx(x2,scale, flav),0.0)
            val=flavfac*prefac*topfac*self.propagator(x1,x2)*(xfx3*xfx4+xfx1*xfx2)*fstate
            if x1>1.0 or x2>1.0:
                val= 0.0

            sum=val+sum
        #average over spins and colours
        sum=sum/36.0
        return sifac*sum
    def propagator(self,x1,x2):
        shat=(x1*x2)*((self.sqrts)**2)
        width=self.wcalc.Width(self.mzp)*self.mzp/100.0
#        prop = ((self.sqrts)**2)/ ((shat-self.mzp**2)**2+shat*(width**2))
        prop = shat/ ((shat-self.mzp**2)**2+shat*(width**2))
        prop = prop
        return prop
    def gfun(self,x):
        from math import log
        return -log(self.sqrts/x)
    def hfun(self,x):
        from math import log
        return log(self.sqrts/x)
    def my_int(self, minm, maxm):
        from scipy.integrate import dblquad,quad        
        mstep=(maxm-minm)/self.nsteps
        mass=minm+mstep/2.0
        sum=0.0
        while mass<maxm:
            val=quad(self.dSdMdy,self.gfun(mass),self.hfun(mass),args=(mass,),epsrel=1e-04)
            sum=sum+2*val[0]*mstep/mass
#            print mass, mstep, val, sum
            mass=mass+mstep
        return sum
    def Sigma(self,mzp=1000.0):
        self.mzp=mzp
  #      sum= dblquad(self.dSdMdy,2*self.tcpar.mt,self.sqrts ,self.gfun, self.hfun, epsrel=1e-03 )
        sum= self.my_int(2*self.tcpar.mt,self.sqrts )
        return sum
        
