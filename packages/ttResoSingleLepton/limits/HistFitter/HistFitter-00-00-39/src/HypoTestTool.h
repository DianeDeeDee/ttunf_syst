// vim: ts=4:sw=4
/* -*- mode: c++ -*- */
// Standard tutorial macro for performing an inverted  hypothesis test 
//
// This macro will perform a scan of tehe p-values for computing the limit
// 

#ifndef HypoTestTool_hh
#define HypoTestTool_hh

#include <string>

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestCalculatorGeneric.h"
#include "TMsgLogger.h"

namespace RooStats {
    class HypoTestInverterResult;
    class HypoTestResult;
}
class RooWorkspace;


// internal class to run the inverter and more

namespace RooStats { 

    class HypoTestTool{

        public:
            HypoTestTool();
            ~HypoTestTool() { this->ResetHypoTestInverter(); this->ResetHypoTestCalculator(); };

            HypoTestInverterResult* RunHypoTestInverter(RooWorkspace * w, 
                        const char * modelSBName, const char * modelBName, 
                        const char * dataName,
                        int type,  int testStatType, 
                        bool useCLs, int npoints, double poimin, double poimax, int ntoys, 
                        bool useNumberCounting = false, 
                        const char * nuisPriorName = 0);

            HypoTestResult* RunHypoTest(RooWorkspace * w, bool doUL=true, 
                        const char * modelSBName="ModelConfig", const char * modelBName="", 
                        const char * dataName="obsData",
                        int type=0,  int testStatType=3, 
                        int ntoys=1000, 
                        bool useNumberCounting = false, 
                        const char * nuisPriorName = 0);

            void AnalyzeResult( HypoTestInverterResult * r,
                        int calculatorType,
                        int testStatType, 
                        bool useCLs,  
                        int npoints, 
                        const char* outfilePrefix = "",
                        const char* outfiletype = ".png" ); ///,const char * fileNameBase = 0 );

            void SetParameter(const char * name, const char * value);
            void SetParameter(const char * name, bool value);
            void SetParameter(const char * name, int value);
            void SetParameter(const char * name, double value);

            inline HypoTestCalculatorGeneric* GetCalculator() { return m_hc; }
            inline HypoTestInverter* GetInverter() { return m_calc; }

        private:
            bool SetupHypoTestInverter(RooWorkspace * w, 
                    const char * modelSBName="ModelConfig", const char * modelBName="", 
                    const char * dataName="obsData",
                    int type=0,  int testStatType=3, 
                    bool useCLs=true, 
                    int npoints=6, double poimin=0, double poimax=5, int ntoys=1000, 
                    bool useNumberCounting = false, 
                    const char * nuisPriorName = 0);

            bool SetupHypoTestCalculator(RooWorkspace * w, bool doUL=true, 
                    const char * modelSBName="ModelConfig", const char *modelBName="", 
                    const char * dataName="obsData",
                    int type=0,  int testStatType=3, 
                    int ntoys=1000, 
                    bool useNumberCounting = false, 
                    const char * nuisPriorName = 0);

            inline void ResetHypoTestCalculator() { if (m_hc!=0) { delete m_hc; m_hc=0; } };
            inline void ResetHypoTestInverter() { if (m_calc!=0) { delete m_calc; m_calc=0; } };

            HypoTestCalculatorGeneric* m_hc;
            HypoTestInverter* m_calc;

            bool mPlotHypoTestResult;
            bool mWriteResult;
            bool mOptimize;
            bool mUseVectorStore;
            bool mGenerateBinned;
            bool mUseProof;
            bool mRebuild;
            int     mNWorkers;
            int     mNToyToRebuild;
            int     mPrintLevel;
            int     mInitialFit; 
            int     mRandomSeed; 
            double  mNToysRatio;
            double  mMaxPoi;
            std::string mMassValue;
            std::string mMinimizerType;  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
            TString     mResultFileName; 

            bool mNoSystematics;
            double mConfLevel;
            TMsgLogger m_logger;

    };

} // end namespace RooStats


#endif