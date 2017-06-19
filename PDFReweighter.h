#ifndef PDFReweighter_h
#define PDFReweighter_h

#include <vector>
#include <string>

class PDFReweighter{
 public:
  static PDFReweighter* Instance(); 
  void SetInputPDFSet(std::string input);
  void SetOutputPDFSet(std::string output);
  void SetErrorSets(std::vector<int> errsets);
  double GetReweightingCentral(float x1, float x2, float q, float pid1, float pid2);
  double GetReweightingCentral(float x1, float x2, float q, float pid1, float pid2,float xfx1, float xfx2);
  std::vector<double> GetReweightingErrorSets(float x1, float x2, float q, float pid1, float pid2,float xfx1, float xfx2);
  std::vector<double> GetReweightingErrorSets(float x1, float x2, float q, float pid1, float pid2);
  void PrintInfo();
  void Initialize();

 private:
  PDFReweighter(){
    m_inputSet="None";
    m_outputSet="MSTW2008nlo68cl.LHgrid";
    std::vector<int> m_outputErrSet;
    m_input_id=0;
    m_output_id=0;
  };
  ~PDFReweighter(){};
  
  int m_input_id;
  int m_output_id;
  std::vector <int> m_used_ids;
  std::vector <int> m_errset_ids;
  std::string m_inputSet;
  std::string m_outputSet;
  std::vector<int> m_outputErrSets;
  
};

#endif
