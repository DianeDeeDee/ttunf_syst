#include "PDFReweighter.h"
#include "LHAPDF/LHAPDF.h"
#include <iostream>


/*

PDF reweighter class - James Ferrando 

uses LHAPDF to return reweighting factors on the fly

*/


PDFReweighter* PDFReweighter::Instance()
{
  static PDFReweighter inst;
  return &inst;
}

void PDFReweighter::Initialize(){
  // can either call this or just set InputSet and OutputSet
}

void PDFReweighter::SetInputPDFSet(std::string  input="None" )
{
  if (m_input_id==0){
    m_inputSet=input;
    if (input!="None"){
      int id= 0 ;
      id = m_used_ids.size()+1;
      m_used_ids.push_back(id);
      m_input_id=id;
      LHAPDF::initPDFSetByName(id,input);
      LHAPDF::usePDFMember(id,0);
    }
  } else{
    std::cout << "PDFReweighter::SetInputPDFSet: You may only set the input PDF once, sorry" << std::endl;
  }
}
void PDFReweighter::SetOutputPDFSet(std::string output="MSTW2008nlo68cl.LHgrid")
{
  if (m_output_id==0){
    m_outputSet=output;
    int id= 0 ;
    id = m_used_ids.size()+1;
    m_used_ids.push_back(id);
    m_output_id=id;
    LHAPDF::initPDFSetByName(id,output);
    LHAPDF::usePDFMember(id,0);
  } else{
    std::cout << "PDFReweighter::SetOutputPDFSet: You may only set the output PDF once, sorry" << std::endl;
  }
}

void PDFReweighter::SetErrorSets(std::vector<int> sets)
{
  if (m_errset_ids.size()==0){
    m_outputErrSets=sets;
    for (unsigned int i=0; i<sets.size();i++){
      int id= 0 ;
      id = m_used_ids.size()+1;
      m_used_ids.push_back(id);
      m_errset_ids.push_back(id);
      LHAPDF::initPDFSetByName(id,m_outputSet);
      LHAPDF::usePDFMember(id,m_outputErrSets.at(i));
      //    std::cout<< "Setting error set: " << sets.at(i) << std::endl; 
      if (id>3) std::cout << "PDFReweighter::SetErrorSets: WARNING: > 3 sets may not work, depending on lhapdf configuration " << std::endl;
    }   
  } else {
    std::cout << "PDFReweighter::SetErrorSets: You may only set the output Error sets once, sorry" << std::endl;
  }
  //  std::cout << "Number of Error Sets: " << m_errset_ids.size() << "" << std::endl;
  
}


double PDFReweighter::GetReweightingCentral(float x1, float x2, float q, float pid1, float pid2)
{
  double my_pdf_weight=1.0;
  
  if (m_inputSet!="None"){
    int my_pid1=pid1;
    int my_pid2=pid2;
    if (my_pid1==21){
      my_pid1=0;
    }
    if (my_pid2==21){
      my_pid2=0;
    }
    double ow1=LHAPDF::xfx(m_input_id,x1,q,my_pid1);
    double ow2=LHAPDF::xfx(m_input_id,x2,q,my_pid2);
    double nw1=LHAPDF::xfx(m_output_id,x1,q,my_pid1);
    double nw2=LHAPDF::xfx(m_output_id,x2,q,my_pid2);
    my_pdf_weight=nw1*nw2/(ow1*ow2);
    if (!(nw1>0.0) || !(nw2 >0.0) || !(ow1>0.0) || !(ow2 >0.0)){
      // don't reweight into negative xfx region - unphysical!
      my_pdf_weight=0.0;
    }
  } else{
    std::cout << "PDFReweighter::GetReweightingCentral: Can't reweight with this method unless you define an input set" << std::endl;
  }
  return my_pdf_weight;
}


  double PDFReweighter::GetReweightingCentral(float x1, float x2, float q, float pid1, float pid2,float xfx1, float xfx2)
{
  double my_pdf_weight=1.0;
  int my_pid1=pid1;
  int my_pid2=pid2;
  if (my_pid1==21){
    my_pid1=0;
  }
  if (my_pid2==21){
    my_pid2=0;
  }
  double nw1=LHAPDF::xfx(m_output_id,x1,q,my_pid1);
  double nw2=LHAPDF::xfx(m_output_id,x2,q,my_pid2);
  my_pdf_weight=nw1*nw2/(xfx1*xfx2);
  if (!(nw1>0.0) || !(nw2 >0.0) || !(xfx1>0.0) || !(xfx2 >0.0)){
    // don't reweight into negative xfx region - unphysical!
    my_pdf_weight=0.0;
  }
  return my_pdf_weight;
}


std::vector<double> PDFReweighter::GetReweightingErrorSets(float x1, float x2, float q, float pid1, float pid2,float xfx1, float xfx2)
{
  std::vector<double> my_weights;
  //  std::cout << "Getting reweightingsets " << m_outputErrSets.size() << std::endl;
  for (unsigned int i=0;i<m_outputErrSets.size();i++){
    double my_pdf_weight=1.0;
    int my_pid1=pid1;
    int my_pid2=pid2;
    if (my_pid1==21){
      my_pid1=0;
    }
    if (my_pid2==21){
      my_pid2=0;
    }
    double nw1=LHAPDF::xfx(m_errset_ids.at(i),x1,q,my_pid1);
    double nw2=LHAPDF::xfx(m_errset_ids.at(i),x2,q,my_pid2);
    my_pdf_weight=nw1*nw2/(xfx1*xfx2);
    if (!(nw1>0.0) || !(nw2>0.0)|| !(xfx1>0.0) || !(xfx2 >0.0)){
      // don't reweight into negative xfx region - unphysical!
      my_pdf_weight=0.0;
    }
    my_weights.push_back(my_pdf_weight);
  }

  return my_weights;
}


std::vector<double> PDFReweighter::GetReweightingErrorSets(float x1, float x2, float q, float pid1, float pid2)
{
  std::vector<double> my_weights;
  if (m_inputSet!="None"){
    double my_pdf_weight=1.0;
    int my_pid1=pid1;
    int my_pid2=pid2;
    if (my_pid1==21){
      my_pid1=0;
    }
    if (my_pid2==21){
      my_pid2=0;
    }
    double ow1=LHAPDF::xfx(1,x1,q,my_pid1);
    double ow2=LHAPDF::xfx(1,x2,q,my_pid2);
    for (unsigned int i=0;i<m_outputErrSets.size();i++){
      double nw1=LHAPDF::xfx(m_errset_ids.at(i),x1,q,my_pid1);
      double nw2=LHAPDF::xfx(m_errset_ids.at(i),x2,q,my_pid2);
      my_pdf_weight=nw1*nw2/(ow1*ow2);
      if (!(nw1>0.0) || !(nw2>0.0) || !(ow1>0.0) || !(ow2 >0.0)){
	// don't reweight into negative xfx region - unphysical!
	my_pdf_weight=0.0;
      }
      my_weights.push_back(my_pdf_weight);
    }
  } else{
    for (unsigned int i=0;i<m_outputErrSets.size();i++){
      my_weights.push_back(1.0);
      std::cout << "PDFReweighter::GetReweightingErrorSets: Can't reweight with this method unless you define an input set" << std::endl;
    }
  }
  
  
  return my_weights;
}


void PDFReweighter::PrintInfo()
{
  std::cout << "---------------------------" << std::endl;
  std::cout << "|   PDF Reweighter Class   |" << std::endl;
  std::cout << "|    v1.0  - J.Ferrando    |" << std::endl;
  std::cout << "----------------------------" << std::endl;
  
  std::cout <<"Internal Info:" << std::endl;
  std::cout <<"Input PDFSet:"  << m_inputSet << std::endl;
  std::cout <<"Output PDFSet:" << m_outputSet << std::endl;
  std::cout <<"Output PDF Error Sets: "  ;
  for (int i=0; i< m_outputErrSets.size();i++){
    std::cout << m_outputErrSets.at(i) << " ";
  }
  std::cout <<	   std::endl;

  std::cout << "----------------------------" << std::endl;
}

