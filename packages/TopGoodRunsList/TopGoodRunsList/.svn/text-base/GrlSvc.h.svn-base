#ifndef GRLSVC_H
#define GRLSVC_H

#include <string>
#include <iostream>

#include "GoodRunsLists/TGoodRunsList.h"

/** @class GrlSvc
 **    
 ** @author W. H. Bell <W.Bell@cern.ch>
 **
 ** @brief Singleton class providing access to a GoodRunsLists (GRL) XML file.
 **
 ** The GoodRunsLists file is read once when the unique instance is created.
 ** This happens on the first call to GrlSvc::svc(std::string grlFileName),
 ** which also return a pointer to the instance. Subsequently it is sufficient
 ** to call GrlSvc::svc() without any arguments.
 */
class GrlSvc {
 public:
  /** Return a pointer to the unique instance, initialize if necessary.
	 *  @param grlFileName path and name of GRL XML file.
	 *  @return pointer to GrlSvc instance, NULL on failure.
	 */
  static GrlSvc *svc(std::string grlFileName);

  /** Return a pointer to the unique instance, or NULL if not initialized.
	 *  @return pointer to GrlSvc instance. */
  static GrlSvc *svc() { return s_instance; }


  /** Provide access to the TGoodRunsList object held by the instance.
	 *  @return pointer to the initialized TGoodRunsList object. */
  const Root::TGoodRunsList* grl() const { return &m_grl; }

  
 private:

  /** Singleton pattern: Private constructor. */
  GrlSvc(std::string grlFileName = "");
  
  /** Singleton pattern: Private destructor. */
  ~GrlSvc();
  
  /** Singleton pattern: Copy constructor is private and undefined. */
  GrlSvc(const GrlSvc &);

  /** Singleton pattern: Assignment operator is private and undefined. */
  GrlSvc &operator=(const GrlSvc &);
  
  /** A pointer to carry the address of the unique instance. */
  static GrlSvc *s_instance;
  
  
  /** A member function to setup the GoodRunsList reader. */
  int initialize(void);

  /** The GoodRunsLists XML file name. */
  std::string m_grlFileName;

  /** A GoodRunsList object to hold the merged list. */
  Root::TGoodRunsList m_grl;
};

#endif


