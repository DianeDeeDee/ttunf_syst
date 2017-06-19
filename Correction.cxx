#include "Correction.h"
#include "Tools.h"

Tools * Correction::globalTools = 0;

Correction::Correction()
  : tools(0) {
  tools = globalTools;
}

Correction::~Correction() {
}

