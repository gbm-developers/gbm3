//  GBM by Greg Ridgeway  Copyright (C) 2003

#include "gbmexcept.h"
#include "node_factory.h"

CNodeFactory::CNodeFactory()
{
}


CNodeFactory::~CNodeFactory()
{
    #ifdef NOISY_DEBUG
    Rprintf("destructing node factory\n");
    #endif
}


void CNodeFactory::NodeFactoryInitialize(unsigned long cDepth) {
  unsigned long i = 0;
  for(i=0; i<NODEFACTORY_NODGBM_RESERVE; i++) {
    TerminalStack.push(&(aBlockTerminal[i]));
    ContinuousStack.push(&(aBlockContinuous[i]));
    CategoricalStack.push(&(aBlockCategorical[i]));
  }
}


CNodeTerminal* CNodeFactory::GetNewNodeTerminal() {
 
  if(TerminalStack.empty()) {
    throw GBM::out_of_nodes();
  }
  
  CNodeTerminal* res = TerminalStack.top();
  TerminalStack.pop();
  res->reset();
  return res;
}


CNodeContinuous* CNodeFactory::GetNewNodeContinuous() {
  if(ContinuousStack.empty()) {
    throw GBM::out_of_nodes();
  }

  CNodeContinuous* res = ContinuousStack.top();
  ContinuousStack.pop();
  
  res->reset();
  
  return res;
}


CNodeCategorical* CNodeFactory::GetNewNodeCategorical() {
  
  if (CategoricalStack.empty()) {
    throw GBM::out_of_nodes();
  }
     
  CNodeCategorical* res = CategoricalStack.top();
  CategoricalStack.pop();
  
  res->reset();
  
  return res;
}
  
void CNodeFactory::RecycleNode(CNodeTerminal *pNode) {
  if(pNode) {
    TerminalStack.push(pNode);
  }
}

void CNodeFactory::RecycleNode(CNodeContinuous *pNode) {
  if(pNode) {
    if(pNode->pLeftNode) pNode->pLeftNode->RecycleSelf(this);
    if(pNode->pRightNode) pNode->pRightNode->RecycleSelf(this);
    if(pNode->pMissingNode) pNode->pMissingNode->RecycleSelf(this);
    ContinuousStack.push(pNode);
  }
}
  
void CNodeFactory::RecycleNode(CNodeCategorical *pNode) {
  if (pNode) {
    if(pNode->pLeftNode) pNode->pLeftNode->RecycleSelf(this);
    if(pNode->pRightNode) pNode->pRightNode->RecycleSelf(this);
    if(pNode->pMissingNode) pNode->pMissingNode->RecycleSelf(this);
    pNode->aiLeftCategory.resize(0);
    CategoricalStack.push(pNode);
  }
}
