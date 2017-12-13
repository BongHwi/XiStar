AliXiStarpp *AddTaskXiStar(bool MCcase=kFALSE, bool AODcase=kFALSE, int CutList=0) {
 
  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }
  
  //____________________________________________//
  // Create tasks
  AliXiStarpp *XiStarTask = new AliXiStarpp("XiStarTask", AODcase, MCcase, CutList);
  if(!XiStarTask) exit(-1);
  mgr->AddTask(XiStarTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF.outputXiStarAnalysis.root";
  
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("MyList", TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
   
  //mgr->ConnectInput(physSelTask,0,cinput);
  mgr->ConnectInput(XiStarTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(XiStarTask, 1, coutput);
  
  
  // Return the task pointer
  return XiStarTask;
}
