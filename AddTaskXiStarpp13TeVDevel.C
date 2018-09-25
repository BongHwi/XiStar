AliXiStarpp13TeVDevel *AddTaskXiStarpp13TeVDevel(bool MCcase=kFALSE, int CutList=0, bool DevelopmentMode=kFALSE, bool HMTrigger=kFALSE, bool PIDOption=kFALSE, bool SetSystematic=kTRUE) {

  //===========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskBF", "No analysis manager to connect to.");
    return NULL;
  }

  //____________________________________________//
  // Create tasks
  AliXiStarpp13TeVDevel *XiStarTask = new AliXiStarpp13TeVDevel("XiStarTask", CutList);
  if(!XiStarTask) exit(-1);
  XiStarTask->SetDevelSetup(DevelopmentMode);
  XiStarTask->SetMCSetup(MCcase);
  XiStarTask->SetHMTSetup(HMTrigger);
  XiStarTask->SetPIDSetup(PIDOption);
  XiStarTask->SetSystematicSetup(SetSystematic);

  if(DevelopmentMode){
    std::cout << "MC Check: " << MCcase << std::endl;
    std::cout << "HMTrigger Check: " << HMTrigger << std::endl;
    std::cout << "PID Check: " << PIDOption << std::endl;
    std::cout << "Systemtaic Check: " << SetSystematic << std::endl;  
  }

  mgr->AddTask(XiStarTask);


  // Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputFileName = AliAnalysisManager::GetCommonFileName();
  outputFileName += ":PWGLF.outputXiStarAnalysis.root";
  AliAnalysisDataContainer *coutXiStar = mgr->CreateContainer(Form("MyList_%d", HMTrigger), TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName.Data());
  mgr->ConnectInput(XiStarTask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(XiStarTask, 1, coutXiStar);


  // Return the task pointer
  return XiStarTask;
}