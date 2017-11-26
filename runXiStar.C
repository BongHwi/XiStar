void runXiStar(const char *dataset = "test1.list") {
    int Nevents=100;
    bool batchmode=kFALSE;  
    bool MCcase=kFALSE; 
    bool AODcase=kFALSE; 
    const char* collectionfile="collection.xml";
    const char* outfilename="MyOutput.root";
    
    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase.so");
    gSystem->Load("libESD.so");
    gSystem->Load("libAOD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");
    
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -g");
    
    
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

    
    // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");
    
    TChain* chain;
    
    if(batchmode==kFALSE) {// local pc mode; gets AOD or ESD files from my machine
        if(MCcase==kFALSE) {
            if(AODcase==kTRUE) {chain = CreateChain("filelistAOD.txt",AODcase);}//My Personal list
            else {chain = CreateChain(dataset,AODcase);}//My Personal list
        }else{
            if(AODcase==kTRUE) chain = CreateChain("filelistAOD_MC.txt",AODcase);
            else chain = CreateChain("filelistESD_MC.txt",AODcase);//Simulation list
        }
    }else {// alien mode
        gSystem->Load("libNetx.so") ;
        gSystem->Load("libgapiUI.so");
        gSystem->Load("libRAliEn.so");
        
        
        printf("*** Connect to AliEn ***\n");
        TGrid::Connect("alien://");
        
        gROOT->LoadMacro("CreateAlienHandler.C");
        AliAnalysisGrid *alienHandler = CreateAlienHandler(MCcase);
        if (!alienHandler) return;
        mgr->SetGridHandler(alienHandler);
    }
    
    if(AODcase) {
        AliAODInputHandler *aodH = new AliAODInputHandler();
        mgr->SetInputEventHandler(aodH);
    }
    else {
        AliVEventHandler* esdH = new AliESDInputHandler();
        ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
        mgr->SetInputEventHandler(esdH);
        if(MCcase){
            AliMCEventHandler *mcHandler  = new AliMCEventHandler();
            mgr->SetMCtruthEventHandler(mcHandler);
        }
    }
    
    
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(MCcase); //! PID response object
    
    // Multiplicity selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *MultSlection = AddTaskMultSelection();

    // centrality selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
    AliCentralitySelectionTask *taskCentrality = AddTaskCentrality(kFALSE);
    
    if (MCcase) {
        taskCentrality->SetMCInput();
        taskCentrality->DontUseCleaning(); // for injected MC
    }

    // event selection
    if(!AODcase){
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(MCcase);
    }
        if(!physSelTask) { Printf("no physSelTask"); return; }
    

    
    //____________________________________________//
    // Create tasks
    gROOT->LoadMacro("AliXiStarppEventCollection.cxx+g");
    gROOT->LoadMacro("AliXiStarpp.cxx+g");
    
    AliXiStarpp *myTask = new AliXiStarpp("MyTask", AODcase, MCcase);
    
     if(!myTask) exit(-1);
    mgr->AddTask(myTask);
    
    
    
    // Create containers for input/output
    
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput = mgr->CreateContainer("MyList", TList::Class(),AliAnalysisManager::kOutputContainer,outfilename);
   
    mgr->ConnectInput(physSelTask,0,cinput);

    //____________________________________________//
    mgr->ConnectInput(myTask, 0, cinput);
    
    
    mgr->ConnectOutput(myTask, 1, coutput);
    
    
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
    
    if(batchmode==kTRUE) mgr->StartAnalysis("grid");
    else mgr->StartAnalysis("local",chain,Nevents);
    
}

//______________________________________________________________________________
TChain *CreateChainFromCollection(const char* xmlfile, bool AODcase)
{
    TString *treename;
    if(AODcase) treename = new TString("aodTree");
    else treename = new TString("esdTree");
    
    // Create a chain from an alien collection.
    TAlienCollection * myCollection  = TAlienCollection::Open(xmlfile);
    
    if (!myCollection) {
        ::Error("CreateChainSingle", "Cannot create an AliEn collection from %s", xmlfile) ;
        return NULL ;
    }
    
    TChain* chain = new TChain(treename->Data());
    myCollection->Reset() ;
    while ( myCollection->Next() ) chain->Add(myCollection->GetTURL("")) ;
    chain->ls();
    return chain;
}

//______________________________________________________________________________
TChain *CreateChain(const char *fileName, bool AODcase)
{
    TString *treename;
    if(AODcase) treename = new TString("aodTree");
    else treename = new TString("esdTree");
    
    TChain* chainGood = new TChain(treename->Data());
    TChain* chainNull = new TChain(treename->Data());
    
    int counter=0;
    char line[500] = new char*;
    TString *name;
    
    ifstream inputstream(fileName);
    if(!inputstream) {cout<<"input file not found"<<endl; return chainNull;}
    
    while(!inputstream.eof()) {
        inputstream.getline(line,500);
        name = new TString(line);
        
        if(AODcase) name->Append("/AliAOD.root");
        else name->Append("/AliESDs.root");
        
        if(!inputstream.eof()) chainGood->Add(name->Data());
        
        counter++;
    }
    
    inputstream.close();
    
    return chainGood;
}


