#if !defined (__CINT__) || defined (__CLING__)
	#include "AliAnalysisAlien.h"
	#include "AliAnalysisManager.h"
	#include "AliAODInputHandler.h"
	#include "AliESDInputHandler.h"
	#include "AliXiStarpp.h"
	#include <TMacro.h>
	R__LOAD_LIBRARY(CreateESDChain.C)
	R__LOAD_LIBRARY($ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C)
	#include "AliPhysicsSelectionTask.h"
#endif

Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,bool MCcase=kFALSE);
void runXiStar(const char *dataset = "test1.list") {
    int Nevents=100;
    bool batchmode=kTRUE;  
    bool MCcase=kFALSE; 
    bool AODcase=kFALSE; 
    bool gridTest=kTRUE; // kTRUE for the test mode, kFALSE for the Alien grid mode
    const char* collectionfile="collection.xml";
    const char* outfilename="MyOutput.root";
    const char* working_directory="pp13TeV_LHC15f";
    /*
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
    
    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY -g");
    */
    #if !defined (__CINT__) || defined (__CLING__)
	    gInterpreter->ProcessLine(".include $ROOTSYS/include");
	    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
	    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
    #else
	    gROOT->ProcessLine(".include $ROOTSYS/include");
    	gROOT->ProcessLine(".include $ALICE_ROOT/include");
    	gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
    #endif
    
    // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");

    // Event Handler
    if(AODcase) { // AOD Case
        AliAODInputHandler *aodH = new AliAODInputHandler();
        mgr->SetInputEventHandler(aodH);
    }
    else {  // ESD Case
        AliVEventHandler* esdH = new AliESDInputHandler();
        ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
        mgr->SetInputEventHandler(esdH);
        if(MCcase){
            AliMCEventHandler *mcHandler  = new AliMCEventHandler();
            mgr->SetMCtruthEventHandler(mcHandler);
        }
    }

    // compile the class and load the add task macro
    #if !defined (__CINT__) || defined (__CLING__)
    	gInterpreter->LoadMacro("AliXiStarppEventCollection.cxx+g");
    	gInterpreter->LoadMacro("AliXiStarpp.cxx+g");
    	AliXiStarpp *myTask = reinterpret_cast<AliXiStarpp*>(gInterpreter->ExecuteMacro("AddTaskXiStar.C"));
	#else
    	gROOT->LoadMacro("AliXiStarppEventCollection.cxx+g");
    	gROOT->LoadMacro("AliXiStarpp.cxx+g");
    	gROOT->LoadMacro("AddTaskXiStar.C");
    	AliXiStarpp *myTask = AddTaskXiStar("MyTask", AODcase, MCcase);
	#endif

	if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

	#if !defined (__CINT__) || defined (__CLING__)
    //gInterpreter->ProcessLine(".include $ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    //gInterpreter->LoadMacro("CreateESDChain.C");
    if(!AODcase){
        ////gInterpreter->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        //AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C+(MCcase)"));
        TMacro = physel("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(physel.Exec(MCcase));
        if(!physSelTask) { Printf("no physSelTask"); return; }
    }
    #else
    gROOT->ProcessLine(".L $ALICE_PHYSICS/PWG/EMCAL/macros/CreateESDChain.C");
    // event selection
    if(!AODcase){
        gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
        AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(MCcase);
        if(!physSelTask) { Printf("no physSelTask"); return; }
    }
    #endif

    if(batchmode==kFALSE) {// local pc mode; gets AOD or ESD files from my machine
    	TChain *chain = new TChain("ESDTree");
    	
        if(MCcase==kFALSE) {
            //if(AODcase==kTRUE) {chain = CreateChain("filelistAOD.txt",AODcase);}//My Personal list
            //else {chain = CreateChain(dataset,AODcase);}//My Personal list
            chain = CreateESDChain(dataset); // for KIAF use
        }else{ // for aod case, need to update
            //if(AODcase==kTRUE) chain = CreateChain("filelistAOD_MC.txt",AODcase);
            //else chain = CreateChain("filelistESD_MC.txt",AODcase);//Simulation list
        }
        mgr->StartAnalysis("local",chain,Nevents);

    }else {// alien mode
    	/*
        gSystem->Load("libNetx.so") ;
        gSystem->Load("libgapiUI.so");
        gSystem->Load("libRAliEn.so");
        */
        AliAnalysisAlien *plugin = new AliAnalysisAlien();
        
        gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
        plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        plugin->SetAnalysisSource("AliXiStarppEventCollection.cxx AliXiStarpp.cxx");
        plugin->SetAdditionalLibs("AliXiStarppEventCollection.h AliXiStarppEventCollection.cxx AliXiStarpp.h AliXiStarpp.cxx");

        //plugin->SetAliROOTVersion("v5-06-15");
    	//plugin->SetAliPhysicsVersion("v5-06-15-01"); //work for kMB untill June 23 2016.  //for data
    	plugin->SetAliPhysicsVersion("vAN-20171030-1");
	    plugin->SetAPIVersion("V1.1x");
	    // Data
		if(!MCcase) plugin->SetRunPrefix("000");   // real or MC data
     	TString *Production=new TString("15f");
		Int_t totruns=0;
	    // Declare input data to be processed.
	    if(!MCcase){// Real data
	        if(Production->Contains("10h")){// ESDs
	            plugin->SetGridDataDir("/alice/data/2010/LHC10h");
	            plugin->SetDataPattern("ESDs/pass2/*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC10h");
	        }
	       
	        if(Production->Contains("11h")){// ESDs
	            plugin->SetGridDataDir("/alice/data/2011/LHC11h_2");
	            plugin->SetDataPattern("ESDs/pass2/*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC11h_2");
	        }
	        if(Production->Contains("10d")){// ESDs
	            plugin->SetGridDataDir("/alice/data/2010/LHC10d");
	            plugin->SetDataPattern("ESDs/pass2/*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC10d");
	        }
	        if(Production->Contains("11a")){// ESDs
	            plugin->SetGridDataDir("/alice/data/2011/LHC11a");
	            plugin->SetDataPattern("ESDs/pass4_without_SDD/*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC11a");
	        }
	        if(Production->Contains("15f")){// ESDs
	            plugin->SetGridDataDir("/alice/data/2015/LHC15f");
	            plugin->SetDataPattern("pass2/*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC15f");
	        }
	    }else {// MC data
	        if(Production->Contains("11a")){// AODs
	            plugin->SetGridDataDir("/alice/sim/2012/LHC12f1b"); //1a : Pythia8, 1b : Phojet
	            plugin->SetDataPattern("*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC11a",MCcase);
	        }
	        if(Production->Contains("10h")){// ESDs
	            plugin->SetGridDataDir("/alice/sim/LHC11a10a_bis");
	            plugin->SetDataPattern("*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC10h",MCcase);
	        }
	        
	        if(Production->Contains("10d")){// David's Xi MC data set
	            plugin->SetGridDataDir("/alice/sim/LHC11a6a");
	            plugin->SetDataPattern("*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC10d",MCcase);
	        }
	               if(Production->Contains("10e")){// AODs
	            plugin->SetGridDataDir("/alice/sim/LHC10e20");
	            plugin->SetDataPattern("*AliESDs.root");
	            totruns += AddGoodRuns(plugin,"LHC10e",MCcase);
	        }
	    }

	    plugin->SetSplitMaxInputFileNumber(40);
	    plugin->SetExecutable("myTask.sh");
	    plugin->SetTTL(10000);
	    plugin->SetJDLName("runCode.jdl");
	    plugin->SetOutputToRunNo(kTRUE);
	    plugin->SetKeepLogs(kTRUE);
	    plugin->SetMaxMergeStages(1);
	    plugin->SetMergeViaJDL(kTRUE);

	    plugin->SetGridWorkingDir(working_directory);
	    plugin->SetGridOutputDir("output");

	    plugin->SetOverwriteMode(kTRUE);
    	plugin->SetUser("blim");

    	mgr->SetGridHandler(plugin);
		if(gridTest) {
            // speficy on how many files you want to run
            plugin->SetNtestFiles(1);
            // and launch the analysis
            plugin->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            plugin->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
    }

    
    
    /*
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
	*/
    
}
Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,bool MCcase=kFALSE) {
    //
    // Adds good runs from the Monalisa Run Condition Table
    //Int_t SetRunNumber = 6;

    
    Int_t nruns=0,ngoodruns=0;
    Int_t runlist[100];
    std::fill_n(runlist, 100, 0);

     if(lhcPeriod=="LHC15f") {
         
         if(MCcase){
         //mc list
         nruns=16;
         Int_t runlist_2[16]={146746, 146747, 146748, 146801, 146802, 146803, 146804, 146805, 146806, 146807, 146817, 146824, 146856, 146858, 146859, 146860};
         }else{
             //data list
            //nruns=48; // 48
            //Int_t runlist_2[48]={226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225587, 225586, 225579, 225578, 225576, 225322, 225314, 225313, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};
            nruns=1;
            Int_t runlist_2[1]={226500};
            for (unsigned i = 0; i < nruns; i++){
            	runlist[i] = runlist_2[i];
            }
         }
         
        // nruns=55; //pass1
        // Int_t runlist[55]={146860, 146859, 146858, 146856, 146824, 146817, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146287, 146282, 146277, 146273, 146272, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146148, 146147, 146141, 146099, 146079, 146072, 146071, 146027, 146026, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145292, 145290, 145289, 145288}; //146402, 146369,146292
    
         
         //167988 // there are only 10 events, no ESD files
     for(unsigned int k=0;k<nruns;k++){
     	if(MCcase) { TString *MCRunName = new TString();
     		*MCRunName += runlist[k];
     		plugin->AddRunNumber(MCRunName->Data());
     	}
     else plugin->AddRunNumber(runlist[k]);
     
     ngoodruns++;
     }
     plugin->SetNrunsPerMaster(1);
     }
    return ngoodruns; 
}