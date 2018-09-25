#if !defined (__CINT__) || defined (__CLING__)
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliXiStarpp13TeVDevel.h"
#endif

Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,Bool_t MCcase=kFALSE);
TChain *CreateChain(const char *fileName);

void runXiStar6(const char *mode = "local") { // local/test/full/terminate
    int Nevents=100;

    Bool_t MCcase=kTRUE;
    Bool_t HMTrigger=kFALSE; // kTRUE for High Multiplicity trigger mode.
    Bool_t PIDOption=kTRUE;
    Bool_t SetSystematic = kFALSE; //
    

    Int_t CutList = 0;
    Bool_t DevelopmentMode=kTRUE;
    TString *Production=new TString("16k");
    const char* working_directory="pp13TeV_LHC16k_MC3";
    const char *dataset = "test1.list";

    gSystem->Load("libTree.so");
    gSystem->Load("libGeom.so");
    gSystem->Load("libVMC.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libSTEERBase. so");
    gSystem->Load("libESD.so");
    gSystem->Load("libANALYSIS.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libANALYSISalice.so");

    gSystem->ListLibraries();

    gSystem->SetIncludePath("-I. -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_ROOT/STEER -I$ALICE_ROOT/ANALYSIS -I$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY -g");
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT 6 MODE
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
    gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
#else
    // ROOT 5 MODE
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
    gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
#endif

    // Make the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("My Manager","My Manager");

    // Event Handler
    AliVEventHandler* esdH = new AliESDInputHandler();
    ((AliESDInputHandler *) esdH)->SetReadFriends(kFALSE);
    mgr->SetInputEventHandler(esdH);
    if(MCcase) {
        AliMCEventHandler *mcHandler  = new AliMCEventHandler();
        mgr->SetMCtruthEventHandler(mcHandler);
    }

    // compile the class and load the add task macro
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT 6 MODE
    //
    // Physics Selection
    AliPhysicsSelectionTask *physSelTask = reinterpret_cast<AliPhysicsSelectionTask*>(gInterpreter->ExecuteMacro(Form("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C(%d)",MCcase)));
    if(!physSelTask) {
        Printf("no physSelTask");
        return;
    }
    // Multiplicity selection
    AliMultSelectionTask *MultSlection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ExecuteMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"));
    if(!MultSlection) {
        Printf("no MultSlection");
        return;
    }
    // PID response
    AliAnalysisTask *fPIDResponse = reinterpret_cast<AliXiStarpp13TeVDevel*>(gInterpreter->ExecuteMacro(Form("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C(%d)",MCcase)));
    if(!fPIDResponse) {
        Printf("no fPIDResponse");
        return;
    }
    //gROOT->LoadMacro("$ALICE_PHYSICS/PWGLF/STRANGENESS/Cascades/Run2/macros/AddTaskWeakDecayVertexer.C");
    //AliAnalysisTaskWeakDecayVertexer *taskWDV = AddTaskWeakDecayVertexer();
    //taskWDV->SetRunV0Vertexer(kTRUE);
    //taskWDV->SetRunCascadeVertexer(kFALSE);
    //taskWDV -> SetDoImprovedCascadeVertexFinding(kTRUE);
    //taskWDV -> SetDoImprovedDCAV0DauPropagation(kTRUE);
    //taskWDV -> SetDoImprovedCascadePosition(kTRUE);
    //taskWDV->SetOnlyCombineMCTrue(kTRUE);
    //taskWDV -> SetCentralityInterval(0,0.5);

    gInterpreter->LoadMacro("AliXiStarpp13TeVDevelEventCollection.cxx+g");
    gInterpreter->LoadMacro("AliXiStarpp13TeVDevel.cxx+g");
    // AddTask
    AliXiStarpp13TeVDevel *myTask = reinterpret_cast<AliXiStarpp13TeVDevel*>(gInterpreter->ExecuteMacro(Form("AddTaskXiStarpp13TeVDevel.C(%d,%i,%d,%d,%d,%d)",MCcase,CutList,DevelopmentMode,HMTrigger,PIDOption,SetSystematic)));
    //AliXiStarpp13TeVDevel *myTask = reinterpret_cast<AliXiStarpp13TeVDevel*>(gInterpreter->ExecuteMacro(Form("AddTaskXiStarpp13TeV.C(%d,%d,%i,%d)",MCcase,CutList,DevelopmentMode)));
#else
    // ROOT 5 MODE
    //
    // Physics Selection
    
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(MCcase);
    if(!physSelTask) {
        Printf("no physSelTask");
        return;
    }

    // Multiplicity selection
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask *MultSlection = AddTaskMultSelection();

    // PID Response
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTask *fPIDResponse = AddTaskPIDResponse(MCcase); //! PID response object

    gROOT->LoadMacro("AliXiStarpp13TeVDevelEventCollection.cxx+g");
    gROOT->LoadMacro("AliXiStarpp13TeVDevel.cxx+g");
    // Add Task
    gROOT->LoadMacro("AddTaskXiStarpp13TeVDevel.C");
    AliXiStarpp13TeVDevel *myTask = AddTaskXiStarpp13TeVDevel(MCcase,CutList,DevelopmentMode,HMTrigger,PIDOption,SetSystematic);
#endif
    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(10);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);
    if(strcmp(mode,"local")==0) {// local pc mode; gets AOD or ESD files from my machine
        TChain *chain = new TChain("ESDTree");
        chain = CreateChain(dataset); // for KIAF use
        mgr->StartAnalysis("local",chain);
    } else {// alien mode
        gSystem->Load("libNetx.so") ;
        gSystem->Load("libgapiUI.so");
        gSystem->Load("libRAliEn.so");
        AliAnalysisAlien *plugin = new AliAnalysisAlien();

        //gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
        plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

        plugin->SetAnalysisSource("AliXiStarpp13TeVDevelEventCollection.cxx AliXiStarpp13TeVDevel.cxx");
        plugin->SetAdditionalLibs("AliXiStarpp13TeVDevelEventCollection.h AliXiStarpp13TeVDevelEventCollection.cxx AliXiStarpp13TeVDevel.h AliXiStarpp13TeVDevel.cxx");
        plugin->SetAliPhysicsVersion("vAN-20180909-1");
        plugin->SetAPIVersion("V1.1x");
        // Data
        if(!MCcase) plugin->SetRunPrefix("000");   // real or MC data
        Int_t totruns=0;
        // Declare input data to be processed.
        if(!MCcase) { // Real data
            if(Production->Contains("15f")) { // ESDs
                plugin->SetGridDataDir("/alice/data/2015/LHC15f");
                plugin->SetDataPattern("pass2/*AliESDs.root");
                totruns += AddGoodRuns(plugin,"LHC15f");
            }
            if(Production->Contains("16k")) { // ESDs
                plugin->SetGridDataDir("/alice/data/2016/LHC16k");
                plugin->SetDataPattern("pass2/*AliESDs.root");
                totruns += AddGoodRuns(plugin,"LHC16k");
            }
            if(Production->Contains("17k")) { // AODs
                plugin->SetGridDataDir("/alice/data/2017/LHC17k");
                plugin->SetDataPattern("/pass1/*AliESDs.root");
                totruns += AddGoodRuns(plugin,"LHC17k");
            }
        } else {// MC data
            if(Production->Contains("16k")) { // AODs
                plugin->SetGridDataDir("/alice/sim/2018/LHC18c6b2");
                //plugin->SetGridDataDir("/alice/sim/2017/LHC17d20a1");
                plugin->SetDataPattern("*AliESDs.root");
                totruns += AddGoodRuns(plugin,"LHC16k", MCcase);
            }
        }

        plugin->SetSplitMaxInputFileNumber(2000);
        plugin->SetExecutable("myTask.sh");
        plugin->SetTTL(10000);
        plugin->SetJDLName("myTask.jdl");
        //plugin->SetOutputToRunNo(kTRUE);
        plugin->SetKeepLogs(kTRUE);
        plugin->SetMaxMergeStages(3);
        plugin->SetMaxMergeFiles(100);
        plugin->SetMergeViaJDL(kTRUE);
        plugin->SetCheckCopy(kFALSE);

        plugin->SetGridWorkingDir(working_directory);
        plugin->SetGridOutputDir("output");

        plugin->SetOverwriteMode(kTRUE);
        plugin->SetUser("blim");

        mgr->SetGridHandler(plugin);
        
        // speficy on how many files you want to run
        if(strcmp(mode,"test")==0)plugin->SetNtestFiles(1);

        plugin->SetRunMode(mode);
        mgr->StartAnalysis("grid");
        
    }

}
Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,Bool_t MCcase) {
    //
    // Adds good runs from the Monalisa Run Condition Table
    //Int_t SetRunNumber = 6;


    Int_t nruns=0,ngoodruns=0;
    Int_t runlist[200];
#if !defined (__CINT__) || defined (__CLING__)
    // ROOT 6 MODE
    std::fill_n(runlist, 200, 0);
#else
    // ROOT 5 MODE
#endif

    if(lhcPeriod=="LHC15f") {

        if(MCcase) {
            //mc list
            nruns=16;
            Int_t runlist_2[16]= {146746, 146747, 146748, 146801, 146802, 146803, 146804, 146805, 146806, 146807, 146817, 146824, 146856, 146858, 146859, 146860};
        } else {
            //data list
            nruns=54;
            Int_t runlist_2[54] ={225000, 225011, 225016, 225026, 225031, 225035, 225037, 225041, 225043, 225050, 225051, 225052, 225106, 225305, 225307, 225309, 225313, 225314, 225322, 225576, 225578, 225579, 225586, 225587, 225707, 225708, 225709, 225710, 225716, 225717, 225719, 225753, 225757, 225762, 225763, 225766, 225768, 226062, 226170, 226220, 226225, 226444, 226445, 226452, 226466, 226468, 226472, 226476, 226483, 226495, 226500, 225310, 225315, 225705 };
            for (unsigned i = 0; i < nruns; i++) {
                runlist[i] = runlist_2[i];
            }
        }
        for(unsigned int k=0; k<nruns; k++) {
            if(MCcase) {
                TString *MCRunName = new TString();
                *MCRunName += runlist[k];
                plugin->AddRunNumber(MCRunName->Data());
            }
            else plugin->AddRunNumber(runlist[k]);

            ngoodruns++;
        }
        plugin->SetNrunsPerMaster(1);
    }
    if(lhcPeriod=="LHC16k") {
        if(MCcase) {
            //mc list
            nruns=1; // 2nd
            Int_t runlist_2[1] = {257605};
            //Int_t runlist_2[194]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
            for (unsigned i = 0; i < nruns; i++) {
                runlist[i] = runlist_2[i];
            }
        } else {
            //data list
            nruns=1; 
            Int_t runlist_2[1] = {257605};
            //Int_t runlist_2[4] = {257605, 257604, 257601, 257595};
            //Int_t runlist_2[104] ={258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014};
            //Int_t runlist_2[194]= {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
            // Total 194 runs
            //(52) 258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 
            
            //(52) 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 

            //(25) 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474
            
            //--merge done above--

            //(27) 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257144, 257141, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257082, 257080, 
            //--run Done--above---
            // (38) 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504
            for (unsigned i = 0; i < nruns; i++) {
                runlist[i] = runlist_2[i];
            }
        }
        for(unsigned int k=0; k<nruns; k++) {
            if(MCcase) {
                TString *MCRunName = new TString();
                *MCRunName += runlist[k];
                plugin->AddRunNumber(MCRunName->Data());
            }
            else plugin->AddRunNumber(runlist[k]);

            ngoodruns++;
        }
        plugin->SetNrunsPerMaster(1);
    }
    if(lhcPeriod=="LHC17k") {

        if(MCcase) {
            //mc list
            nruns=16;
            Int_t runlist_2[16]= {146746, 146747, 146748, 146801, 146802, 146803, 146804, 146805, 146806, 146807, 146817, 146824, 146856, 146858, 146859, 146860};
        } else {
            //data list
            nruns=106; 
            Int_t runlist_2[106]={276508, 276507, 276506, 276462, 276439, 276438, 276437, 276435, 276429, 276351, 276348, 276302, 276297, 276294, 276292, 276290, 276259, 276257, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276140, 276135, 276104, 276102, 276099, 276098, 276097, 275847, 275664, 275661, 275650, 275648, 275647, 275624, 275623, 275622, 275621, 275617, 275612, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275456, 275453, 275452, 275448, 275443, 275406, 275404, 275401, 275372, 275369, 275361, 275360, 275333, 275332, 275328, 275326, 275324, 275322, 275314, 275283, 275247, 275246, 275245, 275239, 275188, 275184, 275180, 275177, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275068, 275067, 274979, 274978, 274886, 274882, 274878, 274877, 274822, 274821, 274815, 274806, 274803, 274802, 274801, 274708, 274690};
            //(total 106) 276508, 276507, 276506, 276462, 276439, 276438, 276437, 276435, 276429, 276351, 276348, 276302, 276297, 276294, 276292, 276290, 276259, 276257, 276230, 276205, 276178, 276177, 276170, 276169, 276166, 276145, 276140, 276135, 276104, 276102, 276099, 276098, 276097, 275847, 275664, 275661, 275650, 275648, 275647, 275624, 275623, 275622, 275621, 275617, 275612, 275559, 275558, 275515, 275472, 275471, 275467, 275459, 275457, 275456, 275453, 275452, 275448, 275443, 275406, 275404, 275401, 275372, 275369, 275361, 275360, 275333, 275332, 275328, 275326, 275324, 275322, 275314, 275283, 275247, 275246, 275245, 275239, 275188, 275184, 275180, 275177, 275174, 275173, 275151, 275150, 275149, 275076, 275075, 275073, 275068, 275067, 274979, 274978, 274886, 274882, 274878, 274877, 274822, 274821, 274815, 274806, 274803, 274802, 274801, 274708, 274690 
        
            for (unsigned i = 0; i < nruns; i++) {
                runlist[i] = runlist_2[i];
            }
        }
        for(unsigned int k=0; k<nruns; k++) {
            if(MCcase) {
                TString *MCRunName = new TString();
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
//______________________________________________________________________________
TChain *CreateChain(const char *fileName)
{
    TString *treename;
    treename = new TString("esdTree");

    TChain* chainGood = new TChain(treename->Data());
    TChain* chainNull = new TChain(treename->Data());

    int counter=0;
    char line[500] = "";
    TString *name;

    ifstream inputstream(fileName);
    if(!inputstream) {
        cout<<"input file not found"<<endl;
        return chainNull;
    }

    while(!inputstream.eof()) {
        inputstream.getline(line,500);
        name = new TString(line);

        name->Append("/AliESDs.root");

        if(!inputstream.eof()) chainGood->Add(name->Data());

        counter++;
    }

    inputstream.close();

    return chainGood;
}
