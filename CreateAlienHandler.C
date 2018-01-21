AliAnalysisGrid* CreateAlienHandler(bool MCcase)
{
    // Check if user has a valid token, otherwise make one. This has limitations.
    // One can always follow the standard procedure of calling alien-token-init then
    //   source /tmp/gclient_env_$UID in the current shell.
    //   if (!AliAnalysisGrid::CreateToken()) return NULL;
    
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    
    //  plugin->SetUser("Your alien username here");
    plugin->SetUser("blim");
    
    // 2016 Jan 6 : gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
    // gSystem->Setenv("alien_CLOSE_SE","ALICE::Trieste::SE");
    // gSystem->Setenv("alien_CLOSE_SE","working_disk_SE");
    plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    // Overwrite all generated files, datasets and output results from a previous session
    plugin->SetOverwriteMode(kTRUE);
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode("full");  // VERY IMPORTANT
    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    //    plugin->SetROOTVersion("v5-34-05");
    //    plugin->SetAliROOTVersion("v5-04-80-AN");
    
 //   plugin->SetROOTVersion("v5-34-08-6");
 //   plugin->SetAliROOTVersion("vAN-20141030");
    
    
   // plugin->SetAliROOTVersion("v5-06-19");
   // plugin->SetAliPhysicsVersion("vAN-20150601");
    
    //plugin->SetAliROOTVersion("v5-06-15");
    //plugin->SetAliPhysicsVersion("v5-06-15-01"); //work for kMB untill June 23 2016.  //for data
    plugin->SetAliPhysicsVersion("vAN-20171030-1");

    
    if(!MCcase) plugin->SetRunPrefix("000");   // real or MC data
    
   
     TString *Production=new TString("16k");
   
    
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
        if(Production->Contains("16k")){// ESDs
            plugin->SetGridDataDir("/alice/data/2016/LHC16k");
            plugin->SetDataPattern("pass1/*AliESDs.root");
            totruns += AddGoodRuns(plugin,"LHC16k");
        }
    }else {// MC data
        if(Production->Contains("15f")){// AODs
            plugin->SetGridDataDir("/alice/sim/2015/LHC18a2c"); //1a : Pythia8, 1b : Phojet
            plugin->SetDataPattern("*AliESDs.root");
            totruns += AddGoodRuns(plugin,"LHC15f",MCcase);
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
    
    
    
    // Define alien work directory where all files will be copied. Relative to alien $HOME.

  //  if(!MCcase) plugin->SetGridWorkingDir("pp276TeV_LHC11a_withoutSDD-24Runs-EM10-2017-woPID");
    if(!MCcase) plugin->SetGridWorkingDir("pp13TeV_LHC16kf_goodrun_20180120");

 //   else plugin->SetGridWorkingDir("pp276TeV_LHC11a_withoutSDD-LHC12f1a-Pythia8-TPCPID");
    else plugin->SetGridWorkingDir("pp13TeV_LHC16kf_goodrun_20180118");

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
    //plugin->SetMaxMergeFiles(50);
    plugin->SetOneStageMerging(kFALSE);
    plugin->SetMaxMergeStages(2);
    plugin->SetMergeViaJDL();
    //plugin->SetUseSubmitPolicy();
    
    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliXiStarpp13TeVEventCollection.cxx AliXiStarpp.cxx");
    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliXiStarpp13TeVEventCollection.h AliXiStarppEventCollection.cxx AliXiStarpp.h AliXiStarpp.cxx");
    
    // No need for output file names. Procedure is automatic.
    plugin->SetDefaultOutputs(kTRUE);
    //plugin->SetOutputFiles("MyOutput.root");
    // No need define the files to be archived. Note that this is handled automatically by the plugin.
    //plugin->SetOutputArchive("log_archive.zip:stdout,stderr");
    // Set a name for the generated analysis macro (default MyAnalysis.C) Make this unique !
    plugin->SetAnalysisMacro("GeneratedMacro.C");
    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore). The optimum for an analysis
    // is correlated with the run time - count few hours TTL per job, not minutes !
    plugin->SetSplitMaxInputFileNumber(100); //50 // 2015Feb
    // Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
    
    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(95); //97
    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(40000);
    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName("runCode.jdl");
    // Optionally modify job price (default 1)
    plugin->SetPrice(1);
    // Optionally modify split mode (default 'se')
    plugin->SetSplitMode("se");
    
    return plugin;
    
}
//___________________________________________________________________-
Int_t AddGoodRuns(AliAnalysisAlien* plugin,TString lhcPeriod,bool MCcase=kFALSE) {
    //
    // Adds good runs from the Monalisa Run Condition Table
    //Int_t SetRunNumber = 6;

    
    Int_t nruns=0,ngoodruns=0;
     if(lhcPeriod=="LHC16k") {
         
         if(MCcase){
         //mc list
         nruns=48;
         Int_t runlist[16]={226500, 226495, 226483, 226476, 226472, 226468, 226466, 226452, 226445, 226444, 226225, 226220, 226170, 226062, 225768, 225766, 225763, 225762, 225757, 225753, 225719, 225717, 225716, 225710, 225709, 225708, 225707, 225587, 225586, 225579, 225578, 225576, 225322, 225314, 225313, 225309, 225307, 225305, 225106, 225052, 225051, 225050, 225043, 225041, 225037, 225035, 225031, 225026};
         }else{
             //data list
             nruns=202; // 51
	     Int_t runlist[202]={258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258048, 258045, 258042, 258041, 258039, 258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, 257979, 257963, 257960, 257958, 257957, 257939, 257937, 257936, 257893, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682, 257644, 257642, 257636, 257635, 257632, 257630, 257606, 257605, 257604, 257601, 257595, 257594, 257592, 257590, 257588, 257587, 257566, 257562, 257561, 257560, 257541, 257540, 257539, 257537, 257531, 257530, 257492, 257491, 257490, 257488, 257487, 257474, 257468, 257457, 257433, 257364, 257358, 257330, 257322, 257320, 257318, 257260, 257224, 257209, 257206, 257204, 257145, 257144, 257142, 257141, 257140, 257139, 257138, 257137, 257136, 257100, 257095, 257092, 257086, 257084, 257083, 257082, 257080, 257077, 257028, 257026, 257021, 257012, 257011, 256944, 256942, 256941, 256697, 256695, 256694, 256692, 256691, 256684, 256681, 256677, 256676, 256658, 256620, 256619, 256592, 256591, 256589, 256567, 256565, 256564, 256562, 256561, 256560, 256557, 256556, 256554, 256552, 256514, 256512, 256510, 256506, 256504};
         }
         
        // nruns=55; //pass1
        // Int_t runlist[55]={146860, 146859, 146858, 146856, 146824, 146817, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146287, 146282, 146277, 146273, 146272, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146148, 146147, 146141, 146099, 146079, 146072, 146071, 146027, 146026, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145292, 145290, 145289, 145288}; //146402, 146369,146292
    
         
         //167988 // there are only 10 events, no ESD files
     for(Int_t k=0;k<nruns;k++){
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

