/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//
//  This class is used to reconstruct the neutral Xi(1530) resonance in pp 13TeV System.
//  This class essentially combines charged Xi candidates from the Xi Vertexer
//  with primary charged pions.
//  Multiplicity Information added(Bong-Hwi)
//
//  Original author: Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//  Modified by: Jihye Song (jihye.song@cern.ch)
//  Last Modified by: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
//
//  Last Modified Date: 2018/09/25
//
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <math.h>
#include "TChain.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjString.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"


#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliVertex.h"
#include "AliMultSelection.h"
//#include "AliOADBCentrality.h"
//#include "AliOADBPhysicsSelection.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
//#include "AliTriggerAnalysis.h"

#include "AliESDcascade.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"

#include "AliXiStarpp13TeVDevel.h"

#define PI 3.1415927


// Author: Bong-Hwi Lim

ClassImp(AliXiStarpp13TeVDevel)

//________________________________________________________________________
AliXiStarpp13TeVDevel::AliXiStarpp13TeVDevel():
    AliAnalysisTaskSE(),
    fname(0),
    fESD(0x0),
    fOutputList(0x0),
    fTrackCut(0x0),
    fPIDResponse(0x0),

    fEC(0x0),
    fEvt(0x0),
    fDevelopeMode(0),
    fHMTrigger(0),
    fPIDOption(0),
    fSetSystematic(0),

    fTempStruct(0x0),
    fZvertexBins(0),
    fEventsToMix(0),
    fMultBins(0),
    fMCcase(0),
    fEventCounter(0),
    fEventNumber(0),

    fMaxDecayLength(0),
    fMassWindow(0),
    fTrueMassPr(0),
    fTrueMassPi(0),
    fTrueMassK(0),
    fTrueMassLam(0),
    fTrueMassXi(0),

    fESDTrack4(0x0),
    fXiTrack(0x0),
    fCutList(0),
    htotalEvent(0),
    htriggered_INELg0_tracklet(0),
    htriggered_CINT7_tracklet(0),
    htriggered_CINT7_VOM(0),
    htriggered_AliMult_VOM(0),

    fVertexDist1(0),
    fVertexDist3(0),
    fDCADist(0),
    fMultDist3d(0),
    fMultDist1(0),
    fMultDist2(0),
    fMultDist3(0),
    fMultDist4(0),
    fMultDist5(0),
    fMultDist_pp(0),
    hEventSelecInfo(0x0),
    hNumberOfEvent(0),

    fPtDist(0),
    fPhiDist(0),
    fEtaDist(0),
    fXiStarYDist(0),
    fQAXiStarYDist(0),

    fXiStarYDistMC(0),
    fXiYDistMC1(0),
    fXiStarYDistMC1(0),
    fXiYDistMCout(0),
    fXiStarYDistMCout(0),
    fCutEvents(0),

    fTPCNcls_p(0),
    fTPCNcls_pi1(0),
    fTPCNcls_pi2(0),
    fTPCNcls_pi3(0),

    fQATPCNcls_p(0),
    fQATPCNcls_pi1(0),
    fQATPCNcls_pi2(0),
    fQATPCNcls_pi3(0),

    fQATPCNcls_p_L(0),
    fQATPCNcls_pi1_L(0),
    fQATPCNcls_pi2_L(0),
    fQATPCNcls_pi3_L(0),

    fQATPCNcls_p_T(0),
    fQATPCNcls_pi1_T(0),
    fQATPCNcls_pi2_T(0),
    fQATPCNcls_pi3_T(0),

    fDCADist_p(0),
    fDCADist_pi1(0),
    fDCADist_pi2(0),

    fQADCADist_p(0),
    fQADCADist_pi1(0),
    fQADCADist_pi2(0),

    fQADCADist_p_L(0),
    fQADCADist_pi1_L(0),
    fQADCADist_pi2_L(0),

    fQADCADist_p_T(0),
    fQADCADist_pi1_T(0),
    fQADCADist_pi2_T(0),

    fDCADist_lambda(0),
    fDCADist_3rd_pi(0),
    fDCADist_pi_p(0),
    fDCADist_pi_lambda(0),

    fQADCADist_lambda(0),
    fQADCADist_3rd_pi(0),
    fQADCADist_pi_p(0),
    fQADCADist_pi_lambda(0),

    fQADCADist_lambda_L(0),
    fQADCADist_3rd_pi_L(0),
    fQADCADist_pi_p_L(0),
    fQADCADist_pi_lambda_L(0),

    fQADCADist_lambda_T(0),
    fQADCADist_3rd_pi_T(0),
    fQADCADist_pi_p_T(0),
    fQADCADist_pi_lambda_T(0),

    fCosPA_lambda(0),
    fCosPA_Xi(0),
    fQACosPA_lambda(0),
    fQACosPA_Xi(0),
    fQACosPA_lambda_L(0),
    fQACosPA_Xi_L(0),
    fQACosPA_lambda_T(0),
    fQACosPA_Xi_T(0),
    hXiInvMass(0),
    hQAXiInvMass(0),

    hTPCPID(0),
    hTPCPIDpi(0),
    hTPCPIDk(0),
    hTPCPIDp(0),


    hdEdxProton(0),
    hdEdxPion1(0),
    hdEdxPion2(0),

    hdEdxProtonAfter(0),
    hdEdxPion1After(0),
    hdEdxPion2After(0),

    hNSig3rdPion(0),
    hQANSig3rdPion(0),

    hTPCNSigProton(0),
    hTPCNSigPion1(0),
    hTPCNSigPion2(0),

    hQATPCNSigProton(0),
    hQATPCNSigPion1(0),
    hQATPCNSigPion2(0),

    hMCinputTotalXiStar1(0),
    hMCinputTotalXiStarbar1(0),
    hMCinputTotalXi1(0),
    hMCinputTotalXibar1(0),

    hMCinputTotalXiStar3(0),
    hMCinputTotalXiStarbar3(0),
    hMCinputTotalXi3(0),
    hMCinputTotalXibar3(0)

{
    // Default Constructor
    for (Int_t i = 0; i < 21; i++) {
        fCovMatrix[i] = -99999.;
        if (i < 12) fMultLimits[i] = 0;
    }
    for (Int_t i = 0; i < kNCuts; i++) {
        fDecayParameters[i] = 0;
        for (Int_t j = 0; j < kNCutVariations; j++) {
            fCutValues[j][i] = 0;
        }
    }
    //
    for (Int_t cv = 0; cv < kNCutVariations; cv++) {
        CutVar[cv].fXi = 0x0;
        CutVar[cv].fXibar = 0x0;
        CutVar[cv].fXiMinusPiPlus = 0x0;
        CutVar[cv].fXiMinusPiMinus = 0x0;
        CutVar[cv].fXiPlusPiPlus = 0x0;
        CutVar[cv].fXiPlusPiMinus = 0x0;
        //
        CutVar[cv].fXiMinusPiPlusbkg = 0x0;
        CutVar[cv].fXiMinusPiMinusbkg = 0x0;
        CutVar[cv].fXiPlusPiPlusbkg = 0x0;
        CutVar[cv].fXiPlusPiMinusbkg = 0x0;
        //
        CutVar[cv].fMCrecXi = 0x0;
        CutVar[cv].fMCrecXibar = 0x0;
        CutVar[cv].fMCrecXiMinusPiPlus = 0x0;
        CutVar[cv].fMCrecXiPlusPiMinus = 0x0;
    }

}
//________________________________________________________________________
AliXiStarpp13TeVDevel::AliXiStarpp13TeVDevel(const char *name, Int_t CutListOption)
    : AliAnalysisTaskSE(name),
      fname(name),
      fESD(0x0),
      fOutputList(0x0),
      fTrackCut(0x0),
      fPIDResponse(0x0),
      fEC(0x0),
      fEvt(0x0),
      fDevelopeMode(kFALSE),
      fHMTrigger(kFALSE),
      fPIDOption(kTRUE),
      fTempStruct(0x0),
      fZvertexBins(0),
      fEventsToMix(0),
      fMultBins(0),
      fMCcase(kFALSE),
      fEventCounter(0),
      fEventNumber(0),
      fMaxDecayLength(0),
      fMassWindow(0),
      fTrueMassPr(0),
      fTrueMassPi(0),
      fTrueMassK(0),
      fTrueMassLam(0),
      fTrueMassXi(0),
      fSetSystematic(kTRUE),

      fESDTrack4(0x0),
      fXiTrack(0x0),
      fCutList(CutListOption),
      htotalEvent(0),
      htriggered_INELg0_tracklet(0),
      htriggered_CINT7_tracklet(0),
      htriggered_CINT7_VOM(0),
      htriggered_AliMult_VOM(0),

      fVertexDist1(0),
      fVertexDist3(0),
      fDCADist(0),
      fMultDist3d(0),
      fMultDist1(0),
      fMultDist2(0),
      fMultDist3(0),
      fMultDist4(0),
      fMultDist5(0),
      fMultDist_pp(0),
      hEventSelecInfo(0x0),
      hNumberOfEvent(0),

      fPtDist(0),
      fPhiDist(0),
      fEtaDist(0),
      fXiStarYDist(0),
      fQAXiStarYDist(0),

      fXiStarYDistMC(0),
      fXiYDistMC1(0),
      fXiStarYDistMC1(0),
      fXiYDistMCout(0),
      fXiStarYDistMCout(0),
      fCutEvents(0),

      fTPCNcls_p(0),
      fTPCNcls_pi1(0),
      fTPCNcls_pi2(0),
      fTPCNcls_pi3(0),

      fQATPCNcls_p(0),
      fQATPCNcls_pi1(0),
      fQATPCNcls_pi2(0),
      fQATPCNcls_pi3(0),

      fQATPCNcls_p_L(0),
      fQATPCNcls_pi1_L(0),
      fQATPCNcls_pi2_L(0),
      fQATPCNcls_pi3_L(0),

      fQATPCNcls_p_T(0),
      fQATPCNcls_pi1_T(0),
      fQATPCNcls_pi2_T(0),
      fQATPCNcls_pi3_T(0),

      fDCADist_p(0),
      fDCADist_pi1(0),
      fDCADist_pi2(0),

      fQADCADist_p(0),
      fQADCADist_pi1(0),
      fQADCADist_pi2(0),

      fQADCADist_p_L(0),
      fQADCADist_pi1_L(0),
      fQADCADist_pi2_L(0),

      fQADCADist_p_T(0),
      fQADCADist_pi1_T(0),
      fQADCADist_pi2_T(0),

      fDCADist_lambda(0),
      fDCADist_3rd_pi(0),
      fDCADist_pi_p(0),
      fDCADist_pi_lambda(0),

      fQADCADist_lambda(0),
      fQADCADist_3rd_pi(0),
      fQADCADist_pi_p(0),
      fQADCADist_pi_lambda(0),

      fQADCADist_lambda_L(0),
      fQADCADist_3rd_pi_L(0),
      fQADCADist_pi_p_L(0),
      fQADCADist_pi_lambda_L(0),

      fQADCADist_lambda_T(0),
      fQADCADist_3rd_pi_T(0),
      fQADCADist_pi_p_T(0),
      fQADCADist_pi_lambda_T(0),

      fCosPA_lambda(0),
      fCosPA_Xi(0),
      fQACosPA_lambda(0),
      fQACosPA_Xi(0),
      fQACosPA_lambda_L(0),
      fQACosPA_Xi_L(0),
      fQACosPA_lambda_T(0),
      fQACosPA_Xi_T(0),
      hXiInvMass(0),
      hQAXiInvMass(0),

      hTPCPID(0),
      hTPCPIDpi(0),
      hTPCPIDk(0),
      hTPCPIDp(0),


      hdEdxProton(0),
      hdEdxPion1(0),
      hdEdxPion2(0),

      hdEdxProtonAfter(0),
      hdEdxPion1After(0),
      hdEdxPion2After(0),

      hNSig3rdPion(0),
      hQANSig3rdPion(0),

      hTPCNSigProton(0),
      hTPCNSigPion1(0),
      hTPCNSigPion2(0),

      hQATPCNSigProton(0),
      hQATPCNSigPion1(0),
      hQATPCNSigPion2(0),

      hMCinputTotalXiStar1(0),
      hMCinputTotalXiStarbar1(0),
      hMCinputTotalXi1(0),
      hMCinputTotalXibar1(0),

      hMCinputTotalXiStar3(0),
      hMCinputTotalXiStarbar3(0),
      hMCinputTotalXi3(0),
      hMCinputTotalXibar3(0)

{
    // Main Constructor
    for (Int_t i = 0; i < 21; i++) {
        fCovMatrix[i] = -99999.;
        if (i < 12) fMultLimits[i] = 0;
    }
    for (Int_t i = 0; i < kNCuts; i++) {
        fDecayParameters[i] = 0;
        for (Int_t j = 0; j < kNCutVariations; j++) {
            fCutValues[j][i] = 0;
        }
    }
    //
    for (Int_t cv = 0; cv < kNCutVariations; cv++) {
        CutVar[cv].fXi = 0x0;
        CutVar[cv].fXibar = 0x0;
        CutVar[cv].fXiMinusPiPlus = 0x0;
        CutVar[cv].fXiMinusPiMinus = 0x0;
        CutVar[cv].fXiPlusPiPlus = 0x0;
        CutVar[cv].fXiPlusPiMinus = 0x0;
        //
        CutVar[cv].fXiMinusPiPlusbkg = 0x0;
        CutVar[cv].fXiMinusPiMinusbkg = 0x0;
        CutVar[cv].fXiPlusPiPlusbkg = 0x0;
        CutVar[cv].fXiPlusPiMinusbkg = 0x0;
        //
        CutVar[cv].fMCrecXi = 0x0;
        CutVar[cv].fMCrecXibar = 0x0;
        CutVar[cv].fMCrecXiMinusPiPlus = 0x0;
        CutVar[cv].fMCrecXiPlusPiMinus = 0x0;
    }

    // Define output for Tree
    DefineInput (0, TChain::Class()); //jisong added

    // Define output slots here
    // Output slot #1
//   DefineOutput(0, TTree::Class());
    DefineOutput(1, TList::Class());

}
//________________________________________________________________________
AliXiStarpp13TeVDevel::AliXiStarpp13TeVDevel(const AliXiStarpp13TeVDevel &obj )
    : AliAnalysisTaskSE(obj.fname),
      fname(obj.fname),
      fESD(obj.fESD),
      fOutputList(obj.fOutputList),
      fTrackCut(obj.fTrackCut),
      fPIDResponse(obj.fPIDResponse),
      fEC(obj.fEC),
      fEvt(obj.fEvt),
      fTempStruct(obj.fTempStruct),
      fZvertexBins(obj.fZvertexBins),
      fEventsToMix(obj.fEventsToMix),
      fMultBins(obj.fMultBins),
      fMultLimits(),
      fMCcase(obj.fMCcase),
      fEventCounter(obj.fEventCounter),
      fEventNumber(obj.fEventNumber),
      fMaxDecayLength(obj.fMaxDecayLength),
      fMassWindow(obj.fMassWindow),
      fTrueMassPr(obj.fTrueMassPr),
      fTrueMassPi(obj.fTrueMassPi),
      fTrueMassK(obj.fTrueMassK),
      fTrueMassLam(obj.fTrueMassLam),
      fTrueMassXi(obj.fTrueMassXi),

      fESDTrack4(obj.fESDTrack4),
      fXiTrack(obj.fXiTrack),
      fCutList(obj.fCutList)

{
    // Copy constructor
    for (Int_t i = 0; i < 21; i++) {
        fCovMatrix[i] = obj.fCovMatrix[i];
        if (i < 12) fMultLimits[i] = obj.fMultLimits[i];
    }
    for (Int_t i = 0; i < kNCuts; i++) {
        fDecayParameters[i] = obj.fDecayParameters[i];
        for (Int_t j = 0; j < kNCutVariations; j++) {
            fCutValues[j][i] = obj.fCutValues[j][i];
        }
    }

}
//________________________________________________________________________

AliXiStarpp13TeVDevel &AliXiStarpp13TeVDevel::operator=(const AliXiStarpp13TeVDevel &obj )
{
    // Assignment operator
    if (this == &obj )
        return *this;

    fname = obj.fname;
    fESD = obj.fESD;
    fOutputList = obj.fOutputList;
    fTrackCut = obj.fTrackCut;
    fPIDResponse = obj.fPIDResponse;

    fEC = obj.fEC;
    fEvt = obj.fEvt;
    fTempStruct = obj.fTempStruct;
    fZvertexBins = obj.fZvertexBins;
    fEventsToMix = obj.fEventsToMix;
    fMultBins = obj.fMultBins;
    for (Int_t i = 0; i < 12; i++) {
        fMultLimits[i] = obj.fMultLimits[i];
    }
    fMCcase = obj.fMCcase;
    fEventCounter = obj.fEventCounter;
    fEventNumber = obj.fEventNumber;
    fMaxDecayLength = obj.fMaxDecayLength;
    fMassWindow = obj.fMassWindow;
    for (Int_t i = 0; i < 21; i++) {
        fCovMatrix[i] = obj.fCovMatrix[i];
    }

    fTrueMassPr = obj.fTrueMassPr;
    fTrueMassPi = obj.fTrueMassPi;
    fTrueMassK = obj.fTrueMassK;
    fTrueMassLam = obj.fTrueMassLam;
    fTrueMassXi = obj.fTrueMassXi;
    fESDTrack4 = obj.fESDTrack4;
    fXiTrack = obj.fXiTrack;
    fCutList = obj.fCutList;


    for (Int_t i = 0; i < kNCuts; i++) {
        fDecayParameters[i] = obj.fDecayParameters[i];
        for (Int_t j = 0; j < kNCutVariations; j++) {
            fCutValues[j][i] = obj.fCutValues[j][i];
        }
    }


    return (*this);
}
//________________________________________________________________________
AliXiStarpp13TeVDevel::~AliXiStarpp13TeVDevel()
{
    // Destructor

    if (fESD) delete fESD;
    if (fOutputList) delete fOutputList;
    if (fTrackCut) delete fTrackCut;
    if (fPIDResponse) delete fPIDResponse;

    if (fEC) delete fEC;
    if (fEvt) delete fEvt;
    if (fTempStruct) delete fTempStruct;
    if (fESDTrack4) delete fESDTrack4;
    if (fXiTrack) delete fXiTrack;

    for (Int_t cv = 0; cv < kNCutVariations; cv++) {
        if (CutVar[cv].fXi) delete CutVar[cv].fXi;
        if (CutVar[cv].fXibar) delete CutVar[cv].fXibar;
        if (CutVar[cv].fXiMinusPiPlus) delete CutVar[cv].fXiMinusPiPlus;
        if (CutVar[cv].fXiMinusPiMinus) delete CutVar[cv].fXiMinusPiMinus;
        if (CutVar[cv].fXiPlusPiPlus) delete CutVar[cv].fXiPlusPiPlus;
        if (CutVar[cv].fXiPlusPiMinus) delete CutVar[cv].fXiPlusPiMinus;
        //
        if (CutVar[cv].fXiMinusPiPlusbkg) delete CutVar[cv].fXiMinusPiPlusbkg;
        if (CutVar[cv].fXiMinusPiMinusbkg) delete CutVar[cv].fXiMinusPiMinusbkg;
        if (CutVar[cv].fXiPlusPiPlusbkg) delete CutVar[cv].fXiPlusPiPlusbkg;
        if (CutVar[cv].fXiPlusPiMinusbkg) delete CutVar[cv].fXiPlusPiMinusbkg;
        //
        if (CutVar[cv].fMCrecXi) delete CutVar[cv].fMCrecXi;
        if (CutVar[cv].fMCrecXibar) delete CutVar[cv].fMCrecXibar;
        if (CutVar[cv].fMCrecXiMinusPiPlus) delete CutVar[cv].fMCrecXiMinusPiPlus;
        if (CutVar[cv].fMCrecXiPlusPiMinus) delete CutVar[cv].fMCrecXiPlusPiMinus;
    }

}
//________________________________________________________________________
void AliXiStarpp13TeVDevel::XiStarInit()
{
    //
    //Inits cuts and analysis settings
    //
    fEventCounter = 0; // event counter initialization
    if (fDevelopeMode)std::cout << "AliXiStarpp13TeVDevel XiStarInit() call" << std::endl;
    if (fDevelopeMode)std::cout << "MC Mode?: " << fMCcase << std::endl;

    ///////////////////////////////////////////////
    // Track Cuts for ESD analysis
    fTrackCut = new AliESDtrackCuts();
    fTrackCut->SetPtRange(.15, 1000);
    fTrackCut->SetAcceptKinkDaughters(kFALSE);
    //fTrackCut->SetMinNClustersTPC(50);
    fTrackCut->SetRequireTPCRefit(kTRUE);
    fTrackCut->SetMaxChi2PerClusterTPC(4); //From Enrico

    ////////////////////////////////////////////////

    fZvertexBins = 20;
    fMultBins = 11;// This must also be set in AliXiStarpp13TeVDevel.h
    if (fMCcase) fEventsToMix = 0;
    else fEventsToMix = 40; // original 40 jisong

    // multiplicity edges for event mixing bins
    fMultLimits[0] = 0, fMultLimits[1] = 5, fMultLimits[2] = 10, fMultLimits[3] = 15, fMultLimits[4] = 20, fMultLimits[5] = 25;
    fMultLimits[6] = 30, fMultLimits[7] = 35, fMultLimits[8] = 40, fMultLimits[9] = 45, fMultLimits[10] = 50, fMultLimits[11] = 150;


    fEC = new AliXiStarpp13TeVDevelEventCollection **[fZvertexBins];
    for (unsigned short i = 0; i < fZvertexBins; i++) {

        fEC[i] = new AliXiStarpp13TeVDevelEventCollection *[fMultBins];

        for (unsigned short j = 0; j < fMultBins; j++) {

            fEC[i][j] = new AliXiStarpp13TeVDevelEventCollection(fEventsToMix + 1);
        }
    }

    fTempStruct = new AliXiStarpp13TeVDevelTrackStruct[kNbinsM * 8];
    fESDTrack4 = new AliESDtrack();
    fXiTrack = new AliESDtrack();


    fMaxDecayLength = 100.;
    fMassWindow = 0.007;

    /////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////
    // DecayParameters Key (number represents array index)
    // NclustersTPC: 0=proton, 1=pion first, 2=pion second, 3=pion third
    // DCAVtx: 4=proton, 5=pion first, 6=pion second, 7=lambda, 8=pion third
    // 9 = DCA proton-pion
    // 10 = DCA Lambda-pion
    // 11 = Rxy Lambda
    // 12 = Rxy Xi
    // 13 = Cos PA Lambda
    // 14 = Cos PA Xi

    // Set Standard Reconstruction cut values
    fCutValues[0][0] = 50;
    fCutValues[0][1] = 50;
    fCutValues[0][2] = 50;
    fCutValues[0][3] = 50; // for 2010 cut (origin: 70)
    fCutValues[0][4] = 0.04;
    fCutValues[0][5] = 0.04;
    fCutValues[0][6] = 0.05;
    fCutValues[0][7] = 0.07;
    fCutValues[0][8] = 2.0;
    fCutValues[0][9] = 1.6;
    fCutValues[0][10] = 1.6;
    fCutValues[0][11] = 0.97;
    fCutValues[0][12] = 0.97;


    for (int cv = 1; cv < kNCutVariations; cv++) {
        for (int ct = 0; ct < kNCuts; ct++) {
            fCutValues[cv][ct] = fCutValues[0][ct];
        }
    }

    //systematic variation// Loose
    fCutValues[1][0] = 45;
    fCutValues[1][1] = 45;
    fCutValues[1][2] = 45;
    fCutValues[1][3] = 45;// 60 -> 45
    fCutValues[2][4] = 0.03;
    fCutValues[3][5] = 0.03;
    fCutValues[4][6] = 0.04;
    fCutValues[5][7] = 0.06;
    fCutValues[6][8] = 2.1;
    fCutValues[7][9] = 1.7;
    fCutValues[8][10] = 1.7;
    fCutValues[9][11] = 0.95;
    fCutValues[10][12] = 0.965;

    //systematic variation// tight
    fCutValues[11][0] = 55;
    fCutValues[11][1] = 55;
    fCutValues[11][2] = 55;
    fCutValues[11][3] = 55;// 70 -> 55
    fCutValues[12][4] = 0.104;
    fCutValues[13][5] = 0.104;
    fCutValues[14][6] = 0.08;
    fCutValues[15][7] = 0.1;
    fCutValues[16][8] = 1.0;
    fCutValues[17][9] = 0.94;
    fCutValues[18][10] = 1.41;
    fCutValues[19][11] = 0.99;
    fCutValues[20][12] = 0.985;


    // PDG mass values
    fTrueMassPr = .93827, fTrueMassPi = .13957, fTrueMassK = .493677, fTrueMassLam = 1.11568, fTrueMassXi = 1.32171;

    // The following CovMatrix is set so that PropogateToDCA() ignores track errors. Only used to propagate Xi to third pion for XiStar reconstruction
    for (Int_t i = 0; i < 21; i++) fCovMatrix[i] = 0;
    fCovMatrix[0] = 1, fCovMatrix[2] = 1, fCovMatrix[5] = 1, fCovMatrix[9] = 1, fCovMatrix[14] = 1, fCovMatrix[20] = 1;


}
//________________________________________________________________________
void AliXiStarpp13TeVDevel::UserCreateOutputObjects()
{
    XiStarInit();
    // XiStarInit();// Initialize settings original

    // Create histograms
    fOutputList = new TList();
    fOutputList->SetOwner();
    if (fDevelopeMode)std::cout << "MC Mode?(UserCreate): " << fMCcase << std::endl;
    if (fMCcase) {
        // For trigger efficiency
        htotalEvent = new TH1F("htotalEvent", "Number of Events", 10, 0, 10);
        htotalEvent->GetXaxis()->SetBinLabel(1, "All Events");
        htotalEvent->GetXaxis()->SetBinLabel(2, "IsINELg0");
        htotalEvent->GetXaxis()->SetBinLabel(3, "tracklet in |Eta|<1");
        htotalEvent->GetXaxis()->SetBinLabel(4, "CINT7 triggered");
        htotalEvent->GetXaxis()->SetBinLabel(5, "AliMultiSelection");
        fOutputList->Add(htotalEvent);

        htriggered_INELg0_tracklet = new TH1F("htriggered_INELg0_tracklet", "Triggered event in INEL>0", 11, 0, 11);
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(1, "0 to Inf (MB)");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(2, "0 to 5");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(3, "6 to 10");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(4, "11 to 15");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(5, "16 to 20");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(6, "21 to 25");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(7, "26 to 30");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(8, "31 to 35");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(9, "36 to 40");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(10, "41 to 50");
        htriggered_INELg0_tracklet->GetXaxis()->SetBinLabel(11, "51 more");
        fOutputList->Add(htriggered_INELg0_tracklet);

        htriggered_CINT7_tracklet = new TH1F("htriggered_CINT7_tracklet", "Triggered event in INEL>0", 11, 0, 11);
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(1, "0 to Inf (MB)");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(2, "0 to 5");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(3, "6 to 10");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(4, "11 to 15");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(5, "16 to 20");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(6, "21 to 25");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(7, "26 to 30");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(8, "31 to 35");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(9, "36 to 40");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(10, "41 to 50");
        htriggered_CINT7_tracklet->GetXaxis()->SetBinLabel(11, "51 more");
        fOutputList->Add(htriggered_CINT7_tracklet);

        htriggered_CINT7_VOM = new TH1F("htriggered_CINT7_VOM", "Triggered event in INEL>0", 11, 0, 11);
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(1, "0.0 - 100.0% (MB)");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(2, "70.0 - 100.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(3, "50.0 - 70.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(4, "40.0 - 50.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(5, "30.0 - 40.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(6, "20.0 - 30.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(7, "15.0 - 20.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(8, "10.0 - 15.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(9, "5.0 - 10.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(10, "1.0 - 5.0%");
        htriggered_CINT7_VOM->GetXaxis()->SetBinLabel(11, "0.0 - 1.0%");
        fOutputList->Add(htriggered_CINT7_VOM);

        htriggered_AliMult_VOM = new TH1F("htriggered_AliMult_VOM", "Triggered event in INEL>0", 11, 0, 11);
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(1, "0.0 - 100.0% (MB)");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(2, "70.0 - 100.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(3, "50.0 - 70.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(4, "40.0 - 50.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(5, "30.0 - 40.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(6, "20.0 - 30.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(7, "15.0 - 20.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(8, "10.0 - 15.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(9, "5.0 - 10.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(10, "1.0 - 5.0%");
        htriggered_AliMult_VOM->GetXaxis()->SetBinLabel(11, "0.0 - 1.0%");
        fOutputList->Add(htriggered_AliMult_VOM);

        fMultDist_MCpp = new TH1F("fMultDist_MCpp", "Multiplicity Distribution of PP", 200, 0, 200);
        fMultDist_MCpp->GetXaxis()->SetTitle("Multiplicity Percentile");
        fOutputList->Add(fMultDist_MCpp);

        fMultDist_MCpp_selected = new TH1F("fMultDist_MCpp_selected", "Multiplicity Distribution of PP", 200, 0, 200);
        fMultDist_MCpp_selected->GetXaxis()->SetTitle("Multiplicity Percentile");
        fOutputList->Add(fMultDist_MCpp_selected);
        
    }

    fVertexDist1 = new TH3F("fVertexDist1", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 600, -30, 30);
    fVertexDist1->GetXaxis()->SetTitle("X Vertex (cm)");
    fVertexDist1->GetYaxis()->SetTitle("Y Vertex (cm)");
    fVertexDist1->GetZaxis()->SetTitle("Z Vertex (cm)");
    fOutputList->Add(fVertexDist1);

    fVertexDist3 = new TH3F("fVertexDist3", "Vertex Distribution", 20, -1, 1, 20, -1, 1, 600, -30, 30);
    fVertexDist3->GetXaxis()->SetTitle("X Vertex (cm)");
    fVertexDist3->GetYaxis()->SetTitle("Y Vertex (cm)");
    fVertexDist3->GetZaxis()->SetTitle("Z Vertex (cm)");
    fOutputList->Add(fVertexDist3);

    fDCADist = new TH2F("fDCADist", "DCA distribution", kNbinsM, -.5, kNbinsM * 8 - .5, 100, 0, 10);
    fOutputList->Add(fDCADist);


    fMultDist3d = new TH3F("fMultDist3d", "Multiplicity Distribution", kNbinsM, -.5, kNbinsM - .5, kNbinsM, -.5, kNbinsM - .5, kNbinsM, -.5, kNbinsM - .5);
    fMultDist3d->GetXaxis()->SetTitle("Multiplicity");
    fMultDist3d->GetYaxis()->SetTitle("Positive Multiplicity");
    fMultDist3d->GetZaxis()->SetTitle("Negative Multiplicity");
    fMultDist3d->SetMarkerStyle(kFullCircle);
    fOutputList->Add(fMultDist3d);


    fMultDist1 = new TH1F("fMultDist1", "Multiplicity Distribution", kNbinsM, 0, kNbinsM * 8);
    fMultDist1->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist1);

    fMultDist2 = new TH1F("fMultDist2", "Multiplicity Distribution After Vz selection", kNbinsM, -.5, kNbinsM * 8 - .5);
    fMultDist2->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist2);

    fMultDist3 = new TH1F("fMultDist3", "Multiplicity Distribution After reject pile up from SPD", kNbinsM, -.5, kNbinsM * 8 - .5);
    fMultDist3->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist3);

    fMultDist4 = new TH1F("fMultDist4", "Multiplicity Distribution : Primary NContributor > 1 ", kNbinsM, -.5, kNbinsM * 8 - .5);
    fMultDist4->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist4);

    fMultDist5 = new TH1F("fMultDist5", "Multiplicity Distribution of Xitrack", kNbinsM * 5, -.5, kNbinsM * 5 - .5);
    fMultDist5->GetXaxis()->SetTitle("Multiplicity");
    fOutputList->Add(fMultDist5);

    if (!fHMTrigger) {
        fMultDist_pp = new TH1F("fMultDist_pp", "Multiplicity Distribution of PP", 200, 0, 200);
        fMultDist_pp->GetXaxis()->SetTitle("Multiplicity Percentile");
        fOutputList->Add(fMultDist_pp);
    }
    else {
        fMultDist_pp = new TH1F("fMultDist_pp", "Multiplicity Distribution of PP", 200, 0, 2);
        fMultDist_pp->GetXaxis()->SetTitle("Multiplicity Percentile");
        fOutputList->Add(fMultDist_pp);
    }

    hEventSelecInfo = new TH1F("hEventSelecInfo", "hEventSelecInfo", 10, 0, 10);
    hEventSelecInfo->GetXaxis()->SetBinLabel(2, "kMB");
    hEventSelecInfo->GetXaxis()->SetBinLabel(3, "kINT7");
    hEventSelecInfo->GetXaxis()->SetBinLabel(4, "kHighMultV0");
    hEventSelecInfo->GetXaxis()->SetBinLabel(8, "kAny");
    hEventSelecInfo->GetXaxis()->SetBinLabel(9, "kAndMB");
    fOutputList->Add(hEventSelecInfo);

    hNumberOfEvent = new TH1F("hNumberOfEvent", "hNumberOfEvent", 10, 0, 10);
    fOutputList->Add(hNumberOfEvent);
    hNumberOfEvent->GetXaxis()->SetBinLabel(1, "Original");
    hNumberOfEvent->GetXaxis()->SetBinLabel(2, "AfterMulti");
    hNumberOfEvent->GetXaxis()->SetBinLabel(3, "0-100Cut");
    // fPtEtaDist = new TH3F("fPtEtaDist","PtEtaDist",2,-1.1,1.1, 300,0,3., 28,-1.4,1.4);
    // fOutputList->Add(fPtEtaDist);

    // fPhiPtDist = new TH3F("fPhiPtDist","PhiPtDist",2,-1.1,1.1, 120,0,2*PI, 300,0,3.);
    // fOutputList->Add(fPhiPtDist);


    fPtDist = new TH1F("fPtDist", "fPtDist", 90, 0, 9);
    fPtDist->GetXaxis()->SetTitle("fPtDist");
    fOutputList->Add(fPtDist);


    fPhiDist = new TH1F("fPhiDist", "fPhiDist", 70, 0, 7);
    fPhiDist->GetXaxis()->SetTitle("fPhiDist");
    fOutputList->Add(fPhiDist);


    fEtaDist = new TH1F("fEtaDist", "fEtaDist", 60, -3, 3);
    fEtaDist->GetXaxis()->SetTitle("fEtaDist");
    fOutputList->Add(fEtaDist);

    fXiStarYDist = new TH1F("fXiStarYDist", "fXiStarYDist", 400, -2, 2);
    fXiStarYDist->GetXaxis()->SetTitle("fXiStarYDist");
    fOutputList->Add(fXiStarYDist);

    fQAXiStarYDist = new TH1F("fQAXiStarYDist", "fQAXiStarYDist", 400, -2, 2);
    fQAXiStarYDist->GetXaxis()->SetTitle("fQAXiStarYDist");
    fOutputList->Add(fQAXiStarYDist);


    // MC rapidity

    fXiStarYDistMC = new TH1F("fXiStarYDistMC", "fXiStarYDist in MC", 400, -2, 2);
    fXiStarYDistMC->GetXaxis()->SetTitle("fXiStarYDistMC");
    fOutputList->Add(fXiStarYDistMC);


    fXiYDistMC1 = new TH1F("fXiYDistMC1", "fXiYDistMC1 in MC", 400, -2, 2);
    fXiYDistMC1->GetXaxis()->SetTitle("fXiYDistMC1");
    fOutputList->Add(fXiYDistMC1);

    fXiStarYDistMC1 = new TH1F("fXiStarYDistMC1", "fXiStarYDist in MC", 400, -2, 2);
    fXiStarYDistMC1->GetXaxis()->SetTitle("fXiStarYDistMC1");
    fOutputList->Add(fXiStarYDistMC1);


    fXiYDistMCout = new TH1F("fXiYDistMCout", "fXiYDistMC in MC output", 400, -2, 2);
    fXiYDistMCout->GetXaxis()->SetTitle("fXiYDistMCout");
    fOutputList->Add(fXiYDistMCout);


    fXiStarYDistMCout = new TH1F("fXiStarYDistMCout", "fXiYDistMC in MC output", 400, -2, 2);
    fXiStarYDistMCout->GetXaxis()->SetTitle("fXiStarYDistMCout");
    fOutputList->Add(fXiStarYDistMCout);




    fCutEvents = new TH1F("fCutEvents", "fCutEvents", 16, 0, 16);
    fOutputList->Add(fCutEvents);

    fTPCNcls_p = new TH1F("fTPCNcls_p", "TPC Number of cluster proton", 200, 0, 200); //par 0
    fOutputList->Add(fTPCNcls_p),
                fTPCNcls_pi1 = new TH1F("fTPCNcls_pi1", "TPC Number of cluster 1st pion", 200, 0, 200); //par 1
    fOutputList->Add(fTPCNcls_pi1),
                fTPCNcls_pi2 = new TH1F("fTPCNcls_pi2", "TPC Number of cluster 2nd pion", 200, 0, 200); //par 2
    fOutputList->Add(fTPCNcls_pi2),
                fTPCNcls_pi3 = new TH1F("fTPCNcls_pi3", "TPC Number of cluster 3rd pion", 200, 0, 200); //par 3
    fOutputList->Add(fTPCNcls_pi3);


    fQATPCNcls_p = new TH1F("fQATPCNcls_p", " After cut : TPC Number of cluster proton", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_p),
                fQATPCNcls_pi1 = new TH1F("fQATPCNcls_pi1", "  After cut : TPC Number of cluster 1st pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi1),
                fQATPCNcls_pi2 = new TH1F("fQATPCNcls_pi2", "  After cut : TPC Number of cluster 2nd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi2),
                fQATPCNcls_pi3 = new TH1F("fQATPCNcls_pi3", "  After cut : TPC Number of cluster 3rd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi3);


    fQATPCNcls_p_L = new TH1F("fQATPCNcls_p_L", " After loose cut : TPC Number of cluster proton", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_p_L),
                fQATPCNcls_pi1_L = new TH1F("fQATPCNcls_pi1_L", "  After loose cut : TPC Number of cluster 1st pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi1_L),
                fQATPCNcls_pi2_L = new TH1F("fQATPCNcls_pi2_L", "  After loose cut : TPC Number of cluster 2nd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi2_L),
                fQATPCNcls_pi3_L = new TH1F("fQATPCNcls_pi3_L", "  After loose cut : TPC Number of cluster 3rd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi3_L);


    fQATPCNcls_p_T = new TH1F("fQATPCNcls_p_T", " After tight cut : TPC Number of cluster proton", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_p_T),
                fQATPCNcls_pi1_T = new TH1F("fQATPCNcls_pi1_T", "  After tight cut : TPC Number of cluster 1st pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi1_T),
                fQATPCNcls_pi2_T = new TH1F("fQATPCNcls_pi2_T", "  After tight cut : TPC Number of cluster 2nd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi2_T),
                fQATPCNcls_pi3_T = new TH1F("fQATPCNcls_pi3_T", "  After tight cut : TPC Number of cluster 3rd pion", 200, 0, 200);
    fOutputList->Add(fQATPCNcls_pi3_T);





    fDCADist_p = new TH1F("fDCADist_p", "DCA distribution proton", 200, 0, 0.5); //par 4
    fOutputList->Add(fDCADist_p),
                fDCADist_pi1 = new TH1F("fDCADist_pi1", "DCA distribution 1st pion", 200, 0, 0.5); //par 5
    fOutputList->Add(fDCADist_pi1),
                fDCADist_pi2 = new TH1F("fDCADist_pi2", "DCA distribution 2nd pion", 200, 0, 0.5); //par 6
    fOutputList->Add(fDCADist_pi2);

    fQADCADist_p = new TH1F("fQADCADist_p", "After cut :DCA distribution proton", 200, 0, 0.5); //par 4
    fOutputList->Add(fQADCADist_p),
                fQADCADist_pi1 = new TH1F("fQADCADist_pi1", "After cut :DCA distribution 1st pion", 200, 0, 0.5); //par 5
    fOutputList->Add(fQADCADist_pi1),
                fQADCADist_pi2 = new TH1F("fQADCADist_pi2", "After cut :DCA distribution 2nd pion", 200, 0, 0.5); //par 6
    fOutputList->Add(fQADCADist_pi2);

    fQADCADist_p_L = new TH1F("fQADCADist_p_L", "After loose cut :DCA distribution proton", 200, 0, 0.5); //par 4
    fOutputList->Add(fQADCADist_p_L),
                fQADCADist_pi1_L = new TH1F("fQADCADist_pi1_L", "After loose cut :DCA distribution 1st pion", 200, 0, 0.5); //par 5
    fOutputList->Add(fQADCADist_pi1_L),
                fQADCADist_pi2_L = new TH1F("fQADCADist_pi2_L", "After loose cut :DCA distribution 2nd pion", 200, 0, 0.5); //par 6
    fOutputList->Add(fQADCADist_pi2_L);

    fQADCADist_p_T = new TH1F("fQADCADist_p_T", "After tight cut :DCA distribution proton", 200, 0, 0.5); //par 4
    fOutputList->Add(fQADCADist_p_T),
                fQADCADist_pi1_T = new TH1F("fQADCADist_pi1_T", "After tight cut :DCA distribution 1st pion", 200, 0, 0.5); //par 5
    fOutputList->Add(fQADCADist_pi1_T),
                fQADCADist_pi2_T = new TH1F("fQADCADist_pi2_T", "After tight cut :DCA distribution 2nd pion", 200, 0, 0.5); //par 6
    fOutputList->Add(fQADCADist_pi2_T);



    fDCADist_lambda = new TH1F("fDCADist_lambda", "DCA distribution Lambda", 200, 0, 0.5);
    fOutputList->Add(fDCADist_lambda),
                fDCADist_3rd_pi = new TH1F("fDCADist_3rd_pi", "DCA distribution 3rd pion", 300, 0, 3);
    fOutputList->Add(fDCADist_3rd_pi),
                fDCADist_pi_p = new TH1F("fDCADist_pi_p", "DCA distribution Pion-Proton", 300, 0, 3);
    fOutputList->Add(fDCADist_pi_p),
                fDCADist_pi_lambda = new TH1F("fDCADist_pi_lambda", "DCA distribution Pion-Lambda", 300, 0, 3);
    fOutputList->Add(fDCADist_pi_lambda);


    fQADCADist_lambda = new TH1F("fQADCADist_lambda", "After cut :DCA distribution Lambda", 200, 0, 0.5); //7
    fOutputList->Add(fQADCADist_lambda),
                fQADCADist_3rd_pi = new TH1F("fQADCADist_3rd_pi", "After cut :DCA distribution 3rd pion", 300, 0, 3); //8
    fOutputList->Add(fQADCADist_3rd_pi),
                fQADCADist_pi_p = new TH1F("fQADCADist_pi_p", "After cut :DCA distribution Pion-Proton", 300, 0, 3); //9
    fOutputList->Add(fQADCADist_pi_p),
                fQADCADist_pi_lambda = new TH1F("fQADCADist_pi_lambda", "After cut :DCA distribution Pion-Lambda", 300, 0, 3); //10
    fOutputList->Add(fQADCADist_pi_lambda);


    fQADCADist_lambda_L = new TH1F("fQADCADist_lambda_L", "After loose cut :DCA distribution Lambda", 200, 0, 0.5); //7
    fOutputList->Add(fQADCADist_lambda_L),
                fQADCADist_3rd_pi_L = new TH1F("fQADCADist_3rd_pi_L", "After loose cut :DCA distribution 3rd pion", 300, 0, 3); //8
    fOutputList->Add(fQADCADist_3rd_pi_L),
                fQADCADist_pi_p_L = new TH1F("fQADCADist_pi_p_L", "After loose cut :DCA distribution Pion-Proton", 300, 0, 3); //9
    fOutputList->Add(fQADCADist_pi_p_L),
                fQADCADist_pi_lambda_L = new TH1F("fQADCADist_pi_lambda_L", "After loose cut :DCA distribution Pion-Lambda", 300, 0, 3); //10
    fOutputList->Add(fQADCADist_pi_lambda_L);

    fQADCADist_lambda_T = new TH1F("fQADCADist_lambda_T", "After tight cut :DCA distribution Lambda", 200, 0, 0.5); //7
    fOutputList->Add(fQADCADist_lambda_T),
                fQADCADist_3rd_pi_T = new TH1F("fQADCADist_3rd_pi_T", "After tight cut :DCA distribution 3rd pion", 300, 0, 3); //8
    fOutputList->Add(fQADCADist_3rd_pi_T),
                fQADCADist_pi_p_T = new TH1F("fQADCADist_pi_p_T", "After tight cut :DCA distribution Pion-Proton", 300, 0, 3); //9
    fOutputList->Add(fQADCADist_pi_p_T),
                fQADCADist_pi_lambda_T = new TH1F("fQADCADist_pi_lambda_T", "After tight cut :DCA distribution Pion-Lambda", 300, 0, 3); //10
    fOutputList->Add(fQADCADist_pi_lambda_T);

    fCosPA_lambda = new TH1F("fCosPA_lambda", "Cosine pointing angle Lambda", 150, 0.85, 1.0);
    fOutputList->Add(fCosPA_lambda),
                fCosPA_Xi = new TH1F("fCosPA_Xi", "Cosine pointing angle Xi", 125, 0.975, 1.0);
    fOutputList->Add(fCosPA_Xi);


    fQACosPA_lambda = new TH1F("fQACosPA_lambda", "After cut :Cosine pointing angle Lambda", 150, 0.85, 1.0);
    fOutputList->Add(fQACosPA_lambda),
                fQACosPA_Xi = new TH1F("fQACosPA_Xi", "After cut :Cosine pointing angle Xi", 125, 0.975, 1.0);
    fOutputList->Add(fQACosPA_Xi);

    fQACosPA_lambda_L = new TH1F("fQACosPA_lambda_L", "After loose cut :Cosine pointing angle Lambda", 150, 0.85, 1.0);
    fOutputList->Add(fQACosPA_lambda_L),
                fQACosPA_Xi_L = new TH1F("fQACosPA_Xi_L", "After loose cut :Cosine pointing angle Xi", 125, 0.975, 1.0);
    fOutputList->Add(fQACosPA_Xi_L);


    fQACosPA_lambda_T = new TH1F("fQACosPA_lambda_T", "After tight cut :Cosine pointing angle Lambda", 150, 0.85, 1.0);
    fOutputList->Add(fQACosPA_lambda_T),
                fQACosPA_Xi_T = new TH1F("fQACosPA_Xi_T", "After tight cut :Cosine pointing angle Xi", 125, 0.975, 1.0);
    fOutputList->Add(fQACosPA_Xi_T);






    hXiInvMass = new TH1F("hXiInvMass", "Xi invariant mass distribution : cent 0 - 10", 200, 1.2, 1.4);
    fOutputList->Add(hXiInvMass);

    hQAXiInvMass = new TH1F("hQAXiInvMass", "Xi invariant mass distribution after mass window selection : cent 0 - 10", 200, 1.2, 1.4);
    fOutputList->Add(hQAXiInvMass);




    hTPCPID = new TH2F("TPCPID", "PID via TPC", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hTPCPID),
                hTPCPIDpi = new TH2F("hTPCPIDpi", "PID pion", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hTPCPIDpi),
                hTPCPIDk = new TH2F("hTPCPIDk", "PID kaon", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hTPCPIDk),
                hTPCPIDp = new TH2F("hTPCPIDp", "PID proton", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hTPCPIDp);




    hdEdxProton = new TH2F("hdEdxProton", "Xi PID p", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxProton),
                hdEdxPion1 = new TH2F("hdEdxPion1", "Xi PID pi", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxPion1),
                hdEdxPion2 = new TH2F("hdEdxPion2", "Xi PID pi_b", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxPion2);

    hdEdxProtonAfter = new TH2F("hdEdxProtonAfter", "Xi PID p", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxProtonAfter),
                hdEdxPion1After = new TH2F("hdEdxPion1After", "Xi PID pi", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxPion1After),
                hdEdxPion2After = new TH2F("hdEdxPion2After", "Xi PID pi_b", 1000, 0, 20, 1000, 0, 200);
    fOutputList->Add(hdEdxPion2After);




    hNSig3rdPion = new TH2F("hNSig3rdPion", "nSigma 3rd pion", 1000, 0, 20, 1000, -5, 5);
    fOutputList->Add(hNSig3rdPion);
    hQANSig3rdPion = new TH2F("hQANSig3rdPion", "nSigma 3rd pion : QA", 1000, 0, 20, 1000, -5, 5);
    fOutputList->Add(hQANSig3rdPion);

    hTPCNSigProton = new TH1F("hTPCNSigProton", "nSigma proton", 100, -5, 5);
    fOutputList->Add(hTPCNSigProton);
    hTPCNSigPion1 = new TH1F("hTPCNSigPion1", "nSigma 1st pion", 100, -5, 5);
    fOutputList->Add(hTPCNSigPion1);
    hTPCNSigPion2 = new TH1F("hTPCNSigPion2", "nSigma 2nd pion", 100, -5, 5);
    fOutputList->Add(hTPCNSigPion2);

    hQATPCNSigProton = new TH1F("hQATPCNSigProton", "nSigma proton : QA", 100, -5, 5);
    fOutputList->Add(hQATPCNSigProton);
    hQATPCNSigPion1 = new TH1F("hQATPCNSigPion1", "nSigma 1st pion : QA", 100, -5, 5);
    fOutputList->Add(hQATPCNSigPion1);
    hQATPCNSigPion2 = new TH1F("hQATPCNSigPion2", "nSigma 2nd pion : QA", 100, -5, 5);
    fOutputList->Add(hQATPCNSigPion2);


    for (Int_t cv = 0; cv < kNCutVariations; cv++) {
        if (!fSetSystematic && cv > 0) continue;
        if (cv == 0) {
            TString *nameXi = new TString("fXi_");
            TString *nameXibar = new TString("fXibar_");
            *nameXi += cv;
            *nameXibar += cv;
            CutVar[cv].fXi = new TH3F(nameXi->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
            fOutputList->Add(CutVar[cv].fXi);
            CutVar[cv].fXibar = new TH3F(nameXibar->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
            fOutputList->Add(CutVar[cv].fXibar);
            //
            TString *nameMCrecXi = new TString("fMCrecXi_");
            TString *nameMCrecXibar = new TString("fMCrecXibar_");
            *nameMCrecXi += cv;
            *nameMCrecXibar += cv;
            CutVar[cv].fMCrecXi = new TH3F(nameMCrecXi->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
            CutVar[cv].fMCrecXibar = new TH3F(nameMCrecXibar->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
            fOutputList->Add(CutVar[cv].fMCrecXi);
            fOutputList->Add(CutVar[cv].fMCrecXibar);
        }
        //
        TString *nameXiMinusPiPlus = new TString("fXiMinusPiPlus_");
        TString *nameXiMinusPiMinus = new TString("fXiMinusPiMinus_");
        TString *nameXiPlusPiPlus = new TString("fXiPlusPiPlus_");
        TString *nameXiPlusPiMinus = new TString("fXiPlusPiMinus_");
        TString *nameXiMinusPiPlusbkg = new TString("fXiMinusPiPlusbkg_");
        TString *nameXiMinusPiMinusbkg = new TString("fXiMinusPiMinusbkg_");
        TString *nameXiPlusPiPlusbkg = new TString("fXiPlusPiPlusbkg_");
        TString *nameXiPlusPiMinusbkg = new TString("fXiPlusPiMinusbkg_");
        *nameXiMinusPiPlus += cv;
        *nameXiMinusPiMinus += cv;
        *nameXiPlusPiPlus += cv;
        *nameXiPlusPiMinus += cv;
        *nameXiMinusPiPlusbkg += cv;
        *nameXiMinusPiMinusbkg += cv;
        *nameXiPlusPiPlusbkg += cv;
        *nameXiPlusPiMinusbkg += cv;
        // Change the Y Info to multiplicity -2,2, 200 -> 0,100,100
        CutVar[cv].fXiMinusPiPlus  = new TH3F(nameXiMinusPiPlus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiMinusPiMinus = new TH3F(nameXiMinusPiMinus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiPlusPiPlus   = new TH3F(nameXiPlusPiPlus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiPlusPiMinus  = new TH3F(nameXiPlusPiMinus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiMinusPiPlusbkg  = new TH3F(nameXiMinusPiPlusbkg->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiMinusPiMinusbkg = new TH3F(nameXiMinusPiMinusbkg->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiPlusPiPlusbkg   = new TH3F(nameXiPlusPiPlusbkg->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fXiPlusPiMinusbkg  = new TH3F(nameXiPlusPiMinusbkg->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);

        fOutputList->Add(CutVar[cv].fXiMinusPiPlus);
        fOutputList->Add(CutVar[cv].fXiMinusPiMinus);
        fOutputList->Add(CutVar[cv].fXiPlusPiPlus);
        fOutputList->Add(CutVar[cv].fXiPlusPiMinus);
        fOutputList->Add(CutVar[cv].fXiMinusPiPlusbkg);
        fOutputList->Add(CutVar[cv].fXiMinusPiMinusbkg);
        fOutputList->Add(CutVar[cv].fXiPlusPiPlusbkg);
        fOutputList->Add(CutVar[cv].fXiPlusPiMinusbkg);
        //


        TString *nameMCrecXiMinusPiPlus = new TString("fMCrecXiMinusPiPlus_");
        TString *nameMCrecXiPlusPiMinus = new TString("fMCrecXiPlusPiMinus_");
        *nameMCrecXiMinusPiPlus += cv;
        *nameMCrecXiPlusPiMinus += cv;
        CutVar[cv].fMCrecXiMinusPiPlus  = new TH3F(nameMCrecXiMinusPiPlus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        CutVar[cv].fMCrecXiPlusPiMinus  = new TH3F(nameMCrecXiPlusPiMinus->Data(), "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        fOutputList->Add(CutVar[cv].fMCrecXiMinusPiPlus);
        fOutputList->Add(CutVar[cv].fMCrecXiPlusPiMinus);
        //

    }




    //////////////////////
    // MC input histos
    //
    if(fMCcase){    
        hMCinputTotalXiStar1 = new TH3F("hMCinputTotalXiStar1", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        hMCinputTotalXiStarbar1 = new TH3F("hMCinputTotalXiStarbar1", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        fOutputList->Add(hMCinputTotalXiStar1);
        fOutputList->Add(hMCinputTotalXiStarbar1);

        hMCinputTotalXi1 = new TH3F("hMCinputTotalXi1", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
        hMCinputTotalXibar1 = new TH3F("hMCinputTotalXibar1", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
        fOutputList->Add(hMCinputTotalXi1);
        fOutputList->Add(hMCinputTotalXibar1);

        //

        hMCinputTotalXiStar3 = new TH3F("hMCinputTotalXiStar3", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        hMCinputTotalXiStarbar3 = new TH3F("hMCinputTotalXiStarbar3", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.4, 1.6);
        fOutputList->Add(hMCinputTotalXiStar3);
        fOutputList->Add(hMCinputTotalXiStarbar3);

        hMCinputTotalXi3 = new TH3F("hMCinputTotalXi3", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
        hMCinputTotalXibar3 = new TH3F("hMCinputTotalXibar3", "Invariant Mass Distribution", 100, 0, 10, 100, 0, 100, 200, 1.2, 1.4);
        fOutputList->Add(hMCinputTotalXi3);
        fOutputList->Add(hMCinputTotalXibar3);
    }
    //


    /// V0A and V0C information

    /*   hV0AC = new TH2F("hV0AC","V0A and V0C timing",500,-20,30,800,-30,50);
       hV0AC->GetXaxis()->SetTitle("V0A-V0C(ns)");
       hV0AC->GetYaxis()->SetTitle("V0A+V0C(ns)");
       fOutputList->Add(hV0AC);
       hV0ACkMB = new TH2F("hV0ACkMB","V0A and V0C timing kMB select",500,-20,30,800,-30,50);
       hV0ACkMB->GetXaxis()->SetTitle("V0A-V0C(ns)");
       hV0ACkMB->GetYaxis()->SetTitle("V0A+V0C(ns)");
       fOutputList->Add(hV0ACkMB);
       hV0Info = new TH1F("hV0Info","V0 Fired information",5,0,5);
       fOutputList->Add(hV0Info);
    */

    ///////////////////////////////////
    PostData(1, fOutputList);


}

//________________________________________________________________________
void AliXiStarpp13TeVDevel::Exec(Option_t *)
{
    // Main loop
    // Called for each event
    if (fDevelopeMode)std::cout << "===========  Event # " << fEventCounter + 1 << "  ===========" << std::endl;
    fEventCounter++;


    fESD = dynamic_cast<AliESDEvent*> (InputEvent());
    if (!fESD) {
        Printf("ERROR: fESD not available");
        return;
    }

    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if(!man){
        Printf("ERROR: AliAnalysisManager not available");
        return;
    }
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if(!inputHandler){
        Printf("ERROR: AliInputEventHandler not available");
        return;
    }

    //
    AliMCEvent  *mcEvent        = 0x0;
    AliStack    *mcstack        = 0x0;
    if (fDevelopeMode)std::cout << "TEST" << std::endl;
    //for mc study
    Bool_t IsINELg0 = kFALSE;
    Bool_t isSelectedkINT7 = kFALSE;

    if (fMCcase) {
        AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
        if(eventHandler){
            AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>(eventHandler);
            if(mcEventHandler) mcEvent = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent();
        }
        if (!mcEvent) {
            if (fDevelopeMode)std::cout << "Return: No mcEvent" << std::endl;
            return;
        }
        mcstack = mcEvent->Stack();
        if (!mcstack) {
            if (fDevelopeMode)std::cout << "Return: No mcStack" << std::endl;
            return;
        }

        htotalEvent->Fill(0); // Total N of event

        // for INELg0 check.
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            AliMCParticle *mcInputTrack = (AliMCParticle*)mcEvent->GetTrack(it);
            if (mcInputTrack->IsPhysicalPrimary() && TMath::Abs(mcInputTrack->Charge()) && TMath::Abs(mcInputTrack->Eta()) < 1 ) IsINELg0 = true;
        }
        if (IsINELg0) {
            if (fDevelopeMode)std::cout << "It's INELg0 event!" << std::endl;
            htotalEvent->Fill(1); // Total N of INELg0 event

            const AliMultiplicity* mult = fESD->GetMultiplicity();
            Int_t fSpdT_origin = mult->GetNumberOfTracklets();
            if (fDevelopeMode)std::cout << "# of total tracklet : " << fSpdT_origin << std::endl;

            Int_t fNSpdT = 0;
            for (Int_t i = 0; i < fSpdT_origin; ++i) if (TMath::Abs(mult->GetEta(i)) < 0.8) fNSpdT++;
            if (fDevelopeMode)std::cout << "# of passed tracklet : " << fNSpdT << std::endl;

            // |Eta| < 1
            if (fNSpdT > 0) {
                htotalEvent->Fill(2); // Total N of triggered event.

                if (fNSpdT > 0) htriggered_INELg0_tracklet->Fill(0); // INEL>0

                if (fNSpdT > 51) htriggered_INELg0_tracklet->Fill(10); // INEL>0
                else if (fNSpdT > 41) htriggered_INELg0_tracklet->Fill(9); // INEL>0
                else if (fNSpdT > 36) htriggered_INELg0_tracklet->Fill(8); // INEL>0
                else if (fNSpdT > 31) htriggered_INELg0_tracklet->Fill(7); // INEL>0
                else if (fNSpdT > 26) htriggered_INELg0_tracklet->Fill(6); // INEL>0
                else if (fNSpdT > 21) htriggered_INELg0_tracklet->Fill(5); // INEL>0
                else if (fNSpdT > 16) htriggered_INELg0_tracklet->Fill(4); // INEL>0
                else if (fNSpdT > 11) htriggered_INELg0_tracklet->Fill(3); // INEL>0
                else if (fNSpdT > 6) htriggered_INELg0_tracklet->Fill(2); // INEL>0
                else if (fNSpdT > 0) htriggered_INELg0_tracklet->Fill(1); // INEL>0

                // CINT7 Triggered event.
                isSelectedkINT7 = inputHandler->IsEventSelected() & AliVEvent::kINT7;
                if (isSelectedkINT7) {

                    htotalEvent->Fill(3); // Total N of triggered event.

                    if (fNSpdT > 0) htriggered_CINT7_tracklet->Fill(0); // INEL>0 in CINT7

                    if (fNSpdT > 51) htriggered_CINT7_tracklet->Fill(10); // INEL>0 in CINT7
                    else if (fNSpdT > 41) htriggered_CINT7_tracklet->Fill(9); // INEL>0 in CINT7
                    else if (fNSpdT > 36) htriggered_CINT7_tracklet->Fill(8); // INEL>0 in CINT7
                    else if (fNSpdT > 31) htriggered_CINT7_tracklet->Fill(7); // INEL>0 in CINT7
                    else if (fNSpdT > 26) htriggered_CINT7_tracklet->Fill(6); // INEL>0 in CINT7
                    else if (fNSpdT > 21) htriggered_CINT7_tracklet->Fill(5); // INEL>0 in CINT7
                    else if (fNSpdT > 16) htriggered_CINT7_tracklet->Fill(4); // INEL>0 in CINT7
                    else if (fNSpdT > 11) htriggered_CINT7_tracklet->Fill(3); // INEL>0 in CINT7
                    else if (fNSpdT > 6) htriggered_CINT7_tracklet->Fill(2); // INEL>0 in CINT7
                    else if (fNSpdT > 0) htriggered_CINT7_tracklet->Fill(1); // INEL>0 in CINT7


                    // AliMultSelection
                    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
                    Float_t lPerc = 300; //nonsense

                    //Quality check                                                                          // it should be same with 1.
                    lPerc = MultSelection->GetMultiplicityPercentile("V0M");
                    fMultDist_MCpp->Fill(lPerc);

                    if (lPerc < 101) htriggered_CINT7_VOM->Fill(0); // INEL>0 in CINT7

                    if (lPerc > 70) htriggered_CINT7_VOM->Fill(1); // INEL>0 in CINT7
                    else if (lPerc > 50) htriggered_CINT7_VOM->Fill(2); // INEL>0 in CINT7
                    else if (lPerc > 40) htriggered_CINT7_VOM->Fill(3); // INEL>0 in CINT7
                    else if (lPerc > 30) htriggered_CINT7_VOM->Fill(4); // INEL>0 in CINT7
                    else if (lPerc > 20) htriggered_CINT7_VOM->Fill(5); // INEL>0 in CINT7
                    else if (lPerc > 15) htriggered_CINT7_VOM->Fill(6); // INEL>0 in CINT7
                    else if (lPerc > 10) htriggered_CINT7_VOM->Fill(7); // INEL>0 in CINT7
                    else if (lPerc > 5) htriggered_CINT7_VOM->Fill(8); // INEL>0 in CINT7
                    else if (lPerc > 1) htriggered_CINT7_VOM->Fill(9); // INEL>0 in CINT7
                    else if (lPerc > 0) htriggered_CINT7_VOM->Fill(10); // INEL>0 in CINT7

                    if (MultSelection->IsEventSelected()) {
                        htotalEvent->Fill(4); // Total N of Multi selected event
                        fMultDist_MCpp_selected->Fill(lPerc);

                        if (lPerc < 101) htriggered_AliMult_VOM->Fill(0); // INEL>0 in CINT7

                        if (lPerc > 70) htriggered_AliMult_VOM->Fill(1); // INEL>0 in CINT7
                        else if (lPerc > 50) htriggered_AliMult_VOM->Fill(2); // INEL>0 in CINT7
                        else if (lPerc > 40) htriggered_AliMult_VOM->Fill(3); // INEL>0 in CINT7
                        else if (lPerc > 30) htriggered_AliMult_VOM->Fill(4); // INEL>0 in CINT7
                        else if (lPerc > 20) htriggered_AliMult_VOM->Fill(5); // INEL>0 in CINT7
                        else if (lPerc > 15) htriggered_AliMult_VOM->Fill(6); // INEL>0 in CINT7
                        else if (lPerc > 10) htriggered_AliMult_VOM->Fill(7); // INEL>0 in CINT7
                        else if (lPerc > 5) htriggered_AliMult_VOM->Fill(8); // INEL>0 in CINT7
                        else if (lPerc > 1) htriggered_AliMult_VOM->Fill(9); // INEL>0 in CINT7
                        else if (lPerc > 0) htriggered_AliMult_VOM->Fill(10); // INEL>0 in CINT7

                    }// IsEventSelected in AliMultSelection
                }//eta < 1
            }//CINT7
        }//INEL>0
    }
    
    // check  : events are selected by physics selection class
    UInt_t fSelectMask = fInputHandler->IsEventSelected();
    
    isSelectedkINT7 = fSelectMask & AliVEvent::kINT7;
    if(isSelectedkINT7) ((TH1F*)hEventSelecInfo)->Fill(2);
    
    Bool_t isSelectedkHighMultV0 = fSelectMask & AliVEvent::kHighMultV0;
    if(isSelectedkHighMultV0) ((TH1F*)hEventSelecInfo)->Fill(3);
    

    if (fHMTrigger) {
        if (!isSelectedkHighMultV0) {
            if (fDevelopeMode)std::cout << "Event Rejected: No kHighMultV0 trigger" << std::endl;
            return;
        }
    }
    else {
        if (!isSelectedkINT7) {
            if (fDevelopeMode)std::cout << "Event Rejected: No kINT7 trigger" << std::endl;
            return;
        }
    }

    if (fDevelopeMode)std::cout << "TEST3" << std::endl;
    ///////////////////////////////////////////////////////////
    const AliESDVertex *PrimaryVertexESD;

    // ---- AliPIDResponse ---- //
    fPIDResponse = inputHandler->GetPIDResponse();
    double nSigTPCPID = 3.0;


    // TClonesArray *mcArray       = 0x0;
    TParticle   *MCLamD1esd     = 0x0;
    TParticle   *MCLamD2esd     = 0x0;
    TParticle   *MCLamesd       = 0x0;
    TParticle   *MCXiD2esd      = 0x0;
    TParticle   *MCXiesd        = 0x0;
    TParticle   *MCXiStarD2esd  = 0x0;
    TParticle   *MCXiStaresd    = 0x0;

    Double_t px1, py1, pz1, px2, py2, pz2;
    Double_t p1sq, p2sq, e1, e2, angle;
    Double_t dca3d;
    Float_t dca2[2];
    Double_t xiVtx[3];//, xiStarVtx[3];
    Double_t xiP[3], xiStarP[3];
    Double_t xiStarMom;
    Double_t xiMass, xiStarMass;
    Double_t xiPt, xiStarPt;
    Double_t xiY, xiStarY, MCxiY, MCxiStarY, MCxiYout, MCxiStarYout ;
    Double_t xiCharge;
    Double_t decayLengthXY;
    Double_t pDaughter1[3];
    Double_t pDaughter2[3];
    Double_t xDaughter1[3];
    Double_t xDaughter2[3];
    //
    Double_t bField = 0;
    UInt_t status = 0;
    Int_t positiveTracks = 0, negativeTracks = 0;
    Int_t myTracks = 0;
    Int_t myxiTracks = 0;
    Int_t myMCTracks = 0;
    //
    Double_t primaryVtx[3] = {0};
    Int_t mBin = 0;
    Int_t zBin = 0;
    Double_t zStep = 2 * 10 / Double_t(fZvertexBins), zStart = -10.;

    //
    Bool_t mcXiFilled = kFALSE; // So that mctracks are never used uninitialized

    fMultDist1->Fill(fESD->GetNumberOfTracks());
    PrimaryVertexESD = fESD->GetPrimaryVertex();
    if (!PrimaryVertexESD) return;
    hNumberOfEvent->Fill(0);
    fEventNumber = ((((ULong64_t)fESD->GetPeriodNumber()) << 36 ) | (((ULong64_t)fESD->GetOrbitNumber()) << 12 ) | ((ULong64_t)fESD->GetBunchCrossNumber()));

    primaryVtx[0] = PrimaryVertexESD->GetX();
    primaryVtx[1] = PrimaryVertexESD->GetY();
    primaryVtx[2] = PrimaryVertexESD->GetZ();
    fVertexDist1->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);

    if (fMCcase) {
        /////////////////////////////////////////////////
        // Lam mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrack = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrack) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }
            //fTreeVariableMCinputTotalSongPID = -1;
            //fTreeVariableEventNumber4 = GetMCEventNumber();
            if (mcInputTrack->GetPdgCode() != +kXiCode && mcInputTrack->GetPdgCode() != -kXiCode && mcInputTrack->GetPdgCode() != +kXiStarCode && mcInputTrack->GetPdgCode() != -kXiStarCode) continue;
            myMCTracks++;
            MCxiStarY = mcInputTrack->Y();
            fXiStarYDistMC->Fill(MCxiStarY);
            // Xi
            if (mcInputTrack->GetPdgCode() == +kXiCode) {
                hMCinputTotalXi1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            if (mcInputTrack->GetPdgCode() == -kXiCode) {
                hMCinputTotalXibar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            // XiStar
            if (mcInputTrack->GetPdgCode() == +kXiStarCode) {
                hMCinputTotalXiStar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
            if (mcInputTrack->GetPdgCode() == -kXiStarCode) {
                hMCinputTotalXiStarbar1->Fill(mcInputTrack->Pt(), mcInputTrack->Y(), mcInputTrack->GetCalcMass());
            }
        }
    }


    // Vertex systematic study default : 10 , loose : 11 , tight : 9 (cm)
    if (fabs(primaryVtx[2]) > 10.) return; // Z-Vertex Cut
    hNumberOfEvent->Fill(1);
    fMultDist2->Fill(fESD->GetNumberOfTracks());

    /* // Disabled due to this is implimented in AliMultSelection::IsEventSelected()
    // INEL > 0 Check
    Int_t nINEL=AliESDtrackCuts::GetReferenceMultiplicity(fESD,AliESDtrackCuts::kTracklets,1.,0.);
    if(nINEL<1){
        if(fDevelopeMode)std::cout << "Reject: INEL < 1" << std::endl;;
        return;
    }*/


    fMultDist3->Fill(fESD->GetNumberOfTracks());
    fVertexDist3->Fill(primaryVtx[0], primaryVtx[1], primaryVtx[2]);


    // Before the AliMulti
    hNumberOfEvent->Fill(2);

    // IncompleteDAQ Check
    if (fESD->IsIncompleteDAQ()) {
        if (fDevelopeMode)std::cout << "Reject: IsIncompleteDAQ" << std::endl;;
        return;
    }

    // Muliplicity Selection
    AliMultSelection *MultSelection = (AliMultSelection*) fESD->FindListObject("MultSelection");
    Float_t lPerc = 300; //nonsense

    if (MultSelection) {
        if (!(MultSelection->IsEventSelected())) {
            AliInfo("This event is not selected: AliMultSelection");
            lPerc += MultSelection->GetMultiplicityPercentile("V0M");

        }
        lPerc = MultSelection->GetMultiplicityPercentile("V0M");
    }
    else {
        //If this happens, re-check if AliMultSelectionTask ran before your task!
        AliInfo("Didn't find MultSelection!");
    }
    if (fDevelopeMode)std::cout << "Multiplicity: " << lPerc << std::endl;

    fMultDist_pp->Fill(lPerc);

    // After the AliMulti
    hNumberOfEvent->Fill(3);

    // multiplicity
    /* // Disabled due to this is implimented in
    if(PrimaryVertexESD->GetNContributors() >= 1) fMultDist4->Fill(fESD->GetNumberOfTracks());
    if(PrimaryVertexESD->GetNContributors() < 1) return; // Enrico cut // 2018.08.15 -> looks like INEL cut...? right? (blim)
    hNumberOfEvent->Fill(5);
    */

    if (fDevelopeMode)std::cout << "There are " << fESD->GetNumberOfTracks() << " tracks in this event" << std::endl;;

    bField = fESD->GetMagneticField();

    // Track loop
    for (Int_t i = 0; i < fESD->GetNumberOfTracks(); i++) {
        AliESDtrack* esdtrack = fESD->GetTrack(i);
        if (!esdtrack) continue;
        status = esdtrack->GetStatus();

        if (!fTrackCut->AcceptTrack(esdtrack)) continue;

        Bool_t goodMomentum = esdtrack->GetPxPyPz(fTempStruct[myTracks].fP);
        if (!goodMomentum) continue;
        esdtrack->GetXYZ( fTempStruct[myTracks].fX);
        //=========checking PID =========//
        //// *** TPC *** ////
        Float_t fTPCPIDmom = esdtrack->GetTPCmomentum();
        Float_t sigTPC = esdtrack->GetTPCsignal();
        Float_t nsigpi = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
        Float_t nsigk = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kKaon));
        Float_t nsigpr = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kProton));

        hTPCPID->Fill(fTPCPIDmom, sigTPC);
        if (nsigpi < 3.) hTPCPIDpi->Fill(fTPCPIDmom, sigTPC);
        if (nsigk < 3.) hTPCPIDk->Fill(fTPCPIDmom, sigTPC);
        if (nsigpr < 3.) hTPCPIDp->Fill(fTPCPIDmom, sigTPC);

        hNSig3rdPion->Fill(fTPCPIDmom, fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));

        //=========selecting 3rd pion using PID=========//
        /* TPC OFF */   //        if(nsigpi>nSigTPCPID) continue; // last update // 20140715 // TPC

        hQANSig3rdPion->Fill(fTPCPIDmom, fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));

        esdtrack->GetCovarianceXYZPxPyPz( fTempStruct[myTracks].fCov);
        //esdtrack->GetImpactParameters(dca2, cov);
        dca2[0] = sqrt( pow(fTempStruct[myTracks].fX[0] - primaryVtx[0], 2) + pow(fTempStruct[myTracks].fX[1] - primaryVtx[1], 2));
        dca2[1] = sqrt( pow(fTempStruct[myTracks].fX[2] - primaryVtx[2], 2));
        dca3d = sqrt( pow(dca2[0], 2) + pow(dca2[1], 2));

        fDCADist->Fill(fESD->GetNumberOfTracks(), dca3d);
        // fPhiPtDist->Fill(esdtrack->Charge(), esdtrack->Phi(), esdtrack->Pt());
        // fPtEtaDist->Fill(esdtrack->Charge(), esdtrack->Pt(), esdtrack->Eta());

        fPtDist->Fill(esdtrack->Pt());
        fPhiDist->Fill(esdtrack->Phi());
        fEtaDist->Fill(esdtrack->Eta());

        fTempStruct[myTracks].fStatus = status;
        fTempStruct[myTracks].fID = esdtrack->GetID();
        fTempStruct[myTracks].fLabel = esdtrack->GetLabel();
        fTempStruct[myTracks].fPhi = atan2(fTempStruct[myTracks].fP[1], fTempStruct[myTracks].fP[0]);
        if (fTempStruct[myTracks].fPhi < 0) fTempStruct[myTracks].fPhi += 2 * PI;
        fTempStruct[myTracks].fPt = sqrt(pow(fTempStruct[myTracks].fP[0], 2) + pow(fTempStruct[myTracks].fP[1], 2));
        fTempStruct[myTracks].fMom = sqrt( pow(fTempStruct[myTracks].fPt, 2) + pow(fTempStruct[myTracks].fP[2], 2) );
        fTempStruct[myTracks].fEta = esdtrack->Eta();
        fTempStruct[myTracks].fCharge = esdtrack->Charge();
        fTempStruct[myTracks].fDCAXY = dca2[0];
        fTempStruct[myTracks].fDCAZ = dca2[1];
        fTempStruct[myTracks].fDCA = dca3d;
        fTempStruct[myTracks].fNSigmaPi = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kPion));
        fTempStruct[myTracks].fNSigmaK = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kKaon));
        fTempStruct[myTracks].fNSigmaPr = fabs(fPIDResponse->NumberOfSigmasTPC(esdtrack, AliPID::kProton));
        fTempStruct[myTracks].fNclusTPC = esdtrack->GetTPCNcls();

        if (esdtrack->Charge() > 0) positiveTracks++;
        else negativeTracks++;

        //if(fTempStruct[myTracks].fNclusTPC < 50) continue;  //60 to 50
        //   if(dca2[1]>3) continue;
        //   if(dca2[0]>3) continue;
        myTracks++;
    }
    hNumberOfEvent->Fill(6);
    if (myTracks >= 1) {
        fMultDist5->Fill(myTracks);
        fMultDist3d->Fill(positiveTracks + negativeTracks, positiveTracks, negativeTracks);
    }

    if (fDevelopeMode)std::cout << "There are " << myTracks << "  myTracks" << std::endl;

    // set Z Vertex bin
    for (Int_t i = 0; i < fZvertexBins; i++) {
        if (fDevelopeMode)std::cout << "fZvertexBins: " << fZvertexBins << std::endl;
        if (fDevelopeMode)std::cout << "loop: " << i << std::endl;
        if ( (primaryVtx[2] > zStart + i * zStep) && (primaryVtx[2] < zStart + (i + 1)*zStep) ) {
            zBin = i;
            break;
        }
    }

    // set Multiplicity bin
    for (Int_t i = 0; i < fMultBins; i++) {
        if (fDevelopeMode)std::cout << "Multi bin: " << fMultBins << std::endl;
        if (fDevelopeMode)std::cout << "loop: " << i << std::endl;
        if ( ( myTracks > fMultLimits[i]) && ( myTracks <= fMultLimits[i + 1]) ) {
            mBin = i;
            break;
        }
    }
    if (fDevelopeMode)std::cout << "zBin : " << zBin << "mBin: " << mBin << std::endl;

    if (fDevelopeMode)std::cout << "01" << std::endl;

    ////////////////////////////////////
    // Add event to buffer if > 0 tracks
    if (myTracks > 0) {
        fEC[zBin][mBin]->FIFOShift();
        (fEvt) = fEC[zBin][mBin]->fEvtStr;
        (fEvt)->fNTracks = myTracks;
        for (Int_t i = 0; i < myTracks; i++) (fEvt)->fTracks[i] = fTempStruct[i];
    }



    if (fMCcase) { // get Input MC information for ESD case

        /////////////////////////////////////////////////
        // Xi mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrackXi = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrackXi) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }


            //if(!mcstack->IsPhysicalPrimary(it)) continue;
            if (abs(mcInputTrackXi->GetPdgCode()) != kXiCode) continue;

            fXiYDistMC1->Fill(mcInputTrackXi->Y());


            if (mcInputTrackXi->GetPdgCode() == +kXiCode) {
                hMCinputTotalXi3->Fill(mcInputTrackXi->Pt(),  mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());
            }
            else {
                hMCinputTotalXibar3->Fill(mcInputTrackXi->Pt(),  mcInputTrackXi->Y(), mcInputTrackXi->GetCalcMass());

            }



        }


        /////////////////////////////////////////////////
        // XiStar mc input
        /////////////////////////////////////////////////
        for (Int_t it = 0; it < mcstack->GetNprimary(); it++) {
            TParticle *mcInputTrackXiStar = (TParticle*)mcstack->Particle(it);
            if (!mcInputTrackXiStar) {
                Error("UserExec", "Could not receive track %d", it);
                continue;
            }
            if (abs(mcInputTrackXiStar->GetPdgCode()) != kXiStarCode) continue;


            fXiStarYDistMC1->Fill(mcInputTrackXiStar->Y());

            if (mcInputTrackXiStar->GetPdgCode() == +kXiStarCode) {
                hMCinputTotalXiStar3->Fill(mcInputTrackXiStar->Pt(), lPerc, mcInputTrackXiStar->GetCalcMass());
            }
            else {
                hMCinputTotalXiStarbar3->Fill(mcInputTrackXiStar->Pt(), lPerc, mcInputTrackXiStar->GetCalcMass());

            }


        }
    }
    hNumberOfEvent->Fill(7);
    ////////////////////////////////////////////////
    // Reconstruction
    if (fDevelopeMode)std::cout << "02" << std::endl;
    for (Int_t i = 0; i < fESD->GetNumberOfCascades(); i++) {

        AliESDcascade *Xicandidate = fESD->GetCascade(i);


        if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetNindex())) continue;
        if (TMath::Abs( Xicandidate->GetPindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;
        if (TMath::Abs( Xicandidate->GetNindex()) == TMath::Abs( Xicandidate->GetBindex())) continue;



        AliESDtrack *pTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetPindex()));
        AliESDtrack *nTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetNindex()));
        AliESDtrack *bTrackXi   = fESD->GetTrack(TMath::Abs( Xicandidate->GetBindex()));

        // Standard track QA cuts
        if (!fTrackCut->AcceptTrack(pTrackXi)) continue;
        if (!fTrackCut->AcceptTrack(nTrackXi)) continue;
        if (!fTrackCut->AcceptTrack(bTrackXi)) continue;



        //////////////////////
        // DecayParameters Key (number represents array index)
        // NclustersTPC: 0=proton, 1=pion first, 2=pion second, 3=pion third
        // DCAVtx: 4=proton, 5=pion first, 6=pion second, 7=lambda, 8=pion third
        // 9 = DCA proton-pion
        // 10 = DCA Lambda-pion
        // 11 = Cos PA Lambda
        // 12 = Cos PA Xi

        //myxiTracks++;

        fDecayParameters[2] = bTrackXi->GetTPCNcls();
        fTPCNcls_pi2->Fill(fDecayParameters[2]);

        Double_t fTPCNSigProton = 10;
        Double_t fTPCNSigPion1 = 10;
        Double_t fTPCNSigPion2 = 10;


        Double_t fTPCPIDMomXi[3] = { -10, -10, -10};
        Double_t fNSigTPCXi[3] = { -10, -10, -10};




        if (Xicandidate->Charge() == -1) {
            fDecayParameters[0] = pTrackXi->GetTPCNcls();
            fDecayParameters[1] = nTrackXi->GetTPCNcls();
            fDecayParameters[4] = fabs(pTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx proton
            fDecayParameters[5] = fabs(nTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion first

            fTPCNSigProton = fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton);
            fTPCNSigPion1 = fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion);
            fTPCNSigPion2 = fPIDResponse->NumberOfSigmasTPC(bTrackXi, AliPID::kPion);

            fTPCPIDMomXi[0] = pTrackXi->GetTPCmomentum();
            fNSigTPCXi[0] = pTrackXi->GetTPCsignal();

            fTPCPIDMomXi[1] = nTrackXi->GetTPCmomentum();
            fNSigTPCXi[1] = nTrackXi->GetTPCsignal();

            fTPCPIDMomXi[2] = bTrackXi->GetTPCmomentum();
            fNSigTPCXi[2] = bTrackXi->GetTPCsignal();

        } else {
            fDecayParameters[0] = nTrackXi->GetTPCNcls();
            fDecayParameters[1] = pTrackXi->GetTPCNcls();
            fDecayParameters[4] = fabs(nTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx proton
            fDecayParameters[5] = fabs(pTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion first

            fTPCNSigProton = fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton);
            fTPCNSigPion1 = fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion);
            fTPCNSigPion2 = fPIDResponse->NumberOfSigmasTPC(bTrackXi, AliPID::kPion);


            fTPCPIDMomXi[0] = nTrackXi->GetTPCmomentum();
            fNSigTPCXi[0] = nTrackXi->GetTPCsignal();

            fTPCPIDMomXi[1] = pTrackXi->GetTPCmomentum();
            fNSigTPCXi[1] = pTrackXi->GetTPCsignal();

            fTPCPIDMomXi[2] = bTrackXi->GetTPCmomentum();
            fNSigTPCXi[2] = bTrackXi->GetTPCsignal();
        }

        fTPCNcls_p->Fill(fDecayParameters[0]);
        fTPCNcls_pi1->Fill(fDecayParameters[1]);
        fDCADist_p->Fill(fDecayParameters[4]);
        fDCADist_pi1->Fill(fDecayParameters[5]);

        hTPCNSigProton->Fill(fTPCNSigProton);
        hTPCNSigPion1->Fill(fTPCNSigPion1);
        hTPCNSigPion2->Fill(fTPCNSigPion2);

        hdEdxProton->Fill(fTPCPIDMomXi[0], fNSigTPCXi[0]);
        hdEdxPion1->Fill(fTPCPIDMomXi[1], fNSigTPCXi[1]);
        hdEdxPion2->Fill(fTPCPIDMomXi[2], fNSigTPCXi[2]);


        if (fTPCNSigProton > -nSigTPCPID && fTPCNSigProton < nSigTPCPID)hdEdxProtonAfter->Fill(fTPCPIDMomXi[0], fNSigTPCXi[0]);
        if (fTPCNSigPion1 > -nSigTPCPID && fTPCNSigPion1 < nSigTPCPID) hdEdxPion1After->Fill(fTPCPIDMomXi[1], fNSigTPCXi[1]);
        if (fTPCNSigPion2 > -nSigTPCPID && fTPCNSigPion2 < nSigTPCPID)hdEdxPion2After->Fill(fTPCPIDMomXi[2], fNSigTPCXi[2]);


        // PID Cuts
        if (fPIDOption && abs(fTPCNSigProton) > nSigTPCPID) continue; // PID for proton
        if (fPIDOption && abs(fTPCNSigPion1) > nSigTPCPID) continue; // PID for 1st pion
        if (fPIDOption && abs(fTPCNSigPion2) > nSigTPCPID) continue; // PID for 2nd pion

        hQATPCNSigProton->Fill(fTPCNSigProton);
        hQATPCNSigPion1->Fill(fTPCNSigPion1);
        hQATPCNSigPion2->Fill(fTPCNSigPion2);


        fDecayParameters[6] = fabs(bTrackXi->GetD(primaryVtx[0], primaryVtx[1], bField)); // DCA Vtx pion second
        fDCADist_pi2->Fill(fDecayParameters[6]);
        fDecayParameters[7] = fabs(Xicandidate->GetD(primaryVtx[0], primaryVtx[1], primaryVtx[2])); // DCA Vtx Lambda
        fDCADist_lambda->Fill(fDecayParameters[7]);
        fDecayParameters[9] = fabs(Xicandidate->GetDcaV0Daughters());// DCA proton-pion
        fDCADist_pi_p->Fill(fDecayParameters[9]);
        fDecayParameters[10] = fabs(Xicandidate->GetDcaXiDaughters());// DCA Lambda-pion
        fDCADist_pi_lambda->Fill(fDecayParameters[10]);

        Double_t tempX[3] = {0};
        Xicandidate->GetXYZ(tempX[0], tempX[1], tempX[2]);

        //    fDecayParameters[11] = sqrt( pow(tempX[0],2) + pow(tempX[1],2));// Rxy Lambda
        //    fRXY_lambda->Fill(fDecayParameters[11]);

        fDecayParameters[11] = Xicandidate->GetV0CosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]); // Cos PA Lambda
        fCosPA_lambda->Fill(fDecayParameters[11]);

        fDecayParameters[12] = Xicandidate->GetCascadeCosineOfPointingAngle(primaryVtx[0], primaryVtx[1], primaryVtx[2]); // Cos PA Xi
        fCosPA_Xi->Fill(fDecayParameters[12]);

        decayLengthXY = sqrt( pow(xiVtx[0] - primaryVtx[0], 2) + pow(xiVtx[1] - primaryVtx[1], 2) );
        //     fDecayParameters[12] = decayLengthXY;// Rxy Xi
        //     fRXY_Xi->Fill(fDecayParameters[12]);


        xiP[0] = Xicandidate->Px();
        xiP[1] = Xicandidate->Py();
        xiP[2] = Xicandidate->Pz();
        xiVtx[0] = Xicandidate->Xv();
        xiVtx[1] = Xicandidate->Yv();
        xiVtx[2] = Xicandidate->Zv();
        xiPt = Xicandidate->Pt();
        xiY = Xicandidate->RapXi();
        xiMass = Xicandidate->M();
        xiCharge = Xicandidate->Charge();

        myxiTracks++;



        if (sqrt( pow(tempX[0], 2) + pow(tempX[1], 2) ) > fMaxDecayLength) continue;
        if (decayLengthXY > fMaxDecayLength) continue;

        Bool_t StandardXi = kTRUE;
        if (StandardXi)fCutEvents->Fill(1, 1);

        if (fDecayParameters[0] < fCutValues[0][0]) StandardXi = kFALSE; // Nclus proton
        if (StandardXi)fCutEvents->Fill(2, 1);

        if (fDecayParameters[1] < fCutValues[0][1]) StandardXi = kFALSE; // Nclus pion first
        if (StandardXi)fCutEvents->Fill(3, 1);

        if (fDecayParameters[2] < fCutValues[0][2]) StandardXi = kFALSE; // Nclus pion second
        if (StandardXi)fCutEvents->Fill(4, 1);

        //
        if (fDecayParameters[4] < fCutValues[0][4]) StandardXi = kFALSE; // DCAVtx proton
        if (StandardXi)fCutEvents->Fill(5, 1);

        if (fDecayParameters[5] < fCutValues[0][5]) StandardXi = kFALSE; // DCAVtx pion first
        if (StandardXi)fCutEvents->Fill(6, 1);

        if (fDecayParameters[6] < fCutValues[0][6]) StandardXi = kFALSE; // DCAVtx pion second
        if (StandardXi)fCutEvents->Fill(7, 1);

        if (fDecayParameters[7] < fCutValues[0][7]) StandardXi = kFALSE; // DCAVtx Lambda
        if (StandardXi)fCutEvents->Fill(8, 1);

        //
        if (fDecayParameters[9] > fCutValues[0][9]) StandardXi = kFALSE; // DCAV proton-pion
        if (StandardXi)fCutEvents->Fill(9, 1);

        if (fDecayParameters[10] > fCutValues[0][10]) StandardXi = kFALSE; // DCAV Lambda-pion
        if (StandardXi)fCutEvents->Fill(10, 1);

        //
        // if(fDecayParameters[11] < fCutValues[0][11]) StandardXi=kFALSE;// Rxy Lambda
        // if(StandardXi)fCutEvents->Fill(11,1);

        // if(fDecayParameters[12] < fCutValues[0][12]) StandardXi=kFALSE;// Rxy Xi
        // if(StandardXi)fCutEvents->Fill(12,1);

        //
        if (fDecayParameters[11] < fCutValues[0][11]) StandardXi = kFALSE; // Cos PA Lambda
        if (StandardXi)fCutEvents->Fill(11, 1);

        if (fDecayParameters[12] < fCutValues[0][12]) StandardXi = kFALSE; // Cos PA Xi
        if (StandardXi)fCutEvents->Fill(12, 1);

        if (StandardXi) {
            if (xiCharge == -1) {
                CutVar[0].fXi->Fill(xiPt, xiY, xiMass);
                hXiInvMass->Fill(xiMass);
            }
            else {
                CutVar[0].fXibar->Fill(xiPt, xiY, xiMass);
                hXiInvMass->Fill(xiMass);
            }
        }



        if (fDevelopeMode)std::cout << "001" << std::endl;
        // MC associaton
        mcXiFilled = kFALSE;
        if (fMCcase ) {

            MCXiD2esd = (TParticle*)mcstack->Particle(abs(bTrackXi->GetLabel()));

            if (abs(MCXiD2esd->GetPdgCode()) == kPionCode) {

                MCLamD1esd = (TParticle*)mcstack->Particle(abs(pTrackXi->GetLabel()));
                MCLamD2esd = (TParticle*)mcstack->Particle(abs(nTrackXi->GetLabel()));

                if (MCLamD1esd->GetMother(0) == MCLamD2esd->GetMother(0)) {
                    if (abs(MCLamD1esd->GetPdgCode()) == kProtonCode || abs(MCLamD2esd->GetPdgCode()) == kProtonCode) {
                        if (abs(MCLamD1esd->GetPdgCode()) == kPionCode || abs(MCLamD2esd->GetPdgCode()) == kPionCode) {

                            MCLamesd = (TParticle*)mcstack->Particle(abs(MCLamD1esd->GetMother(0)));
                            if (abs(MCLamesd->GetPdgCode()) == kLambdaCode) {

                                if (MCLamesd->GetMother(0) == MCXiD2esd->GetMother(0)) {
                                    MCXiesd = (TParticle*)mcstack->Particle(abs(MCLamesd->GetMother(0)));
                                    if (abs(MCXiesd->GetPdgCode()) == kXiCode) {
                                        mcXiFilled = kTRUE;

                                        if (StandardXi) {

                                            fXiYDistMCout->Fill(xiY);

                                            if (Xicandidate->Charge() == -1) {

                                                CutVar[0].fMCrecXi->Fill(xiPt, xiY, xiMass);
                                            } else {
                                                CutVar[0].fMCrecXibar->Fill(xiPt, xiY, xiMass);
                                            }



                                        }

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }// MC association


        if (fabs(xiMass - fTrueMassXi) > fMassWindow) continue;

        if (StandardXi)hQAXiInvMass->Fill(xiMass);

        fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));

        if (fDevelopeMode)std::cout << "002" << std::endl;

        //////////////////////////////////////////////////////////
        // Reconstruct Xi(1530)
        for (Int_t EN = 0; EN < fEventsToMix + 1; EN++) { // Event buffer loop

            for (Int_t l = 0; l < (fEvt + EN)->fNTracks; l++) { // Present(EN=0) and Past(EN from 1 to fEventsToMix) event track loop

                if (EN == 0) {
                    if ((fEvt + EN)->fTracks[l].fID == pTrackXi->GetID()) continue;
                    if ((fEvt + EN)->fTracks[l].fID == nTrackXi->GetID()) continue;
                    if ((fEvt + EN)->fTracks[l].fID == bTrackXi->GetID()) continue;
                }

                fXiTrack->Set(xiVtx, xiP, fCovMatrix, Short_t(xiCharge));

                if (!fESDTrack4) continue;
                fESDTrack4->Set((fEvt + EN)->fTracks[l].fX, (fEvt + EN)->fTracks[l].fP, (fEvt + EN)->fTracks[l].fCov, (fEvt + EN)->fTracks[l].fCharge);

                fDecayParameters[8] = (fEvt + EN)->fTracks[l].fDCAXY; // DCA Vtx pion third
                fDCADist_3rd_pi->Fill(fDecayParameters[8]);



                if ((fEvt + EN)->fTracks[l].fDCAZ > 2) continue;
                if ( (((fEvt + EN)->fTracks[l].fStatus)&AliESDtrack::kITSrefit) == 0) continue; // Require itsrefit
                // no Chi^2 cut applied for ESDs.  Info not available in my track structure.


                if (fabs((fEvt + EN)->fTracks[l].fEta) > 0.8) continue;

                fDecayParameters[3] = (fEvt + EN)->fTracks[l].fNclusTPC;
                fTPCNcls_pi3->Fill(fDecayParameters[3]);

                AliVertex *XiStarVtx = new AliVertex((fEvt + EN)->fTracks[l].fX, 0, 0);
                //fESDTrack4->PropagateToDCA(fXiTrack, bField);// Propagate tracks to dca, both tracks are budged
                if (!(fXiTrack->PropagateToDCA(XiStarVtx, bField, 3))) continue; // Propagate tracks to dca, version which assumes fESDTrack4 is already primary
                /////////////
                fXiTrack->GetPxPyPz(pDaughter1);
                fXiTrack->GetXYZ(xDaughter1);
                fESDTrack4->GetPxPyPz(pDaughter2);
                fESDTrack4->GetXYZ(xDaughter2);
                //////////////////////////



                //xiStarVtx[0] = (xDaughter1[0]+xDaughter2[0])/2.;
                //xiStarVtx[1] = (xDaughter1[1]+xDaughter2[1])/2.;
                //xiStarVtx[2] = (xDaughter1[2]+xDaughter2[2])/2.;
                //decayLength = sqrt(pow(xiStarVtx[0]-primaryVtx[0],2)+pow(xiStarVtx[1]-primaryVtx[1],2)+pow(xiStarVtx[2]-primaryVtx[2],2));

                px1 = pDaughter1[0];
                py1 = pDaughter1[1];
                pz1 = pDaughter1[2];
                px2 = pDaughter2[0];
                py2 = pDaughter2[1];
                pz2 = pDaughter2[2];

                p1sq = px1 * px1 + py1 * py1 + pz1 * pz1;
                p2sq = px2 * px2 + py2 * py2 + pz2 * pz2;
                if (p1sq <= 0 || p2sq <= 0) continue;

                e1 = sqrt(p1sq + fTrueMassXi * fTrueMassXi);
                e2 = sqrt(p2sq + fTrueMassPi * fTrueMassPi);
                angle = px1 * px2 + py1 * py2 + pz1 * pz2;
                xiStarMass = fTrueMassXi * fTrueMassXi + fTrueMassPi * fTrueMassPi + 2.*e1 * e2 - 2.*angle;
                if (xiStarMass < 0.) xiStarMass = 1.e-8;
                xiStarMass = sqrt(xiStarMass);


                xiStarP[0] = px1 + px2;
                xiStarP[1] = py1 + py2;
                xiStarP[2] = pz1 + pz2;
                xiStarMom = sqrt(pow(xiStarP[0], 2) + pow(xiStarP[1], 2) + pow(xiStarP[2], 2));
                if (xiStarMom == 0) continue; // So one of the following lines doesnt break
                xiStarPt = sqrt(xiStarP[0] * xiStarP[0] + xiStarP[1] * xiStarP[1]);


                xiStarY = .5 * log( ((e1 + e2) + xiStarP[2]) / ((e1 + e2) - xiStarP[2]));
                fXiStarYDist->Fill(xiStarY);
                if (xiStarY < -0.5 || xiStarY > 0.5) continue; // here selection of rapidity for pPb analysis

                fQAXiStarYDist->Fill(xiStarY);

                for (int cv = 0; cv < kNCutVariations; cv++) {
                    if (!fSetSystematic && cv > 0) continue;

                    if (fDecayParameters[0] < fCutValues[cv][0]) continue; // Nclus proton
                    if (fDecayParameters[1] < fCutValues[cv][1]) continue; // Nclus pion first
                    if (fDecayParameters[2] < fCutValues[cv][2]) continue; // Nclus pion second
                    if (fDecayParameters[3] < fCutValues[cv][3]) continue; // Nclus pion third
                    //
                    if (fDecayParameters[4] < fCutValues[cv][4]) continue; // DCAVtx proton
                    if (fDecayParameters[5] < fCutValues[cv][5]) continue; // DCAVtx pion first
                    if (fDecayParameters[6] < fCutValues[cv][6]) continue; // DCAVtx pion second
                    if (fDecayParameters[7] < fCutValues[cv][7]) continue; // DCAVtx Lambda
                    //if(fDecayParameters[8] > fCutValues[cv][8]) continue; // DCAVtx pion third
                    if (fDecayParameters[8] > (0.0105 + 0.035 / pow((fEvt + EN)->fTracks[l].fPt, 1.1))) continue; // DCAVtx pion third
                    //0.0182 + 0.035/pow((fEvt+EN)->fTracks[l].fPt,1.1 (2010 cut)

                    //
                    if (fDecayParameters[9] > fCutValues[cv][9]) continue; // DCAV proton-pion
                    if (fDecayParameters[10] > fCutValues[cv][10]) continue; // DCAV Lambda-pion
                    //
                    if (fDecayParameters[11] < fCutValues[cv][11]) continue; // Cos PA Lambda
                    if (fDecayParameters[12] < fCutValues[cv][12]) continue; // Cos PA Xi

                    if (EN == 0 && cv == 0) { // cut QA plot for default cut
                        fQATPCNcls_p->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3->Fill(fDecayParameters[3]);

                        fQADCADist_p->Fill(fDecayParameters[4]);
                        fQADCADist_pi1->Fill(fDecayParameters[5]);
                        fQADCADist_pi2->Fill(fDecayParameters[6]);

                        fQADCADist_lambda->Fill(fDecayParameters[7]);
                        fQADCADist_3rd_pi->Fill(fDecayParameters[8]);
                        fQADCADist_pi_p->Fill(fDecayParameters[9]);
                        fQADCADist_pi_lambda->Fill(fDecayParameters[10]);
                        fQACosPA_lambda->Fill(fDecayParameters[11]);
                        fQACosPA_Xi->Fill(fDecayParameters[12]);
                    }

                    if (EN == 0 && cv == 1) {
                        fQATPCNcls_p_L->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1_L->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2_L->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3_L->Fill(fDecayParameters[3]);
                    }
                    if (EN == 0) {

                        if (cv == 2)fQADCADist_p_L->Fill(fDecayParameters[4]);
                        if (cv == 3)fQADCADist_pi1_L->Fill(fDecayParameters[5]);
                        if (cv == 4)fQADCADist_pi2_L->Fill(fDecayParameters[6]);
                        if (cv == 5)fQADCADist_lambda_L->Fill(fDecayParameters[7]);
                        if (cv == 6)fQADCADist_3rd_pi_L->Fill(fDecayParameters[8]);
                        if (cv == 7)fQADCADist_pi_p_L->Fill(fDecayParameters[9]);
                        if (cv == 8)fQADCADist_pi_lambda_L->Fill(fDecayParameters[10]);
                        if (cv == 9)fQACosPA_lambda_L->Fill(fDecayParameters[11]);
                        if (cv == 10)fQACosPA_Xi_L->Fill(fDecayParameters[12]);

                    }

                    if (EN == 0 && cv == 11) {
                        fQATPCNcls_p_T->Fill(fDecayParameters[0]);
                        fQATPCNcls_pi1_T->Fill(fDecayParameters[1]);
                        fQATPCNcls_pi2_T->Fill(fDecayParameters[2]);
                        fQATPCNcls_pi3_T->Fill(fDecayParameters[3]);
                    }
                    if (EN == 0) {

                        if (cv == 12)fQADCADist_p_T->Fill(fDecayParameters[4]);
                        if (cv == 13)fQADCADist_pi1_T->Fill(fDecayParameters[5]);
                        if (cv == 14)fQADCADist_pi2_T->Fill(fDecayParameters[6]);
                        if (cv == 15)fQADCADist_lambda_T->Fill(fDecayParameters[7]);
                        if (cv == 16)fQADCADist_3rd_pi_T->Fill(fDecayParameters[8]);
                        if (cv == 17)fQADCADist_pi_p_T->Fill(fDecayParameters[9]);
                        if (cv == 18)fQADCADist_pi_lambda_T->Fill(fDecayParameters[10]);
                        if (cv == 19)fQACosPA_lambda_T->Fill(fDecayParameters[11]);
                        if (cv == 20)fQACosPA_Xi_T->Fill(fDecayParameters[12]);


                    }


                    if (EN == 0) {
                        if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiMinusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) {
                            // xiStarPt, xiStarY, xiStarMass -> xiStarPt, lPerc, xiStarMass
                            CutVar[cv].fXiMinusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                        }
                        else if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) {
                            CutVar[cv].fXiPlusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);
                        }
                        else CutVar[cv].fXiPlusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                    } else {
                        if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiMinusPiMinusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) CutVar[cv].fXiMinusPiPlusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fXiPlusPiMinusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                        else CutVar[cv].fXiPlusPiPlusbkg->Fill(xiStarPt, lPerc, xiStarMass);
                    }



                    // MC associaton ESD
                    if (fMCcase && mcXiFilled && EN == 0) { // ESD MC's
                        MCXiStarD2esd = (TParticle*)mcstack->Particle(abs((fEvt)->fTracks[l].fLabel));

                        if (abs(MCXiStarD2esd->GetPdgCode()) == kPionCode) {
                            if (MCXiesd->GetMother(0) == MCXiStarD2esd->GetMother(0)) {

                                MCXiStaresd = (TParticle*)mcstack->Particle(abs(MCXiesd->GetMother(0)));
                                if (abs(MCXiStaresd->GetPdgCode()) == kXiStarCode) {

                                    fXiStarYDistMCout->Fill(xiStarY);


                                    if (fXiTrack->Charge() == -1 &&  fESDTrack4->Charge() == +1) CutVar[cv].fMCrecXiMinusPiPlus->Fill(xiStarPt, lPerc, xiStarMass);
                                    if (fXiTrack->Charge() == +1 &&  fESDTrack4->Charge() == -1) CutVar[cv].fMCrecXiPlusPiMinus->Fill(xiStarPt, lPerc, xiStarMass);

                                }
                            }
                        }
                    }


                }// Cut Variation loop
            }// 3rd pion loop


        }// Event mixing loop







    }// Xi loop



    // Post output data.
    PostData(1, fOutputList);

}
//________________________________________________________________________
void AliXiStarpp13TeVDevel::Terminate(Option_t *)
{
    //if(fDevelopeMode)std::cout<<"Done"<<std::endl;
}
//________________________________________________________________________
/*
Double_t AliXiStarpp13TeVDevel::LinearPropagateToDCA(AliESDtrack *v, AliESDtrack *t, Double_t b) {// Adapted from AliCascadeVertexer.cxx
    //--------------------------------------------------------------------
    // This function returns the DCA between the V0 and the track
    //--------------------------------------------------------------------
    Double_t alpha=t->GetAlpha(), cs1=TMath::Cos(alpha), sn1=TMath::Sin(alpha);
    Double_t r[3]; t->GetXYZ(r);
    Double_t x1=r[0], y1=r[1], z1=r[2];
    Double_t p[3]; t->GetPxPyPz(p);
    Double_t px1=p[0], py1=p[1], pz1=p[2];
    Double_t x2[3]={0};
    Double_t p2[3]={0};
    Double_t vx2,vy2,vz2;     // position and momentum of V0
    Double_t px2,py2,pz2;
    v->GetXYZ(x2);
    v->GetPxPyPz(p2);
    vx2=x2[0], vy2=x2[1], vz2=x2[2];
    px2=p2[0], py2=p2[1], pz2=p2[2];
    // calculation dca
    Double_t dd= Det(vx2-x1,vy2-y1,vz2-z1,px1,py1,pz1,px2,py2,pz2);
    Double_t ax= Det(py1,pz1,py2,pz2);
    Double_t ay=-Det(px1,pz1,px2,pz2);
    Double_t az= Det(px1,py1,px2,py2);
    Double_t dca=TMath::Abs(dd)/TMath::Sqrt(ax*ax + ay*ay + az*az);
    //points of the DCA
    Double_t t1 = Det(vx2-x1,vy2-y1,vz2-z1,px2,py2,pz2,ax,ay,az)/
    Det(px1,py1,pz1,px2,py2,pz2,ax,ay,az);
    x1 += px1*t1; y1 += py1*t1; //z1 += pz1*t1;
    //propagate track to the points of DCA
    x1=x1*cs1 + y1*sn1;
    if (!t->PropagateTo(x1,b)) {
        Error("PropagateToDCA","Propagation failed !");
        return 1.e+33;
    }
    return dca;
}
*/

//________________________________________________________________________
Double_t AliXiStarpp13TeVDevel::Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const {
// Taken from AliCascadeVertexer
    //--------------------------------------------------------------------
    // This function calculates locally a 2x2 determinant
    //--------------------------------------------------------------------
    return a00 * a11 - a01 * a10;
}
//________________________________________________________________________
Double_t AliXiStarpp13TeVDevel::Det(Double_t a00, Double_t a01, Double_t a02,
                                    Double_t a10, Double_t a11, Double_t a12,
                                    Double_t a20, Double_t a21, Double_t a22) const {
    // Taken from AliCascadeVertexer
    //--------------------------------------------------------------------
    // This function calculates locally a 3x3 determinant
    //--------------------------------------------------------------------
    return  a00 * Det(a11, a12, a21, a22) - a01 * Det(a10, a12, a20, a22) + a02 * Det(a10, a11, a20, a21);
}