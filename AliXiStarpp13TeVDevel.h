#ifndef ALIXISTARPP13TEVDEVEL_H
#define ALIXISTARPP13TEVDEVEL_H
//
// Class AliXiStarpp13TeVDevel
//
// AliXiStarpp13TeVDevel
// author:
//  (Original Code) Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//  (1st Modification) Jihye Song (jihye.song@cern.ch)
//  (last Modification) Bong-Hwi Lim (bong-hwi.lim@cern.ch)


class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TProfile;
class TTree;

class AliESDEvent;
class AliESDtrackCuts;
class AliESDpid;

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODPid.h"
#include "AliESDpid.h"
#include "AliXiStarpp13TeVDevelEventCollection.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "AliVertex.h"


class AliXiStarpp13TeVDevel : public AliAnalysisTaskSE {
public:

    AliXiStarpp13TeVDevel();
    AliXiStarpp13TeVDevel(const char *name, Bool_t AODdecision, Bool_t MCdecision, Int_t CutListOption = 0, Bool_t DevelopmentMode = kFALSE, Bool_t fHMTrigger = kFALSE, Bool_t fPIDOption = kFALSE, Bool_t SetSystematic = kTRUE);

    virtual ~AliXiStarpp13TeVDevel();
    AliXiStarpp13TeVDevel(const AliXiStarpp13TeVDevel &obj );
    AliXiStarpp13TeVDevel &operator=(const AliXiStarpp13TeVDevel &obj );

    enum {
        kNbinsM              = 200, // mult bins for certain histograms //300
        kXiStarCode          = 3324,// Xi(1530)^0 MC code
        kXiCode              = 3312,// Xi- MC code
        kLambdaCode          = 3122,// Lambda MC code
        kProtonCode          = 2212,// Proton+ MC code
        kPionCode            = 211,// Pion+ MC code
        kNCutVariations      = 21,// number of cut variations // 13
        kNCuts               = 13// number of cut types //15
    };

    //=================================================================================//
    //generated Histograms//
    //=================================================================================//

    struct St_CutType {
        TH3F *fXi; //!
        TH3F *fXibar; //!
        //
        TH3F *fXiMinusPiPlus; //!
        TH3F *fXiMinusPiMinus; //!
        TH3F *fXiPlusPiPlus; //!
        TH3F *fXiPlusPiMinus; //!

        TH3F *fXiMinusPiPlusbkg; //!
        TH3F *fXiMinusPiMinusbkg; //!
        TH3F *fXiPlusPiPlusbkg; //!
        TH3F *fXiPlusPiMinusbkg; //!
        //
        TH3F *fMCrecXi; //!
        TH3F *fMCrecXibar; //!

        TH3F *fMCrecXiMinusPiPlus; //!
        TH3F *fMCrecXiPlusPiMinus; //!

    };
    struct St_CutType CutVar[kNCutVariations]; //!

private:

    virtual void   UserCreateOutputObjects();
    virtual void   Exec(Option_t *option);
    virtual void   Terminate(Option_t *);
    void XiStarInit();// initialization of fixed values
    Double_t LinearPropagateToDCA(AliESDtrack*, AliESDtrack*, Double_t);// for linear propagation
    Double_t Det(Double_t, Double_t, Double_t, Double_t) const;// for linear propagation
    Double_t Det(Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t) const;// for linear propagation
    ULong64_t GetMCEventNumber();

    TH1F *htotalEvent; //! for trigger efficiency study
    TH1F *htriggered_INELg0_tracklet; //! for trigger efficiency study
    TH1F *htriggered_CINT7_tracklet; //! for trigger efficiency study
    TH1F *htriggered_CINT7_VOM; //! for trigger efficiency study
    TH1F *htriggered_AliMult_VOM; //! for trigger efficiency study
    TH1F *fMultDist_MCpp; //!
    TH1F *fMultDist_MCpp_selected; //!


    TH3F *fVertexDist1; //!
    TH3F *fVertexDist3; //!
    TH2F *fDCADist; //!
    TH3F *fMultDist3d; //!
    TH1F *fMultDist1; //!
    TH1F *fMultDist2; //!
    TH1F *fMultDist3; //!
    TH1F *fMultDist4; //!
    TH1F *fMultDist5; //!
    TH1F *fMultDist_pp; //!
    TH1F *hEventSelecInfo; //!
    TH1F *hNumberOfEvent; //!
    TH1F *a; //!

    TH1F *fPtDist; //!
    TH1F *fPhiDist; //!
    TH1F *fEtaDist; //!
    TH1F *fXiStarYDist; //!
    TH1F *fQAXiStarYDist; //!

    TH1F *fXiStarYDistMC; //!
    TH1F *fXiYDistMC1; //!
    TH1F *fXiStarYDistMC1; //!
    TH1F *fXiYDistMCout; //!
    TH1F *fXiStarYDistMCout; //!
    TH1F *fCutEvents; //!

    TH1F *fTPCNcls_p; //!
    TH1F *fTPCNcls_pi1; //!
    TH1F *fTPCNcls_pi2; //!
    TH1F *fTPCNcls_pi3; //!

    TH1F *fQATPCNcls_p; //!
    TH1F *fQATPCNcls_pi1; //!
    TH1F *fQATPCNcls_pi2; //!
    TH1F *fQATPCNcls_pi3; //!

    TH1F *fQATPCNcls_p_L; //!
    TH1F *fQATPCNcls_pi1_L; //!
    TH1F *fQATPCNcls_pi2_L; //!
    TH1F *fQATPCNcls_pi3_L; //!

    TH1F *fQATPCNcls_p_T; //!
    TH1F *fQATPCNcls_pi1_T; //!
    TH1F *fQATPCNcls_pi2_T; //!
    TH1F *fQATPCNcls_pi3_T; //!

    TH1F *fDCADist_p; //!
    TH1F *fDCADist_pi1; //!
    TH1F *fDCADist_pi2; //!

    TH1F *fQADCADist_p; //!
    TH1F *fQADCADist_pi1; //!
    TH1F *fQADCADist_pi2; //!

    TH1F *fQADCADist_p_L; //!
    TH1F *fQADCADist_pi1_L; //!
    TH1F *fQADCADist_pi2_L; //!

    TH1F *fQADCADist_p_T; //!
    TH1F *fQADCADist_pi1_T; //!
    TH1F *fQADCADist_pi2_T; //!

    TH1F *fDCADist_lambda; //!
    TH1F *fDCADist_3rd_pi; //!
    TH1F *fDCADist_pi_p; //!
    TH1F *fDCADist_pi_lambda; //!

    TH1F *fQADCADist_lambda; //!
    TH1F *fQADCADist_3rd_pi; //!
    TH1F *fQADCADist_pi_p; //!
    TH1F *fQADCADist_pi_lambda; //!

    TH1F *fQADCADist_lambda_L; //!
    TH1F *fQADCADist_3rd_pi_L; //!
    TH1F *fQADCADist_pi_p_L; //!
    TH1F *fQADCADist_pi_lambda_L; //!

    TH1F *fQADCADist_lambda_T; //!
    TH1F *fQADCADist_3rd_pi_T; //!
    TH1F *fQADCADist_pi_p_T; //!
    TH1F *fQADCADist_pi_lambda_T; //!

    TH1F *fCosPA_lambda; //!
    TH1F *fCosPA_Xi; //!

    TH1F *fQACosPA_lambda; //!
    TH1F *fQACosPA_Xi; //!

    TH1F *fQACosPA_lambda_L; //!
    TH1F *fQACosPA_Xi_L; //!

    TH1F *fQACosPA_lambda_T; //!
    TH1F *fQACosPA_Xi_T; //!

    TH1F *hXiInvMass; //!
    TH1F *hQAXiInvMass; //!

    TH2F *hTPCPID; //!
    TH2F *hTPCPIDpi; //!
    TH2F *hTPCPIDk; //!
    TH2F *hTPCPIDp; //!


    TH2F *hdEdxProton; //!
    TH2F *hdEdxPion1; //!
    TH2F *hdEdxPion2; //!

    TH2F *hdEdxProtonAfter; //!
    TH2F *hdEdxPion1After; //!
    TH2F *hdEdxPion2After; //!

    TH2F *hNSig3rdPion; //!
    TH2F *hQANSig3rdPion; //!

    TH1F *hTPCNSigProton; //!
    TH1F *hTPCNSigPion1; //!
    TH1F *hTPCNSigPion2; //!

    TH1F *hQATPCNSigProton; //!
    TH1F *hQATPCNSigPion1; //!
    TH1F *hQATPCNSigPion2; //!

    TH3F *hMCinputTotalXiStar1; //!
    TH3F *hMCinputTotalXiStarbar1; //!
    TH3F *hMCinputTotalXi1; //!
    TH3F *hMCinputTotalXibar1; //!

    TH3F *hMCinputTotalXiStar3; //!
    TH3F *hMCinputTotalXiStarbar3; //!
    TH3F *hMCinputTotalXi3; //!
    TH3F *hMCinputTotalXibar3; //!


    const char* fname; //! // name of class
    // AliInputEventHandler *fEventHandler;                              //  for ESDs or AODs
    AliESDEvent            *fESD; //!    // ESD object
    TList                  *fOutputList; //! Compact Output list
    AliESDtrackCuts        *fTrackCut; //! ESD track cuts
    AliPIDResponse         *fPIDResponse; //! PID object

    AliXiStarpp13TeVDevelEventCollection ***fEC; //! The event collection
    AliXiStarpp13TeVDevelEventStruct *fEvt; //! The current event type
    AliXiStarpp13TeVDevelTrackStruct *fTempStruct; //! A temporary track storage.  Eventually put into fEvt

    //

    Int_t fZvertexBins; //! // number of Z-vertex bins for event-mixing
    Int_t fEventsToMix; //! // number of maximum events to mix
    Int_t fMultBins; //! number of multiplicity bins for event-mixing
    Int_t fMultLimits[11 + 1]; //! the multiplicity edges of the mult bins
    Bool_t fMCcase; //! switch for MC data or real data
    Bool_t fAODcase; //! switch for AODs or ESDs
    Bool_t fDevelopeMode; //!
    Bool_t fHMTrigger; //!
    Bool_t fPIDOption; //!
    Bool_t fSetSystematic; //!
    Int_t fEventCounter; //! The event counter
    ULong64_t fEventNumber; //! calcuate event number

    // cut list data members
    Float_t fMaxDecayLength; //! max decay length
    Float_t fMassWindow; //! Mass window of acceptance for Lambda and Xi candidates

    Double_t fCovMatrix[21]; //! Covarience matrix of track
    Double_t fTrueMassPr, fTrueMassPi, fTrueMassK, fTrueMassLam, fTrueMassXi; //! The PDG mass values
    Bool_t IsTPC  (AliESDtrack *track);






    AliESDtrack* fESDTrack4; //! esdtrack for XiStar's daughter pion
    AliESDtrack* fXiTrack; //! esdtrack for XiStar's daughter Xi

    Int_t fCutList; //! Cut List option (mean values or systematic variations)

    Float_t fDecayParameters[kNCuts]; //! array of reconstruction kinematics
    Float_t fCutValues[kNCutVariations][kNCuts]; //! array of reconstruction kinematics

    ClassDef(AliXiStarpp13TeVDevel, 20);
};

#endif