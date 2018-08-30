#ifndef ALIXISTARPP13TEVDEVELEVENTCOLLECTION_H
#define ALIXISTARPP13TEVDEVELEVENTCOLLECTION_H
//
// Class AliXiStarpp13TeVDevelEventCollection, AliXiStarpp13TeVDevelTrackStruct, AliXiStarpp13TeVDevelEventStruct
//
// AliXiStarpp13TeVDevelEventCollection, AliXiStarpp13TeVDevelTrackStruct, AliXiStarpp13TeVDevelEventStruct
// author:
//  (Original Code) Dhevan Gangadharan (dhevan.raja.gangadharan@cern.ch)
//  (1st Modification) Jihye Song (jihye.song@cern.ch)
//  (last Modification) Bong-Hwi Lim (bong-hwi.lim@cern.ch)

#include <iostream>
#include <string>
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TBits.h"
#include "TObject.h"
#include "TVector2.h"
#include "AliESDtrack.h"

using namespace std;


class AliXiStarpp13TeVDevelTrackStruct{
 public:

  AliXiStarpp13TeVDevelTrackStruct();
  virtual ~AliXiStarpp13TeVDevelTrackStruct();
  AliXiStarpp13TeVDevelTrackStruct(const AliXiStarpp13TeVDevelTrackStruct &obj ); 
  AliXiStarpp13TeVDevelTrackStruct &operator=(const AliXiStarpp13TeVDevelTrackStruct &obj );

  UInt_t fStatus;// track status
  UInt_t fFilterMap;// filter map for AOD filterbits
  Int_t fID;// track id
  Double_t fPhi;// track phi angle
  Double_t fPt;// track pt
  Float_t fMom;// track full momentum
  Double_t fP[3];// track 3d momentum
  Int_t fCharge;// track charge
  Double_t fEta;// track eta
  Double_t fMass;// track accepted mass
  Double_t fDCAXY;// track dca to PV in xy
  Double_t fDCAZ;// track dca to PV in z
  Double_t fDCA;// track full dca
  Double_t fX[3];// track x position
  Double_t fCov[21];// track covariance matrix
  Float_t fNSigmaPi;// track Nsigma pion
  Float_t fNSigmaK;// track Nsigma kaon
  Float_t fNSigmaPr;// track Nsigma proton
  Int_t fLabel;// track label for MC studies
  UShort_t fNclusTPC;// TPC N clusters

  ClassDef(AliXiStarpp13TeVDevelTrackStruct, 1);
};

class AliXiStarpp13TeVDevelEventStruct{
 public:

  AliXiStarpp13TeVDevelEventStruct();
  virtual ~AliXiStarpp13TeVDevelEventStruct();
  AliXiStarpp13TeVDevelEventStruct(const AliXiStarpp13TeVDevelEventStruct &obj ); 
  AliXiStarpp13TeVDevelEventStruct &operator=(const AliXiStarpp13TeVDevelEventStruct &obj );

  Int_t fNTracks;// Events track count
  AliXiStarpp13TeVDevelTrackStruct *fTracks;// Events track structure

  ClassDef(AliXiStarpp13TeVDevelEventStruct, 1);
};

class AliXiStarpp13TeVDevelEventCollection {
 public:
  
  AliXiStarpp13TeVDevelEventCollection();
  AliXiStarpp13TeVDevelEventCollection(Short_t);
  virtual ~AliXiStarpp13TeVDevelEventCollection();
  AliXiStarpp13TeVDevelEventCollection(const AliXiStarpp13TeVDevelEventCollection &obj ); 
  AliXiStarpp13TeVDevelEventCollection &operator=(const AliXiStarpp13TeVDevelEventCollection &obj );
  
  Short_t fFIFO; //Size of the Event Storage buffer. FIFO = first-in-first-out
  AliXiStarpp13TeVDevelEventStruct *fEvtStr;// Event structure

  void FIFOShift();// remove event at end of buffer and add the new one
  void SetBuffSize(Short_t a){fFIFO = a;}// set size of event buffer (Nevents max to mix)
          
  ClassDef(AliXiStarpp13TeVDevelEventCollection, 1);
};
#endif

















