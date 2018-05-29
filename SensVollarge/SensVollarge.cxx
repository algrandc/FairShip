// Timing Detector
// 26/01/2017
// Alexander.Korzenev@cern.ch

#include "SensVollarge.h"
#include "SensVollargePoint.h"

#include "FairVolume.h"
#include "FairGeoVolume.h"
#include "FairGeoNode.h"
#include "FairRootManager.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoMedia.h"
#include "FairGeoBuilder.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "ShipDetectorList.h"
#include "ShipStack.h"

#include "TClonesArray.h"
#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoTrd2.h"
#include "TGeoCompositeShape.h"
#include "TGeoTube.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TMath.h" 
#include "TParticle.h" 
#include "TVector3.h"

#include <iostream>
#include <sstream>
using std::cout;
using std::endl;


SensVollarge::SensVollarge()
  : FairDetector("SensVollarge", kTRUE, kSensVollarge),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fzPos2(-3488),
    fxSize2(250),
    fySize2(300),
    fzSize2(5),
    //
    fDetector(0),
    //
    fSensVollargePointCollection(new TClonesArray("SensVollargePoint"))
{
  /*NBars = fNCol * fNRow;
  if(fNCol>1) fxOv = (fxBar*fNCol - fxSize) / (double)(fNCol-1); else fxOv = 0;
  if(fNRow>1) fyOv = (fyBar*fNRow - fySize) / (double)(fNRow-1); else fyOv = 0;*/
  
  FairDetector::Initialize();
}



SensVollarge::SensVollarge(const char* name, Bool_t active) //***************************** CHANGE PARAMS
  : FairDetector(name, active, kSensVollarge),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    fzPos2(-3488),
    fxSize2(250),
    fySize2(300),
    fzSize2(5),
    fxCenter2(0),
    fyCenter2(0),
    //
    fDetector(0),
    //
    fSensVollargePointCollection(new TClonesArray("SensVollargePoint"))
{
  FairDetector::Initialize();
}


SensVollarge::~SensVollarge()
{
  if (fSensVollargePointCollection) {
    fSensVollargePointCollection->Delete();
    delete fSensVollargePointCollection;
  }
}



Int_t SensVollarge::InitMedium(const char* name)
{
  
   static FairGeoLoader *geoLoad=FairGeoLoader::Instance();
   static FairGeoInterface *geoFace=geoLoad->getGeoInterface();
   static FairGeoMedia *media=geoFace->getMedia();
   static FairGeoBuilder *geoBuild=geoLoad->getGeoBuilder();

   FairGeoMedium *ShipMedium=media->getMedium(name);

   if (!ShipMedium)
   {
     Fatal("InitMedium","Material %s not defined in media file.", name);
     return -1111;
   }
   TGeoMedium* medium=gGeoManager->GetMedium(name);
   if (medium!=NULL)
     return ShipMedium->getMediumIndex();

   return geoBuild->createMedium(ShipMedium);
  
  return 0;
}



Bool_t  SensVollarge::ProcessHits(FairVolume* vol)
{
  /** This method is called from the MC stepping */
  //Set parameters at entrance of volume. Reset ELoss.
  if ( gMC->IsTrackEntering() ) {
    fELoss  = 0.;
    fTime   = gMC->TrackTime() * 1.0e09;
    fLength = gMC->TrackLength();
    gMC->TrackPosition(fPos);
    gMC->TrackMomentum(fMom);
  }

  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();

  // Create vetoPoint at exit of active volume
  if ( gMC->IsTrackExiting()    ||
       gMC->IsTrackStop()       ||
       gMC->IsTrackDisappeared()   ) {
    if (fELoss == 0. ) { return kFALSE; }

    fTrackID  = gMC->GetStack()->GetCurrentTrackNumber();

    fVolumeID= vol->getMCid();
    Int_t uniqueId;
    gMC->CurrentVolID(uniqueId);
    Int_t detID=0;
    gMC->CurrentVolID(detID);
    if (fVolumeID == detID){
       return kTRUE;
    }
    gGeoManager->PrintOverlaps();
    TParticle* p = gMC->GetStack()->GetCurrentTrack();
    Int_t pdgCode = p->GetPdgCode(); //****************************** ADD HERE PYTHIA CODE
    TLorentzVector Pos;
    gMC->TrackPosition(Pos);
    TLorentzVector Mom;
    gMC->TrackMomentum(Mom);
    Double_t xmean = (fPos.X()+Pos.X())/2. ;
    Double_t ymean = (fPos.Y()+Pos.Y())/2. ;
    Double_t zmean = (fPos.Z()+Pos.Z())/2. ;



    AddHit(fTrackID, uniqueId, TVector3(xmean, ymean,  zmean),
           TVector3(fMom.Px(), fMom.Py(), fMom.Pz()), fTime, fLength,
           fELoss,pdgCode,TVector3(Pos.X(), Pos.Y(), Pos.Z()),
	   TVector3(Mom.Px(), Mom.Py(), Mom.Pz()) );

    // Increment number of Em det points in TParticle
    ShipStack* stack = (ShipStack*) gMC->GetStack();
    stack->AddPoint(kSensVollarge);
  }
  
  return kTRUE;
}



void SensVollarge::EndOfEvent()
{
  fSensVollargePointCollection->Clear();
}



void SensVollarge::Register()
{

  /** This will create a branch in the output tree called
      SensVollargePoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("SensVollargePoint", "SensVollarge",
                                        fSensVollargePointCollection, kTRUE);
}



TClonesArray* SensVollarge::GetCollection(Int_t iColl) const
{
  if (iColl == 0) { return fSensVollargePointCollection; }
  else { return NULL; }
}



void SensVollarge::Reset()
{
  fSensVollargePointCollection->Clear();
}



void SensVollarge::ConstructGeometry()
{
  TGeoVolume *top = gGeoManager->GetTopVolume();
  
  InitMedium("Lead");
  TGeoMedium *Lead =gGeoManager->GetMedium("Lead");

  //InitMedium("SensVacuum");
  //TGeoMedium *Vac =gGeoManager->GetMedium("SensVacuum");
  InitMedium("vacuum");
  TGeoMedium *Vac =gGeoManager->GetMedium("vacuum");  
  ///////////////////////////////////////////////////////

  fDetector = new TGeoVolumeAssembly("Large Sensitive Volume");





  TGeoVolume* layeredlarge = new TGeoVolumeAssembly("layeredlarge");

//-----------------------------------BOX SHAPE -35--------------------------------------------

  TGeoVolume *plate = gGeoManager->MakeBox("SensVollarge", Vac, fxSize2/2.,fySize2/2.,fzSize2/2.);
  plate->SetLineColor(kRed);
  AddSensitiveVolume(plate);


 /* Int_t NPlates=2;

  for(Int_t n=0; n<NPlates+1; n++)
    {
      layeredlarge->AddNode(plate, n, new TGeoTranslation(0,0,0.05 -ftotallength/2. + n*(fspacing+0.1) ));
    }*/
      layeredlarge->AddNode(plate, 1, new TGeoTranslation(0,0,-0.1 -ftotallength/2. ));
      layeredlarge->AddNode(plate, 2, new TGeoTranslation(0,0,0.1 +ftotallength/2.  ));

 /* Int_t NPlates=1+ftotallength/(fspacing);

  for(Int_t n=0; n<NPlates+1; n++)
    {
      layeredlarge->AddNode(plate, n, new TGeoTranslation(0,0,-ftotallength/2. + n*(fspacing) ));
    }
*/
    fDetector->AddNode(layeredlarge, 1, new TGeoTranslation( 0,0,fzPos2) );


 
  top->AddNode(fDetector, 1, new TGeoTranslation(0,0,0));

  ///////////////////////////////////////////////////////

  return;
}




SensVollargePoint* SensVollarge::AddHit(Int_t trackID, Int_t detID,
			TVector3 pos, TVector3 mom,
			Double_t time, Double_t length,
			Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom)
{
  TClonesArray& clref = *fSensVollargePointCollection;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) SensVollargePoint(trackID, detID, pos, mom,
		         time, length, eLoss, pdgCode,Lpos,Lmom);
}





ClassImp(SensVollarge)
