// Timing Detector
// 26/01/2017
// Alexander.Korzenev@cern.ch

#include "EmShield.h"
#include "EmShieldPoint.h"

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


EmShield::EmShield()
  : FairDetector("EmShield", kTRUE, kEmShield),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    //
    fzPos1(-2505), //-2505
    fxSize1(150),
    fySize1(200),
    fzSize1(5),
    //fxCenter1(0),
    //fyCenter1(0),
    fzPos2(-3488),
    fxSize2(250),
    fySize2(300),
    fzSize2(5),
    //ftotallength(22),
    //fxCenter2(0),
    //fyCenter2(0),
    fXtrap1(),
    ftrapregion1(),
    //
    fDetector(0),
    //
    fEmShieldPointCollection(new TClonesArray("EmShieldPoint"))
{
  /*NBars = fNCol * fNRow;
  if(fNCol>1) fxOv = (fxBar*fNCol - fxSize) / (double)(fNCol-1); else fxOv = 0;
  if(fNRow>1) fyOv = (fyBar*fNRow - fySize) / (double)(fNRow-1); else fyOv = 0;*/
  
  FairDetector::Initialize();
}



EmShield::EmShield(const char* name, Bool_t active) //***************************** CHANGE PARAMS
  : FairDetector(name, active, kEmShield),
    fTrackID(-1),
    fVolumeID(-1),
    fPos(),
    fMom(),
    fTime(-1.),
    fLength(-1.),
    fELoss(-1),
    //
    fzPos1(-2505),
    fxSize1(150),
    fySize1(200),
    fzSize1(5),
    //ftotallength(22),
    //fxCenter1(0),
    //fyCenter1(0),
    fzPos2(-3488),
    fxSize2(250),
    fySize2(300),
    fzSize2(5),
    fxCenter2(0),
    fyCenter2(0),
    fXtrap1(),
    ftrapregion1(),
    //
    fDetector(0),
    //
    fEmShieldPointCollection(new TClonesArray("EmShieldPoint"))
{
  FairDetector::Initialize();
}


EmShield::~EmShield()
{
  if (fEmShieldPointCollection) {
    fEmShieldPointCollection->Delete();
    delete fEmShieldPointCollection;
  }
}



Int_t EmShield::InitMedium(const char* name)
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



Bool_t  EmShield::ProcessHits(FairVolume* vol)
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
    Int_t pdgCode = p->GetPdgCode(); 
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
    stack->AddPoint(kEmShield);
  }
  
  return kTRUE;
}



void EmShield::EndOfEvent()
{
  fEmShieldPointCollection->Clear();
}



void EmShield::Register()
{

  /** This will create a branch in the output tree called
      EmShieldPoint, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("EmShieldPoint", "EmShield",
                                        fEmShieldPointCollection, kTRUE);
}



TClonesArray* EmShield::GetCollection(Int_t iColl) const
{
  if (iColl == 0) { return fEmShieldPointCollection; }
  else { return NULL; }
}



void EmShield::Reset()
{
  fEmShieldPointCollection->Clear();
}



void EmShield::ConstructGeometry()
{
  TGeoVolume *top = gGeoManager->GetTopVolume();
  
  InitMedium("Lead");
  TGeoMedium *Lead =gGeoManager->GetMedium("Lead");

  InitMedium("SensVacuum");
  TGeoMedium *Vac =gGeoManager->GetMedium("SensVacuum");
  
  ///////////////////////////////////////////////////////

  fDetector = new TGeoVolumeAssembly("Electromagnetic Shield");

//-------------------------------------HOURGLASS shape at -35---------------------------------------

  TGeoVolume *plateup = gGeoManager->MakeBox("EmShieldbox", Lead, fxSize1/2.,(fySize1-ftrapregion1)/4.,fzSize1/2.); 
  plateup->SetLineColor(kYellow);
  AddSensitiveVolume(plateup);
  fDetector->AddNode(plateup, 1, new TGeoTranslation( 0.,(fySize1+ftrapregion1)/4.,fzPos1) ); //upper box
  fDetector->AddNode(plateup, 1, new TGeoTranslation( 0.,(fySize1+ftrapregion1)/(-4.),fzPos1) ); //lower box


  TGeoRotation rotup;
  TGeoTranslation transup;
  rotup.RotateX(-90);
  transup.SetTranslation( 0.,ftrapregion1/4.,fzPos1 );
  TGeoCombiTrans combup(transup,rotup);
  TGeoHMatrix *hmatup= new TGeoHMatrix(combup);
  
  TGeoVolume *trapup = gGeoManager->MakeTrd2("EmShieldtrapup", Lead, fXtrap1/2. ,fxSize1/2.,fzSize1/2.,fzSize1/2., ftrapregion1/4.);
  trapup->SetLineColor(kYellow);
  AddSensitiveVolume(trapup);
  fDetector->AddNode(trapup, 1, hmatup ); //upper trapezoid


  TGeoRotation rotdown;
  TGeoTranslation transdown;
  rotdown.RotateX(-90);
  transdown.SetTranslation( 0.,ftrapregion1/(-4.),fzPos1);
  TGeoCombiTrans combdown(transdown,rotdown);
  TGeoHMatrix *hmatdown= new TGeoHMatrix(combdown);

  TGeoVolume *trapdown = gGeoManager->MakeTrd2("EmShieldtrapdown", Lead, fxSize1/2. ,fXtrap1/2.,fzSize1/2.,fzSize1/2., ftrapregion1/4.);
  trapdown->SetLineColor(kYellow);
  AddSensitiveVolume(trapdown);
  fDetector->AddNode(trapdown, 1, hmatdown ); //lower trapezoid








//-------------------------------------HOURGLASS shape at -35---------------------------------------
 /*   TGeoVolume* hourglass = new TGeoVolumeAssembly("hourglass");
  TGeoVolume* layered = new TGeoVolumeAssembly("layered");



  TGeoVolume *plateup = gGeoManager->MakeBox("EmShieldbox", Lead, fxSize1/2.,(fySize1-ftrapregion1)/4.,fzSize1/2.); 
  plateup->SetLineColor(kYellow);
  AddSensitiveVolume(plateup);
  hourglass->AddNode(plateup, 1, new TGeoTranslation( 0.,(fySize1+ftrapregion1)/4.,fzPos1) ); //upper box
  hourglass->AddNode(plateup, 1, new TGeoTranslation( 0.,(fySize1+ftrapregion1)/(-4.),fzPos1) ); //lower box


  TGeoRotation rotup;
  TGeoTranslation transup;
  rotup.RotateX(-90);
  transup.SetTranslation( 0.,ftrapregion1/4.,fzPos1 );
  TGeoCombiTrans combup(transup,rotup);
  TGeoHMatrix *hmatup= new TGeoHMatrix(combup);
  
  TGeoVolume *trapup = gGeoManager->MakeTrd2("EmShieldtrapup", Lead, fXtrap1/2. ,fxSize1/2.,fzSize1/2.,fzSize1/2., ftrapregion1/4.);
  trapup->SetLineColor(kYellow);
  AddSensitiveVolume(trapup);
  hourglass->AddNode(trapup, 1, hmatup ); //upper trapezoid


  TGeoRotation rotdown;
  TGeoTranslation transdown;
  rotdown.RotateX(-90);
  transdown.SetTranslation( 0.,ftrapregion1/(-4.),fzPos1);
  TGeoCombiTrans combdown(transdown,rotdown);
  TGeoHMatrix *hmatdown= new TGeoHMatrix(combdown);

  TGeoVolume *trapdown = gGeoManager->MakeTrd2("EmShieldtrapdown", Lead, fxSize1/2. ,fXtrap1/2.,fzSize1/2.,fzSize1/2., ftrapregion1/4.);
  trapdown->SetLineColor(kYellow);
  AddSensitiveVolume(trapdown);
  hourglass->AddNode(trapdown, 1, hmatdown ); //lower trapezoid
//-----------------------------------MULTIPLE LAYERS------------------------------------------ 
 
  Int_t NPlates=ftotallength/(fzSize1+0.1);//+ 1mm to take into account the gap with the active layers

 for(Int_t n=0; n<NPlates+1; n++)
    {
      layered->AddNode(hourglass, n, new TGeoTranslation(0,0,-ftotallength/2. +0.05 + (2*n+1)*(0.05 +fzSize1/2.)));
    }

    fDetector->AddNode(layered, 1, new TGeoTranslation( 0,0,0) );
*/

 /* Int_t NPlates=ftotallength/(fzSize1);

  for(Int_t n=0; n<NPlates+1; n++)
    {
      layered->AddNode(hourglass, n, new TGeoTranslation(0,0,-ftotallength/2.  + (2*n+1)*(fzSize1/2.)));
    }

    fDetector->AddNode(layered, 1, new TGeoTranslation( 0,0,0) );

*/

//-----------------------------------BOX SHAPE -25--------------------------------------------

 /* TGeoVolume *plate = gGeoManager->MakeBox("EmShield", Lead, fxSize2,fySize2,fzSize2);
  plate->SetLineColor(kYellow);
  AddSensitiveVolume(plate);
  fDetector->AddNode(plate, 1, new TGeoTranslation( fxCenter2,fyCenter2,fzPos2) ); //-25*/


 
  top->AddNode(fDetector, 1, new TGeoTranslation(0,0,0));

  ///////////////////////////////////////////////////////

  return;
}




EmShieldPoint* EmShield::AddHit(Int_t trackID, Int_t detID,
			TVector3 pos, TVector3 mom,
			Double_t time, Double_t length,
			Double_t eLoss, Int_t pdgCode,TVector3 Lpos, TVector3 Lmom)
{
  TClonesArray& clref = *fEmShieldPointCollection;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) EmShieldPoint(trackID, detID, pos, mom,
		         time, length, eLoss, pdgCode,Lpos,Lmom);
}





ClassImp(EmShield)
