#ifndef SENSVOLLARGEHIT_H
#define SENSVOLLARGEHIT_H 1
#include "FairVolume.h"
#include "ShipHit.h"
#include "SensVollargePoint.h"
#include "TObject.h"
#include "TGeoShape.h"
#include "TGeoPhysicalNode.h"


class SensVollargeHit : public ShipHit
{
  public:

    /** Default constructor **/
    SensVollargeHit();

    /** Constructor with arguments
     *@param detID    Detector ID
     *@param digi      digitized/measured ADC
     *@param flag      True/False, false if below threshold
     **/
    SensVollargeHit(Int_t detID, Float_t adc);

    /** Destructor **/
    virtual ~SensVollargeHit();

    /** Accessors **/
    Double_t GetX();
    Double_t GetY();
    Double_t GetZ();
    TVector3 GetXYZ();
    TGeoNode* GetNode();
    /** Modifier **/
    void SetEloss(Double_t val){fdigi=val;}
    void SetTDC(Double_t val){ft=val;}

    /** Output to screen **/

    virtual void Print(Int_t detID) const;
    Float_t GetADC() const {return fdigi;}
    Float_t GetTDC() const {return ft;}
    Double_t GetEloss() {return fdigi;}
    void setInvalid() {flag = false;}
    void setIsValid() {flag = true;}
    bool isValid() const {return flag;}
  private:
    Double_t ft;
    SensVollargeHit(const SensVollargeHit& point);
    SensVollargeHit operator=(const SensVollargeHit& point);

    Float_t flag;   ///< flag

    ClassDef(SensVollargeHit,1);

};

#endif
