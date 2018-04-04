#ifndef vetoHIT_H
#define vetoHIT_H 1
#include "FairVolume.h"
#include "ShipHit.h"
#include "vetoPoint.h"
#include "TObject.h"
#include "TGeoShape.h"
#include "TGeoPhysicalNode.h"

class vetoHit : public ShipHit
{
  public:

    /** Default constructor **/
    vetoHit();

    /** Constructor with arguments
     *@param detID    Detector ID
     *@param digi      digitized/measured ADC
     *@param flag      True/False, false if below threshold (5MeV)
     *@param flag1     True/False, false if below threshold (25MeV)
     *@param flag2     True/False, false if below threshold (45MeV)
     *@param flag3     True/False, false if below threshold (65MeV)
     *@param flag4     True/False, false if below threshold (0MeV)

     **/
    vetoHit(Int_t detID, Float_t adc);
    /** Destructor **/
    virtual ~vetoHit();

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
    Float_t adc() const {return fdigi;}
    Float_t tdc() const {return ft;}
    Double_t GetEloss() {return fdigi;}
    void setInvalid() {flag = false;}
    void setIsValid() {flag = true;}
    void setInvalid1() {flag1 = false;}
    void setIsValid1() {flag1 = true;}
    void setInvalid2() {flag2 = false;} 
    void setIsValid2() {flag2 = true;}
    void setInvalid3() {flag3 = false;}
    void setIsValid3() {flag3 = true;}
    void setInvalid4() {flag4 = false;}
    void setIsValid4() {flag4 = true;}


    bool isValid() const {return flag;}
    bool isValid1() const {return flag1;}
    bool isValid2() const {return flag2;}
    bool isValid3() const {return flag3;}
    bool isValid4() const {return flag4;} 
  private:
    Double_t ft;
    vetoHit(const vetoHit& point);
    vetoHit operator=(const vetoHit& point);

    Float_t flag;   ///< flag
    Float_t flag1;   ///< flag
    Float_t flag2;
    Float_t flag3;    
    Float_t flag4;    
    
    ClassDef(vetoHit,1);

};

#endif
