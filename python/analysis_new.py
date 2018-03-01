import ROOT as r
import math
import os,sys,getopt
import rootUtils as ut
import shipunit as u
inputFile  = None
geoFile    = None
dy         = None
import shipVeto
import rootUtils as ut
import shipunit as u
from ShipGeoConfig import ConfigRegistry
import shipRoot_conf
shipRoot_conf.configure()
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3d




try:
        opts, args = getopt.getopt(sys.argv[1:], "n:f:g:A:Y:i", ["nEvents=","geoFile="])
except getopt.GetoptError:
        # print help information and exit:
        print ' enter file name'
        sys.exit()
for o, a in opts:
        if o in ("-f",):
            inputFile = a
        if o in ("-g", "--geoFile",):
            geoFile = a
        if o in ("-Y",):
            dy = float(a)


countSBT=0
countSVT=0
countUVT=0
countRPC=0

FidCutEv=0
DocaCutEv=0
hdocaIPCut=0
IPCutEv=0#define histograms
SBTCutEv=0
selectedWtihoutSBT=0
selectedEv=0



hcounmu=r.TH1D('hcountmu','Number of muons hitting spectrometer from 0 track', 100,0,100)
hcountSBT=r.TH1D('hcountSBT','Number of events which are rejected as a background using the SBT', 20000,0,20000)
hcountUVT=r.TH1D('hcountUVT','Number of events which are rejected as a background using the UVT', 20000,0,20000)
hcountRPC=r.TH1D('hcountRPC','Number of events which are rejected as a background using the RPC', 20000,0,20000)
hcountSVT=r.TH1D('hcountSVT','Number of events which are rejected as a background using the SVT', 20000,0,20000)

hselectedHNL=r.TH1D('hselectedHNL', 'The number of events selected as HNL ', 5000,0.,5000)
id_DIS=r.TH1D('id_DIS ','The Id of the particle produced in the DIS ;Id ',5000,-5000,5000)
hvtxrec=r.TH2D('hvtxrec ','The reconstructed vertex of all HNL candidates ;Z [cm];X[cm]',800,-2000.,2000.,12000,-30000,30000)
hdoca=r.TH1D('hdoca','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hip=r.TH1D('hip','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
hdocaFidCut=r.TH1D('hdocaFidCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidCut=r.TH1D('hdocaNoFidCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hipFidCut5cm=r.TH1D('hipFidCut5cm','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
hipNoFidCut5cm=r.TH1D('hipNoFidCut5cm','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
hipDocaCut=r.TH1D('hipDocaCut','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
hipNoDocaCut=r.TH1D('hipNoDocaCut','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
hdocaIPCut=r.TH1D('hdocaIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIPCut=r.TH1D('hdocaNoIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaIP250Cut=r.TH1D('hdocaIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIP250Cut=r.TH1D('hdocaNoIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
mass_ip=r.TH2D('mass_ip',";mass[GeV/c];IP[cm]",100,0,10,500,0.,4000);
notmass_ip=r.TH2D('notmass_ip',";mass[GeV/c];IP[cm]",100,0,10,500,0.,4000);
hdocaSBTcut=r.TH1D('hdocaSBTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNotSBTcut=r.TH1D('hdocaNotSBTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaFidDocaCut=r.TH1D('hdocaFidDocaCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidDocaCut=r.TH1D('hdocaNoFidDocaCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaFidIPCut=r.TH1D('hdocaFidIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidIPCut=r.TH1D('hdocaNoFidIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaFidSBTCut=r.TH1D('hdocaFidSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidSBTCut=r.TH1D('hdocaNoFidSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaFidVetoCut=r.TH1D('hdocaFidVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidVetoCut=r.TH1D('hdocaNoFidVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaDocaIPCut=r.TH1D('hdocaDocaIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoDocaIPCut=r.TH1D('hdocaNoDocaIPCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaDocaVetoCut=r.TH1D('hdocaDocaVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoDocaVetoCut=r.TH1D('hdocaNoDocaVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaDocaSBTCut=r.TH1D('hdocaDocaSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocanNoDocaSBTCut=r.TH1D('hdocaNoDocaSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaFidIP250Cut=r.TH1D('hdocaFidIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoFidIP250Cut=r.TH1D('hdocaNoFidIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaDocaIP250Cut=r.TH1D('hdocaDocaIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoDocaIP250Cut=r.TH1D('hdocaNoDocaIP250Cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)


hdocaIPVetoCut=r.TH1D('hdocaIPVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIPVetoCut=r.TH1D('hdocaNoIPVetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaIPSBTCut=r.TH1D('hdocaIPSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIPSBTCut=r.TH1D('hdocaNoIPSBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaIP250VetoCut=r.TH1D('hdocaIP250VetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIP250VetoCut=r.TH1D('hdocaNoIP250VetoCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaIP250SBTCut=r.TH1D('hdocaIP250SBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
hdocaNoIP250SBTCut=r.TH1D('hdocaNoIP250SBTCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)



hmomentumCandidate=r.TH1D('hmomentumCandidate', 'The momentum of the HNL candidate;P[GeV]', 8000,0.,400.)
hHNLevent_2track=r.TH2D('hHNLevent_2track', 'The number of NHL candidates for event vs number of tracks per event',300,0.,300.,20,0.,20.) 




hdocaUVTcut=r.TH1D('hdocaUVTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaSVTcut=r.TH1D('hdocaSVTcut', 'Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)




hmass_rec=r.TH1D('hmass_rec','The invarian mass of the HNL candidate ;M[GeV]; Nentries',100,0.,10)


mass_ip_Noveto=r.TH2D('mass_ip_Noveto',";mass[GeV/c];IP[cm]",100,0,10,500,0.,4000);




vetoDets = {"UVT":(False,0),"SVT":(False,0),"SBT":(False,0),"RPC":(False,0),"TRA":(False,0)}

hdoca2cut=r.TH1D('hdoca2cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hIP5cm=r.TH1D('hIP5cm', 'IP of HNL candidate ;IP[cm];Entries ',200,0.,1000)
hdoca3cut=r.TH1D('hdoca3cut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)

hdocaIP10=r.TH1D('hdocaIP10', 'The doca of the HNL event;DOCA[cm];Nentries',200,0.,1000)
docaNo5cmCut=r.TH1D('docaNo5cmCut','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
ipNo5cmCut=r.TH1D('ipNo5cmCut', 'IP of HNL candidate ;IP[cm];Entries ',200,0.,1000)




import shipVeto
f = r.TFile(inputFile)
sTree = f.Get("cbmsim")
fgeo=r.TFile(geoFile)
fGeo=fgeo.FAIRGeom
ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py",Yheight = 10, tankDesign = 5 , muShieldDesign = 7, nuTauTargetDesign=1)
startDecayVol=ShipGeo.vetoStation.z+20.*u.cm
trackST1=ShipGeo.TrackStation1.z-20*u.cm

print "Start of the decay Vessel=", startDecayVol

	
#define all functions that you need for the analysis
from array import array

def isInFiducial(Z):
   if Z > trackST1 : return False
   if Z < startDecayVol: return False
   return True 


def WeightsAnalysis(P):

        sigmaTot=0	
	fcross = r.TFile.Open('/afs/cern.ch/work/p/pvenkova/12collaborationMeeting/muDIScrossSec.root')
	for x in fcross.GetListOfKeys():
	       sigma = x.ReadObj()
	       if not sigma.GetTitle() == 'PythiaCross' : break  #'AnalyticCross'
	       sig=sigma.Eval(P)
	       #sigmaTot=sigmaTot+sig
	       #print "the moemntum=", p, "sig=",sig, "sigmaTot=", sigmaTot

	return sig




def ImpactParameter(point,tPos,tMom):
  t = 0
  if hasattr(tMom,'P'): 
  	P = tMom.P()

  else:                 
  	P = tMom.Mag()
  for i in range(3):   
  	t += tMom(i)/P*(point(i)-tPos(i))
  dist = 0
  for i in range(3):  
  	dist += (point(i)-tPos(i)-t*tMom(i)/P)**2
  dist = r.TMath.Sqrt(dist)
  return dist




# using Thomas file 1e-30
def prob2int(cross,weight,muWeight):
	prob2int=(cross*1e-27*6.022e+23*weight)
	return prob2int

def Estimate_phi_theta_muon(px,py,pz):
	
    	P=math.sqrt(px*px+py*py+pz*pz)
    	theta=r.TMath.ACos(pz / P)
    	phi=r.TMath.ATan2(py, px)
	return theta,phi

def WeightPZM(z,p,muWeight):
	weightZ=0
	weightP=0
	weightM=0
    	if muWeight< 3000: 
    		weightM=1.67524814479
    	else : 
    		weightM=0.788198666261


    	if z/100 < -23: 
    		weightZ=6.961
    	else : 
    		weightZ=0.3798
	
    	if P < 12: weightP=0.0885969
    	if P>12 and P< 250: weightP=0.964929
    	if P>250: weightP=37.945
	
	#print "The momentum is =", p, "The z position is =", z,"the muweight is=",muWeight, "The weightZPM=", weightP, weightZ, weightM ,"The weihghjtZP=",weightP*weightZ*weightM
	return weightP*weightZ*weightM

def weightDIS(Nint,phi,theta,scale,wPZ):
	w=0
	if theta<0.09 and phi< 0.53:
		w=0.611
	if theta>0.09 and phi< 0.53:
		w=0
	if theta<0.09 and (0.53<phi<1.06):
	 	w=1.5123
	if theta>0.09 and (0.53<phi<1.06):
		w=25.0949
	if theta <0.09 and (1.06<phi<1.6):
		w=0.713
	if theta >0.09 and (1.06<phi<1.6):
		w=224.8

	#print "The phi,theta are=", phi,theta,"The weight is=", w, "the probabiluty =", Nint, "the scale factor is=", scale
	return Nint*w*scale*wPZ


def scaleFactor(muWeight):
	scale=muWeight*(float(20)/267512)
	return scale



V0dict={130:'KL',310:'KS',3122:'Lambda',321:'K+',22:'gamma'}




veto=shipVeto.Task(sTree)
ev_dic={}
ev=0
ev_dic1=[]
ev_det={}
detDigi={}
numHitSeg=1
detDigiNocut={}
numHitSegNoCut=1
NintSimulated=0



for nb in sTree:
    muWeight=sTree.MCTrack[1].GetWeight()
    weight=sTree.MCTrack[3].GetWeight()
    cross=sTree.MCTrack[0].GetWeight()
    scale=scaleFactor(muWeight)
    prob=prob2int(cross,weight,muWeight)
    P=sTree.MCTrack[0].GetP()
    Z=sTree.MCTrack[0].GetStartZ()
    pz=sTree.MCTrack[0].GetPz() 
    px=sTree.MCTrack[0].GetPx()
    py=sTree.MCTrack[0].GetPy()
    theta=r.TMath.ACos(pz / P)
    phi=r.TMath.ATan2(py, px)
    if abs(phi)<1.57: 
    	phi=abs(phi)
    phi=abs(1.57- abs(phi))
    wPZ=WeightPZM(Z,P,muWeight)
    weightNoMagnet=weightDIS(prob,phi,theta,scale,wPZ)
    
    #print "the prob=", prob, "The mu Weight=", muWeight, "scale factor is=", scale, "The total weight is=", weightNoMagnet

    #print "============================================================================================================================"

    nHNLevent=sTree.Particles.GetEntries()
    numTrack=len(nb.goodTracks)
    if numTrack!=2:continue
    #if sTree.Particles.GetEntries()!=1:continue
    for candidate in sTree.Particles:
        indexCandidate=sTree.Particles.index(candidate)
 	vtx = r.TVector3()
	momentumHNL=r.TLorentzVector()
	candidate.GetVertex(vtx)
	candidate.GetMomentum(momentumHNL)
	mass_rec=momentumHNL.Mag2()
	hmass_rec.Fill(mass_rec,weightNoMagnet)
	hvtxrec.Fill(vtx.Z(),vtx.X(),weightNoMagnet)
	doca=candidate.GetDoca()
	hdoca.Fill(doca,weightNoMagnet)
	vtarget=r.TVector3(0,0,ShipGeo.target.z0)
	ip=ImpactParameter(vtarget,vtx,momentumHNL)
	hip.Fill(ip,weightNoMagnet)
	hmomentumCandidate.Fill(momentumHNL.P(),weightNoMagnet)
	distToWall=veto.fiducialCheckSignal(indexCandidate)
    	vetoDets['SBT'] = veto.SBT_decision()
    	vetoDets['SVT'] = veto.SVT_decision()
    	vetoDets['UVT'] = veto.UVT_decision()
    	vetoDets['RPC'] = veto.RPC_decision() 
    	if vetoDets['SBT'][2] != 0:countSBT+=1
    	if vetoDets['SVT'][2] != 0:countSVT+=1
    	if vetoDets['UVT'][2] != 0:countUVT+=1
    	if vetoDets['RPC'][2] != 0:countRPC+=1

	
	if not ( distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True):
		hipFidCut5cm.Fill(ip,weightNoMagnet)
		continue 
	#find N_A,N_B and NotN_A , 
	if distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True:

		FidCutEv=FidCutEv+1
		hdocaFidCut.Fill(doca,weightNoMagnet)
	else :
		hipNoFidCut5cm.Fill(ip,weightNoMagnet)
		hdocaNoFidCut.Fill(doca,weightNoMagnet)
			
	if doca<1:
		hipDocaCut.Fill(ip,weightNoMagnet)
		DocaCutEv=DocaCutEv+1
	else:
		hipNoDocaCut.Fill(ip,weightNoMagnet)

	
	if ip<10:
		hdocaIPCut.Fill(doca,weightNoMagnet)
		IPCutEv=IPCutEv+1
	else:
		hdocaNoIPCut.Fill(doca,weightNoMagnet)

	if ip<250:
		hdocaIP250Cut.Fill(doca,weightNoMagnet)
	else:
				
		hdocaNoIP250Cut.Fill(doca,weightNoMagnet)
	
        #if vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=0 or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False:

	if vetoDets['SBT'][0] or  vetoDets['UVT'][0] or vetoDets['SVT'][0] or vetoDets['RPC'][0]:
		mass_ip.Fill(mass_rec,ip,weightNoMagnet)
	else:
	       	
		notmass_ip.Fill(mass_rec,ip,weightNoMagnet)
	
	if vetoDets['SBT'][0]!=True:
		hdocaSBTcut.Fill(doca,weightNoMagnet)
		SBTCutEv=SBTCutEv+1
	else:
		
		hdocaNotSBTcut.Fill(doca,weightNoMagnet)

	#now investigare , N_{NotA and NotB} for pair of cuts
        #==========================FiDucial====================================================================================================
	if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True) and doca<1:
		hdocaFidDocaCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoFidDocaCut.Fill(doca,weightNoMagnet)


	
	if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True) and ip<10:
		hdocaFidIPCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoFidIPCut.Fill(doca,weightNoMagnet)
	
	if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True) and vetoDets['SBT'][0]!=True:
		hdocaFidSBTCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoFidSBTCut.Fill(doca,weightNoMagnet)

	
	if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True) and (vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=False  or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False)  :
		hdocaFidVetoCut.Fill(doca,weightNoMagnet)
		
	else:
		
		hdocaNoFidVetoCut.Fill(doca,weightNoMagnet)
	#===========================DOCA==========================================================================================================
	if doca<1 and ip<10:
		hdocaDocaIPCut.Fill(doca,weightNoMagnet)
	else
		
		hdocaNoDocaIPCut.Fill(doca,weightNoMagnet)
		
	if doca<1 and (vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=0 or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False)  :
		hdocaDocaVetoCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoDocaVetoCut.Fill(doca,weightNoMagnet)
	
	if doca<1 and vetoDets['SBT'][0]!=False:
		hdocaDocaSBTCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoDocaVetoCut.Fill(doca,weightNoMagnet)

        #=============================IP=============================================================================================================

        
	if ip<10 and (vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=0 or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False)  :
		hdocaIPVetoCut.Fill(doca,weightNoMagnet)
	else:
		

		hdocaNoIPVetoCut.Fill(doca,weightNoMagnet)
	
	if ip<10 and vetoDets['SBT'][0]!=False:
		hdocaIPSBTCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoIPSBTCut.Fill(doca,weightNoMagnet)

	#=============================IP250========================================================================================================
	if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True) and ip<250:
		hdocaFidIP250Cut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoFidIP250Cut.Fill(doca,weightNoMagnet)

	
	if doca<1 and ip<250:
		hdocaDocaIP250Cut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoDocaIP250Cut.Fill(doca,weightNoMagnet)
	if ip<250 and vetoDets['SBT'][0]!=False and vetoDets['UVT'][0]==False and  vetoDets['SVT'][0]==False and vetoDets['RPC'][0]==False  :
		hdocaIP250VetoCut.Fill(doca,weightNoMagnet)
	else:
		

		hdocaNoIP250VetoCut.Fill(doca,weightNoMagnet)
	
	if ip<250 and vetoDets['SBT'][0]!=False:
		hdocaIP250SBTCut.Fill(doca,weightNoMagnet)
	else:
		
		hdocaNoIP250SBTCut.Fill(doca,weightNoMagnet)

	




	

	'''if isInFiducial(vtx.Z())==True:

		hipFidCut5cm.Fill(prob*scale)
                hipFidCut.Fill(ip,prob*scale)
                FidCutEv=FidCutEv+1
                hdocaFidCut.Fill(doca,prob*scale)
	if vetoDets['UVT'][0]!=False:hdocaUVTcut.Fill(doca,prob*scale)
        if vetoDets['SBT'][0]!=False or vetoDets['UVT'][0]!=0 or vetoDets['SVT'][0]!=False or vetoDets['RPC'][0]!=False:
		mass_ip.Fill(mass_rec,ip,prob*scale)
	if vetoDets['SVT'][0]!=False:hdocaSVTcut.Fill(doca,prob*scale)
	if ip>250:continue
	hdoca3cut.Fill(doca,prob*muWeight*0.001)
	selectedWtihoutSBT=selectedWtihoutSBT+1
        if vetoDets['SBT'][0] or vetoDets['UVT'][0] or vetoDets['SVT'][0] or vetoDets['RPC'][0]:continue
	selectedEv=selectedEv+1
	hselectedHNL.Fill(selectedEv)
	print "Found HNL=", selectedEv'''
hcountSBT.Fill(countSBT)
hcountUVT.Fill(countUVT)
hcountRPC.Fill(countRPC)
hcountSVT.Fill(countSVT)
outfile=r.TFile("muonDISstudy.root", "RECREATE")
hmass_rec.Write()
hdoca.Write()
hvtxrec.Write()
hip.Write()
hipFidCut5cm.Write()
hdocaFidCut.Write()
hipNoFidCut5cm.Write()
hdocaNoFidCut.Write()
hipDocaCut.Write()
hipNoDocaCut.Write()
hdocaIPCut.Write()
hdocaIP250Cut.Write()
hdocaNoIP250Cut.Write()
hdocaNoIPCut.Write()
mass_ip.Write()
notmass_ip.Write()
hdocaSBTcut.Write()
hdocaNotSBTcut.Write()
hdocaFidDocaCut.Write()
hdocaNoFidDocaCut.Write()
hdocaFidIPCut.Write()
hdocaNoFidIPCut.Write()
hdocaFidSBTCut.Write()
hdocaNoFidSBTCut.Write()
hdocaFidVetoCut.Write()
hdocaDocaIPCut.Write()
hdocaNoDocaIPCut.Write()
hdocaNoDocaVetoCut.Write()
hdocaNoFidVetoCut.Write()
hdocaDocaSBTCut.Write()
hdocanNoDocaSBTCut.Write()
hdocaFidIP250Cut.Write()
hdocaNoFidIP250Cut.Write()
hdocaDocaIP250Cut.Write()
hdocaNoDocaIP250Cut.Write()
hdocaIPVetoCut.Write()
hdocaNoIPVetoCut.Write()
hdocaIPSBTCut.Write()
hdocaNoIPSBTCut.Write()
hdocaIP250VetoCut.Write()
hdocaIP250SBTCut.Write()
hdocaNoIP250SBTCut.Write()
hdocaNoIP250VetoCut.Write()



