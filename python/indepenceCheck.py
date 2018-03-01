import ROOT as r
import math
import rootUtils as ut
import argparse
import shipVeto
from ShipGeoConfig import ConfigRegistry
import shipunit as u
from numpy import ndarray
import numpy as np

r.PyConfig.IgnoreCommandLineOptions = True
import gc

def IndependenceCheck():
	
	def prob2int(cross,weight):
		#multiple the cross section in mbarn
		prob2int=(cross*1e-27*6.022e+23*weight)
		return prob2int

	
	#DEFINE THE WEIGHT OF MuonWeight,Z,P,RoL
	def WeightPZM(z,p,muWeight,wrl,theta,phi):
		weightZ=0
		weightP=0
		weightM=0
		weightRoL=0
		weightThetaPhi=0
		#weight for Muon Number 
		if muWeight< 3000: 
			weightM=1.67524814479
		else : 
			weightM=0.788198666261

		#weight for the z position
		if z/100 < -23: 
			weightZ=11.0107
		else : 
			weightZ=0.37196

		#weight for the Momentum
		if P < 46: weightP=0.97696#1.01869
		if P>46 and P< 270: weightP=0.385#0.3022
		if P>270 and P< 282: weightP=989.61#906.7063
		if P>283: weightP=0	
		if wrl>0 and wrl<4000:weightRoL= 0.8585#1.0655 #0.75/0.67324
		if wrl>4001 and wrl<6400:weightRoL=1.305#1.3272 #0.15/0.212989
		if wrl>6401 and wrl<10000:weightRoL=2.27#0.879 #0.1/0.114
		if wrl>10000:weightRoL=0
		
		#weight for the Theta,Phi
		if theta<0.09 and phi< 0.53:
			weightThetaPhi=0.7396
		if theta>0.09 and phi< 0.53:
			weightThetaPhi=0
		if theta<0.09 and (0.53<phi<1.06):
			weightThetaPhi=23.8957
		if theta>0.09 and (0.53<phi<1.06):
			weightThetaPhi=22.1094
		if theta <0.09 and (1.06<phi<1.6):
			weightThetaPhi=0.6655
		if theta >0.09 and (1.06<phi<1.6):
			weightThetaPhi=152.854

		return weightM*weightZ*weightP*weightRoL*weightThetaPhi


	def isInFiducial(Z):
  		if Z > trackST1 : return False
   		if Z < startDecayVol: return False
   		return True 

	
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


	
	parser= argparse.ArgumentParser(description='Script to check Independency of the cuts')
    	parser.add_argument(
		'-g',
        	'--geofile',
		 help='''Geometry file to use. ''')
	
    	parser.add_argument(
		'-f',
        	'--inputFile',
		 help='''Input file to use. ''')
	
	
    	parser.add_argument(
        	'-o',
        	'--outputfile',
        	default='checkIndpency.root')
    	args = parser.parse_args()

    	g = r.TFile.Open(args.geofile, 'read')
    	sGeo = g.FAIRGeom
	f = r.TFile.Open(args.inputFile, 'read')
    	o = r.TFile.Open(args.outputfile, 'recreate')
	
	#define the histograms
   	h={}
	
	ut.bookHist(h,'hmass_rec','The invarian mass of the HNL candidate ;M[GeV]; Nentries',100,0.,10)
	ut.bookHist(h,'hdoca','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hip','The Impact Parameter  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNdoca','Distance of closest aproach for events inside the Fid Volume   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotdoca','Distance of closest aproach for events inside the Fid Volume not passing the DOCA cut  ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNip10','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<10)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNnotip10','The Impact Parameter for events inside the Fiducial Volume not  passing the ip cut(ip<10)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNip250','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<250)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNnotip250','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<250)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNip30','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<30)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNnotip30','The Impact Parameter for events inside the Fiducial Volume not  passing the ip cut(ip<30)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNveto', 'The Impact Parameter for events inside the Fiducial Volume passing the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNnotveto','The Impact Parameter for events inside the Fiducial Volume not  passing  the veto cut  ;IP[cm]; Nentries',500,0.,4000) 
	ut.bookHist(h,'hNsbt','The Impact Parameter for events inside the Fiducial Volume  passing  the SBT cut  ;IP[cm]; Nentries',500,0.,4000) 
	ut.bookHist(h,'hNnotsbt','The Impact Parameter for events inside the Fiducial Volume not  passing  the SBT cut  ;IP[cm]; Nentries',500,0.,4000) 
























	
	ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py",Yheight = 10, tankDesign = 5 , muShieldDesign = 7, nuTauTargetDesign=1)
	startDecayVol=ShipGeo.vetoStation.z+20.*u.cm
	trackST1=ShipGeo.TrackStation1.z-20*u.cm
	
	ch = r.TChain('cbmsim')
	ch.Add(args.inputFile)
	nev = ch.GetEntries()
	o.cd()
	print "The number of events is=", nev

	veto=shipVeto.Task(ch)
	#cout the number of events inside Fid Volume
	NdocaFidV=0
	NnotNdocaFidV=0
	Nip10FidV=0
	Nnotip10FidV=0
	Nip30FidV=0
	Nnotip30FidV=0
	Nip250FidV=0
	Nnotip250FidV=0
	NvetoFidV=0
	NnotvetoFidV=0
	NsbtFidV=0
	NnotsbtFidV=0
	#cout the number of events outside the Fid Volume
	NdocaNotFidV=0
	NnotNdocaNotFidV=0
	Nip10NotFidV=0
	Nnotip10NotFidV=0
	Nip30NotFidV=0
	Nnotip30NotFidV=0
	Nip250NotFidV=0
	Nnotip250NotFidV=0
	NvetoNotFidV=0
	NnotvetoNotFidV=0
	NsbtNotFidV=0
	NnotsbtNotFidV=0

	





	for event in ch:

		
		#the initial muon having the weight to have full statistics of the full spill
		muWeight=event.MCTrack[1].GetWeight()
		#the cross section for DIS(Pythia[mbarn] 
		cross=event.MCTrack[0].GetWeight()
		# the rhoL of the initial muon 
    		weight=event.MCTrack[3].GetWeight()
		# the X,Y,Z momentum components, Z position of the initial muon need it to estimate theta, phi to get the weight 
    		Z=event.MCTrack[0].GetStartZ()
    		X=event.MCTrack[0].GetStartX()
    		Y=event.MCTrack[0].GetStartY()
		#the probability for DIS 
		prob=prob2int(cross,weight)
		P=event.MCTrack[0].GetP()
    		pz=event.MCTrack[0].GetPz() 
    		px=event.MCTrack[0].GetPx()
    		py=event.MCTrack[0].GetPy()
    		theta=r.TMath.ACos(pz / P)
   		phi=r.TMath.ATan2(py, px)
		# Weight to account difference beetwen Magnets ON and Magnets OFF
		wPZMThPhi= WeightPZM(Z,P,muWeight,weight,theta,phi)
		totalWeight=prob*muWeight*wPZMThPhi
		numTrack=len(event.goodTracks)
    		if numTrack!=2:continue
		
		for candidate in event.Particles:

        		#define 3,4 Momentum Vector,vtarget
 			vtx = r.TVector3()
			momentumHNL=r.TLorentzVector()
			vtarget=r.TVector3(0,0,ShipGeo.target.z0)
			indexCandidate=event.Particles.index(candidate)
			candidate.GetVertex(vtx)
			candidate.GetMomentum(momentumHNL)
			mass_rec=momentumHNL.Mag2()
			doca=candidate.GetDoca()
			ip=ImpactParameter(vtarget,vtx,momentumHNL)
			distToWall=veto.fiducialCheckSignal(indexCandidate)
			
			h['hmass_rec'].Fill(mass_rec,totalWeight)
			h['hdoca'].Fill(doca,totalWeight)		
			h['hip'].Fill(ip,totalWeight)
			
			#Use veto systems 
			vetoSBT = veto.SBT_decision()
			vetoSVT = veto.SVT_decision()
			vetoUVT = veto.UVT_decision()
			vetoRPC = veto.RPC_decision() 
			#devide the background to two type inside the the Fiducial Volume and outside 
			if (distToWall>5*u.cm and distToWall!=0 and isInFiducial(vtx.Z())==True):

			#n_cut, n_notcut for background having vertex inside the fiducial volume	
				if doca<1:
					h['hNdoca'].Fill(doca,totalWeight)
				else:
					h['hNnotdoca'].Fill(doca,totalWeight)
					


				if ip<10:
					h['hNip10'].Fill(ip,totalWeight)
				else:
					h['hNnotip10'].Fill(ip,totalWeight)


				if ip<250:

					h['hNip250'].Fill(ip,totalWeight)
				
				else:
					h['hNnotip250'].Fill(ip,totalWeight)
						
			
				if ip<30:
					h['hNip30'].Fill(ip,totalWeight)
				else:
					h['hNnotip30'].Fill(ip,totalWeight)
					
				if vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					h['hNveto'].Fill(ip,totalWeight)
				else:
					h['hNnotveto'].Fill(ip,totalWeight)
					
			
				if vetoSBT == False:
					h['hNsbt'].Fill(doca,totalWeight)
				else:
					h['hNnotsbt'].Fill(doca,totalWeight)
					






			else:
				print "THHHHD"








	for key in h:
		h[key].Write()
	o.Close()


if __name__ == '__main__':
    IndependenceCheck()
	
IndependenceCheck()



