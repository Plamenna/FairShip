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
        	'-oIn',
        	'--outputfile',
        	default='checkIndpencyInsideFidV.root')
	
    	parser.add_argument(
        	'-oOut',
        	'--outputfile1',
        	default='checkIndpencyOutsideFidV.root')



    	args = parser.parse_args()

    	g = r.TFile.Open(args.geofile, 'read')
    	sGeo = g.FAIRGeom
	f = r.TFile.Open(args.inputFile, 'read')
    	oIn = r.TFile.Open(args.outputfile, 'recreate')
    	oOut = r.TFile.Open(args.outputfile1, 'recreate')
	
	#define the histograms for background INSIDE the Fiducial Volume
   	h={}
	
	ut.bookHist(h,'hmass_rec','The invarian mass of the HNL candidate ;M[GeV]; Nentries',100,0.,10)
	ut.bookHist(h,'hdoca','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hdocaInside','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
	
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
	#histograms to for pair of Cuts
	ut.bookHist(h,'hNdocaIP10','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotdocaIP10','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hNip10Veto','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<10) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNnotip10Veto','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<10) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNdocaVeto','Distance of closest aproach for events inside the Fid Volume passing veto and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hnotNdocaVeto','Distance of closest aproach for events inside the Fid Volume not passing veto and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hNip10SBT','Distance of closest aproach for events inside the Fid Volume passing ip  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotip10SBT','Distance of closest aproach for events inside the Fid Volume not passing ip  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hNdocaSBT','Distance of closest aproach for events inside the Fid Volume passing doca  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)

	
	ut.bookHist(h,'hNnotdocaSBT','Distance of closest aproach for events inside the Fid Volume not passing doca  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNdocaIP30','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotdocaIP30','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hNip30Veto','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNnotip30Veto','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNip30SBT','Distance of closest aproach for events inside the Fid Volume passing ip  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotip30SBT','Distance of closest aproach for events inside the Fid Volume not passing ip  and sbt cut   ;DOCA[cm]; Nentries',200,0.,1000)
	
	ut.bookHist(h,'hNdocaIP250','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNnotdocaIP250','Distance of closest aproach for events inside the Fid Volume passing IP and DOCA cut   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hNip250Veto','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNnotip250Veto','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(h,'hNip250SBT','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h,'hNnotip250SBT','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<30) and the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	
	ut.bookHist(h, 'hdocaIP', "The IP vs DOCA of all reconstructed HNL candidates; DOCA[cm]; IP[cm]",200,0.,1000,500,0.,4000)
	#HISTOGRAMS TO CROSS CHECK A|b
	ut.bookHist(h,'hIP10ifdoca','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hIP30ifdoca','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hIP250ifdoca','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hdocaifIP10','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hdocaifIP30','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'hdocaifIP250','Distance of closest aproach ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(h,'Dist2wall', 'Distance to the closest boundary;dist[cm]', 10,0.,1)





	
	#define the histograms for background OUTSIDE the Fiducial Volume
   	hOut={}
	ut.bookHist(hOut,'hdocaOutside','Distance of closest aproach  ;DOCA[cm]; Nentries',200,0.,1000)
	
	
	ut.bookHist(hOut,'hNdoca','Distance of closest aproach for events inside the Fid Volume   ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(hOut,'hNnotdoca','Distance of closest aproach for events inside the Fid Volume not passing the DOCA cut  ;DOCA[cm]; Nentries',200,0.,1000)
	ut.bookHist(hOut,'hNip10','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<10)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNnotip10','The Impact Parameter for events inside the Fiducial Volume not  passing the ip cut(ip<10)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNip250','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<250)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNnotip250','The Impact Parameter for events inside the Fiducial Volume not passing the ip cut(ip<250)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNip30','The Impact Parameter for events inside the Fiducial Volume passing the ip cut(ip<30)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNnotip30','The Impact Parameter for events inside the Fiducial Volume not  passing the ip cut(ip<30)  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNveto', 'The Impact Parameter for events inside the Fiducial Volume passing the veto cut  ;IP[cm]; Nentries',500,0.,4000)
	ut.bookHist(hOut,'hNnotveto','The Impact Parameter for events inside the Fiducial Volume not  passing  the veto cut  ;IP[cm]; Nentries',500,0.,4000) 
	ut.bookHist(hOut,'hNsbt','The Impact Parameter for events inside the Fiducial Volume  passing  the SBT cut  ;IP[cm]; Nentries',500,0.,4000) 
	ut.bookHist(hOut,'hNnotsbt','The Impact Parameter for events inside the Fiducial Volume not  passing  the SBT cut  ;IP[cm]; Nentries',500,0.,4000) 
	























	
	ShipGeo = ConfigRegistry.loadpy("$FAIRSHIP/geometry/geometry_config.py",Yheight = 10, tankDesign = 5 , muShieldDesign = 7, nuTauTargetDesign=1)
	startDecayVol=ShipGeo.vetoStation.z+20.*u.cm
	trackST1=ShipGeo.TrackStation1.z-20*u.cm
	
	ch = r.TChain('cbmsim')
	ch.Add(args.inputFile)
	nev = ch.GetEntries()
	oIn.cd()
	#oOut.cd()

	print "The number of events is=", nev

	veto=shipVeto.Task(ch)
	#cout the number of events inside Fid Volume
	nIn=0
	NdocaFidV=0
	NnotdocaFidV=0
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
	#count events inside pairs 
        NdocaIP10FidV=0	
	NnotdocaIP10FidV=0
	Nip10VetoFidV=0
	Nnotip10VetoFidV=0
	NdocaVetoFidV=0
	NnotdocaVetoFidV=0
	Nip10SBTFidV=0
	Nnotip10SBTFidV=0
	NdocaSBTFidV=0
	NnotdocaSBTFidV=0
	NdocaIP30FidV=0
	NnotdocaIP30FidV=0
	Nip30VetoFidV=0
	Nnotip30VetoFidV=0
	Nip30SBTFidV=0
	Nnotip30SBTFidV=0
	NdocaIP250FidV=0
	NnotdocaIP250FidV=0
        Nip250VetoFidV=0	
	Nnotip250VetoFidV=0
        Nip250SBTFidV=0	
	Nnotip250SBTFidV=0

	#cout the number of events outside the Fid Volume
	nOut=0
	NdocaNotFidV=0
	NnotdocaNotFidV=0
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
	Nip250VetoFidV=0
	Nnotip250NotFidV=0

	





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
    		if abs(phi)<1.57: 
    			phi=abs(phi)
    		phi=abs(1.57- abs(phi))
		wPZMThPhi= WeightPZM(Z,P,muWeight,weight,theta,phi)
		totalWeight=prob*muWeight*wPZMThPhi
		numTrack=len(event.goodTracks)
    		if numTrack!=2:continue
		
		for candidate in event.Particles:
			#h['hdocaIP'].Fill(doca,ip, totalWeight)
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
			if distToWall>5*u.cm and isInFiducial(vtx.Z())==True:  
				h['Dist2wall'].Fill(distToWall)
		 	
		 	#if distToWall==0  and isInFiducial(vtx.Z())==True:
			#if distToWall>5*u.cm  and isInFiducial(vtx.Z())==True:
				h['hdocaInside'].Fill(doca,totalWeight)
				nIn+=1*totalWeight

			#n_cut, n_notcut for background having vertex inside the fiducial volume	
				if doca<1:
					h['hNdoca'].Fill(doca,totalWeight)
					NdocaFidV+=1*totalWeight
					if ip<10:
						h['hIP10ifdoca'].Fill(doca,totalWeight)
							
					if ip<30: 
						h['hIP30ifdoca'].Fill(doca,totalWeight)
					if ip<250: 
						h['hIP250ifdoca'].Fill(doca,totalWeight)
					



				else:
					h['hNnotdoca'].Fill(doca,totalWeight)
					NnotdocaFidV+=1*totalWeight	
					


				if ip<10:
					h['hNip10'].Fill(ip,totalWeight)
					Nip10FidV+=1*totalWeight
					if doca<1:
						h['hdocaifIP10'].Fill(doca,totalWeight)
				else:
					h['hNnotip10'].Fill(ip,totalWeight)
					Nnotip10FidV+=1*totalWeight
					
				if ip<250:

					h['hNip250'].Fill(ip,totalWeight)
					Nip250FidV+=1*totalWeight
					if doca<1:
						h['hdocaifIP250'].Fill(doca,totalWeight)
					
				


					
				else:
					h['hNnotip250'].Fill(ip,totalWeight)
					Nnotip250FidV+=1*totalWeight
					
			
				if ip<30:
					h['hNip30'].Fill(ip,totalWeight)
					Nip30FidV+=1*totalWeight
					if doca<1:
						h['hdocaifIP30'].Fill(doca,totalWeight)
					


				else:
					h['hNnotip30'].Fill(ip,totalWeight)
					Nnotip30FidV+=1*totalWeight
					

				if vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					h['hNveto'].Fill(ip,totalWeight)
					NvetoFidV+=1*totalWeight
				else:
					h['hNnotveto'].Fill(ip,totalWeight)
					NnotvetoFidV+=1*totalWeight
					
			
				if vetoSBT == False:
					h['hNsbt'].Fill(doca,totalWeight)
					NsbtFidV+=1*totalWeight
				else:
					h['hNnotsbt'].Fill(doca,totalWeight)
					NnotsbtFidV+=1*totalWeight
				#Pairs Of cuts
				if doca<1 and ip< 10:
					h['hNdocaIP10'].Fill(doca, totalWeight)
					NdocaIP10FidV+=1*totalWeight
			
				if doca>1 and ip>10:
				
					h['hNnotdocaIP10'].Fill(doca, totalWeight)
					NnotdocaIP10FidV+=1*totalWeight
					
				if ip<10 and vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					h['hNip10Veto'].Fill(ip,totalWeight)
					Nip10VetoFidV+=1*totalWeight
				if ip>10 and (vetoSBT == True or  vetoSVT == True and vetoUVT == True and vetoRPC == True) :
				
					h['hNnotip10Veto'].Fill(ip,totalWeight)
					Nnotip10VetoFidV+=1*totalWeight
				if doca<1 and  vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False:
					h['hNdocaVeto'].Fill(doca,totalWeight)
					NdocaVetoFidV+=1*totalWeight

				if doca>1 and (vetoSBT == True or  vetoSVT == True and vetoUVT == True and vetoRPC == True) :  
				
					h['hNnotdocaVeto'].Fill(doca,totalWeight)
					NnotdocaVetoFidV+=1*totalWeight
					
				if vetoSBT==False and ip<10:
					h['hNip10SBT'].Fill(ip,totalWeight)
					Nip10SBTFidV+=1*totalWeight
				if vetoSBT==True and ip>10:
				
					h['hNnotip10SBT'].Fill(ip,totalWeight)	
					Nnotip10SBTFidV+=1*totalWeight
					
				if vetoSBT==False and doca<1:
					h['hNdocaSBT'].Fill(doca,totalWeight)
					NdocaSBTFidV+=1*totalWeight
				if vetoSBT==True and doca>1:
				
					h['hNnotdocaSBT'].Fill(doca,totalWeight)	
					NnotdocaSBTFidV+=1*totalWeight
				
				if doca<1 and ip< 30:
					h['hNdocaIP30'].Fill(doca, totalWeight)
					NdocaIP30FidV+=1*totalWeight
				if doca>1 and ip> 30:
					
					h['hNnotdocaIP30'].Fill(doca, totalWeight)
					NnotdocaIP30FidV+=1*totalWeight
				
				if ip<30 and vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					h['hNip30Veto'].Fill(ip,totalWeight)
					Nip30VetoFidV+=1*totalWeight
				if ip>30 and (vetoSBT == True or  vetoSVT == True or vetoUVT == True or vetoRPC == True) :
				
					h['hNnotip30Veto'].Fill(ip,totalWeight)
					Nnotip30VetoFidV+=1*totalWeight
					

				if vetoSBT==False and ip<30:
					h['hNip30SBT'].Fill(ip,totalWeight)
					Nip30SBTFidV+=1*totalWeight
				if vetoSBT==True and ip>30:
				
					h['hNnotip30SBT'].Fill(ip,totalWeight)	
					Nnotip30SBTFidV+=1*totalWeight
				
				
				if doca<1 and ip< 250:
					h['hNdocaIP250'].Fill(doca, totalWeight)
					NdocaIP250FidV+=1*totalWeight
				
				if doca>1 and ip>250:
					
					h['hNnotdocaIP250'].Fill(doca, totalWeight)
					NnotdocaIP250FidV+=1*totalWeight

				
				if ip<250 and vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					h['hNip250Veto'].Fill(ip,totalWeight)
					Nip250VetoFidV+=1*totalWeight
				if ip>250 and (vetoSBT == True or  vetoSVT == True or vetoUVT == True or vetoRPC == True) :
				
					h['hNnotip250Veto'].Fill(ip,totalWeight)
					Nnotip250VetoFidV+=1*totalWeight
				if vetoSBT==False and ip<250:
					h['hNip250SBT'].Fill(ip,totalWeight)
					Nip250SBTFidV+=1*totalWeight
				
				if vetoSBT==True and ip>250:
					h['hNnotip250SBT'].Fill(ip,totalWeight)	
					Nnotip250SBTFidV+=1*totalWeight
				





			else:

			        hOut['hdocaOutside'].Fill(doca)	
				if doca<1:
					hOut['hNdoca'].Fill(doca,totalWeight)
					NdocaNotFidV+=1*totalWeight

				else:
					hOut['hNnotdoca'].Fill(doca,totalWeight)
					NnotdocaNotFidV+=1*totalWeight	
					


				if ip<10:
					hOut['hNip10'].Fill(ip,totalWeight)
					Nip10NotFidV+=1*totalWeight
				else:
					hOut['hNnotip10'].Fill(ip,totalWeight)
					Nnotip10NotFidV+=1*totalWeight
					
				if ip<250:

					hOut['hNip250'].Fill(ip,totalWeight)
					Nip250NotFidV+=1*totalWeight
					
				else:
					hOut['hNnotip250'].Fill(ip,totalWeight)
					Nnotip250NotFidV+=1*totalWeight
					
			
				if ip<30:
					hOut['hNip30'].Fill(ip,totalWeight)
					Nip30NotFidV+=1*totalWeight
				else:
					hOut['hNnotip30'].Fill(ip,totalWeight)
					Nnotip30NotFidV+=1*totalWeight
					

				if vetoSBT == False and  vetoSVT == False and vetoUVT == False and vetoRPC == False :
					hOut['hNveto'].Fill(ip,totalWeight)
					NvetoNotFidV+=1*totalWeight
				else:
					hOut['hNnotveto'].Fill(ip,totalWeight)
					NnotvetoNotFidV+=1*totalWeight
					
			
				if vetoSBT == False:
					hOut['hNsbt'].Fill(doca,totalWeight)
					NsbtNotFidV+=1*totalWeight
				else:
					hOut['hNnotsbt'].Fill(doca,totalWeight)
					NnotsbtNotFidV+=1*totalWeight

	for key in h:
		h[key].Write()
	'''for key1 in hOut:
		hOut[key1].Write()'''
	
	oIn.Close()
	#oOut.Close()

#nameN=['N_doca','N_ip10', 'N_ip30','N_ip250','N_notdoca','N_notip10', 'N_notip30','N_notip250','
#numN=['NdocaFidV','Nip10FidV', 'Nip30FidV','Nip250FidV','NnotdocaFidV','Nnotip10FidV', 'Nnotip30FidV','Nnotip250FidV',
	dic={'NallIn':nIn,'N_doca': NdocaFidV, 'N_ip10':Nip10FidV, 'N_ip30':Nip30FidV,'N_ip250':Nip250FidV,'Nveto':NvetoFidV,'Nsbt':NsbtFidV,'N_notdoca':NnotdocaFidV, 'N_notip10':Nnotip10FidV, 'N_notip30':Nnotip30VetoFidV,'N_notip250':Nnotip250FidV, 'NnotVeto': NnotvetoFidV,'Nnotsbt': NnotsbtFidV,'Ndocaip10':NdocaIP10FidV,'Nnotdocaip10':NnotdocaIP10FidV,'Nip10veto':Nip10VetoFidV,'Nnotip10veto':Nnotip10VetoFidV,'NdocaVeto':NdocaVetoFidV,'NnotdocaVeto':NnotdocaVetoFidV,'Nip10sbt':Nip10SBTFidV,'Nnotip10sbt':Nnotip10SBTFidV,'Ndocasbt':NdocaSBTFidV,'Nnotdocasbt':NnotdocaSBTFidV,'NdocaIP30':NdocaIP30FidV, 'NnotdocaIP30': NnotdocaIP30FidV,'Nip30Veto':Nip30VetoFidV,'Nnotip30Veto':Nnotip30VetoFidV,'Nip30sbt':Nip30SBTFidV,'Nnotip30sbt':Nnotip10SBTFidV,'NdocaIp250':NdocaIP30FidV,'NnotIp250':NnotdocaIP250FidV,'Nip250Veto':Nip250VetoFidV,'Nnotip250Veto':Nnotip250VetoFidV,'Nip250sbt':Nip250SBTFidV,'Nnotip250SBT':Nnotip250SBTFidV}
	print dic.items()



if __name__ == '__main__':
    IndependenceCheck()
	



