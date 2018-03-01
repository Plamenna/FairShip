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


