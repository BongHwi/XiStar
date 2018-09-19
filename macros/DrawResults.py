#Last Working backup


# DrawResults.py
# Draw histograms with specific criteria and fit
# 
# Author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
# 
# How to use:
# python DrawResults.py AnalysisResults.root 1
# parameters:
#   - AnalysisResults.root: Input root file
#   - 1: 1 or 0, 1: save all outputs 0: save only merged one
#
#from ROOT import *
import ROOT
from math import *
import numpy as np
import datetime
import sys

currenttime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
PrintOption = 1 #0 for no print, 1 for print

if(PrintOption == 1): print("Current time: %s"%currenttime)

if(len(sys.argv)>=2): 
	Inputfile01 = sys.argv[1]
	if(PrintOption == 1): print("First input file: %s"%(Inputfile01))

	if(len(sys.argv)>=3): 
		Inputfile02 = sys.argv[2]
		if(PrintOption == 1): print("Second input file: %s"%(Inputfile02))

#===================
# Load the input root file 01
if(PrintOption == 1): print("====================")
if(PrintOption == 1): print("Load file: "+Inputfile01)
f1 = ROOT.TFile(Inputfile01)
mydir1 = f1.Get("PWGLF.outputXiStarAnalysis.root;1")
mydir1.ls()
mylist1 = mydir1.Get("MyList;1")
if(mylist1.FindObject("fMCinputTotalXiStar3").GetEntries()>0): 
	f_MC = f1
	mylist_mc = mylist1
	if(PrintOption == 1): print("First input: MC")
else:
	f_data = f1
	mylist_data = mylist1
	if(PrintOption == 1): print("First input: data")
#===================
# Load the input root file 02
if(len(sys.argv)>=3):
	if(PrintOption == 1): print("====================")
	if(PrintOption == 1): print("Load file: "+Inputfile02)
	f2 = ROOT.TFile(Inputfile02)
	mydir2 = f2.Get("PWGLF.outputXiStarAnalysis.root;1")
	mydir2.ls()
	mylist2 = mydir2.Get("MyList;1")
	if(mylist2.FindObject("fMCinputTotalXiStar3").GetEntries()>0): 
		f_MC = f2
		mylist_mc = mylist2
		if(PrintOption == 1): print("Second input: MC")
	else:
		f_data = f2
		mylist_data = mylist2
		if(PrintOption == 1): print("Second input: data")

def myLevyPtFunc(x, par):
	lMass = 1.5318 #Xi* mass
	ldNdy  = par[0] #dN/dy
	l2pi   = 2*ROOT.TMath.Pi() #2pi
	lTemp = par[1] #Temperature
	lPower = par[2] #power=n
	lBigCoef = ((lPower-1)*(lPower-2)) / (l2pi*lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)))
	lInPower = 1 + (ROOT.TMath.Sqrt(x[0]*x[0]+lMass*lMass)-lMass) / (lPower*lTemp)
	return l2pi * ldNdy * x[0] * lBigCoef * ROOT.TMath.Power(lInPower,(-1)*lPower)

def DrawResults():
	#====================
	# Inputs
	outputFileName = "Analysis_Results_Xi1530_%s.root"%currenttime
	
	pTRange = [0, 0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6, 8.8,15]
	#pTRange = [0,30]
	MultiplicityRage = [0,5,15,30,50,100]
	#MultiplicityRage = [0,15,50,100]
	MultiplicitySet = 0 #0 of MB, 1 for Multi add
	
	massRebin = 1
	AxisRange = [1.48,1.59]
	#NormalizeRange = [1.49,1.51] #for LHC15f
	NormalizeRange = [1.56,1.58] #for LHC16k
	NormalizeRange_L = [1.49,1.51]
	NormalizeRange_R = [1.56,1.58]
	fitRange = [1.52,1.55]
	fitintrange = [1.52,1.55]
	hNames = ["fXiMinusPiPlus_0", "fXiPlusPiMinus_0"] #fXiMinusPiPlus_0: Xi(1530), fXiPlusPiMinus_0: Anti-Xi(1530)
	hNamesXi = ["",""] # Xi histogram
	#hNames_mix = ["fXiMinusPiPlusbkg_0", "fXiPlusPiMinusbkg_0"] #fXiMinusPiPlus_0: Xi(1530), fXiPlusPiMinus_0: Anti-Xi(1530)
	hNames_mix = ["fXiMinusPiMinus_0", "fXiPlusPiPlus_0"] # Like sign #fXiMinusPiPlus_0: Xi(1530), fXiPlusPiMinus_0: Anti-Xi(1530)
	hNames_MCinput = ["fMCinputTotalXiStar3", "fMCinputTotalXiStarbar3"] #generated MC
	hNames_MCinput_prePV = ["fMCinputTotalXiStar1", "fMCinputTotalXiStarbar1"] #generated MC pre PV cut
	hNames_Mcrecon = ["fMCrecXiMinusPiPlus_0", "fMCrecXiPlusPiMinus_0"] #recon MC
	QAhisto = ["fMultDist_pp", "hEventSelecInfo", "fCutEvents","TPCPID","hdEdxProton","hdEdxPion1","hdEdxPion2","hdEdxProtonAfter","hdEdxPion1After","hdEdxPion2After"] #for QA
	Fit_type = 0 #0: sigma fit(MB study) 1: width fit
	FitSigma = 0.0020 #0.0020 for default value
	Save = 0 #0: skip saving pdf   images, 1: save pdf   images
	SaveType = "png" #png, pdf
	#====================
	#====================
	# Default lists
	QAh = []
	QAh_mc = []

	pTCenter=[]
	pTRange_err = []
	MultCenter=[]
	MultRange_err = []

	Xi = []
	XipT = []

	XiStar = []
	XiStar_mix = []
	XiStarpT = []
	XiStarpT_mix = []

	XiStarMCinput = []
	XiStarMCinput_prePV = []
	XiStarMCrecon = []
	XiStarMCinputpT = []
	XiStarMCinputpT_prePV = []
	XiStarMCreconpT = []

	XiStarMCinputpT_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	XiStarMCinputpT_prePV_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	XiStarMCreconpT_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]

	MCEfficiency = []
	MCEfficiency_7TeV = [0,0.006,0.0237,0.0433,0.0705,0.1,0.142,0.158,0.157,0,0]
	MCEfficiency_13TeV_LHC17 = [0, 0.003375927, 0.01433215, 0.03185342, 0.05328012, 0.0841077, 0.114816, 0.1273939, 0.1275806, 0.1091881, 0.0629708]
	MCEfficiency_13TeV_LHC16 = [0, 0.003411181, 0.01415889, 0.03134302, 0.05217555, 0.08471585, 0.1153396, 0.1278187, 0.1275664, 0.1092717, 0.06316746]
	MCEfficiency_13TeV_LHC15 = [0.0001308913, 0.004535666, 0.01680443, 0.03632307, 0.05875776, 0.091464, 0.1245708, 0.1364504, 0.1363351, 0.1163355, 0.06741513]
	MCEfficiency_Error = []
	MCEfficiency_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	MCEfficiency_Error_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]

	integratedNumberforNormalizepT = []
	integratedNumberforNormalizepT_mix = []
	integratedNumberforNormalizepT_mult = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]
	integratedNumberforNormalizepT_mult_mix = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]

	fitChi2ndf=[]
	fitParameter_mean = []
	fitParameter_mean_err = []
	fitParameter_sigma = []
	fitParameter_sigma_err = []

	fitChi2ndf_mult=[[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	fitParameter_mean_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	fitParameter_mean_err_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	fitParameter_sigma_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	fitParameter_sigma_err_mult = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]
	
	XiStarpT_Multi = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]
	XiStarpT_Multi_mix = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]
	
	RawYield = []
	RawYield_err = []
	RawYield_Multi = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]
	RawYield_Multi_err = [[0 for col in range(len(MultiplicityRage)-1)] for row in range(len(pTRange)-1)]
	
	h_Xispectrum_mult = []

	for i in range(0,len(pTRange)-1):
		pTCenter.append(pTRange[i]+(pTRange[i+1]-pTRange[i])/2)
		pTRange_err.append((pTRange[i+1]-pTRange[i])/2)
	for i in range(0,len(MultiplicityRage)-1):
		MultCenter.append(MultiplicityRage[i]+(MultiplicityRage[i+1]-MultiplicityRage[i])/2)
		MultRange_err.append((MultiplicityRage[i+1]-MultiplicityRage[i])/2)

	#===================
	if(PrintOption == 1): print("====================")
	if(PrintOption == 1): print("make file: "+outputFileName)
	fo = ROOT.TFile(outputFileName,"RECREATE");
	fo.cd()
	if(PrintOption == 1): print("====================")
	#===================
	# QA Plots
	c0 = ROOT.TCanvas('can','canvas',1280,720)
	c0.SetBorderMode(1)					 
	c0.cd() 
	try:
		mylist_data
	except NameError:
		print("None data QA for data")
	else:
		for i in range(0,len(QAhisto)):
			QAh.append(mylist_data.FindObject(QAhisto[i]))
			if (QAh[i].GetName() == "fMultDist_pp"):
				QAh[i].SetAxisRange(0,110,"X")
				QAh[i].SetStats(0)
				QAh[i].GetXaxis().SetTitle("Multiplicity percentile (%)")
			if (QAh[i].GetName() in ["TPCPID","hdEdxProton","hdEdxPion1","hdEdxPion2","hdEdxProtonAfter","hdEdxPion1After","hdEdxPion2After"]):
				ROOT.gPad.SetLogz(1)
				QAh[i].SetStats(0)
				QAh[i].SetTitle("")
				QAh[i].SetAxisRange(0.3,20,"X")
				QAh[i].GetXaxis().SetTitle("p (GeV/c)")
				c0.SetLogx()
				QAh[i].GetYaxis().SetTitle("dE/dx (arb. units)")
				ROOT.gStyle.SetPalette(1)
			QAh[i].Write("%s"%QAhisto[i])
			QAh[i].Draw("colz")
			c0.SaveAs("figs/QA/QA_%s.%s"%(QAhisto[i],SaveType))
	try:
		mylist_mc
	except NameError:
		print("None data QA for MC")
	else:
		for i in range(0,len(QAhisto)):
			QAh_mc.append(mylist_mc.FindObject(QAhisto[i]))
			if (QAh_mc[i].GetName() == "fMultDist_pp"):
				QAh_mc[i].SetAxisRange(0,110,"X")
				QAh_mc[i].SetStats(0)
				QAh_mc[i].GetXaxis().SetTitle("Multiplicity percentile (%)")
			QAh_mc[i].Write("%s_mc"%QAhisto[i])
			QAh_mc[i].Draw("colz")
			c0.SaveAs("figs/QA/QA_%s_mc.%s"%(QAhisto[i],SaveType))
		Events_PhysSel=mylist_mc.FindObject("fMultDist1")
		Events_PhysSel.Write("fMultDist1")
		Events_PhysSel.Draw()
		c0.SaveAs("figs/QA/QA_Events_PhysSel_mc.%s"%SaveType)

	#===================
	#ROOT.gStyle.SetOptStat(1)
	try:
		mylist_mc
	except NameError:
		print("No MC input: skip")
	else: #MC Case.
		i = 0
		for name in hNames_MCinput:
			if(PrintOption == 1): print("Load Histogram: "+name)
			XiStarMCinput.append(mylist_mc.FindObject(name))
			XiStarMCinput[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			XiStarMCinput[i].GetYaxis().SetTitle("Multiplicity percentile (%)")
			XiStarMCinput[i].GetZaxis().SetTitle("Mass (GeV/c^{2})")
			XiStarMCinput[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		if(PrintOption == 1): print("====================")

		i = 0
		for name in hNames_MCinput_prePV:
			if(PrintOption == 1): print("Load Histogram: "+name)
			XiStarMCinput_prePV.append(mylist_mc.FindObject(name))
			XiStarMCinput_prePV[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			XiStarMCinput_prePV[i].GetYaxis().SetTitle("Multiplicity percentile (%)")
			XiStarMCinput_prePV[i].GetZaxis().SetTitle("Mass (GeV/c^{2})")
			XiStarMCinput_prePV[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		if(PrintOption == 1): print("====================")
		

		i = 0
		for name in hNames_Mcrecon:
			if(PrintOption == 1): print("Load Histogram: "+name)
			XiStarMCrecon.append(mylist_mc.FindObject(name))
			XiStarMCrecon[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			XiStarMCrecon[i].GetYaxis().SetTitle("Multiplicity percentile (%)")
			XiStarMCrecon[i].GetZaxis().SetTitle("Mass (GeV/c^{2})")
			XiStarMCrecon[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		if(PrintOption == 1): print("====================")

		# Add them each other
		XiStarMCinput[0].Add(XiStarMCinput[1])
		XiStarMCinput_prePV[0].Add(XiStarMCinput_prePV[1])
		XiStarMCrecon[0].Add(XiStarMCrecon[1])

		# Generated Xi(1530)
		temp1 = ROOT.TH1F()
		temp1 = XiStarMCinput[0].Project3D("X")
		temp1.GetXaxis().SetTitle("p_{T} (GeV/c)")
		temp1.SetTitle("Generated Xi(1530)")
		temp1.SetAxisRange(1e-1, 1e6,"Y")
		temp1.Write("Generated Xi1530")
		if(PrintOption == 1): print("Generated Integral: %f"%XiStarMCinput[0].Integral(XiStarMCinput[0].GetXaxis().FindBin(pTRange[0]),XiStarMCinput[0].GetXaxis().FindBin(pTRange[-1]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[0]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[-1]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[-1])))

		# Generated Xi(1530) pre PV cut
		temp3 = ROOT.TH1F()
		temp3 = XiStarMCinput_prePV[0].Project3D("X")
		temp3.GetXaxis().SetTitle("p_{T} (GeV/c)")
		temp3.SetTitle("Generated Xi(1530) pre PV")
		temp3.SetAxisRange(1e-1, 1e6,"Y")
		temp3.Write("Generated Xi1530_prePV")
		if(PrintOption == 1): print("Generated Integral prePV: %f"%temp3.Integral(temp3.FindBin(0),temp3.FindBin(8)))

		# Reconstructed Xi(1530)
		temp2 = ROOT.TH1F()
		temp2 = XiStarMCrecon[0].Project3D("X")
		temp2.GetXaxis().SetTitle("p_{T} (GeV/c)")
		temp2.SetTitle("Reconstructed Xi(1530)")
		temp2.Write("Reconstructed Xi1530")
		if(PrintOption == 1): print("Reconstructed Integral: %f"%XiStarMCrecon[0].Integral(XiStarMCrecon[0].GetXaxis().FindBin(pTRange[0]),XiStarMCrecon[0].GetXaxis().FindBin(pTRange[-1]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[0]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[-1]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[-1])))
		
		### Draw MC QA Histos
		c1 = ROOT.TCanvas('can1','canvas',1280,720)
		c1.SetBorderMode(1)					 
		c1.cd() 
		c1.SetLogy()

		temp1.Draw("E")
		c1.Write("Xi1530_MCQA_NumberWithpT_generated")
		c1.SaveAs("figs/mc/Xi1530_MCQA_NumberWithpT_generated.%s"%SaveType)

		temp2.SetTitle("Reconstructed Xi(1530)")
		temp2.Draw("E")
		c1.Write("Xi1530_MCQA_NumberWithpT_reconstructed")
		c1.SaveAs("figs/mc/Xi1530_MCQA_NumberWithpT_reconstructed.%s"%SaveType)

		temp2.SetLineColor(1)
		temp1.SetTitle("Generated/Reconstructed Xi(1530)")
		temp1.Draw("E")
		temp2.Draw("E SAME")
		c1.Write("Xi1530_MCQA_NumberWithpT_same")
		c1.SaveAs("figs/mc/Xi1530_MCQA_NumberWithpT_same.%s"%SaveType)

		# Get Entries with given pT bin for generated MC
		for i in range(0,len(pTRange)-1):
			if(PrintOption == 1): print("=========MC=========")
			if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
			XiStarMCinputpT.append(XiStarMCinput[0].Integral(XiStarMCinput[0].GetXaxis().FindBin(pTRange[i]),XiStarMCinput[0].GetXaxis().FindBin(pTRange[i+1]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[0]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[-1]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[-1])))
			if(PrintOption == 1): print("Normalizing factor from generated: %f"%XiStarMCinputpT[i])
			XiStarMCreconpT.append(XiStarMCrecon[0].Integral(XiStarMCrecon[0].GetXaxis().FindBin(pTRange[i]),XiStarMCrecon[0].GetXaxis().FindBin(pTRange[i+1]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[0]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[-1]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[-1])))
			if(PrintOption == 1): print("Normalizing factor from reconstructed: %f"%XiStarMCreconpT[i])
			MCEfficiency.append(XiStarMCreconpT[i]/XiStarMCinputpT[i])
			MCEfficiency_Error.append(sqrt(pow(sqrt(XiStarMCreconpT[i])/XiStarMCinputpT[i],2) + pow(sqrt(XiStarMCinputpT[i])*XiStarMCreconpT[i]/pow(XiStarMCinputpT[i],2),2)))
			if(PrintOption == 1): print("Efficiency: %f"%MCEfficiency[i])
			if(PrintOption == 1): print("Efficiency_Error: %f"%MCEfficiency_Error[i])
			for j in range(0,len(MultiplicityRage)-1):
				if(PrintOption == 1): print("=========Mult %d to %d========="%(MultiplicityRage[j],MultiplicityRage[j+1]))
				if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
				XiStarMCinputpT_mult[i][j] = XiStarMCinput[0].Integral(XiStarMCinput[0].GetXaxis().FindBin(pTRange[i]),XiStarMCinput[0].GetXaxis().FindBin(pTRange[i+1]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[j]),XiStarMCinput[0].GetYaxis().FindBin(MultiplicityRage[j+1]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCinput[0].GetZaxis().FindBin(fitintrange[-1]))
				if(PrintOption == 1): print("Normalizing factor from generated: %f"%XiStarMCinputpT_mult[i][j])
				XiStarMCreconpT_mult[i][j] = XiStarMCrecon[0].Integral(XiStarMCrecon[0].GetXaxis().FindBin(pTRange[i]),XiStarMCrecon[0].GetXaxis().FindBin(pTRange[i+1]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[j]),XiStarMCrecon[0].GetYaxis().FindBin(MultiplicityRage[j+1]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[0]),XiStarMCrecon[0].GetZaxis().FindBin(fitintrange[-1]))
				if(PrintOption == 1): print("Normalizing factor from reconstructed: %f"%XiStarMCreconpT_mult[i][j])
				MCEfficiency_mult[i][j] = XiStarMCreconpT_mult[i][j]/XiStarMCinputpT_mult[i][j]
				MCEfficiency_Error_mult = sqrt(pow(sqrt(XiStarMCreconpT_mult[i][j])/XiStarMCinputpT_mult[i][j],2) + pow(sqrt(XiStarMCinputpT_mult[i][j])*XiStarMCreconpT_mult[i][j]/pow(XiStarMCinputpT_mult[i][j],2),2))
				if(PrintOption == 1): print("Efficiency: %f"%MCEfficiency[i])
				if(PrintOption == 1): print("Efficiency_Error: %f"%MCEfficiency_Error[i])
			


		gr_Eff = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(MCEfficiency, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(MCEfficiency_Error, 'd'))
		gr_Eff.SetMarkerStyle(20)
		gr_Eff.SetMinimum(0)
		gr_Eff.SetMaximum(0.4)
		gr_Eff.SetTitle("")
		gr_Eff.GetXaxis().SetTitle("p_{T} (GeV/c)")
		gr_Eff.GetYaxis().SetTitle("Acceptance x Efficiency x BR")
		gr_Eff.Draw("AP")

		gr_Eff_7TeV = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(MCEfficiency_7TeV, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(MCEfficiency_Error, 'd'))
		gr_Eff_7TeV.SetMarkerStyle(20)
		gr_Eff_7TeV.SetMarkerColor(2)
		gr_Eff_7TeV.SetMinimum(0)
		gr_Eff_7TeV.SetMaximum(0.4)
		gr_Eff_7TeV.SetTitle("")
		gr_Eff_7TeV.Draw("P")

		gr_Eff_13TeV_17 = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(MCEfficiency_13TeV_LHC17, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(MCEfficiency_Error, 'd'))
		gr_Eff_13TeV_17.SetMarkerStyle(20)
		gr_Eff_13TeV_17.SetMarkerColor(3)
		gr_Eff_13TeV_17.SetMinimum(0)
		gr_Eff_13TeV_17.SetMaximum(0.4)
		gr_Eff_13TeV_17.SetTitle("")
		gr_Eff_13TeV_17.Draw("P")

		gr_Eff_13TeV_16 = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(MCEfficiency_13TeV_LHC16, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(MCEfficiency_Error, 'd'))
		gr_Eff_13TeV_16.SetMarkerStyle(20)
		gr_Eff_13TeV_16.SetMarkerColor(4)
		gr_Eff_13TeV_16.SetMinimum(0)
		gr_Eff_13TeV_16.SetMaximum(0.4)
		gr_Eff_13TeV_16.SetTitle("")
		#gr_Eff_13TeV_16.Draw("P")

		leg = ROOT.TLegend(.65,.20,.80,.40)
		leg.SetBorderSize(0)
		leg.SetFillColor(0)
		leg.SetFillStyle(0)
		leg.SetTextFont(42)
		leg.SetTextSize(0.035)
		leg.AddEntry(gr_Eff,"13 TeV (LHC18c6a)","P")
		leg.AddEntry(gr_Eff_13TeV_17,"13 TeV (Old, wrong)","P")
		#leg.AddEntry(gr_Eff_13TeV_16,"13 TeV (LHC16)","P")
		leg.AddEntry(gr_Eff_7TeV,"7 TeV (Old)","P")
		leg.Draw()

		c1.Write("Xi1530_MCQA_Efficiency")
		c1.SaveAs("figs/mc/Xi1530_MCQA_Efficiency.%s"%SaveType)
		
	#---------------------------------------------------------------------------------------
	try:
		mylist_data
	except NameError:
		print("No data input: skip")
	else: # Data Case
		i = 0
		for name in hNames:
			if(PrintOption == 1): print("Load Histogram: "+name)
			XiStar.append(mylist_data.FindObject(name))
			XiStar[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			XiStar[i].GetYaxis().SetTitle("Multiplicity percentile (%)")
			XiStar[i].GetZaxis().SetTitle("Mass (GeV/c^{2})")
			XiStar[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			XiStar[i].RebinZ(massRebin)
			i = i+1
		if(PrintOption == 1): print("====================")

		i = 0
		for name in hNames_mix:
			if(PrintOption == 1): print("Load Histogram: "+name)
			XiStar_mix.append(mylist_data.FindObject(name))
			XiStar_mix[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			XiStar_mix[i].GetYaxis().SetTitle("Multiplicity percentile (%)")
			XiStar_mix[i].GetZaxis().SetTitle("Mass (GeV/c^{2})")
			XiStar_mix[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			XiStar_mix[i].RebinZ(massRebin)
			i = i+1
		if(PrintOption == 1): print("====================")	

		if(PrintOption == 1): print("Load Histogram: "+"fXi_0")
		Xi.append(mylist_data.FindObject("fXi_0"))
		Xi[0].GetXaxis().SetTitle("p_{T} (GeV/c)")
		Xi[0].GetYaxis().SetTitle("Multiplicity percentile (%)")
		Xi[0].GetZaxis().SetTitle("Mass (GeV/c^{2})")
		if(PrintOption == 1): print("====================")

		if(PrintOption == 1): print("Load Histogram: "+"fXibar_0")
		Xi.append(mylist_data.FindObject("fXibar_0"))
		Xi[1].GetXaxis().SetTitle("p_{T} (GeV/c)")
		Xi[1].GetYaxis().SetTitle("Multiplicity percentile (%)")
		Xi[1].GetZaxis().SetTitle("Mass (GeV/c^{2})")
		if(PrintOption == 1): print("====================")

		Xi[0].Add(Xi[1])
		Xi[0].Write("Original_Xi_pT_Multi_Mass") #Original 3D histogram
		
		# QA Plots
		c1 = ROOT.TCanvas('can1','canvas',1280,720)
		c1.SetBorderMode(1)					 
		c1.cd() 

		# Add Antiparticle
		if(PrintOption == 1): print("Entry of Xi(1530)       : %d"%XiStar[0].GetEntries())
		if(PrintOption == 1): print("Entry of anti-Xi(1530)  : %d"%XiStar[1].GetEntries())
		XiStar[0].Add(XiStar[1])
		XiStar_mix[0].Add(XiStar_mix[1])
		if(PrintOption == 1): print("Entry of merged Xi(1530): %d"%XiStar[0].GetEntries())

		XiStar[0].Write("Original_Xi1530_pT_Multi_Mass") #Original 3D histogram
		
		#Slice Histogram with given pT bin
		for i in range(0,len(pTRange)-1):
			if(PrintOption == 1): print("========Xi1530==========")
			if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
			XiStarpT.append(XiStar[0].Clone())
			XiStarpT[i].SetAxisRange(pTRange[i],pTRange[i+1],"X")
			XiStarpT[i].SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			XiStarpT[i].Write("Xi1530_pT_%.1f_to_%.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			temp = ROOT.TH1F()
			temp = XiStarpT[i].Project3D("Z")
			temp.SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			#Normalization Factor
			nbin1 = temp.GetXaxis().FindBin(NormalizeRange_L[0])
			nbin2 = temp.GetXaxis().FindBin(NormalizeRange_L[1])
			nbin3 = temp.GetXaxis().FindBin(NormalizeRange_R[0])
			nbin4 = temp.GetXaxis().FindBin(NormalizeRange_R[1])
			#integratedNumberforNormalizepT.append(temp.Integral(nbin1,nbin2)+temp.Integral(nbin3,nbin4)) #BOTH
			#integratedNumberforNormalizepT.append(temp.Integral(nbin3,nbin4)) #R
			integratedNumberforNormalizepT.append(temp.Integral(nbin1,nbin2)) #L
			if (integratedNumberforNormalizepT[i] == 0): integratedNumberforNormalizepT[i] = 1
			if(PrintOption == 1): print("Normalizing factor: %f"%integratedNumberforNormalizepT[i])
			temp.Draw('E')
			temp.Write("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			if(Save): c1.SaveAs("figs/pTBin/sig/Xi1530_pT_%.1f_to_%.1f.%s"%(pTRange[i],pTRange[i+1],SaveType))
			XipT.append(Xi[0].Clone())
			XipT[i].SetAxisRange(pTRange[i],pTRange[i+1],"X")
			XipT[i].SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			XipT[i].Write("Xi1530_pT_%.1f_to_%.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			temp2 = ROOT.TH1F()
			temp2 = XipT[i].Project3D("Z")
			temp2.Draw('E')
			temp2.Write("Xi_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			if(Save): c1.SaveAs("figs/QA/Xi/pTBin/Xi_pT_%.1f_to_%.1f.%s"%(pTRange[i],pTRange[i+1],SaveType))
		
		#Slice Histogram with given pT bin for event mix
		for i in range(0,len(pTRange)-1):
			if(PrintOption == 1): print("=====Event Mixing=====")
			if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
			XiStarpT_mix.append(XiStar_mix[0].Clone())
			XiStarpT_mix[i].SetAxisRange(pTRange[i],pTRange[i+1],"X")
			XiStarpT_mix[i].SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			XiStarpT_mix[i].Write("Xi1530_pT_%.1f_to_%.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			temp = ROOT.TH1F()
			temp = XiStarpT_mix[i].Project3D("Z")
			temp.SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			#Normalization Factor
			nbin1 = temp.GetXaxis().FindBin(NormalizeRange_L[0])
			nbin2 = temp.GetXaxis().FindBin(NormalizeRange_L[1])
			nbin3 = temp.GetXaxis().FindBin(NormalizeRange_R[0])
			nbin4 = temp.GetXaxis().FindBin(NormalizeRange_R[1])
			#integratedNumberforNormalizepT_mix.append(temp.Integral(nbin1,nbin2)+temp.Integral(nbin3,nbin4)) #BOTH
			#integratedNumberforNormalizepT_mix.append(temp.Integral(nbin3,nbin4)) #R
			integratedNumberforNormalizepT_mix.append(temp.Integral(nbin1,nbin2)) #L
			if (integratedNumberforNormalizepT_mix[i] == 0): integratedNumberforNormalizepT_mix[i] = 1
			if(PrintOption == 1): print("Normalizing factor: %f"%integratedNumberforNormalizepT_mix[i])
			temp.Draw('E')
			temp.Write("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			if(Save): c1.SaveAs("figs/pTBin/bkg/Xi1530_pT_%.1f_to_%.1f_EventMix.%s"%(pTRange[i],pTRange[i+1],SaveType))
		
		# Draw Together(pT)
		for i in range(0,len(pTRange)-1):
			### First, Draw together.
			c1.cd()
			temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			temp1.SetLineColor(0)
			temp1.Draw("E")
			temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			temp2.Scale(integratedNumberforNormalizepT[i]/integratedNumberforNormalizepT_mix[i])
			temp2.SetLineColor(1)
			temp2.Draw("E, same")
			temp2.SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			c1.Write("Xi1530_pT_%.1f_to_%.1f_with_bkg"%(pTRange[i],pTRange[i+1]))
			c1.SaveAs("figs/pTBin/same/Xi1530_pT_%.1f_to_%.1f_with_bkg.%s"%(pTRange[i],pTRange[i+1],SaveType))

			### Second, Subtract bkg.
			temp1.Add(temp2,-1)
			temp1.Draw('E')
			temp1.Write("Xi1530_pT_%.1f_to_%.1f_with_bkg_substraction"%(pTRange[i],pTRange[i+1]))
			temp1.SaveAs("figs/pTBin/sub/Xi1530_pT_%.1f_to_%.1f_with_bkg_substraction.%s"%(pTRange[i],pTRange[i+1],SaveType))

			### Third, Fit it.
			# Reference1: $ROOTSYS/tutorials/roofit/rf204_extrangefit.C
			# Reference2: https://github.com/clelange/roofit/blob/master/rf204_extrangefit.py
			if(PrintOption == 1): print("==================================================================================%.1f to %.1f case =================================================================================="%(pTRange[i],pTRange[i+1]))
			x = ROOT.RooRealVar("x", "x", AxisRange[0], AxisRange[1])
			x.setRange("IntegralRange", fitintrange[0], fitintrange[1])
			mean1 = ROOT.RooRealVar("mean1", "mean of gaussians", 1.530,1.525,1.545)
			if(Fit_type == 0): #sigma fit, MB study
				sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 0.003,0.000001,0.01)
				width = ROOT.RooRealVar("width", "width of voigtian", 0.0091)
			if(Fit_type == 1): #width fit
				sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", FitSigma)
				width = ROOT.RooRealVar("width", "width of voigtian", 0.0091,0,0.02)
			signal = ROOT.RooVoigtian("signal", "Signal component 1", x, mean1, width, sigma1)

			a0 = ROOT.RooRealVar("a0","a0",1) ;
			a1 = ROOT.RooRealVar("a1","a1",0,-1,1) ;
			a2 = ROOT.RooRealVar("a2","a2",1,0,10) ;
			background = ROOT.RooPolynomial("p2","p2",x,ROOT.RooArgList(a0,a1,a2),0) ;
			
			# --- Construct signal+background PDF ---
			nsig = ROOT.RooRealVar("nsig","#signal events",2000,0.,10e5) ;
			#nbkg = ROOT.RooRealVar("nbkg","#background events",800,0.,1000) ;
			#model = ROOT.RooAddPdf("model","g+a",ROOT.RooArgList(signal,background),ROOT.RooArgList(nsig,nbkg))
			model = ROOT.RooAddPdf("model","g+a",ROOT.RooArgList(signal),ROOT.RooArgList(nsig))

			data = ROOT.RooDataHist("RooDataHist","dataset with Gauss  Fit",ROOT.RooArgList(x),temp1)

			c1.cd()
			ROOT.gPad.SetTickx(2)
			xframe = x.frame()
			#signal.fitTo(data,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
			model.fitTo(data,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
			#signal.fitTo(data)
			#model.fitTo(data)
 			#model.fitTo(data,ROOT.RooFit.Range(1.52,1.55))

			data.plotOn(xframe)
			#signal.plotOn(xframe,ROOT.RooFit.LineStyle(ROOT.kDashed))
			model.plotOn(xframe,ROOT.RooFit.Range(fitRange[0],fitRange[1]),ROOT.RooFit.LineStyle(ROOT.kDashed))
			#model.plotOn(xframe)
			#model.plotOn(xframe,ROOT.RooFit.Components("p2"),ROOT.RooFit.LineStyle(ROOT.kDashed))
			#model.plotOn(xframe,ROOT.RooFit.Components("p2"),RooFit.LineColor(kGreen))
			#sig.paramOn(xframe, ROOT.RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
			xframe.SetMarkerStyle(1)
			xframe.Draw("E")		
			xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)"%(pTRange[i],pTRange[i+1]))
			c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting(dash).%s"%(pTRange[i],pTRange[i+1],SaveType))
			c1.Clear()
			data.plotOn(xframe)
			model.plotOn(xframe,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
			#model.plotOn(xframe,ROOT.RooFit.LineStyle(ROOT.kDashed))
			#signal.plotOn(xframe,ROOT.RooFit.Components("background"),ROOT.RooFit.LineStyle(ROOT.kDashed))
			#model.plotOn(xframe,ROOT.RooFit.Components("p2"),ROOT.RooFit.LineStyle(ROOT.kDashed))
			#model.plotOn(xframe,ROOT.RooFit.Components("p2"),RooFit.LineColor(kGreen))
			#sig.paramOn(xframe, ROOT.RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
			xframe.SetMarkerStyle(2)
			xframe.Draw("E")		
			xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting"%(pTRange[i],pTRange[i+1]))
			c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting.%s"%(pTRange[i],pTRange[i+1],SaveType))

			#chi2/ndf
			fitChi2ndf.append(xframe.chiSquare())
			if(PrintOption == 1): print("chi2: %.4f"%xframe.chiSquare())  # NEED TO ADD NDF INFO!!
			

			# --- Extract fit parameters ---
			#if(PrintOption == 1): print("Fit Result: Mean: %.3f, Sigma: %.5f"%(mean1.getValV(),sigma1.getValV()))
			fitParameter_mean.append(mean1.getValV())
			fitParameter_mean_err.append(mean1.getError())
			if(Fit_type == 0): #sigma fit, MB study
				fitParameter_sigma.append(sigma1.getValV())
				fitParameter_sigma_err.append(sigma1.getError())
			if(Fit_type == 1):
				fitParameter_sigma.append(width.getValV())
				fitParameter_sigma_err.append(width.getError())
		
			# Get the raw yield
			#fitIntegral = signal.createIntegral(mesArgSet,ROOT.RooFit.NormSet(mesArgSet),ROOT.RooFit.Range("IntegralRange"))
			fitIntegral = model.createIntegral(ROOT.RooArgSet(x),ROOT.RooArgSet(x),"IntegralRange")
			RawYield.append(nsig.getValV()*fitIntegral.getVal())
			fitIntegral_err = fitIntegral.getValV()*nsig.getError()/nsig.getValV()
			RawYield_err.append(fitIntegral_err)

			

		if(MultiplicitySet=="1"):
			#==========================	MULTIPLICITY============================================
			#Slice Histogram with given Multiplicity
			for i in range(0,len(pTRange)-1):
				if(PrintOption == 1): print("====================")
				if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))	
				for j in range(0,len(MultiplicityRage)-1):
					if(PrintOption == 1): print("--Multiplicity Rage: "+str(MultiplicityRage[j])+" to "+str(MultiplicityRage[j+1]))
					XiStarpT_Multi[i][j] = XiStarpT[i].Clone()
					XiStarpT_Multi[i][j].SetAxisRange(MultiplicityRage[j],MultiplicityRage[j+1],"Y")
					temp = ROOT.TH1F()
					temp = XiStarpT_Multi[i][j].Project3D("Z")
					#Normalization Factor
					nbin1 = temp.GetXaxis().FindBin(NormalizeRange_L[0])
					nbin2 = temp.GetXaxis().FindBin(NormalizeRange_L[1])
					nbin3 = temp.GetXaxis().FindBin(NormalizeRange_R[0])
					nbin4 = temp.GetXaxis().FindBin(NormalizeRange_R[1])
					integratedNumberforNormalizepT_mult[i][j] = temp.Integral(nbin1,nbin2)+temp.Integral(nbin3,nbin4) #BOTH
					#integratedNumberforNormalizepT.append(temp.Integral(nbin3,nbin4)) #R
					#integratedNumberforNormalizepT.append(temp.Integral(nbin1,nbin2)) #L
					if (integratedNumberforNormalizepT_mult[i][j] == 0): integratedNumberforNormalizepT_mult[i][j] = 1
					if(PrintOption == 1): print("Normalizing factor: %f"%integratedNumberforNormalizepT_mult[i][j])	
					temp.SetTitle("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp.Draw('E')
					temp.Write("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					if(Save): c1.SaveAs("figs/Multibin/Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

			#Slice Histogram with given pT bin for event mix with given Multiplicity
			for i in range(0,len(pTRange)-1):
				if(PrintOption == 1): print("====================")
				if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))	
				for j in range(0,len(MultiplicityRage)-1):
					if(PrintOption == 1): print("=====Event Mixing=====")
					if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
					XiStarpT_Multi_mix[i][j] = XiStarpT_mix[i].Clone()
					XiStarpT_Multi_mix[i][j].SetAxisRange(pTRange[i],pTRange[i+1],"X")
					XiStarpT_Multi_mix[i][j].SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					XiStarpT_Multi_mix[i][j].Write("Xi1530_pT_%.1f_to_%.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					temp = ROOT.TH1F()
					temp = XiStarpT_Multi_mix[i][j].Project3D("Z")
					temp.SetTitle("pT Range %.1f to %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					#Normalization Factor
					nbin1 = temp.GetXaxis().FindBin(NormalizeRange_L[0])
					nbin2 = temp.GetXaxis().FindBin(NormalizeRange_L[1])
					nbin3 = temp.GetXaxis().FindBin(NormalizeRange_R[0])
					nbin4 = temp.GetXaxis().FindBin(NormalizeRange_R[1])
					integratedNumberforNormalizepT_mult_mix[i][j] = temp.Integral(nbin1,nbin2)+temp.Integral(nbin3,nbin4) #BOTH
					#integratedNumberforNormalizepT_mix.append(temp.Integral(nbin3,nbin4)) #R
					#integratedNumberforNormalizepT_mix.append(temp.Integral(nbin1,nbin2)) #L
					if (integratedNumberforNormalizepT_mult_mix[i][j] == 0): integratedNumberforNormalizepT_mult_mix[i][j] = 1
					if(PrintOption == 1): print("Normalizing factor: %f"%integratedNumberforNormalizepT_mult_mix[i][j])
					temp.SetTitle("Xi1530_pT_%.1f_to_%.1f_Multi_Mix_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp.Draw('E')
					temp.Write("Xi1530_pT_%.1f_to_%.1f_Multi_Mix_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					if(Save): c1.SaveAs("figs/Multibin/Xi1530_pT_%.1f_to_%.1f_Multi_Mix_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
			
			# Draw Together(pT, Mult)
			for i in range(0,len(pTRange)-1):
				if(PrintOption == 1): print("====================")
				if(PrintOption == 1): print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))	
				for j in range(0,len(MultiplicityRage)-1):
					### First, Draw together.
					c1.cd()
					temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp1.SetLineColor(0)
					temp1.Draw("E")
					temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_Mix_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp2.Scale(integratedNumberforNormalizepT_mult[i][j]/integratedNumberforNormalizepT_mult_mix[i][j])
					temp2.SetLineColor(1)
					temp2.Draw("E, same")
					temp2.SetTitle("pT Range %.1f to %.1f GeV/c, Multiplicity %d to %d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					c1.Write("Xi1530_pT_%.1f_to_%.1f_with_bkg_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					c1.SaveAs("figs/pTBin/same/Xi1530_pT_%.1f_to_%.1f_with_bkg_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

					### Second, Subtract bkg.
					temp1.Add(temp2,-1)
					temp1.Draw('E')
					temp1.Write("Xi1530_pT_%.1f_to_%.1f_with_bkg_substraction_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp1.Draw("E")
					c1.SaveAs("figs/pTBin/sub/Xi1530_pT_%.1f_to_%.1f_with_bkg_substraction_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

					### Third, Fit it.
					# Reference1: $ROOTSYS/tutorials/roofit/rf204_extrangefit.C
					# Reference2: https://github.com/clelange/roofit/blob/master/rf204_extrangefit.py
					if(PrintOption == 1): print("==================================================================================%.1f to %.1f case =================================================================================="%(pTRange[i],pTRange[i+1]))
					x = ROOT.RooRealVar("x", "x", AxisRange[0], AxisRange[1])
					mean1 = ROOT.RooRealVar("mean1", "mean of gaussians", 1.530,1.515,1.555)
					if(Fit_type == 0):
						sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", 0.003,0.000001,0.01)
						width = ROOT.RooRealVar("width", "width of voigtian", 0.0091)
					if(Fit_type == 1):
						sigma1 = ROOT.RooRealVar("sigma1", "width of gaussians", FitSigma)
						width = ROOT.RooRealVar("width", "width of voigtian",0.0091,0.000,0.05)
					signal = ROOT.RooVoigtian("signal", "Signal component 1", x, mean1, width, sigma1)

					a0 = ROOT.RooRealVar("a0","a0",1) ;
					a1 = ROOT.RooRealVar("a1","a1",0,-1,1) ;
					a2 = ROOT.RooRealVar("a2","a2",1,0,10) ;
					background = ROOT.RooPolynomial("p2","p2",x,ROOT.RooArgList(a0,a1,a2),0) ;
					
					# --- Construct signal+background PDF ---
					nsig = ROOT.RooRealVar("nsig","#signal events",2000,0.,10e5) ;
					#nbkg = ROOT.RooRealVar("nbkg","#background events",800,0.,1000) ;
					#model = ROOT.RooAddPdf("model","g+a",ROOT.RooArgList(signal,background),ROOT.RooArgList(nsig,nbkg))
					model = ROOT.RooAddPdf("model","g+a",ROOT.RooArgList(signal),ROOT.RooArgList(nsig))

					data = ROOT.RooDataHist("RooDataHist","dataset with Gauss  Fit",ROOT.RooArgList(x),temp1)

					c1.cd()
					ROOT.gPad.SetTickx(2)
					xframe = x.frame()
					#signal.fitTo(data,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
					model.fitTo(data,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
					#signal.fitTo(data)
					#model.fitTo(data)
		 			#model.fitTo(data,ROOT.RooFit.Range(1.52,1.55))

					data.plotOn(xframe)
					#signal.plotOn(xframe,ROOT.RooFit.LineStyle(ROOT.kDashed))
					model.plotOn(xframe,ROOT.RooFit.Range(fitRange[0],fitRange[1]),ROOT.RooFit.LineStyle(ROOT.kDashed))
					#model.plotOn(xframe)
					#model.plotOn(xframe,ROOT.RooFit.Components("p2"),ROOT.RooFit.LineStyle(ROOT.kDashed))
					#model.plotOn(xframe,ROOT.RooFit.Components("p2"),RooFit.LineColor(kGreen))
					#sig.paramOn(xframe, ROOT.RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
					xframe.SetMarkerStyle(1)
					xframe.Draw("E")
					xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
					c1.Clear()
					data.plotOn(xframe)
					#signal.plotOn(xframe)
					model.plotOn(xframe,ROOT.RooFit.Range(fitRange[0],fitRange[1]))
					#model.plotOn(xframe,ROOT.RooFit.LineStyle(ROOT.kDashed))
					#signal.plotOn(xframe,ROOT.RooFit.Components("background"),ROOT.RooFit.LineStyle(ROOT.kDashed))
					#model.plotOn(xframe,ROOT.RooFit.Components("p2"),ROOT.RooFit.LineStyle(ROOT.kDashed))
					#model.plotOn(xframe,ROOT.RooFit.Components("p2"),RooFit.LineColor(kGreen))
					#sig.paramOn(xframe, ROOT.RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
					xframe.SetMarkerStyle(2)
					xframe.Draw("E")		
					xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

								#chi2/ndf
					fitChi2ndf.append(xframe.chiSquare())
					if(PrintOption == 1): print("chi2: %.4f"%xframe.chiSquare())  # NEED TO ADD NDF INFO!!
					

					# --- Extract fit parameters ---
					#if(PrintOption == 1): print("Fit Result: Mean: %.3f, Sigma: %.5f"%(mean1.getValV(),sigma1.getValV()))
					fitParameter_mean_mult[j][i] = mean1.getValV()
					fitParameter_mean_err_mult[j][i] = mean1.getError()
					if(Fit_type == 0):
						fitParameter_sigma_mult[j][i] = sigma1.getValV()
						fitParameter_sigma_err_mult[j][i] = sigma1.getError()
					if(Fit_type == 1):
						fitParameter_sigma_mult[j][i] = width.getValV()
						fitParameter_sigma_err_mult[j][i] = width.getError()

					# Get the raw yield
					#fitIntegral = signal.createIntegral(mesArgSet,ROOT.RooFit.NormSet(mesArgSet),ROOT.RooFit.Range("IntegralRange"))
					fitIntegral = model.createIntegral(ROOT.RooArgSet(x),ROOT.RooArgSet(x),"IntegralRange")
					RawYield_Multi[i][j].append(nsig.getValV()*fitIntegral.getVal())
					fitIntegral_err = fitIntegral.getValV()*nsig.getError()/nsig.getValV()
					RawYield_Multi[i][j].append(fitIntegral_err)
		
		#For Result Presensation
		c2 = ROOT.TCanvas('c2','Invariant Mass Distribution with pT bins',1920,1080)
		c2.Clear()
		if(len(pTRange)<3):
			print("no need to divide")
		elif(len(pTRange)<9):
			c2.Divide(3,3,0.01,0.01)
		else:
			c2.Divide(4,3,0.01,0.01)

		t = ROOT.TLatex()
		t.SetNDC()
		t.SetTextSize(0.075)
		
		# RAW signal with pT 
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT.%s"%SaveType)


		# RAW signal with pT of Xi
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/QA/Invmass_pT_Xi.%s"%SaveType)

		# RAW signal+bkg with pT 
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			temp1.SetLineColor(0)
			temp1.Draw("E")
			temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			temp2.Scale(integratedNumberforNormalizepT[i]/integratedNumberforNormalizepT_mix[i])
			temp2.SetLineColor(1)
			temp2.Draw("E, same")
			temp2.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_with_bkg.%s"%SaveType)
		
		# Fitted plot with pT
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			#t.DrawLatex(0.6,0.8,"\chi^{2} / ndf = %.1f"%fitChi2ndf[i])
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_fitted.%s"%SaveType)

		# Fitted plot with pT(Dashed)
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
			#t.DrawLatex(0.6,0.8,"\chi^{2} / ndf = %.1f"%fitChi2ndf[i])
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_fitted(dash).%s"%SaveType)


		c1.cd()
		# Fit mean value plot with pT
		fitMeanAspT = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(fitParameter_mean, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(fitParameter_mean_err, 'd'))
		fitMeanAspT.SetMarkerStyle(20)
		fitMeanAspT.SetMarkerColor(2)
		fitMeanAspT.SetTitle("")
		fitMeanAspT.GetYaxis().SetTitle("Mass (GeV/c^{2})")
		fitMeanAspT.GetYaxis().SetTitleOffset(1.3)
		fitMeanAspT.GetXaxis().SetTitle("p_{T} (GeV/c)")
		fitMeanAspT.Draw("AP")
		pdg_massline = ROOT.TF1("pdg_massline","1.53178",-1,pTRange[-1]+1)
		pdg_massline.Draw("same")
		
		c1.Write("Xi1530_fit_mean")
		c1.SaveAs("figs/fit_means.%s"%SaveType)

		# Fit sigma value plot with pT
		fitSigmaAspT = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(fitParameter_sigma, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(fitParameter_sigma_err, 'd'))
		fitSigmaAspT.SetMarkerStyle(20)
		fitSigmaAspT.SetMinimum(0);
		if(Fit_type == 0):
			fitSigmaAspT.SetMaximum(0.0050);
			fitSigmaAspT.GetYaxis().SetTitle("#sigma (GeV/c^{2})")
		if(Fit_type == 1):
			fitSigmaAspT.SetMaximum(0.020);
			fitSigmaAspT.GetYaxis().SetTitle("#Gamma (GeV/c^{2})")
		fitSigmaAspT.SetTitle("")
		fitSigmaAspT.GetXaxis().SetTitle("p_{T} (GeV/c)")
		fitSigmaAspT.Draw("AP")
		if(Fit_type == 1):
			pdg_widthline = ROOT.TF1("pdg_widthline","0.0091",-1,pTRange[-1]+1)
			pdg_widthline.Draw("same")
		
		c1.Write("Xi1530_fit_sigma")
		c1.SaveAs("figs/fit_sigma.%s"%SaveType)

		# RAW Yield plot with pT
		RawYieldAspT = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(RawYield, 'd') , np.asarray(pTRange_err, 'd'), np.asarray(RawYield_err, 'd'))
		RawYieldAspT.SetMarkerStyle(20)
		RawYieldAspT.SetMinimum(0);
		RawYieldAspT.GetYaxis().SetTitle("Raw Yield")
		RawYieldAspT.SetTitle("")
		RawYieldAspT.GetXaxis().SetTitle("p_{T} (GeV/c)")
		RawYieldAspT.Draw("AP")
		
		c1.Write("Xi1530_Raw_Yield")
		c1.SaveAs("figs/Raw_Yield.%s"%SaveType)

		if(MultiplicitySet=="1"):
			#For Result Presensation with multiplicity bin
			c3 = ROOT.TCanvas('c3','Invariant Mass Distribution with pT bins',1920,1080)
			if(len(pTRange)<3):
				print("no need to divide")
			elif(len(pTRange)<9):
				c3.Divide(3,3,0.01,0.01)
			else:
				c3.Divide(4,3,0.01,0.01)

			t2 = ROOT.TLatex()
			t2.SetNDC()
			t2.SetTextSize(0.075)
			for j in range(0,len(MultiplicityRage)-1):
				for i in range(0,len(pTRange)-1):
					c3.cd(i+1)
					hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					hist.SetTitle("pT %.1f to %.1f with Multiplicity %d - %d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					hist.Draw('E')
					t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
				c3.Draw()
				c3.SaveAs("figs/Invmass_Multi_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

			# RAW signal+bkg with pT 
			t2 = ROOT.TLatex()
			t2.SetNDC()
			t2.SetTextSize(0.075)
			for j in range(0,len(MultiplicityRage)-1):
				for i in range(0,len(pTRange)-1):
					c3.cd(i+1)
					temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp1.SetLineColor(0)
					temp1.Draw("E")
					temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_Mix_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					temp2.Scale(integratedNumberforNormalizepT_mult[i][j]/integratedNumberforNormalizepT_mult_mix[i][j])
					temp2.SetLineColor(1)
					temp2.Draw("E, same")
					temp2.SetTitle("pT Range %.1f to %.1f GeV/c, Multiplicity %d to %d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
				c3.Draw()
				c3.SaveAs("figs/Invmass_with_bkg_Multi_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

			# Fitted plot with pT
			for j in range(0,len(MultiplicityRage)-1):
				for i in range(0,len(pTRange)-1):
					c3.cd(i+1)
					hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					hist.Draw('E')
					t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					#t.DrawLatex(0.6,0.8,"\chi^{2} / ndf = %.1f"%fitChi2ndf[i])
				c3.Draw()
				c3.SaveAs("figs/Invmass_pT_fitted_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
			
			# Fitted plot with pT(Dashed)
			for j in range(0,len(MultiplicityRage)-1):
				for i in range(0,len(pTRange)-1):
					c3.cd(i+1)
					hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
					hist.SetTitle("pT %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					hist.Draw('E')
					t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f GeV/c"%(pTRange[i],pTRange[i+1]))
					#t.DrawLatex(0.6,0.8,"\chi^{2} / ndf = %.1f"%fitChi2ndf[i])
				c3.Draw()
				c3.SaveAs("figs/Invmass_pT_fitted(dash)_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

			
			c1.cd()
			# Fit mean value plot with pT
			for j in range(0,len(MultiplicityRage)-1):
				fitMeanAspT_mult = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(fitParameter_mean_mult[j], 'd') , np.asarray(pTRange_err, 'd'), np.asarray(fitParameter_mean_err_mult[j], 'd'))
				fitMeanAspT_mult.SetMarkerStyle(20)
				fitMeanAspT_mult.SetMarkerColor(2)
				fitMeanAspT_mult.SetTitle("")
				fitMeanAspT_mult.GetYaxis().SetTitle("Mass (GeV/c^{2})")
				fitMeanAspT_mult.GetYaxis().SetTitleOffset(1.3)
				fitMeanAspT_mult.GetXaxis().SetTitle("p_{T} (GeV/c)")
				fitMeanAspT_mult.Draw("AP")
				pdg_massline = ROOT.TF1("pdg_massline","1.53178",-1,pTRange[-1]+1)
				pdg_massline.Draw("same")
				
				c1.Write("Xi1530_fit_mean_%d_to_%d"%(MultiplicityRage[j],MultiplicityRage[j+1]))
				c1.SaveAs("figs/fit_means_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
			
			# Fit sigma value plot with pT
			for j in range(0,len(MultiplicityRage)-1):
				fitSigmaAspT_mult = ROOT.TGraphErrors(len(pTRange)-1,np.asarray(pTCenter, 'd'),np.asarray(fitParameter_sigma_mult[j], 'd') , np.asarray(pTRange_err, 'd'), np.asarray(fitParameter_sigma_err_mult[j], 'd'))
				fitSigmaAspT_mult.SetMarkerStyle(20)
				fitSigmaAspT_mult.SetMinimum(0);
				if(Fit_type == 0):
					fitSigmaAspT_mult.SetMaximum(0.0050);
					fitSigmaAspT_mult.GetYaxis().SetTitle("#sigma (GeV/c^{2})")
				if(Fit_type == 1):
					fitSigmaAspT_mult.SetMaximum(0.020);
					fitSigmaAspT_mult.GetYaxis().SetTitle("#Gamma (GeV/c^{2})")
				fitSigmaAspT_mult.SetTitle("")
				fitSigmaAspT_mult.GetXaxis().SetTitle("p_{T} (GeV/c)")
				fitSigmaAspT_mult.Draw("AP")
				if(Fit_type == 1):
					pdg_widthline = ROOT.TF1("pdg_widthline","0.0091",-1,pTRange[-1]+1)
					pdg_widthline.Draw("same")
					print("Fit width in %d_to_%d:"%(MultiplicityRage[j],MultiplicityRage[j+1]))
					print(fitParameter_sigma_mult[j])
					print(fitParameter_sigma_err_mult[j])
				
				c1.Write("Xi1530_fit_sigma_%d_to_%d"%(MultiplicityRage[j],MultiplicityRage[j+1]))
				c1.SaveAs("figs/fit_sigma_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
			#Fit width plot with Multiplicity with full pT range
			if(Fit_type == 1):
				fitSigmaAsmult = ROOT.TGraphErrors(len(MultiplicityRage)-1,np.asarray(MultCenter, 'd'),np.asarray(fitParameter_sigma_mult, 'd') , np.asarray(MultRange_err, 'd'), np.asarray(fitParameter_sigma_err_mult, 'd'))
				fitSigmaAsmult.SetMarkerStyle(20)
				fitSigmaAsmult.SetMinimum(0);
				fitSigmaAsmult.SetMaximum(0.020);
				fitSigmaAsmult.GetYaxis().SetTitle("#Gamma (GeV/c^{2})")
				fitSigmaAsmult.SetTitle("")
				fitSigmaAsmult.GetXaxis().SetTitle("Multiplicity percentile (%)")
				fitSigmaAsmult.Draw("AP")
				pdg_widthline = ROOT.TF1("pdg_widthline","0.0091",-1,MultiplicityRage[-1])
				pdg_widthline.Draw("same")
				
				c1.Write("Xi1530_fit_gamma")
				c1.SaveAs("figs/fit_gamma.%s"%SaveType)
		#print(RawYield)
		#print(RawYield_err)
	try:
		mylist_data, mylist_mc
	except NameError:
		print("")
	else: # Using data/mc together
		c0.cd()
		c0.SetLogy()
		
		# Corrected Yield Spectra in MB

		h_Xispectrum = ROOT.TH1D("h_Xispectrum","Xi spectrum",len(pTRange)-1,np.asarray(pTRange, 'd'))
		h_Xispectrum.SetMarkerStyle(20)
		h_Xispectrum.SetMarkerSize(1.0)
		h_Xispectrum.SetTitle("#Xi(1530) spectrum pp #sqrt{s}= 13 TeV")
		h_Xispectrum.GetXaxis().SetTitle("p_{T} (GeV/c)")
		h_Xispectrum.GetYaxis().SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}")
		h_Xispectrum.SetMinimum(3.0e-6)
		h_Xispectrum.SetMaximum(0.2)
		#h_Xispectrum.GetXaxis().SetLimits(0,5.6)

		h_Xispectrum.SetBinContent(1,+100)
		for pT in range(1,len(pTRange)-1):
			if(pT == len(pTRange)-2):
				continue
			#base = MCEfficiency[ptbin]*LostRatio*(raprange[1]-raprange[0])*(2*pt_points_e[ptbin])*Events_PhysSel->GetEntries()/0.852
			#LostRatio -> Calculated the ratio between prePV cut and after. in this case, it was 1.
			#raprange -> -0.5 to 0.5 -> 1
			#pt_points_e
			#print(pT)
			base = MCEfficiency[pT]*(2*pTRange_err[pT])*Events_PhysSel.GetEntries()/0.9
			h_Xispectrum.SetBinContent(pT+1, RawYield[pT]/base)
			#print(RawYield[pT]/base)
			h_Xispectrum.SetBinError(pT+1, RawYield_err[pT]/base)
		h_Xispectrum.Draw()
		c0.SaveAs("figs/Xi1530_spectra.%s"%SaveType)

		#Levy fit

		Levyfit = ROOT.TF1("Levyfit",myLevyPtFunc,0,8.5,3)
		Levyfit.SetParName(0,"dN/dy")
		Levyfit.SetParName(1,"C")
		Levyfit.SetParName(2,"n")
		Levyfit.SetParameter(0,.008)
		Levyfit.SetParameter(1,.3)
		Levyfit.SetParameter(2,15)
		Levyfit.SetParLimits(0,.001,.2)
		Levyfit.SetParLimits(1,.1,1)
		Levyfit.SetParLimits(2,1,500)

		h_Xispectrum.Fit(Levyfit,"IME","",pTRange[1],pTRange[-2])
		Levyfit.Draw("same")
		legend = ROOT.TLegend(0.1,0.7,0.35,0.9,"","brNDC")
		legend.SetBorderSize(1)
		legend.SetTextSize(.04) #small .03; large .036 
		#legend->SetLineColor(0)
		legend.SetFillColor(0)

		legend.AddEntry(h_Xispectrum,"(#Xi(1530) + cc)/2","p");
		legend.AddEntry(Levyfit,"Levy","l");
		c0.SaveAs("figs/Xi1530_spectra_levy.%s"%SaveType)

		#Do same thing in several multiplicities

		for i in range(0,len(MultiplicityRage)-1):
			# Corrected Yield Spectra in Mutli

			h_Xispectrum_mult.append(ROOT.TH1D("h_Xispectrum","Xi spectrum",len(pTRange)-1,np.asarray(pTRange, 'd')))
			h_Xispectrum_mult[i].SetMarkerStyle(20)
			h_Xispectrum_mult[i].SetMarkerSize(1.0)
			h_Xispectrum_mult[i].SetTitle("#Xi(1530) spectrum pp #sqrt{s}= 13 TeV")
			h_Xispectrum_mult[i].GetXaxis().SetTitle("p_{T} (GeV/c)")
			h_Xispectrum_mult[i].GetYaxis().SetTitle("1/N_{E}d^{2}N/dydp_{T} (GeV/c)^{-1}")
			h_Xispectrum_mult[i].SetMinimum(3.0e-6)
			h_Xispectrum_mult[i].SetMaximum(0.2)
			#h_Xispectrum_mult[i].GetXaxis().SetLimits(0,5.6)

			h_Xispectrum_mult[i].SetBinContent(1,+100)
			for pT in range(1,len(pTRange)-1):
				if(pT == len(pTRange)-2):
					continue
				#base = MCEfficiency[ptbin]*LostRatio*(raprange[1]-raprange[0])*(2*pt_points_e[ptbin])*Events_PhysSel->GetEntries()/0.852
				#LostRatio -> Calculated the ratio between prePV cut and after. in this case, it was 1.
				#raprange -> -0.5 to 0.5 -> 1
				#pt_points_e
				#print(pT)
				base = MCEfficiency[pT]*(2*pTRange_err[pT])*Events_PhysSel.GetEntries()/0.9
				h_Xispectrum_mult[i].SetBinContent(pT+1, RawYield[pT]/base)
				#print(RawYield[pT]/base)
				h_Xispectrum_mult[i].SetBinError(pT+1, RawYield_err[pT]/base)
			h_Xispectrum_mult[i].Draw()
			c0.SaveAs("figs/Xi1530_spectra_%d_to_%d.%s"%(MultiplicityRage[i],MultiplicityRage[i+1],SaveType))



if __name__ == "__main__":
		DrawResults()