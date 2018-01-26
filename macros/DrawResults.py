# DrawResults.py
# Draw histograms with specific criteria and fit
# 
# Author: Bong-Hwi Lim (bong-hwi.lim@cern.ch)
# 
# How to use:
# python DrawResults.py AnalysisResults.root 1 MC
# parameters:
#   - AnalysisResults.root: Input root file
#   - 1: 1 or 0, 1: save all outputs 0: save only merged one
#   - MC: MC or blank, MC: make a plot from MC data, or not.
#
from ROOT import *
import numpy as np
import datetime
import sys

currenttime = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
print("Current time: %s"%currenttime)

Inputfile = "AnalysisResults.root"
InputSaveType = 0
isMC = 0
if(len(sys.argv)>=2): 
	Inputfile = sys.argv[1]
	print("input file: %s"%(Inputfile))

	if(len(sys.argv)>=3): 
		InputSaveType = sys.argv[2]
		print("Outputs will be saved")

		if(len(sys.argv)>=4): 
			isMC = sys.argv[3]
			print("MC data mode")

def DrawResults(InputSaveType, Inputfile, isMC):
	#====================
	# Inputs
	inputFileName = Inputfile
	outputFileName = "Analysis_Results_Xi1530_%s.root"%currenttime
	pTRange = [0, 0.8, 1.2, 1.6, 2.0, 2.4, 3.2, 4.0, 4.8, 5.6]
	massRebin = 10
	AxisRange = [1.48,1.59]
	#NormalizeRange = [1.49,1.51] #for LHC15f
	NormalizeRange = [1.56,1.58] #for LHC16k
	fitRange = [1.5,1.55]
	MultiplicityRage = [0,5,15,30,50,100]
	hNames = ["fXiMinusPiPlus_0", "fXiPlusPiMinus_0"] #fXiMinusPiPlus_0: Xi(1530), fXiPlusPiMinus_0: Anti-Xi(1530)
	hNames_mix = ["fXiMinusPiPlusbkg_0", "fXiPlusPiMinusbkg_0"] #fXiMinusPiPlus_0: Xi(1530), fXiPlusPiMinus_0: Anti-Xi(1530)
	hNames_MCinput = ["fMCinputTotalXiStar3", "fMCinputTotalXiStarbar3"] #generated MC
	hNames_Mcrecon = ["fMCrecXiMinusPiPlus_0", "fMCrecXiPlusPiMinus_0"] #recon MC
	Save = InputSaveType #0: skip saving pdf   images, 1: save pdf   images
	SaveType = "png" #png, pdf
	#====================
	#====================
	# Default lists
	XiStar = []
	XiStar_mix = []
	XiStarpT = []
	XiStarpT_mix = []

	XiStarMCinput = []
	XiStarMCrecon = []

	integratedNumberforNormalizepT = []
	integratedNumberforNormalizepT_mix = []
	
	XiStarpT_Multi = [[0 for col in range(len(pTRange)-1)] for row in range(len(MultiplicityRage)-1)]

	NormalizationFactor = [[0 for col in range(len(pTRange)-1)] for row in range(len(hNames)/4)]
	RawYield = [[0 for col in range(len(pTRange)-1)] for row in range(len(hNames)/4)]
	
	#===================
	#===================
	# Load the input root file
	print("====================")
	print("Load file: "+inputFileName)
	f = TFile(inputFileName)
	mydir = f.Get("PWGLF.outputXiStarAnalysis.root;1")
	mylist = mydir.Get("MyList;1")
	print("====================")
	print("Load file: "+outputFileName)
	fo = TFile(outputFileName,"RECREATE");
	fo.cd()
	print("====================")
	#===================
	#===================
	if(isMC): #MC Case.
		i = 0
		for name in hNames_MCinput:
			print("Load Histogram: "+name)
			XiStarMCinput.append(mylist.FindObject(name))
			XiStarMCinput[i].GetXaxis().SetTitle("pT_Range")
			XiStarMCinput[i].GetYaxis().SetTitle("Multiplicity")
			XiStarMCinput[i].GetZaxis().SetTitle("Inv_Mass")
			XiStarMCinput[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		print("====================")

		i = 0
		for name in hNames_Mcrecon:
			print("Load Histogram: "+name)
			XiStarMCrecon.append(mylist.FindObject(name))
			XiStarMCrecon[i].GetXaxis().SetTitle("pT_Range")
			XiStarMCrecon[i].GetYaxis().SetTitle("Multiplicity")
			XiStarMCrecon[i].GetZaxis().SetTitle("Inv_Mass")
			XiStarMCrecon[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		print("====================")

		# Add them each other
		XiStarMCinput[0].Add(XiStarMCinput[1])
		XiStarMCrecon[0].Add(XiStarMCrecon[1])

		### Draw MC QA Histos
		c1 = TCanvas('can','canvas',1280,720)
		c1.SetBorderMode(1)					 
		c1.cd() 

		# Generated Xi(1530)
		temp = TH1F()
		temp = XiStarMCinput[0].Project3D("X")
		temp.GetXaxis().SetTitle("pT_Range")
		temp.SetTitle("Generated Xi(1530)")

		# Reconstructed Xi(1530)
		temp2 = TH1F()
		temp2 = XiStarMCrecon[0].Project3D("X")
		temp2.GetXaxis().SetTitle("pT_Range")
		temp2.SetTitle("Reconstructed Xi(1530)")

		temp.Draw("E")
		c1.SetLogy()
		c1.SaveAs("figs/mc/Xi1530_MCQA_NumberWithpT_generated.pdf")

		temp2.Draw("E")
		c1.SetLogy()
		c1.SaveAs("figs/mc/Xi1530_MCQA_NumberWithpT_reconstructed.pdf")		
		
	else: # Data Case
		i = 0
		for name in hNames:
			print("Load Histogram: "+name)
			XiStar.append(mylist.FindObject(name))
			XiStar[i].GetXaxis().SetTitle("pT_Range")
			XiStar[i].GetYaxis().SetTitle("Multiplicity")
			XiStar[i].GetZaxis().SetTitle("Inv_Mass")
			XiStar[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		print("====================")

		i = 0
		for name in hNames_mix:
			print("Load Histogram: "+name)
			XiStar_mix.append(mylist.FindObject(name))
			XiStar_mix[i].GetXaxis().SetTitle("pT_Range")
			XiStar_mix[i].GetYaxis().SetTitle("Multiplicity")
			XiStar_mix[i].GetZaxis().SetTitle("Inv_Mass")
			XiStar_mix[i].SetAxisRange(AxisRange[0],AxisRange[1],"Z")
			i = i+1
		print("====================")	

		#Make a Canvas
		c1 = TCanvas('can','canvas',1280,720)
		c1.SetBorderMode(1)					 
		c1.cd() 

		# Add Antiparticle
		print("Entry of Xi(1530)       : %d"%XiStar[0].GetEntries())
		print("Entry of anti-Xi(1530)  : %d"%XiStar[1].GetEntries())
		XiStar[0].Add(XiStar[1])
		print("Entry of merged Xi(1530): %d"%XiStar[0].GetEntries())

		XiStar[0].Write("Original_Xi1530_pT_Multi_Mass") #Original 3D histogram

		#Slice Histogram with given pT bin
		for i in range(0,len(pTRange)-1):
			print("====================")
			print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
			XiStarpT.append(XiStar[0].Clone())
			XiStarpT[i].SetAxisRange(pTRange[i],pTRange[i+1],"X")
			XiStarpT[i].SetTitle("pT Range %.1f to %.1f"%(pTRange[i],pTRange[i+1]))
			XiStarpT[i].Write("Xi1530_pT_%.1f_to_%.1f"%(pTRange[i],pTRange[i+1]))
			temp = TH1F()
			temp = XiStarpT[i].Project3D("Z")
			temp.SetTitle("pT Range %.1f to %.1f"%(pTRange[i],pTRange[i+1]))
			#Normalization Factor
			nbin1 = temp.GetXaxis().FindBin(NormalizeRange[0])
			nbin2 = temp.GetXaxis().FindBin(NormalizeRange[1])
			integratedNumberforNormalizepT.append(temp.Integral(nbin1,nbin2))
			print("Normalizing factor: %f"%integratedNumberforNormalizepT[i])
			temp.Draw('E')
			temp.Write("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			if(Save): c1.SaveAs("figs/pTBin/sig/Xi1530_pT_%.1f_to_%.1f.%s"%(pTRange[i],pTRange[i+1],SaveType))

		#Slice Histogram with given pT bin for event mix
		for i in range(0,len(pTRange)-1):
			print("=====Event Mixing=====")
			print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))
			XiStarpT_mix.append(XiStar_mix[0].Clone())
			XiStarpT_mix[i].SetAxisRange(pTRange[i],pTRange[i+1],"X")
			XiStarpT_mix[i].SetTitle("pT Range %.1f to %.1f"%(pTRange[i],pTRange[i+1]))
			XiStarpT_mix[i].Write("Xi1530_pT_%.1f_to_%.1f"%(pTRange[i],pTRange[i+1]))
			temp = TH1F()
			temp = XiStarpT_mix[i].Project3D("Z")
			temp.SetTitle("pT Range %.1f to %.1f"%(pTRange[i],pTRange[i+1]))
			#Normalization Factor
			nbin1 = temp.GetXaxis().FindBin(NormalizeRange[0])
			nbin2 = temp.GetXaxis().FindBin(NormalizeRange[1])
			integratedNumberforNormalizepT_mix.append(temp.Integral(nbin1,nbin2))
			print("Normalizing factor: %f"%integratedNumberforNormalizepT_mix[i])
			temp.Draw('E')
			temp.Write("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			if(Save): c1.SaveAs("figs/pTBin/bkg/Xi1530_pT_%.1f_to_%.1f_EventMix.%s"%(pTRange[i],pTRange[i+1],SaveType))

		# Draw Together(pT)
		for i in range(0,len(pTRange)-1):
			### First, Draw together.
			c1.cd()
			temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			temp1.SetLineColor(kBlack)
			temp1.Draw("E")
			temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			temp2.Scale(integratedNumberforNormalizepT[i]/integratedNumberforNormalizepT_mix[i])
			temp2.SetLineColor(kRed)
			temp2.Draw("E, same")
			temp2.SetTitle("pT Range %.1f to %.1f"%(pTRange[i],pTRange[i+1]))
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
		
			mes = RooRealVar("mes", "m_{Xi} (GeV)", AxisRange[0], AxisRange[1])
			mean1 = RooRealVar("mean1", "mean of gaussians", 1.530,1.48,1.59)
			sigma1 = RooRealVar("sigma1", "width of gaussians", 0.03,0,0.3)
			sigma2 = RooRealVar("sigma2", "width of gaussians", 0.2,0,0.6)
			signal = RooBreitWigner("signal", "Signal component 1", mes, mean1, sigma1)

			a0 = RooRealVar("a0","a0",1) ;
			a1 = RooRealVar("a1","a1",0,-1,1) ;
			a2 = RooRealVar("a2","a2",1,0,10) ;
			background = RooPolynomial("p2","p2",mes,RooArgList(a0,a1,a2),0) ;
			
			# --- Construct signal+background PDF ---
			nsig = RooRealVar("nsig","#signal events",200,0.,1000) ;
			nbkg = RooRealVar("nbkg","#background events",800,0.,1000) ;
			model = RooAddPdf("model","g+a",RooArgList(signal,background),RooArgList(nsig,nbkg))

			data = RooDataHist("RooDataHist","dataset with Gauss  Fit",RooArgList(mes),temp1)

			c1.cd()
			gPad.SetTickx(2)
			xframe = mes.frame()
			model.fitTo(data)
			#modle.fitTo(data,RooFit.Range("fitRange"))

			data.plotOn(xframe)
			#model.plotOn(xframe,RooFit.Components("sig"),RooFit.LineColor(kRed))
			model.plotOn(xframe,RooFit.LineStyle(kDashed))
			model.plotOn(xframe,RooFit.Components("background"),RooFit.LineStyle(kDashed))
			#sig.paramOn(xframe, RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
			xframe.SetMarkerStyle(2)
			xframe.Draw("E")		
			xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)"%(pTRange[i],pTRange[i+1]))
			c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting(dash).%s"%(pTRange[i],pTRange[i+1],SaveType))
			c1.Clear()
			data.plotOn(xframe)
			model.plotOn(xframe)
			model.plotOn(xframe,RooFit.Components("background"),RooFit.LineStyle(kDashed))
			#sig.paramOn(xframe, RooFit.Layout(0.6),RooFit.Format("NEU",RooFit.AutoPrecision(1)))
			xframe.SetMarkerStyle(2)
			xframe.Draw("E")		
			xframe.Write("Xi1530_pT_%.1f_to_%.1f_with_fitting"%(pTRange[i],pTRange[i+1]))
			c1.SaveAs("figs/pTBin/fit/Xi1530_pT_%.1f_to_%.1f_with_fitting.%s"%(pTRange[i],pTRange[i+1],SaveType))
			

		#Slice Histogram with given Multiplicity
		for i in range(0,len(pTRange)-1):
			print("====================")
			print("pT Rage: "+str(pTRange[i])+" to "+str(pTRange[i+1]))	
			for j in range(0,len(MultiplicityRage)-1):
				print("--Multiplicity Rage: "+str(MultiplicityRage[j])+" to "+str(MultiplicityRage[j+1]))
				XiStarpT_Multi[j][i] = XiStarpT[i].Clone()
				XiStarpT_Multi[j][i].SetAxisRange(MultiplicityRage[j],MultiplicityRage[j+1],"Y")
				temp = TH1F()
				temp = XiStarpT_Multi[j][i].Project3D("Z")			
				temp.SetTitle("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
				temp.Draw('E')
				temp.Write("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
				if(Save): c1.SaveAs("figs/Multibin/Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d.%s"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1],SaveType))
		
		#For Result Presensation
		c2 = TCanvas('c2','Invariant Mass Distribution with pT bins',1440,1080)
		c2.Clear()
		c2.Divide(3,3,0.01,0.01)

		t = TLatex()
		t.SetNDC()
		t.SetTextSize(0.075)
		
		# RAW signal with pT 
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT.%s"%SaveType)

		# RAW signal+bkg with pT 
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			temp1 = fo.Get("Xi1530_pT_%.1f_to_%.1f_projected"%(pTRange[i],pTRange[i+1]))
			temp1.SetLineColor(kBlack)
			temp1.Draw("E")
			temp2 = fo.Get("Xi1530_pT_%.1f_to_%.1f_EventMix"%(pTRange[i],pTRange[i+1]))
			temp2.Scale(integratedNumberforNormalizepT[i]/integratedNumberforNormalizepT_mix[i])
			temp2.SetLineColor(kRed)
			temp2.Draw("E, same")
			temp2.SetTitle("pT %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_with_bkg.%s"%SaveType)
		
		# Fitted plot with pT
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_fitted.%s"%SaveType)

		# Fitted plot with pT(Dashed)
		for i in range(0,len(pTRange)-1):
			c2.cd(i+1)
			hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_with_fitting(dash)"%(pTRange[i],pTRange[i+1]))
			hist.SetTitle("pT %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
			hist.Draw('E')
			t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
		c2.Draw()
		c2.SaveAs("figs/Invmass_pT_fitted(dash).%s"%SaveType)

		#For Result Presensation with multiplicity bin
		c3 = TCanvas('c3','Invariant Mass Distribution with pT bins',1440,1080)
		c3.Divide(3,3,0.01,0.01)

		t2 = TLatex()
		t2.SetNDC()
		t2.SetTextSize(0.075)
		for j in range(0,len(MultiplicityRage)-1):
			for i in range(0,len(pTRange)-1):
				c3.cd(i+1)
				hist=fo.Get("Xi1530_pT_%.1f_to_%.1f_Multi_%d_to_%d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
				hist.SetTitle("pT %.1f to %.1f with Multiplicity %d - %d"%(pTRange[i],pTRange[i+1],MultiplicityRage[j],MultiplicityRage[j+1]))
				hist.Draw('E')
				t.DrawLatex(0.15,0.8,"P_{T}: %.1f - %.1f"%(pTRange[i],pTRange[i+1]))
			c3.Draw()
			c3.SaveAs("figs/Invmass_Multi_%d_to_%d.%s"%(MultiplicityRage[j],MultiplicityRage[j+1],SaveType))

if __name__ == "__main__":
		DrawResults(InputSaveType, Inputfile, isMC)