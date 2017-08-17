
import  sys
#sys.argv.append( '-b-' )
import ROOT, os, re, string, csv, math, sys, gc
from ROOT import TCanvas, TTree, TFile, TPaveStats, TProfile, TNtuple, TH1F, TH1D, TH2F, TF1, TGaxis, TPad, TPaveLabel, TPaveText, TLegend, TLatex, THStack, TLine, TMath, TGraph, TGraphErrors, TMultiGraph, TStyle
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad, gStyle
from subprocess import call
from array import array
ROOT.gROOT.SetBatch()
#from tdrStyle import *
#setTDRStyle()
#gROOT.ProcessLine(".L ./tdrstyle.C")
#gROOT.ProcessLine("setTDRStyle()")
gROOT.ProcessLine(".L ./rootlogon.C")
gStyle.SetPalette(1);
gStyle.SetOptStat("mner");
#gStyle.SetOptStat(0);



#Dictoinaries for later use
f1 ={}
f2 = {}
Tree = {}
Histograms = {}

#Do not include additional spaces before or after the strings in any of these arrays!
Files = ["ntuple3TeV.root", "ntuple6TeV.root"]

#These arrays correspond to the variables of interest; be sure that they all have the same size
Data = ["genpart_KKPhoton_pt", "genpart_KKPhoton_mass", "genpart_KKPhoton_eta"]
HistogramTitles = ["KK Photon pt", "KK Photon Mass", "KK Photon eta"]
HistogramNames = [["3TeV pT", "6TeV pt"], ["3TeV Mass", "6TeV Mass"], ["3TeV eta", "6TeV eta"]]
ImageNames = ["Test1", "Test2", "Test3"]
Bins = [80, 300, 80]
X_Low = [0, 2500, 0]
X_High = [150, 6500, 15]
Y_Low = [0, 0, 0]
Y_High = [135, 2000, 120]
X_Label = ["pT (GeV/c)", "Mass (GeV/c^2)", "eta"]
Y_Label = ["Number of Events", "Number of Events", "Number of Events"]
LogPlot = [0, 1, 0]

#These arrays have slightly different sizes; they correspond to [[variable 1 for file 1, variable 1 for file 2, ...], [[variable 2 for file 1, variable 2 for file 2, ...], ...]
LineColors = [[ROOT.kOrange, ROOT.kYellow+1],  [ROOT.kCyan-4, ROOT.kMagenta], [ROOT.kCyan, ROOT.kAzure+1]]
MarkerColors = [[ROOT.kGreen+1, ROOT.kGreen+1],  [ROOT.kGreen+1, ROOT.kGreen+1], [ROOT.kCyan, ROOT.kAzure+1]]
norm = 1

'''
Colors to choose from
ROOT.kOrange, ROOT.kYellow+1,  ROOT.kCyan-4, ROOT.kMagenta, ROOT.kYellow, ROOT.kYellow-4, 
ROOT.kYellow-6, ROOT.kYellow-7, ROOT.kYellow-9,
ROOT.kYellow-10,  ROOT.kOrange-2,  ROOT.kOrange-3, ROOT.kOrange-4, ROOT.kOrange-5,
ROOT.kCyan+1, ROOT.kCyan,ROOT.kCyan-6, ROOT.kCyan-7, ROOT.kCyan-9,
ROOT.kCyan-10, ROOT.kAzure+1, ROOT.kAzure+2,  ROOT.kAzure+6, ROOT.kAzure+6, ROOT.kAzure+8,
ROOT.kMagenta+1,  ROOT.kMagenta-4, ROOT.kMagenta-6, ROOT.kMagenta-7, ROOT.kMagenta-9,
ROOT.kMagenta-10, ROOT.kPink+2, ROOT.kPink+5,  ROOT.kPink+6, ROOT.kPink+7, ROOT.kPink+8,
ROOT.kGreen+1, ROOT.kGreen, ROOT.kGreen-4, ROOT.kGreen-6, ROOT.kGreen-7, ROOT.kGreen-9,
ROOT.kGreen-10, ROOT.kSpring-1, ROOT.kSpring-2,  ROOT.kSpring-3, ROOT.kSpring-4, ROOT.kSpring-5
ROOT.kAzure+2,  ROOT.kAzure+6,  ROOT.kAzure+8,
ROOT.kGreen+1, ROOT.kGreen, ROOT.kGreen-6,
'''


for a in range(0, len(Data)):
  c4 = TCanvas("c4","Signal MC Models", 1000, 800) 
  Legend = TLegend(0.9,0.9,0.7,0.7, "Legend");
  Legend.SetHeader("Legend")
  Legend.SetBorderSize(10)
  Legend.SetTextSize(0.03) 
  Legend.SetTextFont(62)
  Legend.SetFillColor(0);
  for b in range(0, len(Files)):
    f1[b] = TFile(Files[b]);
    f2[b] = f1[b].GetDirectory("GenAnalyzer");
    Tree[b] = TTree();
    f2[b].GetObject("emJetTree", Tree[b])
    Histograms[b] = TH1F(HistogramNames[a][b], HistogramTitles[a], Bins[a], X_Low[a], X_High[a])
    Histograms[b].SetTitle(HistogramTitles[a])
    Histograms[b].SetAxisRange(Y_Low[a], Y_High[a],"Y")
    Histograms[b].SetLineColor(LineColors[a][b]);
    #Histograms[b].SetFillColor(LineColors[a][b])
    Histograms[b].SetMarkerSize(1)
    #Histograms[b].SetMarkerStyle(8)
    Histograms[b].SetMarkerColor(MarkerColors[a][b])
    Histograms[b].GetYaxis().CenterTitle()
    Histograms[b].GetXaxis().SetTitle(X_Label[a])
    Histograms[b].GetYaxis().SetTitle(Y_Label[a])
    Histograms[b].SetLineWidth(4)
    Tree[b].Draw(Data[a] + ">>" + HistogramNames[a][b])
    Legend.AddEntry(Histograms[b], HistogramNames[a][b] , "f") 
    #TPaveStats StatsBox = (TPaveStats*)Histograms[b].GetListOfFunctions().FindObject("stats")
    #StatsBox.SetX1NDC(100)
    #StatsBox.SetX2NDC(200)
    #norm = Histograms[b].GetEntries()
    #Histograms[b].Scale(1/norm)
    #Histograms[b].SetAxisRange(0, 1,"Y")

  for d in range(0, len(Files)):
    if d == 0:
      Histograms[0].DrawNormalized();
      #Histograms[0].DrawNormalized();
    else:
      Histograms[d].DrawNormalized("SAMES");
      #Histograms[d].DrawNormalized("SAMES");
    Legend.Draw("SAME");
  #latex = TLatex()
  #latex.SetNDC()
  #latex.SetTextSize(0.04)
  #latex.SetTextAlign(31) 
  #latex.DrawLatex(0.60, 0.95, HistogramTitles[a]);
  #if (LogPlot[a] == 1):
    #c4.SetLogy() 
  #Histograms[b].SetMaximum(10)
  #Histograms[b].SetMinimum(0.1)
  c4.SaveAs(ImageNames[a] + ".gif")
  c4.SaveAs(ImageNames[a] + ".root")
  del c4

