#!/usr/bin/env python

import sys
import ROOT
import math
import array
import argparse
import os

parser = argparse.ArgumentParser() #basically this section allows us to set command line prompts and inputs, here the files we want to interact with.
parser.add_argument('--input', action='store', default="input.txt") #example, this specifies a new command line command --input [] where instead of brackets we place the input file of our choosing.
parser.add_argument('--inputTree', action='store', default="allev/hftree") # the other 2 bits store parts of those files
parser.add_argument('--background', action='store', default="background.txt")
parser.add_argument('--backgroundTree', action='store', default="allev/hftree")
parser.add_argument('--background2', action='store', default="background2.txt")
parser.add_argument('--backgroundTree2', action='store', default="allev/hftree")
parser.add_argument('--output', action='store', default="hist.root")
parser.add_argument('--lumi',action='store', default=1000.) # 1 fb-1
parser.add_argument('--debug',action='store_true')
args=parser.parse_args()




chain = ROOT.TChain(args.inputTree) #these chains allow us to run over the values in our Ttrees, which is where the data is stored, each file gets a tree, which contains all the information
chain2 = ROOT.TChain(args.backgroundTree) ##create second and third chain analysis for background info
chain3 = ROOT.TChain(args.backgroundTree2)

if ( args.input[-4:] == 'root' ):  #input file loop to be compared
  print ( "Running over single root file:" )
  print ( "   > %s" % args.input )
  chain.Add(args.input) #the meat of the loop, put all the input arguments into a the chain.
else:
  print ( "Running over list of root files:" )
  for line in open(args.input):
    print ("   > " + line.rstrip('\n'))
    chain.Add(line.rstrip('\n'))

if ( args.background[-4:] == 'root' ):   ##iteration for background file 1
  print ( "Running over single root file:" )
  print ( "   > %s" % args.background )
  chain2.Add(args.background)
else:
  print ( "Running over list of root files:" )
  for line in open(args.background):
    print ("   > " + line.rstrip('\n'))
    chain2.Add(line.rstrip('\n'))
    
if ( args.background2[-4:] == 'root' ):   ##iteration for background file 2
  print ( "Running over single root file:" )
  print ( "   > %s" % args.background2 )
  chain3.Add(args.background2)
else:
  print ( "Running over list of root files:" )
  for line in open(args.background2):
    print ("   > " + line.rstrip('\n'))
    chain3.Add(line.rstrip('\n'))
    
numFiles=chain.GetNtrees()  
numFiles2=chain2.GetNtrees()
numFiles3=chain3.GetNtrees()
print ( "Loaded %s chains..." % numFiles )

# Prevent the canvas from displaying
ROOT.gROOT.SetBatch(True)

# a histogram for our output
outfile=ROOT.TFile.Open(args.output,"RECREATE")




# make a histogram
h_MET = ROOT.TH1F("","",100,0,1500)
b_MET1 = ROOT.TH1F("","",100,0,1500)
b_MET2 = ROOT.TH1F("","",100,0,1500)
h_jPT = ROOT.TH1F("","",100,0,1500)
b_jPT1 = ROOT.TH1F("","",100,0,1500)
b_jPT2 = ROOT.TH1F("","",100,0,1500)

h_Mdi = ROOT.TH1F("","",100,0,15000)
b_Mdi1 = ROOT.TH1F("","",100,0,15000)
b_Mdi2 = ROOT.TH1F("","",100,0,15000)
h_Etadi = ROOT.TH1F("","",100,-10,10)
b_Etadi1 = ROOT.TH1F("","",100,-10,10)
b_Etadi2 = ROOT.TH1F("","",100,-10,10)
h_Phidi = ROOT.TH1F("","",100,-10,10)
b_Phidi1 = ROOT.TH1F("","",100,-10,10)
b_Phidi2 = ROOT.TH1F("","",100,-10,10)
h_Ptdi = ROOT.TH1F("","",100,0,1500)
b_Ptdi1 = ROOT.TH1F("","",100,0,1500)
b_Ptdi2 = ROOT.TH1F("","",100,0,1500)

### Loop through all events in chain


###Loops for MET (Missing Transverse Energy) and filling those histograms.
entry = 0
for event in chain:
  entry += 1

  if ( entry != 0 and entry%10000 == 0 ):
    print ("%d events processed" % entry)
    sys.stdout.flush()

  # this is how we know how much this event "counts" when looking at large collections of events.
  weight=event.weight/numFiles
  jet1_tlv = ROOT.TLorentzVector()
  jet1_tlv.SetPtEtaPhiM(event.j1PT,event.j1Eta,event.j1Phi,0)

  jet2_tlv = ROOT.TLorentzVector()
  jet2_tlv.SetPtEtaPhiM(event.j2PT,event.j2Eta,event.j2Phi,0)


  dijet_tlv = jet1_tlv+jet2_tlv

  # missing transverse energy (MET)
  MET=event.MET #to make a new event to analyze, simply add a new variable and call for event.[insert variable here]
  jPT=event.j1PT
  Mdi = dijet_tlv.M() 
  Etadi = dijet_tlv.Eta()
  Phidi = dijet_tlv.Phi()
  Ptdi = dijet_tlv.Pt()

  # require event to have some properties.
  if event.mjj > 100:
    # if the event passes, fill the histogram.
    h_MET.Fill(MET,weight)
    h_jPT.Fill(jPT,weight)
    h_Mdi.Fill(Mdi,weight)
    h_Etadi.Fill(Etadi,weight)
    h_Phidi.Fill(Phidi,weight)
    h_Ptdi.Fill(Ptdi,weight)


entry2 = 0
for event2 in chain2:
  entry2 += 1

  if ( entry2 != 0 and entry2%10000 == 0 ):
    print ("%d events processed" % entry2)
    sys.stdout.flush()

  # this is how we know how much this event "counts" when looking at large collections of events.
  weight2=event2.weight/numFiles2
  jet1_tlv = ROOT.TLorentzVector()
  jet1_tlv.SetPtEtaPhiM(event2.j1PT,event2.j1Eta,event2.j1Phi,0)

  jet2_tlv = ROOT.TLorentzVector()
  jet2_tlv.SetPtEtaPhiM(event2.j2PT,event2.j2Eta,event2.j2Phi,0)


  bdijet_tlv1 = jet1_tlv+jet2_tlv

  # missing transverse energy (MET)
  bMET1=event2.MET
  bjPT1=event2.j1PT
  bMdi1 = bdijet_tlv1.M() 
  bEtadi1 = bdijet_tlv1.Eta()
  bPhidi1 = bdijet_tlv1.Phi()
  bPtdi1 = bdijet_tlv1.Pt()

  # require event to have some properties.
  if event2.mjj > 100:
    # if the event passes, fill the histogram.
    b_MET1.Fill(bMET1,weight2)
    b_jPT1.Fill(bjPT1,weight2)
    b_Mdi1.Fill(bMdi1,weight2)
    b_Etadi1.Fill(bEtadi1,weight2)
    b_Phidi1.Fill(bPhidi1,weight2)
    b_Ptdi1.Fill(bPtdi1,weight2)
    
entry3 = 0
for event3 in chain3:
  entry3 += 1

  if ( entry3 != 0 and entry3%10000 == 0 ):
    print ("%d events processed" % entry3)
    sys.stdout.flush()

  # this is how we know how much this event "counts" when looking at large collections of events.
  weight3=event3.weight/numFiles3
  jet1_tlv = ROOT.TLorentzVector()
  jet1_tlv.SetPtEtaPhiM(event3.j1PT,event3.j1Eta,event3.j1Phi,0)

  jet2_tlv = ROOT.TLorentzVector()
  jet2_tlv.SetPtEtaPhiM(event3.j2PT,event3.j2Eta,event3.j2Phi,0)


  bdijet_tlv2 = jet1_tlv+jet2_tlv

  # missing transverse energy (MET)
  bMET2=event3.MET
  bjPT2=event3.j1PT
  bMdi2 = bdijet_tlv2.M() 
  bEtadi2 = bdijet_tlv2.Eta()
  bPhidi2 = bdijet_tlv2.Phi()
  bPtdi2 = bdijet_tlv2.Pt()

  # require event to have some properties.
  if event3.mjj > 100:
    # if the event passes, fill the histogram.
    b_MET2.Fill(bMET2,weight3)
    b_jPT2.Fill(bjPT2,weight3)
    b_Mdi2.Fill(bMdi2,weight3)
    b_Etadi2.Fill(bEtadi2,weight3)
    b_Phidi2.Fill(bPhidi2,weight3)
    b_Ptdi2.Fill(bPtdi2,weight3)
    




b_MET1.SetStats(0)
b_MET2.SetStats(0)
h_MET.SetStats(0)
b_jPT1.SetStats(0)
b_jPT2.SetStats(0)
h_jPT.SetStats(0)
b_Mdi1.SetStats(0)
b_Mdi2.SetStats(0)
h_Mdi.SetStats(0)
b_Etadi1.SetStats(0)
b_Etadi2.SetStats(0)
h_Etadi.SetStats(0)
b_Phidi1.SetStats(0)
b_Phidi2.SetStats(0)
h_Phidi.SetStats(0)
b_Ptdi1.SetStats(0)
b_Ptdi2.SetStats(0)
h_Ptdi.SetStats(0)


# write the histogram to the output file.
###MET comparison Canvas
c1 = ROOT.TCanvas( 'c1', 'MET SvB') 

METBackgroundsum = b_MET1+b_MET2 #you can add histograms like this
METBackgroundsum.SetLineColor(ROOT.kOrange) #set histogram colors
h_MET.SetLineColor(ROOT.kBlue)
b_MET1.SetLineColor(ROOT.kRed)
b_MET2.SetLineColor(ROOT.kGreen)

METBackgroundsum.SetAxisRange(1,10**7,"Y")
h_MET.SetAxisRange(1,10**7,"Y")
b_MET1.SetAxisRange(1,10**7,"Y")
b_MET2.SetAxisRange(1,10**7,"Y")

METBackgroundsum.Draw("SAME C") #Draw command adds to canvas, "SAME" does not overwrite the canvas
h_MET.Draw("SAME C")
b_MET1.Draw("SAME")
b_MET2.Draw("SAME")

c1.SetLogy(ROOT.kTRUE) #sets graph to logarithmic scale

L1 = ROOT.TLegend(0.65,0.65,0.9,0.9) #creates legend, numbers are fractions of screen.
L1.AddEntry(h_MET,"h_MET","l") #add entry to legend, first entry is histogram to match to, 2nd is title, 3rd is specification, in this case l=line
L1.AddEntry(METBackgroundsum,"Backgroundsum(MET)","l")
L1.AddEntry(b_MET1,"METbackground1","l")
L1.AddEntry(b_MET2,"METbackground2","l")
L1.Draw("SAME") 

c1.SetTitle("MET BvS")

c1.Update() #update command refreshes and sets the current changes
c1.Write() #write command saves to root file

###jPT comparison Canvas
c2 = ROOT.TCanvas( 'c2', 'jPT SvB')

jPTBackgroundsum = b_jPT1+b_jPT2
jPTBackgroundsum.SetLineColor(ROOT.kOrange)
h_jPT.SetLineColor(ROOT.kBlue)
b_jPT1.SetLineColor(ROOT.kRed)
b_jPT2.SetLineColor(ROOT.kGreen)

jPTBackgroundsum.SetAxisRange(1,1500,"X")
h_jPT.SetAxisRange(1,1500,"X")
b_jPT1.SetAxisRange(1,1500,"X")
b_jPT2.SetAxisRange(1,1500,"X")

jPTBackgroundsum.SetAxisRange(1,10**7,"Y")
h_jPT.SetAxisRange(1,10**7,"Y")
b_jPT1.SetAxisRange(1,10**7,"Y")
b_jPT2.SetAxisRange(1,10**7,"Y")

jPTBackgroundsum.Draw("SAME C")
h_jPT.Draw("SAME C")
b_jPT1.Draw("SAME")
b_jPT2.Draw("SAME")

c2.SetLogy(ROOT.kTRUE)

L2 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L2.AddEntry(h_jPT,"h_jPT","l")
L2.AddEntry(jPTBackgroundsum,"Backgroundsum(jPT)","l")
L2.AddEntry(b_jPT1,"jPTbackground1","l")
L2.AddEntry(b_jPT2,"jPTbackground2","l")
L2.Draw("SAME")

c2.Update()
c2.Write()

#TLorentz Vector graphs
##Mdi Comparison Canvas
c3 = ROOT.TCanvas( 'c3', 'M SvB')

MdiBackgroundsum = b_Mdi1+b_Mdi2
MdiBackgroundsum.SetLineColor(ROOT.kOrange)
h_Mdi.SetLineColor(ROOT.kBlue)
b_Mdi1.SetLineColor(ROOT.kRed)
b_Mdi2.SetLineColor(ROOT.kGreen)

MdiBackgroundsum.SetAxisRange(1,15000,"X")

MdiBackgroundsum.Draw("SAME C")
h_Mdi.Draw("SAME C")
b_Mdi1.Draw("SAME")
b_Mdi2.Draw("SAME")

c3.SetLogy(ROOT.kTRUE)

L3 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L3.AddEntry(h_Mdi,"Mdi","l")
L3.AddEntry(MdiBackgroundsum,"Backgroundsum(Mdi)","l")
L3.AddEntry(b_Mdi1,"Mdibackground1","l")
L3.AddEntry(b_Mdi2,"Mdibackground2","l")
L3.Draw("SAME")

c3.Update()
c3.Write()
##Eta Comparison Canvas
c4 = ROOT.TCanvas( 'c4', 'Eta SvB')

EtadiBackgroundsum = b_Etadi1+b_Etadi2
EtadiBackgroundsum.SetLineColor(ROOT.kOrange)
h_Etadi.SetLineColor(ROOT.kBlue)
b_Etadi1.SetLineColor(ROOT.kRed)
b_Etadi2.SetLineColor(ROOT.kGreen)



EtadiBackgroundsum.Draw("SAME C")
h_Etadi.Draw("SAME C")
b_Etadi1.Draw("SAME")
b_Etadi2.Draw("SAME")

c4.SetLogy(ROOT.kTRUE)

L4 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L4.AddEntry(h_Etadi,"Etadi","l")
L4.AddEntry(EtadiBackgroundsum,"Backgroundsum(Etadi)","l")
L4.AddEntry(b_Etadi1,"Etadibackground1","l")
L4.AddEntry(b_Etadi2,"Etadibackground2","l")
L4.Draw("SAME")

c4.Update()
c4.Write()
##PhiDi comparison canvas
c5 = ROOT.TCanvas( 'c5', 'Phidi SvB')

PhidiBackgroundsum = b_Phidi1+b_Phidi2
PhidiBackgroundsum.SetLineColor(ROOT.kOrange)
h_Phidi.SetLineColor(ROOT.kBlue)
b_Phidi1.SetLineColor(ROOT.kRed)
b_Phidi2.SetLineColor(ROOT.kGreen)

Phidiackgroundsum.SetAxisRange(1,10**7,"Y")
h_Phidi.SetAxisRange(1,10**7,"Y")
b_Phidi1.SetAxisRange(1,10**7,"Y")
b_Phidi2.SetAxisRange(1,10**7,"Y")

PhidiBackgroundsum.Draw("SAME C")
h_Phidi.Draw("SAME C")
b_Phidi1.Draw("SAME")
b_Phidi2.Draw("SAME")

c5.SetLogy(ROOT.kTRUE)

L5 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L5.AddEntry(h_Phidi,"Phidi","l")
L5.AddEntry(PhidiBackgroundsum,"Backgroundsum(Phi)","l")
L5.AddEntry(b_Phidi1,"Phibackground1","l")
L5.AddEntry(b_Phidi2,"Phibackground2","l")
L5.Draw("SAME")

c5.Update()
c5.Write()
##P5di Comparison Canvas
c6 = ROOT.TCanvas( 'c6', 'Ptdi SvB')

PtdiBackgroundsum = b_Ptdi1+b_Ptdi2
PtdiBackgroundsum.SetLineColor(ROOT.kOrange)
h_Ptdi.SetLineColor(ROOT.kBlue)
b_Ptdi1.SetLineColor(ROOT.kRed)
b_Ptdi2.SetLineColor(ROOT.kGreen)



PtdiBackgroundsum.Draw("SAME C")
h_Ptdi.Draw("SAME C")
b_Ptdi1.Draw("SAME")
b_Ptdi2.Draw("SAME")

PtdiBackgroundsum.SetAxisRange(1,10**7,"Y")
h_Ptdi.SetAxisRange(1,10**7,"Y")
b_Ptdi1.SetAxisRange(1,10**7,"Y")
b_Ptdi2.SetAxisRange(1,10**7,"Y")

c6.SetLogy(ROOT.kTRUE)

L6 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L6.AddEntry(h_Ptdi,"Ptdi","l")
L6.AddEntry(PtdiBackgroundsum,"Backgroundsum(Ptdi)","l")
L6.AddEntry(b_Ptdi1,"Ptdibackground1","l")
L6.AddEntry(b_Ptdi2,"Ptdibackground2","l")
L6.Draw("SAME")

c6.Update()
c6.Write()

    # if the event passes, fill the histogram.

 




#h_Mdi.SetLineColor(ROOT.kBlue)
#h_Mdi.Draw("SAME")
#h_Etadi.SetLineColor(ROOT.kRed)
#h_Etadi.Draw("SAME")
#h_Phidi.SetLineColor(ROOT.kGreen)
#h_Phidi.Draw("SAME")
#h_Ptdi.SetLineColor(ROOT.kOrange)
#h_Ptdi.Draw("SAME")

#c3.Update()
#c3.Write()











# close the output file.
outfile.Close()
  
print ("Done!")
