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
h_MET = ROOT.TH1F("","",100,0,600)
b_MET1 = ROOT.TH1F("","",100,0,600)
b_MET2 = ROOT.TH1F("","",100,0,600)
h_jPT = ROOT.TH1F("","",100,0,600)
b_jPT1 = ROOT.TH1F("","",100,0,600)
b_jPT2 = ROOT.TH1F("","",100,0,600)
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

  # missing transverse energy (MET)
  MET=event.MET #to make a new event to analyze, simply add a new variable and call for event.[insert variable here]
  jPT=event.j1PT

  # require event to have some properties.
  if event.mjj > 100:
    # if the event passes, fill the histogram.
    h_MET.Fill(MET,weight)
    h_jPT.Fill(jPT,weight)


entry2 = 0
for event2 in chain2:
  entry2 += 1

  if ( entry2 != 0 and entry2%10000 == 0 ):
    print ("%d events processed" % entry2)
    sys.stdout.flush()

  # this is how we know how much this event "counts" when looking at large collections of events.
  weight2=event2.weight/numFiles2

  # missing transverse energy (MET)
  bMET1=event2.MET
  bjPT1=event2.j1PT

  # require event to have some properties.
  if event2.mjj > 100:
    # if the event passes, fill the histogram.
    b_MET1.Fill(bMET1,weight2)
    b_jPT1.Fill(bjPT1,weight2)
    
entry3 = 0
for event3 in chain3:
  entry3 += 1

  if ( entry3 != 0 and entry3%10000 == 0 ):
    print ("%d events processed" % entry3)
    sys.stdout.flush()

  # this is how we know how much this event "counts" when looking at large collections of events.
  weight3=event3.weight/numFiles3

  # missing transverse energy (MET)
  bMET2=event3.MET
  bjPT2=event3.j1PT

  # require event to have some properties.
  if event3.mjj > 100:
    # if the event passes, fill the histogram.
    b_MET2.Fill(bMET2,weight3)
    b_jPT2.Fill(bjPT2,weight3)
    




b_MET1.SetStats(0)
b_MET2.SetStats(0)
h_MET.SetStats(0)
b_jPT1.SetStats(0)
b_jPT2.SetStats(0)
h_jPT.SetStats(0)



# write the histogram to the output file.
###MET comparison Canvas
c1 = ROOT.TCanvas( 'c1', 'MET SvB') 

METBackgroundsum = b_MET1+b_MET2 #you can add histograms like this
METBackgroundsum.SetFillColor(ROOT.kOrange) #set histogram colors
METBackgroundsum.SetLineColor(ROOT.kOrange)
h_MET.SetFillColor(ROOT.kBlue)
h_MET.SetLineColor(ROOT.kBlue)

METBackgroundsum.Draw("SAME C") #Draw command adds to canvas, "SAME" does not overwrite the canvas
h_MET.Draw("SAME C")

c1.SetLogy(ROOT.kTRUE) #sets graph to logarithmic scale

L1 = ROOT.TLegend(0.65,0.65,0.9,0.9) #creates legend, numbers are fractions of screen.
L1.AddEntry(h_MET,"h_MET","l") #add entry to legend, first entry is histogram to match to, 2nd is title, 3rd is specification, in this case l=line
L1.AddEntry(METBackgroundsum,"Backgroundsum(MET)","l")
L1.Draw("SAME") 

c1.SetTitle("MET BvS")

c1.Update() #update command refreshes and sets the current changes
c1.Write() #write command saves to root file

###jPT comparison Canvas
c2 = ROOT.TCanvas( 'c2', 'jPT SvB')

jPTBackgroundsum = b_jPT1+b_jPT2
jPTBackgroundsum.SetFillColor(ROOT.kOrange)
jPTBackgroundsum.SetLineColor(ROOT.kOrange)
h_jPT.SetFillColor(ROOT.kBlue)
h_jPT.SetLineColor(ROOT.kBlue)

jPTBackgroundsum.Draw("SAME C")
h_jPT.Draw("SAME C")

c2.SetLogy(ROOT.kTRUE)

L2 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L2.AddEntry(h_jPT,"h_jPT","l")
L2.AddEntry(jPTBackgroundsum,"Backgroundsum(jPT)","l")
L2.Draw("SAME")

c2.Update()
c2.Write()

###BackgroundOnly (MET)
c3 = ROOT.TCanvas( 'c3', 'METBackground Only')

b_MET1.SetLineColor(ROOT.kRed)
b_MET2.SetLineColor(ROOT.kGreen)

METBackgroundsum.Draw("SAME C")
b_MET1.Draw("SAME")
b_MET2.Draw("SAME")


c3.SetLogy(ROOT.kTRUE)

L3 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L3.AddEntry(METBackgroundsum,"METbackgroundsum","l")
L3.AddEntry(b_MET1,"METbackground1","l")
L3.AddEntry(b_MET2,"METbackground2","l")
L3.Draw("SAME")

c3.Update()
c3.Write()
###BackgroundOnly (jpt)
c4 = ROOT.TCanvas( 'c4', 'jPTBackground Only')

b_jPT1.SetLineColor(ROOT.kRed)
b_jPT2.SetLineColor(ROOT.kGreen)

jPTBackgroundsum.Draw("SAME C")
b_jPT1.Draw("SAME")
b_jPT2.Draw("SAME")

c4.SetLogy(ROOT.kTRUE)

L4 = ROOT.TLegend(0.65,0.65,0.9,0.9)
L4.AddEntry(jPTBackgroundsum,"jPTbackgroundsum","l")
L4.AddEntry(b_jPT1,"jPTbackground1","l")
L4.AddEntry(b_jPT2,"jPTbackground2","l")
L4.Draw("SAME")

c4.Update()
c4.Write()

#TLorentz Vector graphs















# close the output file.
outfile.Close()
  
print ("Done!")
