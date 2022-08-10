import parsecsv
import sequtils
import strutils
import stats
import math
import os
import plotly
import parseopt
import algorithm

var version=0.1
var tme=newSeq[float](0)
var sig=newSeq[float](0)
var y=newSeq[float](0)
var count= 0

var fromF=0.0
var toF=0.0
var ofac=1.0
var alpha=0.05
var res: tuple[freq:seq[float], pn:seq[float], level: float, p: float, peak: float, peakAt: float]
var filename:string
var output=1
var save=0
var plotit=0
var timesvar="t"
var signalvar="y"

# Call possibilities:
#./lomb data.csv --times=t --signal=y --from=0.5 --to=1.5 --ofac=5 --alpha=0.001 -p -s
# Anything starting with - or -- is optional, only filname is required (no - or --signs)
# -p : plot, -s :save 
# Default analysis is of columns t and y

var pc = initOptParser()
while true:
  pc.next()
  case pc.kind
  of cmdEnd: break
  of cmdShortOption, cmdLongOption:
    if pc.val == "":
      if pc.key=="p":
          plotit=1
      if pc.key=="s":
          save=1
      if pc.key=="v":         
         echo ("Version:  0.1 by Thomas Ruf; Vetmeduni Vienna Austria; GPL 3.0")   
         echo("  ")
      if pc.key=="h":
          echo "-p plot"
          echo "-s save results"
          echo "-h help"
          echo "-v Version"
          echo "--from=start fequency"
          echo "--to=maximum fequency"
          echo "--ofac=oversampling factor"
          echo "--alpha=significance threshold"
          echo "--times=name of time column"
          echo "--signal=name of signal column"
          echo ""
  
    else:
      if pc.key=="from": 
          fromF=parseFloat(pc.val)
      if pc.key=="to":
          toF=parseFloat(pc.val)
      if pc.key=="ofac":
          ofac=parseFloat(pc.val)
      if pc.key=="alpha":
          alpha=parseFloat(pc.val)        
      if pc.key=="times":
          timesvar=pc.val
      if pc.key=="signal":
          signalvar=pc.val          
  of cmdArgument:
    filename=pc.key
    


if ofac<1.0:
     ofac=1.0
     
if filename=="":
    quit("No Filename provided, sorry")
if fromF > toF:
    quit ("from > to. Nothing to do.")
    
var p: CsvParser
p.open(filename) 
p.readHeaderRow()
    

while p.readRow():
  #for col in items(p.headers):
      #echo "##", col, ":", p.rowEntry(col), "##"
  count=count+1
  tme.add(parseFloat(p.rowEntry(timesvar)))
  sig.add(parseFloat(p.rowEntry(signalvar)))  
p.close()


var mn=0.0
var sd=0.0
var rs:RunningStat
rs.push(sig)
mn= rs.mean()
sd=rs.standardDeviationS()
y=map(sig, proc(x:float): float=(x-mn)) # scale signal, remove mean

##############################################

proc ggamma (N:float):float=
    (sqrt(2.0 / N) * exp(lgamma(N / 2.0) - lgamma((N - 1.0) / 2.0)))    # helper function,copied from astropy

proc pbaluev (Z:float, fmax: float, tm:seq) :float=  #p value according to Baluev for standarized (0-1) power values
    #code copied from astropy timeseries
    var rs:RunningStat
    var N=0.0
    var Dt=0.0
    var NH=0.0
    var NK=0.0
    var fsingle=0.0
    var Teff=0.0
    var W=0.0
    var tau=0.0
    var m1=0.0
    var m2=0.0
    var p=0.0
    var sm=0.0
    
    N=float(tm.len)
    for a in tm:
        sm=sm+pow(a,2.0)
    m1=sm/(N)
    rs.push(tm)
    m2=rs.mean()
    m2=pow(m2,2.0)
    Dt=m1-m2
    #Dt=mean(tm^2.0)-mean(tm)^2.0
    NH=N-1
    NK=N-3
    fsingle=pow( (1.0 - Z) , (0.5 * NK) )
    Teff = sqrt(4.0 * PI * Dt) # Effective baseline
    W = fmax * Teff
    tau=ggamma(NH) * W * pow( (1.0 - Z) , (0.5 * (NK - 1.0)) )*sqrt(0.5 * NH * Z)
    p= -(exp(-tau)-1.0) + fsingle * exp(-tau)
    return(p)

proc getZ(guess:float, alpha:float,fmax:float, tm:seq):float= # get Significance level for certain alpha
    var fib= @    [ 1,1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946]
    var now=0.0
    var dir =1.0
    var here=guess
    
    for i in countdown(20,0):
        now=pbaluev(here,fmax,tm)
        if now < alpha:
            dir= -1.0
        else:
            dir= 1.0
        here=float(here+float(fib[i])/10000.0*dir) 
        if here<0.0: here=0.0
        if here>1.0: here=1.0

    return(here)


##############################################
proc calc (t:seq, y:seq, fromF:float, toF:float, ofac:float, alpha:float): tuple = # main proc to compute LS
  var x=newSeq[float](0)
  var w=newSeq[float](0)
  var pn=newSeq[float](0)
  var ss=newSeq[float](0)
  var cc=newSeq[float](0)
  var arg=newSeq[float](0)
  var cs=newSeq[float](0)
  var sn=newSeq[float](0)
  var aa=newSeq[float](0)
  
  var n:int
  var nf=1
  var nout=0
  var peak=0.0
  var peakindex=0
  var peakAt=0.0 
  var tspan=0.0
  var frd=0.0
  var step=0.0
  var fmax=0.0
  var top=0.0
  var norm=1.0
  var wi=0.0
  var tau=0.0
  var A=0.0
  var B=0.0
  var C=0.0
  var D=0.0
  var p=0.0
  var level=0.0
  var sm=0.0
  
  n=t.len
  tspan=t[n-1]-t[0]
  frd=1/tspan
  step=1/(tspan*ofac)
  
  if toF==0.0:
       fmax=floor(0.5*float(n)*ofac)*step
  else:
       fmax=toF  

  var freq=newSeq[float](0)
  
  top=frd
  while(top<=fmax):
      freq.add(top)
      top+=step
     

  if fromF>0.0:
         for i in 1..freq.len:
               if freq[i-1]>=fromF:
                        nf=i
                        break
         freq=freq[nf-1..freq.len-1]
 
  nout=freq.len
  x=map(t, proc (a:float): float=a*2.0*PI)
  
  sm=0.0
  for a in y:
      sm=sm+pow(a,2.0)
  norm=1/sm
  
  w=map(freq,proc(a:float): float=a*2.0*PI)
  
  for i in 0..nout-1:
      wi=w[i]
      ss=map(t, proc (a:float): float=sin(a*wi) )
      cc=map(t, proc (a:float): float=cos(a*wi) )
      tau=0.5*arctan2(ss.sum,cc.sum)/wi
      arg=map(t, proc(a:float): float=(a-tau)*wi)
      cs.setLen(0)
      sn.setLen(0)
      for a in arg:
          cs.add(cos(a))
          sn.add(sin(a))
      aa.setLen(0)
      for j in 1.. y.len:
          aa.add(y[j-1]*cs[j-1])
      A=pow(sum(aa),2)
      aa.setLen(0)
      for a in cs:
          aa.add(a*a)
      B=sum(aa)
      aa.setLen(0)
      for j in 1.. y.len:
          aa.add(y[j-1]*sn[j-1])
      C=pow(sum(aa),2)
      aa.setLen(0)
      for a in sn:
          aa.add(a*a)
      D=sum(aa)
      pn.add(A/B+C/D)
      
  for i in 1..pn.len:
      pn[i-1]=pn[i-1]*norm
  peak=max(pn)
  for i in 1..pn.len:
      if pn[i-1]==peak:
          peakindex=i-1
          break
  peakAT=freq[peakindex]
  p=pbaluev(peak, fmax, t)
  level=getZ(0.01,alpha,fmax,t)

  res=(freq,pn,level,p,peak,peakAt)
  return res

################################
res=calc(tme,y, fromF,toF,ofac, alpha)
  
if output>0:
      echo "Peak at    :  ", res.peakAt,"   ",1.0/res.peakAt
      echo "Peak height : ", res.peak
      echo "P-value     : ", res.p
      echo "Sign. level : ", res.level
      echo "alpha       : ",alpha
      echo "N out       : ", res.freq.len

if save==1:   
    let fo=open("lomb_results.csv",fmWrite) 
    fo.writeLine("Variable , Value")
    fo.writeLine("Peak at ,",res.peakAt)
    fo.writeLine("Peak height ,",res.peak)
    fo.writeLine("P-value ,", res.p)
    fo.writeLine("Significance level,", res.level)    
    fo.writeLine("alpha ,",alpha)
    fo.writeLine("N out ,",res.freq.len)

    for line in res.freq:
          fo.writeLine("Frequency ,",line)

    for line in res.pn:
         fo.writeLine("Power , ",line)   
             
    fo.flushFile()
    fo.close()

if output>0:
    echo "done"
    
if plotit==1:
    let
        d = Trace[float](mode: PlotMode.Lines, `type`: PlotType.Scatter, name:"Power")
        d1 = Trace[float](mode: PlotMode.Lines, `type`: PlotType.Scatter, name:"Sig. level")

    d.xs = res.freq
    d.ys = res.pn
    d1.xs = res.freq

    var pp=newSeq[float](0)
    for i in 1..res.pn.len:
        pp.add(res.level)

    d1.ys=pp

    let
            layout = Layout(title: "Lomb-Scargle Periodogram", width: 800, height: 400,
            xaxis: Axis(title:"Frequency"),
            yaxis: Axis(title: "Normalized power"),
            autosize: false)
            
            pl = Plot[float](layout: layout, traces: @[d,d1])
  

    #echo pl.save()
    pl.show()
