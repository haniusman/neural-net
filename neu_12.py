import csv
import random
import math
#import numpy as np

with open("dataX.csv") as csvfile:
    readCSV =csv.reader(csvfile,delimiter=',') 
    dataX = []
    for row in readCSV:
        date=row[0]
        dataX.append(date)
        
#print(dataX[1])


with open("dataY.csv") as csvfile:
    readCSV =csv.reader(csvfile,delimiter=',') 
    dataY = []
    for row in readCSV:
        date=row[0]
        dataY.append(date)
#print(dataY)

with open("dataZ.csv") as csvfile:
    readCSV =csv.reader(csvfile,delimiter=',') 
    dataZ = []
    for row in readCSV:
        date=row[0]
        dataZ.append(date)
#print(dataZ)

with open("dataT.csv") as csvfile:
    readCSV =csv.reader(csvfile,delimiter=',') 
    dataT = []
    for row in readCSV:
        date=row[0]
        dataT.append(date)
#print(dataT)]

def rdm():
    w = round(random.uniform(-1,1), 2)
    return w;

def output(x,y,z,w1,w2,w3):
    op = (x*w1) + (y*w2) + (z*w3)
    op = op*(-1)
    sig = 1/(1+math.exp(op))
    return sig;
    
def output_k(outa, outb, outc, outd,oute,outf, outg, outh, outi, outj, outm, outn, wa, wb, wc, wd,we, wf, wg, wh, wi,wj, wm,wn):
    op = (outa*wa) + (outb*wb) + (outc*wc) + (outd*wd) + (oute*we) + (outf*wf) +(outg*wg) +(outh*wh) + (outi*wi)+(outj*wj)+(outm*wm)+(outn*wn)
    sig = 1/(1+math.exp(-op))
    return sig;


def error(errk, w, out):
    err = (errk * w * out * (1-out))
    return err;

def delta(err, out):
    cw = lr * err * out
    return cw;

def scale(x):
    min = 4800
    max = 13000
    sx = (x - min)/(max-min)
    return sx;

lr = 0.3 #learning rate
#sf = 0.1
pvalue = []


for i in range(0,705):
    dx = float(dataX[i])
    dy = float(dataY[i])
    dz = float(dataZ[i])
    dt = float(dataT[i])

    dx = scale(dx)
    dy = scale(dy)
    dz = scale(dz)
    dt = scale(dt)
    
    #print("Inputs")
    #print (dx,dy,dz,dt)


    t1 = dt-(0.1*dt)
    t2 = dt+(0.1*dt)
    
    w1a = rdm() 
    w2a = rdm()  
    w3a = rdm()
    
    w1b = rdm()  
    w2b = rdm()  
    w3b = rdm()
    
    w1c = rdm()
    w2c = rdm()
    w3c = rdm()
    
    w1d = rdm()
    w2d = rdm()
    w3d = rdm()

    w1e = rdm()
    w2e = rdm()
    w3e = rdm()

    w1f = rdm()
    w2f = rdm()
    w3f = rdm()

    w1g = rdm()
    w2g = rdm()
    w3g = rdm()

    w1h = rdm()
    w2h = rdm()
    w3h = rdm()

    w1i = rdm()
    w2i = rdm()
    w3i = rdm()

    w1j = rdm()
    w2j = rdm()
    w3j = rdm()

    w1m = rdm()
    w2m = rdm()
    w3m = rdm()

    w1n = rdm()
    w2n = rdm()
    w3n = rdm()
    
    wak = rdm()
    wbk = rdm()
    wck = rdm()
    wdk = rdm()
    wek = rdm()
    wfk = rdm()
    wgk = rdm()
    whk = rdm()
    wik = rdm()
    wjk = rdm()
    wmk = rdm()
    wnk = rdm()
    
    oa = output(dx,dy,dz,w1a,w2a,w3a)
    ob = output(dx,dy,dz,w1b,w2b,w3b)
    oc = output(dx,dy,dz,w1c,w2c,w3c)
    od = output(dx,dy,dz,w1d,w2d,w3d)
    oe = output(dx,dy,dz,w1e,w2e,w3e)
    of = output(dx,dy,dz,w1f,w2f,w3f)
    og = output(dx,dy,dz,w1g,w2g,w3g)
    oh = output(dx,dy,dz,w1h,w2h,w3h)
    oi = output(dx,dy,dz,w1i,w2i,w3i)
    oj = output(dx,dy,dz,w1j,w2j,w3j)
    om = output(dx,dy,dz,w1m,w2m,w3m)
    on = output(dx,dy,dz,w1n,w2n,w3n)
    
    
    ok = output_k(oa,ob,oc,od,oe,of,og,oh,oi,oj,om,on,wak,wbk,wck,wdk,wek,wfk,wgk,whk,wik,wjk,wmk,wnk)
    #print("Output_k")
    #print (ok)

    count = 0
    while (ok < t1 or ok > t2):
        #print (count)
        errk = (dt-ok)*ok*(1-ok)
        
        erra = error(errk, wak, oa)
        errb = error(errk, wbk, ob)
        errc = error(errk, wck, oc)
        errd = error(errk, wdk, od)
        erre = error(errk, wek, oe)
        errf = error(errk, wfk, of)
        errg = error(errk, wgk, og)
        errh = error(errk, whk, oh)
        erri = error(errk, wik, oi)
        errj = error(errk, wjk, oj)
        errm = error(errk, wmk, om)
        errn = error(errk, wnk, on)
        
        #print("Error")
        #print(errk,erri,errj)

        cwak = delta(errk, oa)
        cwbk = delta(errk, ob)
        cwck = delta(errk, oc)
        cwdk = delta(errk, od)
        cwek = delta(errk, oe)
        cwfk = delta(errk, of)
        cwgk = delta(errk, og)
        cwhk = delta(errk, oh)
        cwik = delta(errk, oi)
        cwjk = delta(errk, oj)
        cwmk = delta(errk, om)
        cwnk = delta(errk, on)

        
        cw1a = delta(erra, dx)
        cw2a = delta(erra, dy)
        cw3a = delta(erra, dz)
        
        cw1b = delta(errb, dx)
        cw2b = delta(errb, dy)
        cw3b = delta(errb, dz)
        
        cw1c = delta(errc, dx)
        cw2c = delta(errc, dy)
        cw3c = delta(errc, dz)
        
        cw1d = delta(errd, dx)
        cw2d = delta(errd, dy)
        cw3d = delta(errd, dz)

        cw1e = delta(erre, dx)
        cw2e = delta(erre, dy)
        cw3e = delta(erre, dz)

        cw1f = delta(errf, dx)
        cw2f = delta(errf, dy)
        cw3f = delta(errf, dz)

        cw1g = delta(errg, dx)
        cw2g = delta(errg, dy)
        cw3g = delta(errg, dz)

        cw1h = delta(errh, dx)
        cw2h = delta(errh, dy)
        cw3h = delta(errh, dz)

        cw1i = delta(erri, dx)
        cw2i = delta(erri, dy)
        cw3i = delta(erri, dz)

        cw1j = delta(errj, dx)
        cw2j = delta(errj, dy)
        cw3j = delta(errj, dz)

        cw1m = delta(errm, dx)
        cw2m = delta(errm, dy)
        cw3m = delta(errm, dz)

        cw1n = delta(errn, dx)
        cw2n = delta(errn, dy)
        cw3n = delta(errn, dz)

        wak = wak + cwak
        wbk = wbk + cwbk
        wck = wck + cwck
        wdk = wdk + cwdk
        wek = wek + cwek
        wfk = wfk + cwfk
        wgk = wgk + cwgk
        whk = whk + cwhk
        wik = wik + cwik
        wjk = wjk + cwjk
        wmk = wmk + cwmk
        wnk = wnk + cwnk
        
        w1a = w1a + cw1a 
        w2a = w2a + cw2a  
        w3a = w3a + cw3a
        
        w1b = w1b + cw1b  
        w2b = w2b + cw2b
        w3b = w3b + cw3b
        
        w1c = w1c + cw1c
        w2c = w2c + cw2c
        w3c = w2c + cw3c
        
        w1d = w1d + cw1d
        w2d = w2d + cw2d
        w3d = w3d + cw3d

        w1e = w1e + cw1e
        w2e = w2e + cw2e
        w3e = w3e + cw3e

        w1f = w1f + cw1f
        w2f = w2f + cw2f
        w3f = w3f + cw3f

        w1g = w1g + cw1g
        w2g = w2g + cw2g
        w3g = w3g + cw3g

        w1h = w1h + cw1h
        w2h = w2h + cw2h
        w3h = w3h + cw3h

        w1i = w1i + cw1i
        w2i = w2i + cw2i
        w3i = w3i + cw3i

        w1j = w1j + cw1j
        w2j = w2j + cw2j
        w3j = w3j + cw3j

        w1m = w1m + cw1m
        w2m = w2m + cw2m
        w3m = w3m + cw3m

        w1n = w1n + cw1n
        w2n = w2n + cw2n
        w3n = w3n + cw3n


        oa = output(dx,dy,dz,w1a,w2a,w3a)
        ob = output(dx,dy,dz,w1b,w2b,w3b)
        oc = output(dx,dy,dz,w1c,w2c,w3c)
        od = output(dx,dy,dz,w1d,w2d,w3d)
        oe = output(dx,dy,dz,w1e,w2e,w3e)
        of = output(dx,dy,dz,w1f,w2f,w3f)
        og = output(dx,dy,dz,w1g,w2g,w3g)
        oh = output(dx,dy,dz,w1h,w2h,w3h)
        oi = output(dx,dy,dz,w1i,w2i,w3i)
        oj = output(dx,dy,dz,w1j,w2j,w3j)
        om = output(dx,dy,dz,w1m,w2m,w3m)
        on = output(dx,dy,dz,w1n,w2n,w3n)
    
    
        ok = output_k(oa,ob,oc,od,oe,of,og,oh,oi,oj,om,on,wak,wbk,wck,wdk,wek,wfk,wgk,whk,wik,wjk,wmk,wnk)
    
        
        count = count+1

    #print ("while Loop ended")
    ok = ok * (13000-4800)+ 4800
    #print(ok,i)
    
    pvalue.append(ok)


for j in range(0,len(pvalue)):
    print(round(pvalue[j],2),sep='\n')
