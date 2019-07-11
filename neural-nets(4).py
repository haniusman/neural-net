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
    
def output_k(outa, outb, outc, outd, wa, wb, wc, wd):
    op = (outa*wa) + (outb*wb) + (outc*wc) + (outd*wd)
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
    
    wak = rdm()
    wbk = rdm()
    wck = rdm()
    wdk = rdm()
    
    oa = output(dx,dy,dz,w1a,w2a,w3a)
    ob = output(dx,dy,dz,w1b,w2b,w3b)
    oc = output(dx,dy,dz,w1c,w2c,w3c)
    od = output(dx,dy,dz,w1d,w2d,w3d)
    
    ok = output_k(oa,ob,oc,od,wak,wbk,wck,wdk)
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
        #print("Error")
        #print(errk,erri,errj)

        cwak = delta(errk, oa)
        cwbk = delta(errk, ob)
        cwck = delta(errk, oc)
        cwdk = delta(errk, od)             

        
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

        wak = wak + cwak
        wbk = wbk + cwbk
        wck = wck + cwck
        wdk = wdk + cwdk
        
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

        oa = output(dx,dy,dz,w1a,w2a,w3a)
        ob = output(dx,dy,dz,w1b,w2b,w3b)
        oc = output(dx,dy,dz,w1c,w2c,w3c)
        od = output(dx,dy,dz,w1d,w2d,w3d)
    
        ok = output_k(oa,ob,oc,od,wak,wbk,wck,wdk)
        count = count+1

    #print ("while Loop ended")
    ok = ok * (13000-4800)+ 4800
    #print(ok,i)
    
    pvalue.append(ok)


for j in range(0,len(pvalue)):
    print(round(pvalue[j],2),sep='\n')
        










        
       


    

    





    
