import numpy as np
import math
import matplotlib.pyplot as plt
from math import sqrt
from tkinter import *
from PIL import Image, ImageTk
import pylab
from mpl_toolkits.mplot3d import Axes3D

def stop(L,U,x,y):
    maxx = 0
    for i in range(len(x)):
        for j in range(len(y)):
            if abs(U[i,j] - L[i,j]) > maxx:
                maxx = abs(U[i,j] - L[i,j])
    return maxx

def norma(a):
    norm = 0
    for i in range(len(a)):
        norm += a[i]**2
    return sqrt(norm)

def f(x,y):
    return 0

def phi1(y):
    return np.exp(-y)*np.cos(y)

def phi2(y):
    return 0

def phi3(x):
    return np.cos(x)

def phi4(x):
    return 0

def TrueSolve(x, y):
    return np.exp(-y)*np.cos(x)*np.cos(y)

d = 1
a = 0
b = 2
c = 3

alfa1 = 0
alfa2 = 1
beta1 = 0
beta2 = 1
gama1 = 0
gama2 = 1
delta1 = 0
delta2 = 1

lbx = 0
ubx = np.pi/2
nx = 30
hx = (ubx - lbx)/nx


lby = 0
uby = np.pi/2
ny = 30
hy = (uby - lby)/ny


x = np.arange(lbx, ubx + hx, hx)
y = np.arange(lby, uby + hy, hy)

print(x,'\n',y)

eps = 0.001

def print_matr(U,x,y):
    print('            ', end = ' ')
    for j in range(len(x)):
        print('x = ', "%.2f" %x[j], end = '  ')
    print(' ')
    for i in range(len(y)):
        print('y = ',"%.3f" % y[i], end = ' |')
        for j in range(len(x)):
            print("%.5f" % U[j,i], end = '  ')
        print(' ')
    return

def MPI(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = phi3(x[i])/gama2
        U[i,-1] = phi4(x[i])/delta2    
    for j in range(len(y)):
        U[0,j] = phi1(y[j])/alfa2
        U[-1,j] = phi2(y[j])/beta2
    k = 0
    for j in range(1,len(y)-1):
        for i in range(1,len(x)-1):
            U[i,j] = U[0,j] + (x[i] - x[0])*(U[-1,j] - U[0,j])/(x[-1] - x[0])           
    while True:
        k = k+1
        print('k = ', k, end = '    ')
        #L = copy_matr(U,x,y)
        L = U.copy()
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = phi3(x[i])/gama2
            U[i,-1] = phi3(x[i])/delta2    
        for j in range(1, len(y)-1):
            U[0,j] = phi1(y[j])/alfa2
            U[-1,j] = phi2(y[j])/beta2
        
        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = (hx*hx*f(x[i],y[j]) - (L[i+1,j] + L[i-1,j]) - d*hx*hx*(L[i,j+1] + L[i,j-1])/(hy*hy) - a*hx*0.5*(L[i+1,j] - L[i-1,j]) - b*hx*hx*(L[i,j+1] - L[i,j-1])/(2*hy) )/(c*hx*hx - 2*(hy*hy + d*hx*hx)/(hy*hy))

        print(' eps = ',stop(L,U,x,y))
        if stop(L,U,x,y) <= eps:
            break
    return U, k


def ZEI(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = phi3(x[i])/gama2
        U[i,-1] = phi4(x[i])/delta2    
    for j in range(1,len(y)-1):
        U[0,j] = phi1(y[j])/alfa2
        U[-1,j] = phi2(y[j])/beta2
    k = 0
    for j in range(1,len(y)-1):
        for i in range(1,len(x)-1):
            U[i,j] = U[0,j] + (x[i] - x[0])*(U[-1,j] - U[0,j])/(x[-1] - x[0])           
    while True:
        k = k+1
        print('k = ', k, end = '    ')
        #L = copy_matr(U,x,y)
        L = U.copy()
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = phi3(x[i])/gama2
            U[i,-1] = phi3(x[i])/delta2    
        for j in range(1, len(y)-1):
            U[0,j] = phi1(y[j])/alfa2
            U[-1,j] = phi2(y[j])/beta2
            
        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = (hx*hx*f(x[i],y[j]) - (L[i+1,j] + U[i-1,j]) - d*hx*hx*(L[i,j+1] + U[i,j-1])/(hy*hy) - a*hx*0.5*(L[i+1,j] - U[i-1,j]) - b*hx*hx*(L[i,j+1] - U[i,j-1])/(2*hy) )/(c*hx*hx - 2*(hy*hy + d*hx*hx)/(hy*hy))

        print(' eps = ',stop(L,U,x,y))
        if stop(L,U,x,y) <= eps:
            break
    return U, k


def MPIR(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby, w):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U = np.zeros((len(x),len(y)))
    for i in range(len(x)):
        U[i,0] = phi3(x[i])/gama2
        U[i,-1] = phi4(x[i])/delta2    
    for j in range(1,len(y)-1):
        U[0,j] = phi1(y[j])/alfa2
        U[-1,j] = phi2(y[j])/beta2
    k = 0
    for j in range(1,len(y)-1):
        for i in range(1,len(x)-1):
            U[i,j] = U[0,j] + (x[i] - x[0])*(U[-1,j] - U[0,j])/(x[-1] - x[0])           
    while True:
        k = k+1
        print('k = ', k, end = '    ')
        #L = copy_matr(U,x,y)
        L = U.copy()
        U = np.zeros((len(x),len(y)))
        for i in range(len(x)):
            U[i,0] = phi3(x[i])/gama2
            U[i,-1] = phi3(x[i])/delta2    
        for j in range(1, len(y)-1):
            U[0,j] = phi1(y[j])/alfa2
            U[-1,j] = phi2(y[j])/beta2
            
        for j in range(1, len(y)- 1):
            for i in range(1, len(x)- 1):
                U[i,j] = ((hx*hx*f(x[i],y[j]) - (L[i+1,j] + U[i-1,j]) - d*hx*hx*(L[i,j+1] + U[i,j-1])/(hy*hy) - a*hx*0.5*(L[i+1,j] - U[i-1,j]) - b*hx*hx*(L[i,j+1] - U[i,j-1])/(2*hy) )/(c*hx*hx - 2*(hy*hy + d*hx*hx)/(hy*hy)))*w + (1 - w)*L[i,j]
        print(' eps = ',stop(L,U,x,y))
        if stop(L,U,x,y) <= eps:
            break
    return U, k

import pylab
from mpl_toolkits.mplot3d import Axes3D
import numpy

def makeData (lbx,lby,ubx,uby,hx,hy):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    xgrid, ygrid = numpy.meshgrid(x, y)
    zgrid = TrueSolve(xgrid, ygrid)
    return xgrid, ygrid, zgrid

def plotU(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby, w, ky):
    x = np.arange(lbx, ubx + hx, hx)
    y = np.arange(lby, uby + hy, hy)
    U1, k1 = MPI(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby)
    U2, k2 = ZEI(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby)
    U3, k3 = MPIR(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby, w)
    plt.figure(1)
    tit = 'U(x, ' + str(ky) + ')'
    plt.title(tit)
    label1 = 'МПИ (количество итераций: ' + str(k1) + ')'
    label2 = 'метод Зейделя (количество итераций: ' + str(k2) + ' )'
    label3 = 'МПИ с верхней релаксацией (количество итераций: ' + str(k3) + ')'
    plt.plot(x,TrueSolve(x,y[ky]),color = 'red', label = 'аналитическое решение')
    plt.plot(x,U1[:,ky],color = 'green', label = label1)
    plt.plot(x,U2[:,ky],color = 'blue', label = label2)
    plt.plot(x,U3[:,ky],color = 'orange', label = label3)
    plt.legend()
    plt.grid()
    
    
    plt.figure(3)
    plt.title('Погрешность')
    plt.xlim((0, 1.6))
    plt.grid()
    eps = []
    for i in range(len(x)):
        tmp = abs(TrueSolve(x[i], y[ky]) - U1[i,ky])
        eps = np.append(eps, tmp)
    plt.plot(x, eps, color = 'green')

    eps = []
    for i in range(len(x)):
        tmp = abs(TrueSolve(x[i], y[ky]) - U2[i,ky])
        eps = np.append(eps, tmp)
    plt.plot(x, eps, color = 'blue')
    
    eps = []
    for i in range(len(x)):
        tmp = abs(TrueSolve(x[i], y[ky]) - U3[i,ky])
        eps = np.append(eps, tmp)
    plt.plot(x, eps, color = 'orange')

    Solve = np.zeros((len(x), len(y)))
    for j in range(len(y)):
        for i in range(len(x)):
            Solve[i, j] = TrueSolve(x[i], y[j])

    x, y, z = makeData(lbx,lby,ubx,uby,hx,hy)



    fig = pylab.figure()
    axes = Axes3D(fig)
    axes.plot_surface(x, y, U1[:, :], color = 'green', label = label1)
    axes.plot_surface(x, y, U2[:, :],color = 'blue', label = label2)
    axes.plot_surface(x, y, U3[:, :],color = 'orange', label = label3)
    axes.plot_surface(x, y, Solve,color = 'red')
    plt.show()    
    return

root = Tk()
root.title("Лабораторная работа №7")
root["bg"] = "white"
w = root.winfo_screenwidth() # ширина экрана
h = root.winfo_screenheight() # высота экрана
ww = str(int(w/2))
hh = str(int(h-100))
a = ww + 'x' + hh
w = w//2 # середина экрана
h = h//2 
w = w - 200 # смещение от середины
h = h - 200
root.geometry(a + '+0+0'.format(w, h))

label4 = Label(text = "Параметры задачи:", justify = LEFT, font = "Arial 12")
label4.place(x = 10, y = 80)

labela = Label(text = "a = ", justify = LEFT, font = "Arial 12", bg = "white")
labela.place(x = 10, y = 110)
labelb = Label(text = "b = ", justify = LEFT, font = "Arial 12", bg = "white")
labelb.place(x = 10, y = 140)
labelc = Label(text = "c = ", justify = LEFT, font = "Arial 12", bg = "white")
labelc.place(x = 10, y = 170)
labeld = Label(text = "d = ", justify = LEFT, font = "Arial 12", bg = "white")
labeld.place(x = 10, y = 200)
labeleps = Label(text = "eps = ", justify = LEFT, font = "Arial 10", bg = "white")
labeleps.place(x = 10, y = 230)
labelw = Label(text = "w = ", justify = LEFT, font = "Arial 10", bg = "white")
labelw.place(x = 10, y = 260)


entrya = Entry(root, justify = RIGHT)
entrya.place(x = 50, y = 110)

entryb = Entry(root, justify = RIGHT)
entryb.place(x = 50, y = 140)

entryc = Entry(root, justify = RIGHT)
entryc.place(x = 50, y = 170)

entryd = Entry(root, justify = RIGHT)
entryd.place(x = 50, y = 200)

entryeps = Entry(root, justify = RIGHT)
entryeps.place(x = 50, y = 230)

entryw = Entry(root, justify = RIGHT)
entryw.place(x = 50, y = 260)

entrya.insert(0, 0)
entryb.insert(0, 2)
entryc.insert(0, 3)
entryd.insert(0, 1)
entryeps.insert(0, 0.001)
entryw.insert(0, 1.5)

label5 = Label(text = "Параметры конечно-разностной сетки:", justify = LEFT, font = "Arial 12")
label5.place(x = 300, y = 80)
labellx = Label(text = "lx = ", justify = LEFT, font = "Arial 12", bg = "white")
labellx.place(x = 300, y = 110)
labelnx = Label(text = "nx = ", justify = LEFT, font = "Arial 12", bg = "white")
labelnx.place(x = 300, y = 140)
labelly = Label(text = "ly = ", justify = LEFT, font = "Arial 12", bg = "white")
labelly.place(x = 300, y = 170)
labelny = Label(text = "ny = ", justify = LEFT, font = "Arial 12", bg = "white")
labelny.place(x = 300, y = 200)

entrylx = Entry(root, justify = RIGHT)
entrylx.place(x = 350, y = 110)

entrynx = Entry(root, justify = RIGHT)
entrynx.place(x = 350, y = 140)

entryly = Entry(root, justify = RIGHT)
entryly.place(x = 350, y = 170)

entryny = Entry(root, justify = RIGHT)
entryny.place(x = 350, y = 200)

entrylx.insert(0, np.pi/2)
entrynx.insert(0, 40)
entryly.insert(0, np.pi/2)
entryny.insert(0, 40)

labelhx = Label(text = " ", justify = LEFT, font = "Arial 12", bg = "white")
labelhx.place(x = 650, y = 110)
labelhy = Label(text = " ", justify = LEFT, font = "Arial 12", bg = "white")
labelhy.place(x = 650, y = 170)
labelerror = Label(text = " ", justify = LEFT, font = "Arial 12", bg = "white")
labelerror.place(x = 650, y = 140)


def params():
    try:
        lbx = 0
        ubx = float(entrylx.get())
        lby = 0
        uby = float(entryly.get())
        nx = float(entrynx.get())
        ny = float(entryny.get())
        labelhx.config(text = "hx = {}".format((ubx - lbx)/nx))
        labelhy.config(text = "hy = {}".format((uby - lby)/ny))
    except ValueError:
        labelerror.config(text = "Заполните все поля", fg = "red")

but = Button(root, text = "Показать все параметры сетки", command = params, font = "Arial 9")
but.place(x = 650, y = 679)

labelky = Label(text = "ky = ", justify = LEFT, font = "Arial 10", bg = "white")
labelky.place(x = 10, y = 320)

entryky = Entry(root, justify = RIGHT)
entryky.place(x = 50, y = 320)
entryky.insert(0, 1)

ky = float(entryky.get())

d = float(entryd.get())
a = float(entrya.get())
b = float(entryb.get())
c = float(entryc.get())

alfa1 = 0
alfa2 = 1
beta1 = 0
beta2 = 1
gama1 = 0
gama2 = 1
delta1 = 0
delta2 = 1

lbx = 0
ubx = float(entrylx.get())
nx = float(entrynx.get())
hx = (ubx - lbx)/nx


lby = 0
uby = float(entryly.get())
ny = float(entryny.get())
hy = (uby - lby)/ny

eps = float(entryeps.get())
w = float(entryw.get())

        

def solvv():
    try:
        ky = int(entryky.get())
        d = float(entryd.get())
        a = float(entrya.get())
        b = float(entryb.get())
        c = float(entryc.get())

        alfa1 = 0
        alfa2 = 1
        beta1 = 0
        beta2 = 1
        gama1 = 0
        gama2 = 1
        delta1 = 0
        delta2 = 1

        lbx = 0
        ubx = float(entrylx.get())  
        nx = float(entrynx.get())
        hx = (ubx - lbx)/nx

        lby = 0
        uby = float(entryly.get())
        ny = float(entryny.get())
        hy = (uby - lby)/ny

        eps = float(entryeps.get())
        w = float(entryw.get())
        plotU(a, b, c, d, alfa1, alfa2, beta1, beta2, gama1, gama2, delta1, delta2, hx, hy, eps, lbx, lby, ubx, uby, w, ky)
    except ValueError:
        labelerror.config(text = "Заполните все поля", fg = "red")



solv = Button(root, text = "  Решить  ", font = "Arial 14", command = solvv)
solv.place(x = 650, y = 300)

root.mainloop() 
