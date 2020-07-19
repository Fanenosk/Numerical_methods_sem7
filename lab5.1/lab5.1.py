import numpy as np
import math
import matplotlib.pyplot as plt
from math import sqrt
from tkinter import *

def norma(a):
    norm = 0
    for i in range(len(a)):
        norm += a[i]**2
    return sqrt(norm)

def Trioganal(a, b, c, d):
    n = len(d)
    x = np.zeros(n)
    p = [-c[0] / b[0]]
    q = [d[0] / b[0]]
    for i in range(1, n):
        p.append(-c[i] / (b[i] + a[i] * p[i - 1]))
        q.append((d[i] - a[i] * q[i - 1]) / (b[i] + a[i] * p[i - 1]))
    x[-1] = q[-1]
    for i in reversed(range(n-1)):
        x[i] = p[i] * x[i + 1] + q[i]
    return x

def psi(x):
    return  np.cos(x)

def phi0(t, a, b):
    return - np.exp(-a*t) * (np.cos(b*t) + np.sin(b*t))

def phi1(t, a, b):
    return np.exp(-a*t) * (np.cos(b*t) + np.sin(b*t))

def TrueSolve(x, t, a, b):
    return np.exp(-a*t) * np.cos(x + b*t)


def explicit(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    U = np.zeros((len(t), len(x)))
    for j in range(len(x)):
        U[0, j] = psi(x[j])
    for i in range(1, len(t)):
        for j in range(1, len(x) - 1):
            U[i, j] = U[i - 1, j] + (a * tau / h ** 2) * (U[i - 1, j - 1] - 2 * U[i - 1, j] + U[i - 1, j + 1]) + (
                        b * tau * 0.5 / h) * (- U[i - 1, j - 1] + U[i - 1, j + 1])
            U[i, 0] = (h * phi0(t[i], a, b) - alfa * U[i, 1]) / (h * beta - alfa)
            U[i, -1] = (h * phi1(t[i], a, b) + gama * U[i, -2]) / (h * delta + gama)
    return U


def plot_ex(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, k):
    x1 = np.arange(lb, ub + h[0], h[0])
    x2 = np.arange(lb, ub + h[1], h[1])
    x3 = np.arange(lb, ub + h[2], h[2])
    t1 = np.arange(0, T + tau[0], tau[0])
    t2 = np.arange(0, T + tau[1], tau[1])
    t3 = np.arange(0, T + tau[2], tau[2])
    plt.figure(1)
    plt.title('Явная схема, t = ' + str(t3[k]))
    plt.grid()
    plt.plot(x3, TrueSolve(x3, t3[k], a, b), color='red', label='аналитическое решение')
    colors = ['green', 'orange', 'blue']
    U1 = explicit(a, b, alfa, beta, gama, delta, lb, ub, h[0], tau[0], T)
    plt.plot(x1, U1[k, :], color=colors[0], label='n = 10')
    U2= explicit(a, b, alfa, beta, gama, delta, lb, ub, h[1], tau[1], T)
    plt.plot(x2, U2[k, :], color=colors[1], label='n = 20')
    U3 = explicit(a, b, alfa, beta, gama, delta, lb, ub, h[2], tau[2], T)
    plt.plot(x3, U3[k, :], color=colors[2], label='n = 40')
    plt.legend()
    plt.xlim((0, ub))
    plt.figure(2)
    plt.title('Погрешность')
    plt.grid()
    eps = []
    for i in range(len(t1)):
        tmp = TrueSolve(x1, t1[i], a, b) - U1[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t1, eps, color='green', label='n=10')

    eps = []
    for i in range(len(t2)):
        tmp = TrueSolve(x2, t2[i], a, b) - U2[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t2, eps, color='orange', label='n=20')

    eps = []
    for i in range(len(t3)):
        tmp = TrueSolve(x3, t3[i], a, b) - U3[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t3, eps, color='blue', label='n=40')
    plt.legend()
    plt.show()
    return


def implicit(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    U = np.zeros((len(t),len(x)))
    for j in range(len(x)):
        U[0,j] = psi(x[j])
    for i in range(1, len(t)):
        aa = np.zeros(len(x))
        bb = np.zeros(len(x))
        cc = np.zeros(len(x))
        dd = np.zeros(len(x))
        dd[0] = phi0(t[i],a,b)*h
        dd[-1] = phi1(t[i],a,b)*h
        bb[0] = h*beta - alfa
        bb[-1] = h*delta + gama
        cc[0] = alfa
        aa[-1] = - gama
        for j in range(1,len(x)-1):
            aa[j] = (2*a - b*h)*tau*0.5/(h**2)
            bb[j] = - 2*a*tau/(h**2) - 1
            cc[j] = (2*a + h*b)*tau*0.5/(h**2)
            dd[j] = -U[i-1, j]
        xx = Trioganal(aa, bb, cc, dd)
        for j in range(len(x)):
            U[i, j] = xx[j]
    return U


def plot_im(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, k):
    x1 = np.arange(lb, ub + h[0], h[0])
    x2 = np.arange(lb, ub + h[1], h[1])
    x3 = np.arange(lb, ub + h[2], h[2])
    t1 = np.arange(0, T + tau[0], tau[0])
    t2 = np.arange(0, T + tau[1], tau[1])
    t3 = np.arange(0, T + tau[2], tau[2])
    plt.figure(1)
    plt.title('Неявная схема, t = ' + str(t3[k]))
    plt.grid()
    plt.plot(x3, TrueSolve(x3, t3[k], a, b), color='red', label='аналитическое решение')
    colors = ['green', 'orange', 'blue']
    U1 = implicit(a, b, alfa, beta, gama, delta, lb, ub, h[0], tau[0], T)
    plt.plot(x1, U1[k, :], color=colors[0], label='n = 10')
    U2 = implicit(a, b, alfa, beta, gama, delta, lb, ub, h[1], tau[1], T)
    plt.plot(x2, U2[k, :], color=colors[1], label='n = 20')
    U3 = implicit(a, b, alfa, beta, gama, delta, lb, ub, h[2], tau[2], T)
    plt.plot(x3, U3[k, :], color=colors[2], label='n = 40')
    plt.legend()
    plt.xlim((0, ub))
    plt.figure(2)
    plt.title('Погрешность')
    plt.grid()
    eps = []
    for i in range(len(t1)):
        tmp = TrueSolve(x1, t1[i], a, b) - U1[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t1, eps, color='green', label='n=10')

    eps = []
    for i in range(len(t2)):
        tmp = TrueSolve(x2, t2[i], a, b) - U2[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t2, eps, color='orange', label='n=20')

    eps = []
    for i in range(len(t3)):
        tmp = TrueSolve(x3, t3[i], a, b) - U3[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t3, eps, color='blue', label='n=40')
    plt.legend()
    plt.show()
    return


def eximplicit(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, teta):
    x = np.arange(lb, ub + h, h)
    t = np.arange(0, T + tau, tau)
    U = np.zeros((len(t), len(x)))
    for j in range(len(x)):
        U[0, j] = psi(x[j])
    for i in range(1, len(t)):
        aa = np.zeros(len(x))
        bb = np.zeros(len(x))
        cc = np.zeros(len(x))
        dd = np.zeros(len(x))
        dd[0] = phi0(t[i], a, b) * h
        dd[-1] = phi1(t[i], a, b) * h
        bb[0] = h * beta - alfa
        bb[-1] = h * delta + gama
        cc[0] = alfa
        aa[-1] = - gama
        for j in range(1, len(x) - 1):
            aa[j] = (2 * a - b * h) * teta * tau * 0.5 / (h ** 2)
            bb[j] = - 2 * a * teta * tau / (h ** 2) - 1
            cc[j] = (2 * a + h * b) * tau * teta * 0.5 / (h ** 2)
            dd[j] = -U[i - 1, j] - (1 - teta) * (
                        (a * tau / h ** 2) * (U[i - 1, j - 1] - 2 * U[i - 1, j] + U[i - 1, j + 1]) + (
                            b * tau * 0.5 / h) * (- U[i - 1, j - 1] + U[i - 1, j + 1]))
        xx = Trioganal(aa, bb, cc, dd)
        for j in range(len(x)):
            U[i, j] = xx[j]
    return U

def plot_exim(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, teta, k):
    x1 = np.arange(lb, ub + h[0], h[0])
    x2 = np.arange(lb, ub + h[1], h[1])
    x3 = np.arange(lb, ub + h[2], h[2])
    t1 = np.arange(0, T + tau[0], tau[0])
    t2 = np.arange(0, T + tau[1], tau[1])
    t3 = np.arange(0, T + tau[2], tau[2])
    plt.figure(1)
    if teta == 0.5:
        plt.title('Схема Кранка - Николсона, t = ' + str(t3[k]))
    else:
        plt.title('Явно-неявная схема, t = ' + str(t3[k]))
    plt.grid()
    plt.plot(x3, TrueSolve(x3, t3[k], a, b), color='red', label='аналитическое решение')
    colors = ['green', 'orange', 'blue']
    U1 = eximplicit(a, b, alfa, beta, gama, delta, lb, ub, h[0], tau[0], T, teta)
    plt.plot(x1, U1[k, :], color=colors[0], label='n = 10')
    U2 = eximplicit(a, b, alfa, beta, gama, delta, lb, ub, h[1], tau[1], T, teta)
    plt.plot(x2, U2[k, :], color=colors[1], label='n = 20')
    U3 = eximplicit(a, b, alfa, beta, gama, delta, lb, ub, h[2], tau[2], T, teta)
    plt.plot(x3, U3[k, :], color=colors[2], label='n = 40')
    plt.legend()
    plt.xlim((0, ub))
    plt.figure(2)
    plt.title('Погрешность')
    plt.grid()
    eps = []
    for i in range(len(t1)):
        tmp = TrueSolve(x1, t1[i], a, b) - U1[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t1, eps, color='green', label='n=10')

    eps = []
    for i in range(len(t2)):
        tmp = TrueSolve(x2, t2[i], a, b) - U2[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t2, eps, color='orange', label='n=20')

    eps = []
    for i in range(len(t3)):
        tmp = TrueSolve(x3, t3[i], a, b) - U3[i, :]
        eps = np.append(eps, norma(tmp))
    plt.plot(t3, eps, color='blue', label='n=40')
    plt.legend()
    plt.show()
    return


root = Tk()
root.title("Лабораторная работа №5")
root["bg"] = "white"
w = root.winfo_screenwidth()  # ширина экрана
h = root.winfo_screenheight()  # высота экрана
ww = str(int(w / 2))
hh = str(int(h - 250))
a = ww + 'x' + hh
w = w // 2  # середина экрана
h = h // 2
w = w - 200  # смещение от середины
h = h - 200
root.geometry(a + '+0+0'.format(w, h))

label4 = Label(text="Параметры задачи:", justify=LEFT, font="Arial 12")
label4.place(x=10, y=30)
labela = Label(text="a = ", justify=LEFT, font="Arial 12", bg="white")
labela.place(x=10, y=60)
labelb = Label(text="b = ", justify=LEFT, font="Arial 12", bg="white")
labelb.place(x=10, y=90)
labell = Label(text="l = ", justify=LEFT, font="Arial 12", bg="white")
labell.place(x=10, y=120)

entrya = Entry(root, justify=RIGHT)
entrya.place(x=50, y=60)

entryb = Entry(root, justify=RIGHT)
entryb.place(x=50, y=90)

entryl = Entry(root, justify=RIGHT)
entryl.place(x=50, y=120)

entrya.insert(0, 1)
entryb.insert(0, 2)
entryl.insert(0, np.pi)

label5 = Label(text="Параметры конечно-разностной сетки:", justify=LEFT, font="Arial 12")
label5.place(x=300, y=30)
#labeln = Label(text="n = ", justify=LEFT, font="Arial 12", bg="white")
#labeln.place(x=300, y=60)
labelsigm = Label(text="σ = ", justify=LEFT, font="Arial 12", bg="white")
labelsigm.place(x=300, y=90)
labelT = Label(text="T = ", justify=LEFT, font="Arial 12", bg="white")
labelT.place(x=300, y=120)

#entryn = Entry(root, justify=RIGHT)
#entryn.place(x=350, y=60)

entrysigm = Entry(root, justify=RIGHT)
entrysigm.place(x=350, y=90)

entryT = Entry(root, justify=RIGHT)
entryT.place(x=350, y=120)

#entryn.insert(0, 40)
entrysigm.insert(0, 0.1)
entryT.insert(0, 2)

labelh = Label(text=" ", justify=LEFT, font="Arial 12", bg="white")
labelh.place(x=650, y=60)
labeltau = Label(text=" ", justify=LEFT, font="Arial 12", bg="white")
labeltau.place(x=650, y=90)
labelerror = Label(text=" ", justify=LEFT, font="Arial 12", bg="white")
labelerror.place(x=650, y=150)
labelK = Label(text=" ", justify=LEFT, font="Arial 12", bg="white")
labelK.place(x=650, y=120)


def params():
    try:
        lb = 0.0
        n = [10.0, 20.0, 40.0]
        h = np.zeros(3)
        ub = float(entryl.get())
        teta = float(entryteta.get())
        a = float(entrya.get())
        T = float(entryT.get())
        #n = float(entryn.get())
        sigm = float(entrysigm.get())
        labelh.config(text="h = {}".format((ub - lb) / n[0]))
        for i in range(0, 3):
            h[i] = (ub - lb) / n[i]
            tau[i] = sigm * h[i] * h[i] / a ** 2
            K[i] = int(T / tau[i])
        labeltau.config(text="τ = {}".format(sigm * h[0] * h[0] / a))
        labelK.config(text="K = {}".format(int(T * a / (sigm * h[0] * h[0]))))
    except ValueError as e:
        labelerror.config(e)


but = Button(root, text="Показать все параметры сетки", command=params, font="Arial 9")
but.place(x=650, y=527)

label6 = Label(text="Параметр схемы:", justify=LEFT, font="Arial 12")
label6.place(x=10, y=200)
labelteta = Label(text="ϴ = ", justify=LEFT, font="Arial 12", bg="white")
labelteta.place(x=10, y=230)
entryteta = Entry(root, justify=RIGHT)
entryteta.place(x=50, y=230)
entryteta.insert(0, 0)

label7 = Label(text="Отобразить решение на шаге по времени:", justify=LEFT, font="Arial 12")
label7.place(x=300, y=200)
labelk = Label(text="k = ", justify=LEFT, font="Arial 12", bg="white")
labelk.place(x=300, y=230)
entryk = Entry(root, justify=RIGHT)
entryk.place(x=340, y=230)

a = float(entrya.get())
ub = float(entryl.get())
lb = 0.0
n = [10.0, 20.0, 40.0]
h = np.zeros(3)
tau = np.zeros(3)
K = np.zeros(3)
teta = float(entryteta.get())
T = float(entryT.get())
sigm = float(entrysigm.get())
#n = float(entryn.get())
for i in range(0, 2):
    h[i] = (ub - lb) / n[i]
    tau[i] = sigm * h[i] * h[i] / a**2
    K[i] = int(T / tau[i])

entryk.insert(0, 12)


def solver(a, b, alfa, beta, gama, delta, lb, ub, n, K, h, tau, T, sigm, teta, k):
    if teta == 0:
        plot_ex(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, k)
    if teta == 1:
        plot_im(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, k)
    if teta > 0 and teta < 1:
        plot_exim(a, b, alfa, beta, gama, delta, lb, ub, h, tau, T, teta, k)
    return

    plt.figure()
    plt.grid()
    plt.scatter(a, b)
    plt.show()
    return


def solvv():
    try:
        teta = float(entryteta.get())
        a = float(entrya.get())
        b = float(entryb.get())
        ub = float(entryl.get())
        lb = 0
        T = float(entryT.get())
        sigm = float(entrysigm.get())
        #n = float(entryn.get())
        n = [10.0, 20.0, 40.0]
        h = np.zeros(3)
        tau = np.zeros(3)
        K = np.zeros(3)
        for i in range(0, 3):
            h[i] = (ub - lb) / n[i]
            tau[i] = sigm * h[i] * h[i] / a**2
            K[i] = int(T / tau[i])
        k = int(entryk.get())
        alfa = 1
        beta = -1
        gama = 1
        delta = -1
        solver(a, b, alfa, beta, gama, delta, lb, ub, n, K, h, tau, T, sigm, teta, k)
    except ValueError as e:
        print(e)
        labelerror.config(text="Заполните все поля", fg="red")


solv = Button(root, text="  Решить  ", font="Arial 14", command=solvv)
solv.place(x=650, y=200)

root.mainloop()
root.mainloop()