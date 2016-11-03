#!/usr/bin/python3
#Runs with Python 3.2-3.4
__doc__="""
Module clsexexFit
Fits formula 
Ax + B + Cexp(-exp(-D(t-te)))
to the relative elongation of stems.
The data are read into self.DataTuples (time, rate [, error]).
Conversion is done to the basic units whenever applicable.
They do not have to be equidistant.
The callable function is called self.funcexex().
Requires matplotlib, numpy, scipy, sys, os, glob
P. Zajdel, A. Haduch-Sendecka, M. Pietruszka
Katowice 2013-2016
Dec 2014 - updated due to changes in Tcl/Tk
Feb 2016 - updated due to changes in Tcl/Tk, NavigationToolbar is using pack() as default and we need grid()
Mar 2016 - changed exp(exp()) function to avoid over and uderflow
Nov 2016 - added tangent line at te
"""
import os #.path.
import sys
if sys.hexversion >= 0x02060000 and sys.hexversion < 0x03000000: 
    inPy26 = True
else: 
    inPy26 = False

if sys.hexversion >= 0x03000000: 
    inPy3 = True
else: 
    inPy3 = False
    
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg 
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr

#from mpl_toolkits.mplot3d.axes3D import Axes3D
#from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
np.seterr(all = 'warn', over = 'raise', under = 'raise')

import scipy as scp
import scipy.optimize as sci_opti

try:
    import Tkinter as tkTkinter
    from Tkinter import *
except:
    import tkinter as tkTkinter
""" For Python 3.0 """
if inPy3:
    import tkinter.scrolledtext as tkScrolledText
    import tkinter.messagebox as tkmb
    import tkinter.simpledialog as tksimdial
    import tkinter.filedialog as tkFileDialog
elif inPy26 and not inPy3:
    import ScrolledText as tkScrolledText
    import tkMessageBox as tkmb
    import tkSimpleDialog as tksimdial
    import tkFileDialog as tkFileDialog
else:
    print("Unsupported (<2.6) Python version")
    quit()

import glob as glb

class clsFitexex(tkTkinter.Tk):
    #All numerical types, which can be passed to function
    NumericalTypes = (int,float,np.float64,np.float32,np.float)
# numpy.float32 is C compatible
# numpy.float64 is Python compatible
# numpy.float is what most people would use in numpy
# int and float can be passed from python
    def __init__(self, parent):
        tkTkinter.Tk.__init__(self, parent)
        self.parent=parent
        self.grid()
        self.protocol("WM_DELETE_WINDOW", self._quit)
        self.frame = tkTkinter.Frame(self)
        self.toolbar_frame = tkTkinter.Frame(self)
        self.toolbar_frame.grid()
        self.A = np.float64(0.0)
        self.B = np.float64(0.0)
        self.C = np.float64(1.0)
        self.D = np.float64(1.0)
#        self.F = 1        
        self.te = np.float64(0.0)
        self.tmin = np.float64(0.0)
        self.tmax = np.float64(100.0)
        self.tstep = np.float64(1.0)
        self.DebugLevel = 0
        self.createWidgets()
        self.ClearMe()
        self.UpdateValues()

    def ClearMe(self):
        """ 
Method ClearMe
Clears and initializes variables
        """
        #self.strFileName.set("")
        self.filename = ""
        self.A = np.float64(0.0)
        self.B = np.float64(0.0)
        self.C = np.float64(1.0)
        self.D = np.float64(1.0)
#        self.F = 1
        self.te = np.float64(0.0)
        self.tmin = np.float64(0.0)
        self.tmax = np.float64(100.0)
        self.tstep = np.float64(1.0)
        self.DataTuples = []
        self.TimeExper = []
        self.RelElExper = []
        self.ErrorExper = []
        self.Npoints = 0
        # inEstimate informs plotting function if there is initial estimate of parameters
        # we do not want to plot unless function will roughly fit in the window
        # = 0 no estimate
        # = 1 start of estimation. Set by cmdEstimate
        # = 2 waiting for user input. Set by onClick
        # = 3 done estimating. After returning from onClick. Set by cmdEstimate
        self.inEstimate = 0
        # For full fit
        self.t0 = 0
        self.tau1 = 1
        self.tau2 = 1
        self.FiBar = 1
        self.R0 = 0
        self.LowCutOff = 0
        self.LabelX = ""
        self.LabelY = ""
        # end of parameters for full fit
        self.butPlot.config(state="disabled")
        self.butEstim.config(state="disabled")
        self.butLSQ.config(state="disabled")
        self.butfmin.config(state="disabled")
        self.butSaveReport.config(state="disabled")
        self.butA.config(state="disabled")
        self.butB.config(state="disabled")
        self.butC.config(state="disabled")
        self.butD.config(state="disabled")
#        self.butF.config(state="disabled")
        self.butt0.config(state="disabled")
        self.buttmax.config(state="disabled")
        self.buttmin.config(state="disabled")
        self.buttstep.config(state="disabled")        

    def RoundMe(self, tValErr):
        """ 
Method RoundMe
Receives 2-tuple (value,error) as np.float
Returns 2-tuple ("value","error" ) as string 
1. Brings both values to common order
2. Leaves 2 decimal places in error 
        """
        self._tV = tValErr[0]
        self._tVsign = 1
        if self._tV < 0.0:
            self._tVsign = -1
            self._tV = -self._tV
        self._tE = tValErr[1]
        self._tOrder = int(np.log10(self._tE))
        if self._tOrder>=0:
            self._tOrder -=1
        else:
            self._tOrder -=2
        #print(self._tV,self._tE)            
        self._tV = int(self._tV/10**self._tOrder + 0.5)
        self._tE = int(self._tE/10**self._tOrder + 0.9)
        #print(self._tV,self._tE)
        self._tV = self._tV* 10**self._tOrder
        self._tE = self._tE* 10**self._tOrder
        self._tV = self._tV * self._tVsign
        return ("{0:>10g}".format(self._tV),"{0:>10g}".format(self._tE))
            
    def createWidgets(self):
        #Quits app
        self.QUIT = tkTkinter.Button(self, text = "QUIT", fg   = "blue", command =  self._quit)
        #QPlots current data and fir
        self.butPlot = tkTkinter.Button(self,text = "Plot", command = self.PlotMeBut, width=10)
        # opens and reads file
        self.butReadFile = tkTkinter.Button(self, text = "Read data", command = self.cmdReadFile,width=10)
        # user guided estimation of parameters
        self.butEstim = tkTkinter.Button(self, text = "Get estimates", command = self.cmdEstimate,width=10)
        # simples minimalization
        self.butfmin = tkTkinter.Button(self, text = "fmin", command = self.cmdfmin,width=10)
        # Levenberg-Marq. min
        self.butLSQ = tkTkinter.Button(self, text = "LSQ", command = self.cmdLSQ,width=10)
        # save txt and html
        self.butSaveReport = tkTkinter.Button(self,text  = "Save report",command = self.cmdSaveReport,width=10)
        # labels for text fields
        self.lblTlow = tkTkinter.Label(self)
        self.lblTlow["text"] = "0"

        self.lblThigh = tkTkinter.Label(self)
        self.lblThigh["text"] = "300"

        self.lblTstep = tkTkinter.Label(self)
        self.lblTstep["text"] = "1"
        # fir parameters
        self.butA = tkTkinter.Button(self, text = "A",command = self.setA)
        self.butB = tkTkinter.Button(self,text= "B", command = self.setB)
        self.butC = tkTkinter.Button(self,text= "C", command = self.setC)
        self.butD = tkTkinter.Button(self, text= "D", command = self.setD)
#        self.butF = tkTkinter.Button(self, text= "F", command = self.setF)
        self.butt0 = tkTkinter.Button(self, text= "te", command = self.sett0)
        self.buttmin = tkTkinter.Button(self, text = "tmin", command = self.settmin)
        self.buttmax = tkTkinter.Button(self, text = "tmax", command = self.settmax)
        self.buttstep = tkTkinter.Button(self, text = "tstep", command = self.settstep)
        
        self.lblA = tkTkinter.Label(self,width=10, height=1, text=str(self.A))
        self.lblB = tkTkinter.Label(self, width=10, height=1, text=str(self.B))
        self.lblC = tkTkinter.Label(self, width=10, height=1, text=str(self.C))
        self.lblD = tkTkinter.Label(self, width=10, height=1, text = str(self.D))
#       self.lblF = tkTkinter.Label(self, width=10, height=1, text = str(self.F))
        self.lblt0 = tkTkinter.Label(self, width=10, height=1, text = str(self.te))
        self.lbltmin = tkTkinter.Label(self, width=10, height=1, text = str(self.tmin))
        self.lbltmax = tkTkinter.Label(self, width=10, height=1, text = str(self.tmax))
        self.lbltstep = tkTkinter.Label(self, width=10, height=1, text = str(self.tstep))
        # fir quality fields
        self.lblchi2 = tkTkinter.Label(self, width=10, height=1, text = "Chi2")
        self.lblchi2value = tkTkinter.Label(self, width=10, height=1, text = "None")
        self.txtInfo = tkScrolledText.ScrolledText(self,wrap="word", width=40, height = 5)
        
        self.myfig = plt.figure()
        #print( self.myfig.number)
        self.canvas = FigureCanvasTkAgg(self.myfig, master=self)
        #self.canvas.get_tk_widget().grid(column=0, row=2, columnspan=2)
        #print("Pre")
        #self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame)
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.toolbar_frame)
        #print("Post")
        self.ax1 = self.myfig.add_subplot(111)
        self.ax1.grid(True)
        #self.ax1.axhline(0, color='black', lw=2)
#position widgets
        # row1 buttons
        self.c=0
        self.r=0
        self.QUIT.grid(column=self.c,row=self.r, sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butPlot.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butReadFile.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butEstim.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butfmin.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butLSQ.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)
        self.c+=1
        self.butSaveReport.grid(column=self.c,row=self.r,sticky=tkTkinter.E+tkTkinter.W)        

        self.r=1
        self.c=0
        self.butA.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lblA.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.butB.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lblB.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.butC.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lblC.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.butD.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lblD.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
#        self.butF.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
#        self.r+=1
#        self.lblF.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
#        self.r+=1
        self.butt0.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lblt0.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.buttmin.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lbltmin.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.buttmax.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lbltmax.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.buttstep.grid(column=self.c,row=self.r, sticky=tkTkinter.W+tkTkinter.E+tkTkinter.S)
        self.r+=1
        self.lbltstep.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.lblchi2.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        self.r+=1
        self.lblchi2value.grid(column=self.c,row=self.r, sticky=tkTkinter.N)
        #row 2, column 2 
        self.c=1
        self.r=1
        self.toolbar_frame.grid(column=self.c,row=self.r, columnspan=6, sticky=tkTkinter.N+tkTkinter.W+tkTkinter.S)
        self.canvas._tkcanvas.grid(column=self.c,row=2, columnspan=6, rowspan=16,sticky=tkTkinter.N+tkTkinter.E+tkTkinter.W)
        self.txtInfo.grid(column=self.c,row=2+16, columnspan=6, rowspan=3,sticky=tkTkinter.N+tkTkinter.E+tkTkinter.W)
        self.addInfo("Hello!")

    def addInfo(self, info):
        """ 
Method addInfo(str)
Adds information to the status window
        """
        self.txtInfo.insert(tkTkinter.END,info+os.linesep )
        self.txtInfo.update()

        
    def _quit(self):
        """ Clean up windows and quit app """
        #self.newfig.close()
        self.destroy()
        self.quit()
    
    def PlotMeBut(self, full = 0):
        """ Callback from PlotMe button. Allows to change labels of the axes """
        self.LabelX = tksimdial.askstring("X Label", "X Label", initialvalue=self.LabelX)
        self.LabelY = tksimdial.askstring("Y Label", "Y Label", initialvalue=self.LabelY)
        full = tksimdial.askinteger("Plot type", "0 - no tangent\n 1 - with tangent", initialvalue=1)
        self.PlotMe(full)
        
    def PlotMe(self, full = 0):
        """ 
Method PlotMe(full=0)
Plots data and fit 
full=0 (default) restrict to (tmin, tmax)
full=1 with tangent line
"""
        #print("In plot")
        self.myfig.clf()
        self.ax1 = self.myfig.add_subplot(111)
        self.ax1.clear()
        if self.LabelX != "" or self.LabelY != "":
            self.ax1.xlabel = self.LabelX
            self.ax1.ylabel = self.LabelY
        self.mex = []
        self.mey = []
        for self._i in self.DataTuples:
            #if full == 1 and self._i[0]<self.LowCutOff: continue
            if self._i[0] >= self.tmin and self._i[0] <= self.tmax:
                self.mex.append(self._i[0])
                self.mey.append(self._i[1])
        # if time is normalized to seconds, changes format of x-axis
        if self.tmax > 1000000:
            self.ax1.xaxis.set_major_formatter(tkr.FormatStrFormatter('%6.3g'))
        else:
            self.ax1.xaxis.set_major_formatter(tkr.FormatStrFormatter('%.01f'))
        self.ax1.yaxis.set_major_formatter(tkr.FormatStrFormatter('%6.3g'))
        self.ax1.plot(self.mex, self.mey, "ro", label="data")
        self.ax1.plot(self.mex,self.funcexex(self.mex), "b-", label="model", linewidth = 3)
        if self.inEstimate > 1:
            #plot functions during (=1) or after estimation of fit(=2)
            # there is no sense in plotting functions, which are not at least partially comparable to data
            self._tx = np.linspace(self.tmin,self.tmax, 101)
            self.ax1.plot(self._tx,self.funcBack(self._tx), "c-", label="linear",linewidth = 3)
            self.ax1.plot(self._tx,self.funcGrow(self._tx), "m-", label="nonlin",linewidth = 3)
#        if full == 1: 
#            for self._i in self.DataTuples:
#                if self._i[0]<self.LowCutOff: 
#                    continue
#                else:
#                    self.mex.append(self._i[0])
#                    self.mey.append(self._i[1])
#            self.ax1.plot(self._tx,self.funcFull(self._tx), "g-", label="Full",linewidth = 3)
        self.ax1.axvline(x=self.te, linewidth=3, color='y')    
        #if self.LabelX != "" or self.LabelY != "":
        # ******************** UPDATE from here
        self.tanga = self.A + self.C * self.D / np.e        
        self.tangb = self.B + self.C/np.e*(1.0 - self.D * self.te)        
        self.tint = - (self.tangb)/(self.tanga)                
        if full == 1:            
            #Tangent line y = ax + b            
            self._tx = np.linspace(self.tint, 2*self.te-self.tint, 51)            
            self.ax1.plot(self._tx,self.tanga*self._tx+self.tangb, "g--",linewidth = 2)            
            self.ax1.axvline(x=self.tint, linewidth=3, color='g')             
            # End Tangent Line
        self.ax1.set_xlabel(self.LabelX)
        self.ax1.set_ylabel(self.LabelY)
        plt.legend(loc=2)
        self.canvas.show()

    def funcexex(self, *args):
        """ 
Main user-callable function. 
Expects one numerical argument or compound type containing numbers.
Returns numpy.array of numpy.float64 
            Ax + B + Cexp(-exp(-D(t-te)))
        """
        if self.DebugLevel>0: print("in funcexex")
        if len(args)>2:
            if self.DebugLevel>0: print("One argument only: number or iterable")
            return None        
        if type(args[0]) in clsFitexex.NumericalTypes:
            if self.DebugLevel >1: print("Number")
            return self.__funcexex__(args[0])
        result = []
        if self.DebugLevel >2: print(args[0])
        for self._ii,self._tii in enumerate(args):
#from here is can be cut            
            for self._jj in self._tii:
                if self.DebugLevel >2: print("In call Elem ", self._jj, " in arg ", self._ii)
                result.append(self.__funcexex__(self._jj))
# This is split into 3 lines for debug purposes                 
# Otherwise it can be replaced by comprehension list
#        return np.array([self.__funcexex__(self._jj) for self._jj in self._tii ]).astype(np.float64)
        return np.array(result).astype(np.float64)
        
    def __funcexex__(self, gx):
        """ 
This is just a wrapper function, which does basic checks on parameter values.
Returns calls separate functions, which return linear and non-linear parts.
It is not intended to be called by user
It operates on current values of parameters A, B, C, D, te
            Ax + B and Cexp(-exp(-D(t-te)))
            """
            # C must be >0
        if self.C<0.0:
            if self.DebugLevel>1: print("C under 0, correcting")
            self.C = np.float64(0.0000001)
            # D must be >0            
        if self.D<0.0:
            self.D = np.float64(0.0000001)
            if self.DebugLevel>1: print("D under 0, correcting")         
        if self.B<0.0: self.B = np.float64(0.0)
        # A should not be < 0
        # if you are sure that you need A < 0, comment out those lines
        if self.A<0.0:
            self.A = np.float64(0.0)
            if self.DebugLevel>1: print("A under 0, correcting")
        return  self.__funcBack__(gx) + self.__funcGrowParmsCDte__(gx, self.C, self.D, self.te)

    def funcBack(self, *args):
        """ 
User callable function for the linear part. 
Expects one numerical argument or compound type containing numbers.
Returns  numpy.array of numpy.float64 
            Ax + B """
        #Ax + B
        if self.DebugLevel>0: print("in funcback")
        if len(args)>2:
            if self.DebugLevel>0: print("One argument only: number, or compound type")
            return None        
        if type(args[0]) in clsFitexex.NumericalTypes:
            if self.DebugLevel >1: print("Number")
            return self.__funcBack__(args[0])
        result = []
        if self.DebugLevel >2: print(args[0])
        for self._ii,self._tii in enumerate(args):
            for self._jj in self._tii:
                if self.DebugLevel >2: print("In call Elem ", self._jj, " in arg ", self._ii)
                result.append(self.__funcBack__(self._jj))
                pass
# as in the main function. If you do not need debug, replace it with comprehension list 
#        return np.array([self.__funcBack__(self._jj) for self._jj in self._tii ]).astype(np.float64)
        return np.array(result).astype(np.float64)        
        
    def __funcBack__(self, gx):
        """
Returns linear background for a single point.
Expects one numerical argument.
It is not meant to be called directly
        """
        #Ax + B
        return self.A * np.float64(gx) + self.B
        
    def funcGrow(self, *args):
        """ 
User callable function for the non-linear part.
Expects one numerical argument or compound type containing numbers.
Optionally C, D, te
Returns  numpy.array of numpy.float64
            Cexp(-exp(-D(t-te))) """        
        #Cexp(-exp(-D(t-te)))
        if self.DebugLevel>0: print(" in funcGrow")
        if not (len(args)==1 or len(args)==4):
            if self.DebugLevel>0: print("Wrong number of parms")
            return None        
        if type(args[0]) in clsFitexex.NumericalTypes:
            if self.DebugLevel >1: print("Number")
            if len(args)==1: return self.__funcGrowParmsCDte__(args[0], self.C, self.D, self.te)
            if len(args)==4:
                return self.__funcGrowParmsCDte__(args[0],args[1],args[2],args[3])
        result = []
        if self.DebugLevel >2: print(args[0])
        for self._jj in args[0]:
            if self.DebugLevel >2: print("In call Elem ", self._jj, " in arg[0] ")
            if len(args)==1: 
                result.append(self.__funcGrowParmsCDte__(self._jj, self.C, self.D, self.te))
            if len(args)==4: 
                result.append(self.__funcGrowParmsCDte__(self._jj, args[1],args[2],args[3]))

                pass
        return np.array(result).astype(np.float64) 
        
    def __funcGrowParmsCDte__(self, gx, c,d,te):
        """
Returns non-linear contribution for a single point.
Expects one numerical argument and C, D, te
It is not meant to be called directly
        """
        #Cexp(-exp(-D(t-te)))
        gx = np.float64(gx)
        c = np.float64(c)
        d = np.float64(d)
        te = np.float64(te)
        arg = - d *(gx - te)
        if arg > 6.5:
            #The inner exp returns large positive value in gx2
            # then the outer exp gets -gx2
            # and the outer exp returns 0
            #print("Inner arg > 7")
            return np.float64(0.0)
        elif arg < -700.0:
            #The inner exp returns 0
            #it is multiplied by a neg coeff
            #the outer exp will return 1          
            return c * np.float64(1.0) 
         
        try:
        	gx2 = np.exp(arg)
    #		pass
        except FloatingPointError as err:
            print("Should never be here")
            quit()
       # here gx2 should be calculated without erros  
        if gx2 > 700:
            #underflow will occur in outer exp for 710
            # should be covered above
            return np.float64(0.0)
        try:
        	grow = c * np.exp(-gx2)
    #		pass
        except FloatingPointError as err:
            print('Uncatched Overflow or Underflow in __funcgrow__: ' + str(err))
            quit()

        return grow

#    def __funcGrow__(self, gx):
#        """
#Returns non-linear contribution for a single point.
#Expects one numerical argument.
#It is not meant to be called directly
#        """

#        return self.C * np.exp(-np.exp(-self.D*(gx-self.te)))
 
        
    def UpdateValues(self):
        """ Updates info on labels """
        self.lblA.config(text="{0:>12g}".format(self.A))
        self.lblB.config(text="{0:>12g}".format(self.B))
        self.lblC.config(text="{0:>12g}".format(self.C))
        self.lblD.config(text ="{0:>12g}".format(self.D))
        self.lblt0.config(text ="{0:>12g}".format(self.te))
        self.lbltmin.config(text ="{0:>12g}".format(self.tmin))
        self.lbltmax.config(text ="{0:>12g}".format(self.tmax))
        self.lbltstep.config(text ="{0:>12g}".format(self.tstep))
                  
    def setA(self):
        self.A = tksimdial.askfloat("A", "A", initialvalue=self.A)
        self.UpdateValues()
        pass 
   
    def setB(self):
        self.B = tksimdial.askfloat("B", "B", initialvalue=self.B)
        self.UpdateValues()
        pass
        
    def setC(self):
        self.C = tksimdial.askfloat("C", "C", initialvalue=self.C)
        self.UpdateValues()
        pass
        
    def setD(self):
        self.D = tksimdial.askfloat("D", "D", initialvalue=self.D)
        self.UpdateValues()
        pass
    
    def sett0(self):
        self.te = tksimdial.askfloat("te", "te", initialvalue=self.te)
        self.UpdateValues()
        pass 
   
    def settmin(self):
        self.tmin = tksimdial.askfloat("t_min", "t_min", initialvalue=self.tmin)
        self.UpdateValues()
        pass
        
    def settmax(self):
        self.tmax = tksimdial.askfloat("t_max", "t_max", initialvalue=self.tmax)
        self.UpdateValues()
        pass    
        
    def settstep(self):
        self.tstep = tksimdial.askfloat("t_step", "t_step", initialvalue=self.tstep)
        self.UpdateValues()
        pass            
        
    def cmdReadFile(self):
        """ 
Method reads text file in format
time, relative length, error
Data is normalized to get proper scaling 
If no error is given, assume 1
        """
        # Clean if opening a new file
        #if self.filename != "": 
        self.ClearMe()
       #while self.strFileName.get()=="":
            #or not os.access(self.strFileName, os.W_OK):
        self.filename = tkFileDialog.askopenfilename(defaultextension=".dat", initialdir=".", multiple=False)
        #print self.strFileName, os.access(self.strFileName, os.W_OK)
            
        if self.DebugLevel > 1: print("Filename:" + self.filename)
        self.fp = open(self.filename,'r')
        self.lisDataFile = self.fp.readlines()
        self.fp.close()
        self._divisorL = tksimdial.askfloat("Input divisor L", "Divisor that makes the elongation data dimensionless and 1\n e.g. 10000 for micron/cm", initialvalue = str(10000))
        self._multiT = tksimdial.askfloat("Input time factor", "Factor that converts time scale into seconds\n e.g. 60 for if the time in the file is in minutes", initialvalue = str(60))
# initialize tmin and tmax, so we can catch proper values
        self.tmax = -1*self._multiT
        self.tmin = 100*self._multiT
        # process header or until we find [Data]
        for self._line in self.lisDataFile:
            # drop comments
            if self._line[0] == "#": continue
            
            for self._i,self._t in enumerate(self._line.split()):
                if len(self._line.split()) == 0: 
                    if self.DebugLevel >2: print("empty line")
                    continue
                self._error = 1
                if self._i == 0: 
                    #print("Time " + self._t)
                    self._time = self._multiT*float(self._t)
                    if self._time > self.tmax: self.tmax = self._time
                    if self._time < self.tmin: self.tmin = self._time
                    self.TimeExper.append(self._time)
                elif self._i == 1:
                    #print("Elong " + self._t)
                    self._elong = float(self._t)/self._divisorL
                    self.RelElExper.append(self._elong)
                else:
                    self._error = float(self._t)/self._divisorL
                    if self._error == 0: self._error = 0.000001
                    self.ErrorExper.append(self._error)
            self.DataTuples.append((self._time, self._elong, self._error))

        self.Npoints = len(self.DataTuples)
        if len(self.ErrorExper) == 0: self.ErrorExper = np.ones(self.Npoints)
        if self.DebugLevel >0: print("Number of data = ",len(self.DataTuples))
        if self.DebugLevel >2: print(self.DataTuples)
        
        if len(self.DataTuples)<=3:
            tkmb.askokcancel("No data found", "Check your input file")
            return
        self.TimeExper = np.array(self.TimeExper).astype(np.float64)
        self.RelElExper = np.array(self.RelElExper).astype(np.float64)
        self.ErrorExper = np.array(self.ErrorExper).astype(np.float64)
        if self.DebugLevel >2: print(self.TimeExper)
        if self.DebugLevel >2: print(self.RelElExper)
        if self.DebugLevel >2: print(self.ErrorExper)
        #print(self.TimeExper[0],self.RelElExper[0],self.ErrorExper[0] )
        self.UpdateValues()
        self.PlotMe()
        self.butPlot.config(state="normal")
        self.butEstim.config(state="normal")
        self.butA.config(state="normal")
        self.butB.config(state="normal")
        self.butC.config(state="normal")
        self.butD.config(state="normal")        
        self.butt0.config(state="normal")
        self.buttmax.config(state="normal")
        self.buttmin.config(state="normal")
        self.buttstep.config(state="normal")

                
    def cmdEstimate(self):
        """ Estimates initial values of parameters from user visual input """
        if self.inEstimate == 0:
            #starting estimation procedure
            self.B = self.RelElExper[0]
            self.C = self.RelElExper[self.Npoints-1]/2
            self.cid = self.myfig.canvas.mpl_connect('button_press_event', self.onClick)
        #scipy.misc.derivative(func, x0, dx=1.0, n=1, args=(), order=3)
            self.inEstimate = 1
            # inEstimate will be set to 2 in self.onClick
        #self.myfig.canvas.mpl_disconnect(cid)
            self.UpdateValues()
            self.PlotMe()
            tkmb.askokcancel("Press Done when done", "Click on the point with the maximum slope")
            self.butEstim.config(text = "Done")
        elif self.inEstimate == 2:
            # we are here only if user successfully clicks on screen and onClick callback is executed
            # inEstimate is set to 2 in onClick
            #print(self.TimeExper)
            if self.DebugLevel > 0: print("Esti in 2")
            # find x closest to t0
            self._min = 0
            for self._i, self._ti in enumerate(self.TimeExper):
                #print(self.TimeExper[self._i])
                if self.TimeExper[self._i]>self.te: break
                pass
            self._deri = (self.RelElExper[self._i-1]-self.RelElExper[self._i-2])/(self.TimeExper[self._i-1]-self.TimeExper[self._i-2])
            self._deri += (self.RelElExper[self._i]-self.RelElExper[self._i-1])/(self.TimeExper[self._i]-self.TimeExper[self._i-1])
            self._deri += (self.RelElExper[self._i+1]-self.RelElExper[self._i])/(self.TimeExper[self._i+1]-self.TimeExper[self._i])
            self._deri += (self.RelElExper[self._i+2]-self.RelElExper[self._i+1])/(self.TimeExper[self._i+2]-self.TimeExper[self._i+1])    
            self._deri *= 0.25
            self.D = (self._deri - self.A)/(self.RelElExper[self._i]-self.B)
            self.UpdateValues()
            self.PlotMe()
            self.butEstim.config(text = "Estimate")
            self.butfmin.config(state="normal")
            self.inEstimate == 3
            #run it twice to get some idea
            self.cmdfmin()
            self.cmdfmin()
        pass
   
        
    def calc_chi2(self, apar):
        """
Calculates chi2 of the current fit
Expects list of parameters in order A, B, C, D, te
Returns float
Optimized in simplex Nelder-Mead algorithm
ToDo: It needs cleaning to unify simples and LSQ calls
        """
        if self.DebugLevel >1: print("Call in calc_chi2")
        self.A = apar[0]
        self.B = apar[1]
        self.C = apar[2]
        self.D = apar[3]
        self.te = apar[4]       
        self.sqdiff = 0
        for self._i in range(0,self.Npoints,1):
            #print(self._i)
            self._diff = (self.__funcexex__(self.TimeExper[self._i]) - self.RelElExper[self._i])/self.ErrorExper[self._i]
            self.sqdiff += self._diff*self._diff
        if self.DebugLevel >2: print("In optimize sqdiff", self.sqdiff)
        self.lblchi2value.config(text="{0:>12g}".format(self.sqdiff))
        return self.sqdiff

    def funOptiLSQ(self,x,A,B,C,D,te):
        """
This function is a wrapper for LSQ fit
Makes sure that A, C and D do not fall below 0
ToDo: It needs cleaning to unify simples and LSQ calls
        """
        if self.DebugLevel >1: print("Call in OptiLSQ")
        #p0 = A, p1 = B, p2 = C, p3 = D, p4 = te
        if C<0:
            if self.DebugLevel>1: print("C under 0, correcting")
            C = np.float64(0.0000001)
        if D<0:
            D = np.float64(0.0000001)
            if self.DebugLevel>1: print("D under 0, correcting")          
        if self.B<0: self.B = np.float64(0.0)
        if A<0:
            A = np.float64(0.0)
            if self.DebugLevel>1: print("A under 0, correcting")        
        return A*x + B + self.funcGrow(x,C,D,te)
        
    def cmdLSQ(self):
        """
Callback funtion for button LSQ
Fits function using curve_fit from scipy.optimize
Plots fit result
Prepares LSQ report
        """
        #At + B + Cexp(-exp(-D(t-te)))
        # Store old parameters
        self.oldA = self.A
        self.oldB = self.B
        self.oldC = self.C
        self.oldD = self.D      
        self.oldte = self.te
        if self.DebugLevel >1: print("Before LSQ", self.calc_chi2([self.A, self.B, self.C, self.D, self.te]))
        # Main fitting loop using curve_fit
        # we use funOptiLSQ in order to keep physical values of parameters
        # This fragment needs cleaning as I had to use a new wrapper function funOptiLSQ
        # to match calling convention of curve_fit
        # curve_fit calls back funOptiLSQ with x as the first parameter and all other taken from list p0
        self.popt, self.pcov = sci_opti.curve_fit(self.funOptiLSQ,
                                                  self.TimeExper,
                                                  self.RelElExper,
                                                  p0=[self.A, self.B, self.C, self.D, self.te]
                                                  )
        self.A = self.popt[0]
        self.B = self.popt[1]
        self.C = self.popt[2]
        self.D = self.popt[3]
        self.te = self.popt[4]
        # Update chi2
        self.M = self.calc_chi2([self.A, self.B, self.C, self.D, self.te])
        self.UpdateValues()
        if self.DebugLevel >1: print("After LSQ", self.calc_chi2([self.A, self.B, self.C, self.D, self.te]))                                                  

        self.PlotMe(full=1) 
        # Prepare LSQ report
        self.LSQ4Report = []
        self._str = " **************** SOLUTION Least Squares ****************" 
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Fit from file: {}".format(self.filename)
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Fit chi2 with {}-5 degrees of freedom: {}".format(self.Npoints,self.M)
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Normalized chi2(S,GooF assuming you provided correct expeimental errors): {0:>12g}".format(self.M/(self.Npoints-5))
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Formula At + B + Cexp(-exp(-D(t-te)))"
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "Results of Levenberg-Marquardt:"
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "A = {0[0]:>12}, stddev {0[1]:>12} seconds^-1".format( self.RoundMe((self.A,np.sqrt(self.pcov[0][0]) ))) 
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "B = {0[0]:>12}, stddev {0[1]:>12} dimensionless".format( self.RoundMe((self.B,np.sqrt(self.pcov[1][1]) )))  
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "C = {0[0]:>12}, stddev {0[1]:>12} dimensionless".format( self.RoundMe((self.C,np.sqrt(self.pcov[2][2]) ) )) 
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "D = {0[0]:>12}, stddev {0[1]:>12} seconds^-1".format( self.RoundMe((self.D,np.sqrt(self.pcov[3][3]) ) ))
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)        
        self._str = "te = {0[0]:>12}, stddev {0[1]:>12} seconds ".format( self.RoundMe((self.te,np.sqrt(self.pcov[4][4])))  )
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "tint = {0[0]:>12} seconds ".format( self.tint)        
        self.addInfo(self._str)        
        self.LSQ4Report.append(self._str) 
        self._str = "stddev is one standard deviation and measures only statistical uncertainity"
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "In order to get a better estimate of errors compare it to the results from the simplex method"
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Estimated T2: {0[0]:>12} +/- {0[1]:>12} s".format(
            self.RoundMe((1/self.D,1/self.D**2 * np.sqrt(self.pcov[3][3]) ) ))
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)
        self._str = "Estimated (P-Y)n0: {0[0]:>12} +/- {0[1]:>12} MPa".format(
            self.RoundMe((self.A/10**(-6),np.sqrt(self.pcov[0][0])/10**(-6) ) ))
        self.addInfo(self._str)
        self.LSQ4Report.append(self._str)          
        pass
    
#    def LSQCycle(self):
        #At + B + Cexp(-exp(-D(t-t0)))
#        pass
    
    def cmdfmin(self, maxiter=1000):
        """
Performs initial fit using simplex Nelder-Mead algorithm
Maximum number of iterations can be set by giving 1 parameter maxiter (defalut maxiter=1000)
        """
        #At + B + Cexp(-exp(-D(t-t0)))
        self.oldA = self.A
        self.oldB = self.B
        self.oldC = self.C
        self.oldD = self.D       
        self.oldte = self.te
        if self.DebugLevel >1: print("Before optimize", self.calc_chi2([self.A, self.B, self.C, self.D, self.te]))
        sci_opti.fmin(self.calc_chi2, [self.A, self.B, self.C, self.D, self.te], maxiter=maxiter)
        self._lastopti = self.calc_chi2([self.A, self.B, self.C, self.D, self.te])
        if self.DebugLevel >1: print("After optimize", self._lastopti)
        self.UpdateValues()
        self.Simplex4Report = []
        self._str = " *********** SOLUTION SIMPLEX ***************** "
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "Fit from file: {}".format(self.filename)
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "Fit chi2: {:>12g}".format(self._lastopti)
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "Formula At + B + Cexp(-exp(-D(t-t0))) "
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "Results of simplex minimization (no error estimates):"
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "A = {:>12g} seconds^-1".format(self.A) 
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "B = {:>12g} dimensionless".format(self.B) 
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "C = {:>12g} dimensionless".format(self.C) 
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)        
        self._str = "D = {:>12g} seconds^-1".format(self.D) 
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)
        self._str = "te = {:>12g} seconds".format(self.te) 
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)       
        self._str = "tint = {:>12g} seconds".format(self.tint)         
        self.addInfo(self._str)        
        self.Simplex4Report.append(self._str)       
        self._str = "Estimated T2: {0:>12g} s".format(1/self.D )
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)
        self._str = "Estimated (P-Y)*n0: {0:>12g}".format(self.A/10**(-6))
        self.addInfo(self._str)
        self.Simplex4Report.append(self._str)

        self.butLSQ.config(state="normal")
        self.PlotMe(full=1)
        self.butSaveReport.config(state="normal")


        pass
        
    def cmdSaveReport(self):
        """ This method saves fit results as a txt and HTML files\ntogether with one png figure """
        #self._i = self.filename.rfind( '/', 0)
        #self.trunk = self.filename[:self._i+1]
        #self.filen = self.filename[self._i+1:]
        self.filen = os.path.basename(self.filename)
        self.trunk = os.path.dirname(self.filename)
        # cut off the extension
        self._i = self.filen.rfind( '.', 0)
        self.filen = self.filen[:self._i]
        self.figname = os.path.join(self.trunk,"Fig1_"+self.filen + ".png")
        self.htmlname = os.path.join(self.trunk,"Web_"+self.filen + ".html")
        self.reportname = os.path.join(self.trunk,"Report_"+self.filen + ".txt")
        plt.savefig(self.figname, dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format="png",transparent=False, bbox_inches=None, pad_inches=0.1)  
        self.fhtml = open(self.htmlname,'w')
        self.fhtml.write("<HTML><TITLE>Report</TITLE><BODY>\n<CENTER><P>\n")
        self.fp = open(self.reportname,'w')
        for self._i in self.Simplex4Report:
            self.fp.write(self._i + os.linesep)
            self.fhtml.write(self._i + "<BR>"+os.linesep)
        self.fp.write(os.linesep)
        self.fhtml.write("<BR></P><P>"+os.linesep) 
        for self._i in self.LSQ4Report:
            self.fp.write(self._i + os.linesep)
            self.fhtml.write(self._i + "<BR>"+os.linesep)        
        self.fhtml.write("<IMG SRC=" + self.figname +">"+os.linesep)        
        self.fp.close()
        self.fhtml.write("\n</CENTER></P></BODY>\n<HTML>")
        self.fhtml.close()
    
    def onClick(self,event):
        """ 
Callback from click event on Figure
Method captures user on-screen clik in order to estimate function inflection point 
Should only be called during estimation

"""
        #we should only get here if user clicks butEstimate
        #cmdEstimate
        if self.inEstimate == 1:
            self.xpoint0 = event.x
            self.xpoint0 = event.y
            self.te = event.xdata
            self.UpdateValues()
            self.PlotMe()
            self.inEstimate = 2
            # stop capturing
            self.myfig.canvas.mpl_disconnect(self.cid)
            pass
        elif self.inEstimate == 2:
            #should never get here
            self.addInfo("Esti 2. should not get here")
        else:
            if self.DebugLevel > 0: print('button={}, x={}, y={}, xdata={}, ydata={}'.format(event.button, event.x, event.y, event.xdata, event.ydata))

    
if __name__ == "__main__":
    app = clsFitexex(None)
    app.title('Relative growth fit P. Zajdel, A Haduch-Sendecka, M. Pietruszka 2014')
    # based on non-git version 3.0
    #Ver 4.0
    app.mainloop()
