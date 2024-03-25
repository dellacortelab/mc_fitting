import Tkinter, tkFileDialog
import random
import subprocess as sb
import math
try:
    import umat_new as umat
except:
    print 'NO UMAT FOUND, PLEASE COMPILE VIA GUI'
from itertools import combinations 
from Tkinter import N,S,E,W
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from ttk import Frame, Button, Style, Widget
import numpy as np
import time
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


################################################################
# 
#                A basis Monte Carlo Curve Fitting Tool
#                   Version 0.5 for default params
#        Author: Dennis Della Corte (dellacorte86@gmail.com)
#        Affiliation: IFAM Aachen, Forschungszentrum Juelich
#
################################################################

class MC_GUI(Frame):
    '''The GUI class. It can help an unexperienced terminal user to quickly
       load in experiment and to use the most essential monte carlo curve
       fitting capabilities to identify material model properties.
       It will allways be loaded first and the needs to be populated with:
       1. umat
       2. properties from a file
       3. an experiment from a file
        
    '''
    def __init__(self, master):
        '''Initializing method'''
        # Create a container
        self.master = master
        frame = Tkinter.Frame(self.master)
        self.plotframe=frame
        self.frameButtons = Frame(frame, borderwidth=1)
        self.frameParams = Frame(frame, borderwidth=1)
        self.framePlot = Frame(frame, borderwidth=1)
        self.frameCtrl2 = Frame(frame, borderwidth=1)
        self.frameCtrl = Frame(self.frameCtrl2, borderwidth=1)
        self.frameCtrlT = Frame(self.frameCtrl2, borderwidth=1)

        self.plotframe.grid(row=0,column=0,rowspan=3,columnspan=10,sticky=N+S+E+W)
       
        self.frameButtons.grid(row=0,column=0,columnspan=10,sticky=N+S+E+W)
        self.frameParams.grid(row=1,column=0,rowspan=2,columnspan=2,sticky=N+S+E+W)
        self.framePlot.grid(row=1,column=2,rowspan=2,columnspan=6,sticky=N+S+E+W)
        self.frameCtrl2.grid(row=1,column=8,rowspan=2,columnspan=2,sticky=N+S+E+W)
        self.frameCtrl.grid(row=0,column=0,columnspan=2,sticky=N+S+E+W)
        self.frameCtrlT.grid(row=1,column=0,columnspan=2,sticky=N+S+E+W)
        frame.rowconfigure(0, weight=1)
        frame.rowconfigure(1, weight=1)
        for i in range(10):        
            frame.columnconfigure(i, weight=1)
        
        self.frameButtons.rowconfigure(0,weight=1)
        self.frameButtons.columnconfigure(0,weight=1)
        self.frameParams.rowconfigure(0,weight=1)
        self.frameParams.rowconfigure(1,weight=1)
        self.frameParams.columnconfigure(0,weight=1)
        self.frameParams.columnconfigure(1,weight=1)
        self.framePlot.rowconfigure(0,weight=1)
        self.framePlot.rowconfigure(1,weight=1)
        self.framePlot.columnconfigure(0,weight=1)
        self.plotframe.rowconfigure(0,weight=1)
        self.plotframe.rowconfigure(1,weight=1)
        self.plotframe.rowconfigure(2,weight=1)
        for i in range(10):
            self.plotframe.columnconfigure(i,weight=1)
        self.frameCtrl.rowconfigure(0,weight=1)
        self.frameCtrl.columnconfigure(0,weight=1)
        self.frameCtrl.columnconfigure(1,weight=1)
        self.frameCtrl2.columnconfigure(0,weight=1)
        self.frameCtrl2.columnconfigure(1,weight=1)
        #self.frameCtrl2.rowconfigure(0,weight=1)
        for idx in range(30):
            self.frameCtrl2.rowconfigure(idx,weight=1)

        self.frameCtrlT.rowconfigure(0,weight=1)
        self.frameCtrlT.columnconfigure(0,weight=1)
        self.frameCtrlT.columnconfigure(1,weight=1)

        self.initLoad()
        # Create buttons
        self.button_load = Tkinter.Button(self.frameButtons
          ,text="Load Props",command=self.load_props_button)
        self.button_mc = Tkinter.Button(self.frameButtons
          ,text="Start Monte Carlo Search",command=self.monte_carlo_button)
        self.reset_button = Tkinter.Button(self.frameButtons,
                             text="Update Plot", command=self.update_plot_button)
        self.openfile_button = Tkinter.Button(self.frameButtons, 
                             text='Open Experiment', command=self.open_exp_button)
        self.print_button = Tkinter.Button(self.frameButtons, 
                             text='Print Props', command=self.print_current_params)
        self.barplot_button = Tkinter.Button(self.frameButtons, 
                             text='Show Barplot', command=self.barplot)
        self.randexp_button = Tkinter.Button(self.frameButtons, 
                             text='Biaxial Experiment', command=self.randexp)
        self.randexp2_button = Tkinter.Button(self.frameButtons, 
                             text='Uniaxial Experiment', command=self.randexp2)
        self.compile_umat_button = Tkinter.Button(self.frameButtons, 
                             text='Compile Umat', command=self.compile_umat_but)
        self.button_mc.pack(side="left",expand=1)
        self.barplot_button.pack(side="left",expand=1)
        self.reset_button.pack(side="left",expand=1)
        self.randexp2_button.pack(side="left",expand=1)
        self.randexp_button.pack(side="left",expand=1)
        self.compile_umat_button.pack(side="left",expand=1)
        self.button_load.pack(side="left",expand=1)
        self.openfile_button.pack(side='left',expand=1)
        self.print_button.pack(side='left',expand=1)
        
        
        #create labels and entry fields for monte carlo runs
        
        #mc_steps 
        self.mc_steps = Tkinter.StringVar()
        self.mc_steps.set('1000')
        Tkinter.Label(self.frameCtrl,text='MC Steps').grid(\
                        row=0,column=0,sticky=N+E+S+W) 
        self.mc_step_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.mc_steps) 
        self.mc_step_getter.grid(row=0, column=1,sticky=N+E+S+W)
        
        #mc temperature 
        self.temperature = Tkinter.StringVar()
        self.temperature.set('0.05')
        Tkinter.Label(self.frameCtrl,text='MC Temp').grid(\
                        row=1,column=0,sticky=N+E+S+W) 
        self.mc_temp_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.temperature) 
        self.mc_temp_getter.grid(row=1, column=1,sticky=N+E+S+W)
        
        #mc stepsize
        self.mc_stepsize = Tkinter.StringVar()
        self.mc_stepsize.set('0.1')
        Tkinter.Label(self.frameCtrl,text='MC Stepsize').grid(\
                        row=2,column=0,sticky=N+E+S+W) 
        self.mc_stepsize_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.mc_stepsize) 
        self.mc_stepsize_getter.grid(row=2, column=1,sticky=N+E+S+W)

        #annealing steps
        self.anealing_steps = Tkinter.StringVar()
        self.anealing_steps.set('1')
        Tkinter.Label(self.frameCtrl,text='SimAneal Steps').grid(\
                        row=3,column=0,sticky=N+E+S+W) 
        self.mc_astep_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.anealing_steps) 
        self.mc_astep_getter.grid(row=3, column=1,sticky=N+E+S+W)


        #temp multiplier
        self.temperature_multiplier = Tkinter.StringVar()
        self.temperature_multiplier.set('1')
        Tkinter.Label(self.frameCtrl,text='Temp. Multi').grid(\
                        row=4,column=0,sticky=N+E+S+W) 
        self.mc_mtemp_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.temperature_multiplier) 
        self.mc_mtemp_getter.grid(row=4, column=1,sticky=N+E+S+W)

        #step multiplier
        self.stepsize_multiplier = Tkinter.StringVar()
        self.stepsize_multiplier.set('1')
        Tkinter.Label(self.frameCtrl,text='Step Multi').grid(\
                        row=5,column=0,sticky=N+E+S+W) 
        self.mc_mstep_getter = Tkinter.Entry(self.frameCtrl,textvariable =\
                        self.stepsize_multiplier) 
        self.mc_mstep_getter.grid(row=5, column=1,sticky=N+E+S+W)


        #create checkboxes for axis display 
        self.axis_vars = [ Tkinter.IntVar() for i in range(6) ]
        self.axis_button = [ Tkinter.Checkbutton(self.frameCtrlT, text="Component:"+str(i)\
                        , variable=self.axis_vars[i]) for i in range(6) ]
        offset = 5 # after other variables are set this is the legnth of them
        for idx in range(6):
            self.axis_button[idx].grid(row=idx+offset, column=0,sticky=N+E+S+W)
            self.axis_button[idx].select() 
       
        self.status = Tkinter.StringVar() 
        self.status.set('')
        Tkinter.Label(self.frameCtrlT,text='Select for Analysis').grid(\
                        row=offset-1,column=0,sticky=N+E+S+W) 
        Tkinter.Label(self.frameCtrlT,textvariable=\
                self.status).grid(row=14,column=0,sticky=N+E+S+W) 


        #make the main figures
        self.fig = plt.Figure()
        self.ax  = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig,master=self.framePlot)
                
        self.outputfname = './output_props.dat'

        # Add a menu
        menubar = Tkinter.Menu(master)
        master.config(menu=menubar)
        fileMenu = Tkinter.Menu(menubar)
        fileMenu.add_command(label="Exit", command=self.onExit)
        menubar.add_cascade(label="File", menu=fileMenu)


    def compile_umat_but(self):
        global umat
        '''Tricky function that can load in a new umat, if none is specified'''
        umat_fname = self.askopenfilename()
        sb.call('f2py -m umat_new -c '+umat_fname,shell=True)
        try:
            import umat_new as umat
        except:
            print 'ERROR!!! PLEASE CHECK COMPILER FOR HINTS!!!'

    def randexp2(self):
        '''dummy caller'''
        self.randexp(False)    

    def randexp(self, biax = True):
        '''generates randomly an experiment that can be used for fitting
           testing'''
        s = Simulated_Experiments(180)
        if biax:
            strains = s.biaxial(non_iso = True)   #make a random biaxial test
        else:
            strains = s.uniaxial()   #make a random biaxial test
        d = Load_props('./inputs/faser_input.dat')
        u_t = Umat_wrapper(Material_model(d))
        u_t.test_inputs(strains)
        strains, stresses =  u_t.get_strstres()
        self.open_exp(strains, stresses)

    def barplot(self):
        '''Makes the barplot'''
        self.reset_plot()
        keys, ctrs = self.opti.barplot()

        ctrs = np.array(ctrs)

        N = len(keys)
        ind = np.arange(1,2*N+1,2)  # the x locations for the groups
        width = 0.4       # the width of the bars

        ax = self.ax
        rects1 = ax.bar(ind, ctrs, width, color='r')

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Number of Cost Changes')
        ax.set_title('Importance of Parameters for Current Fit')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( keys )

        def autolabel(rects, m ):
            # attach some text labels
            idx = 0
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., m[idx]+(m[idx]/\
                            max(1,abs(m[idx]))*.1), '%.2f'%round(m[idx],2),
                        ha='center', va='bottom')
                idx += 1
        autolabel(rects1, ctrs)
        self.ax.set_ylim([min(ctrs)*1.1,max(ctrs)*1.1])
        plt.setp( ax.xaxis.get_majorticklabels(), rotation=70 )
        self.canvas.show()


    def monte_carlo_button(self):
        '''perform monte carlo simulations, plots most recent calculation 
           afterwards'''
        #get these numbers form the fields
        
        #Get a status parameter up
        self.status.set('MONTE CARLO IS RUNNING!')

        self.opti.set_duration(int(float(self.mc_steps.get())))
        self.opti.set_stepsize(float(self.mc_stepsize.get()))
        self.opti.set_temperature(float(self.temperature.get()))
        self.opti.simulated_annealing(int(float(self.anealing_steps.get())) \
            ,float(self.temperature_multiplier.get())\
            , int(float(self.mc_steps.get())),False,\
            float(self.stepsize_multiplier.get()))
        self.plot_calc()
        self.update_displayed_props()
        self.status.set('')   
 
    def set_axis(self):
        '''set the axis for the current operations, can be checked in boxes'''
        axis = []
        for idx in range(6):
            if self.axis_vars[idx].get() == 1:
                axis.append(idx)
        if len(axis) < 1:
            print 'WARNING: NO COMPONENT SELECTED FOR ANALYSIS!'
        self.axis = axis
        self.opti.set_axis(axis)

    def update_displayed_props(self):
        '''change the displayed params'''
        for key, val in self.matmod.get_sorted_props():
            self.propsInputs[key].set(val) 

    def update_plot_button(self):
        '''make sure that current props are read in'''
        self.update_props()
        self.plot_calc()

    def plot_calc(self):
        '''plots the currently avaibale approximation'''
        self.reset_plot()
        #make sure that current props are loaded
        data = self.umatw.cost_calculations(True, axis = self.axis)
        idx = 0
        for x,y,cost in data:
            self.ax.plot(x,y,'--',label='Fit: '+str(self.axis[idx])+\
            ' Cost: %2.2e' % cost,linewidth = 4)
            idx+=1
        self.plot_ref()        

    def open_exp_button(self):
        '''opens an experimental file as specified by the user and plots it'''
        try:
            self.umatw = Umat_wrapper(self.matmod)
            self.opti = Optimization(self.umatw)
            data = np.loadtxt(self.askopenfilename())
            strains = data[:,0:6]
            stresses = data[:,6:]
            self.opti.set_comparison(strains, stresses)
            #SET which axis to plot
            self.set_axis()
            self.plot_calc()
        except:
            print 'LOAD IN PROPERTIES FIRSTS!'


    def open_exp(self,strains, stresses):
        '''opens an experiment with user defined strains and stresses'''
        self.umatw = Umat_wrapper(self.matmod)
        self.opti = Optimization(self.umatw)
        self.opti.set_comparison(strains, stresses)
        #SET which axis to plot
        self.set_axis()
        self.plot_calc()

    def plot_ref(self):
        '''plot the references into the GUI'''
        data = self.opti.plot_refs()
        idx= 0
        cols = ['black','grey','brown','darkblue','darkcyan','darkgoldenrod'] 
        for x,y in data:
           self.ax.plot(x,y,label='Ref: '+str(self.axis[idx]),linewidth =3,color=cols[idx])
           idx += 1
        self.ax.set_xlabel('Strain')
        self.ax.set_ylabel('Stress')
        self.ax.set_title('Monte Carlo Curve Fitting')
        self.show_plot()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
 

    def get_prop_val(self, key):
        '''get the value of  a certain property'''
        return self.opti.u.m.getv_only(key)

    def set_prop_val(self, key, val):
        '''set the value of a certain key'''
        self.opti.u.m.setv_only(key, val)
        

    def show_properties(self):
        '''after loading props display these on the left of gui'''
        #make sure that properties are StringVars in the gui 
        es = []
        idx = 0
        self.propsInputs = {}
        sorted_props = self.matmod.get_sorted_props()
        dprops = self.matmod.parameters
        for key, val in sorted_props:
            v, t,mi,ma,ctr = dprops[key]
            self.propsInputs[key] = Tkinter.StringVar()
            self.propsInputs[key].set(str(v))
            txt = key + '['+str(mi)+':'+str(ma)+']'
            Tkinter.Label(self.frameParams, text=txt).grid(row=idx,column=0 \
                            ,sticky=N+E+S+W) 
            es.append(Tkinter.Entry(self.frameParams,textvariable =  \
                      self.propsInputs[key]))
            es[-1].grid(row=idx, column=1,sticky=N+E+S+W)
            idx += 1
       #make sure to set all properties acording to the fields before doing any
       #calculation

    def update_props(self):
        '''call this before starting claculation to pass the current params
           into the matmod'''
        for key in self.propsInputs:
            val = self.propsInputs[key].get() 
            self.set_prop_val(key, float(val))

    def load_props_button(self):
        '''Load in the Properties from user selection'''
        props = Load_props(self.askopenfilename())
        d = props.get_start_parameters()
        try:
            for key in d:
                v, t, m, ma, ctr = d[key]
                self.umatw.setv(key,v,t,m,ma,ctr)
            self.update_displayed_props()
            self.update_props()    
            self.matmod = Material_model(props)
            self.show_properties()
            self.plot_calc()
        except:
            self.matmod = Material_model(props)
            self.show_properties()
  
    def load_props(self, fname):
        '''Load in the Properties from specified test file fname'''
        self.matmod = Material_model(Load_props(fname))
        self.show_properties()
 
    def print_current_params(self):
        '''an easy way to get the current set of params to terminal'''
        self.update_props()
        print '\n-------------------------------------'
        print 'Current Parameter Set\n'
        cprops = self.matmod.get_sorted_props() #make output sorted
        with open(self.outputfname,'w') as f: 
            for key,v in cprops: 
                v, t, m, ma, ctr = self.matmod.parameters[key]
                t = t.__name__ 
                print '%12s\t%3.3f\t%8s\t%3.3f\t%3.3f\t%5d' % (key,v,t,m,ma,ctr)
                f.write('%12s\t%3.3f\t%8s\t%3.3f\t%3.3f\t%5d\n' % (key,v,t,m,ma,ctr))
        print '\n-------------------------------------'

        


    def show_plot(self):
        '''just to have a uniform showing'''
        self.ax.set_xlim([.7,1.9])
        #self.ax.set_ylim([0,1100])
        self.ax.legend(loc='upper left',prop={'size':15})
        self.canvas.show()

    def reset_plot(self):
        '''clear the figure before new drawing'''
        self.set_axis()
        self.ax.clear()
        self.show_exp = False
        self.cost_ref = 0




    def onSelect(self, val):
        '''reset values of params'''
        sender = val.widget
        idx = sender.curselection()
        value = sender.get(idx)
        self.var.set(value)
        
 

    def onExit(self):
        '''closes the frame for good'''
        self.master.destroy()

    def initLoad(self):
        '''defines the defaults for the open file dialog'''
        # define options for opening or saving a file
        self.file_opt = options = {}
        options['defaultextension'] = '.dat'
        #options['filetypes'] = [('all files', '.*'), ('text files', '.f')]
        options['initialdir'] = 'C:\\'
        options['initialfile'] = 'myfile.dat'
        options['parent'] = self.master
        options['title'] = 'Open File'

        # defining options for opening a directory
        self.dir_opt = options = {}
        options['initialdir'] = 'C:\\'
        options['mustexist'] = False
        options['parent'] = self.master
        options['title'] = 'Open Directory'


    def askopenfilename(self):
        """Returns the fname of selected file"""
        return tkFileDialog.askopenfilename(**self.file_opt)



#functions for convenience  
#all the stuff from MCCP
def compile_umat():
    '''if you forget, this is how you can compile the umat'''
    sb.call('f2py -m umat -c ElasticAnisotropicStructTensors.f',shell=True)
    print 'Compiled the umat once more'

def lagrangian_strain(dfgrd1):
    '''calculate lagrangian strain from deformation gradient'''
    strain = np.zeros((3,3),float)
    delta = np.zeros((3,3),float)
    delta[0,0] =1 
    delta[1,1] =1 
    delta[2,2] =1 

    for i in range(3):
        for j in range(3):
            for k in range(3):
                strain[i,j] += dfgrd1[k,i] * dfgrd1[k,j]
            strain[i,j] = 1/2. * (strain[i,j] - delta[i,j])
    return strain

def cauchy_green_strain(dfgrd1):
    '''calculate right cauchy green strain from deformation gradient'''
    strain = np.zeros((3,3),float)
    delta = np.zeros((3,3),float)
    delta[0,0] =1 
    delta[1,1] =1 
    delta[2,2] =1 

    for i in range(3):
        for j in range(3):
            for k in range(3):
                strain[i,j] += dfgrd1[k,i] * dfgrd1[k,j]
    return strain

def pk2(cauchy, dfgrd1):
    '''calculate the 2 Piola Kirchhoff Stress fom Cauchy Stress'''
    
    pk2 = np.zeros((3,3),float)
    idfgrd1 = umat.invert(dfgrd1,3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for r in range(3):
                    pk2[i,j] += idfgrd1[i,r] * cauchy[r,k] * idfgrd1[j,k]      
    return pk2

def pk2_voigt(cauchy, dfgrd1):
    '''calculate the 2 Piola Kirchhoff Stress fom Cauchy Stress, return voigt'''
    
    pk2 = np.zeros((3,3),float)
    idfgrd1 = umat.invert(dfgrd1,3)
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for r in range(3):
                    pk2[i,j] += idfgrd1[i,r] * cauchy[r,k] * idfgrd1[j,k]      
    return [pk2[0,0], pk2[1,1], pk2[2,2], pk2[0,1], pk2[0,2], pk2[1,2] ]


def mat6n33a(a6):
      '''! c....transfer 6-dim. Voigt notation to 3x3 matrix'''
      b9 = np.zeros((3,3),float)
      for i in range(3):
        b9[i,i]=a6[i]
      
      b9[0,1]=a6[3]
      b9[1,0]=a6[3]
      b9[1,2]=a6[5]
      b9[2,1]=a6[5]
      b9[2,0]=a6[4]
      b9[0,2]=a6[4]
      return b9


def load_experimental_data(filename):
    '''reads an x y text file'''
    x,y = np.loadtxt(filename,float) 
    return x,y

        
class Material_model():
    '''Material model calss. It requieres and instance of the properties class
       to be initalized. It provides an interface to set and get properties'''
    def __init__(self, props_instance):
        '''define the paramters of the material model'''
        self.parameters = props_instance.get_start_parameters()
        self.mapping = props_instance.get_prop_mapping()
        self.nprops = len(self.parameters)
        self.loaded_data = props_instance

    def get_sorted_props(self):
        '''returns the properties in the order provided in the input props file'''
        return self.loaded_data.get_start_sortet_parameters()
   
    def get_keys(self):
        '''returns all keys in parameter dictionary'''
        return self.parameters.keys()
 
    def reset(self, start_props):
        '''quick way to get back to original constants'''
        self.__init__(start_props)

    def getv(self, key):
        '''gets the value of a certain entry'''
        return self.parameters[key]  

    def getv_only(self, key):
        '''gets the value of a certain entry'''
        return self.parameters[key][0]

    def setv(self, key, value, ty = float, mini = 0, maxi = 1, nr = 0):
        '''sets a sertian value'''
        self.parameters[key] = (value, ty, mini, maxi, nr)
    
    def setv_only(self, key, value):
        '''sets a sertian value'''
        val, ty, mini, maxi, nr = self.getv(key)
        self.parameters[key] = (value, ty, mini, maxi, nr)
    
    def get_umat_input(self):
        '''just returns the sorted list of parameters for the umat call'''
        props = []
        for idx in range(self.nprops):
            props.append( float(self.parameters[ self.mapping[idx] ][0] ))
        return props

class Umat_wrapper():
    '''Class that handles the interfacing with the umat code. Is capable of calling
       fortran code and retrieving stress results from the umat calculaton. It 
       requieres an instance of a material model to be initialized.'''
    def __init__(self, matmod):
        '''init function, creates the instance'''
        self.load_model(matmod)

    def get_props_dic(self):
        '''return the dictionary of properties from the underlying material model'''
        return self.m.parameters

    def setv(self,key, value, ty = float, mini = 0, maxi = 1, nr = 0):
        '''set a value in the underlying material model'''
        self.m.setv(key, value, ty, mini, maxi, nr )

    def getv(self, key):
        '''get a value in the underlying material model'''
        return self.m.getv(key)

    def get_strstres(self):
        '''get the strains and stresses from the loaded experiment'''
        return self.expstrains, self.expstresses

    def plot(self):
        '''after loading plot the data'''
        plt.xlim(1,1.75)
        plt.ylim(0,1100)
        for i in range(2):
            plt.plot(self.expstrains[:,i],self.expstresses[:,i],label='Experiment '+str(i),linewidth=3)

    def load_model(self, matmod):
        '''loads in an instance of the material model parameters'''
        self.m = matmod

    def get_model(self):
        '''retrieves material model'''
        return self.m

    def call_umat(self, dfgrd1, props):
        '''calls the umat with the Gcurrent values in the material model'''
        stress= np.array([ 0. for i in range(6) ] )
        umat.myumat(stress,props,dfgrd1,len(props)) 
        return stress

    def load_inputs(self, inp):
        '''loads an input file of strains and stresses'''
        a = np.loadtxt(inp, float, skiprows= 1)
        self.expstrains  = a[:,:6]
        self.expstresses = a[:,6:]
        self.dfgrds = [ self.convert(strain) for strain in self.expstrains ]

    def exp_inputs(self, strains, stresses):
        '''throw in strains and stresses directly'''
        self.expstrains  = np.array(strains)
        self.expstresses = np.array(stresses)
        self.dfgrds = [ self.convert(strain) for strain in self.expstrains ]

    def test_inputs(self, strains):
        '''throw in dformation gradients for a certain test'''
        self.dfgrds = [ self.convert(strain) for strain in strains ]
        #create experiment stresses and strains based on the default params
        props = self.m.get_umat_input()
        calcstresses = np.array([pk2_voigt(mat6n33a( self.call_umat(df, props) ) ,df) for df in self.dfgrds ] )
        self.expstresses = np.array(calcstresses)
        self.expstrains  = np.array(strains)

    def helper(self, p, *strains):
        '''helper function that solve a system of equations'''
        x11,x12,x13,x22,x23,x33 = p
        strains = strains[0]
        return ( x11**2 +  x12**2  + x13**2  - strains[0], 
                 x12**2 +  x22**2  + x23**2  - strains[1],
                 x13**2 +  x23**2  + x33**2  - strains[2],
                 x11*x12 + x12*x22 + x13*x23 - strains[3],
                 x11*x13 + x12*x23 + x13*x33  -strains[4],
                 x12*x13 + x22*x23 + x23*x33 - strains[5]
                )
            


    def convert(self, strains):
        '''converts strains into deformgradient'''
        #assume that there is no rotation, then F = U, U^2 = C, and C is
        #only diagonal
        x11,x12,x13,x22,x23,x33 =  \
                    fsolve(self.helper, (1, 0, 0, 1, 0, 1), args=strains)
        dfgrad = np.zeros((3,3))
        dfgrad[0][0] = x11 
        dfgrad[0][1] = x12 
        dfgrad[0][2] = x13 
        dfgrad[1][0] = x12 
        dfgrad[1][1] = x22 
        dfgrad[1][2] = x23 
        dfgrad[2][0] = x13 
        dfgrad[2][1] = x23 
        dfgrad[2][2] = x33 
        return dfgrad
            

    def cost_calculations(self, get_plots = False, axis = []):
        '''will calculate the cost for a paramter set with a given input,
           if axis = (0) gets cost of 0 axis, can also be (0,1,2)'''
        props = self.m.get_umat_input()
        #the cost is now the difference between the calculated and the 
        #given stresses
        calcstresses = np.array([pk2_voigt(mat6n33a( self.call_umat(df, props) ) ,df) for df in self.dfgrds ] )
        #here I could specify a certain axis to get only the desired cost function of that axis 
        if len(axis) == 0:    
            print 'WARNING! NO COMPONENT SPECIFIED!!!'''
            cost = 0
            if get_plots:
                return ([0],[0],0)
            
        else:
            cost = 0
            for i in axis:
                cost += np.sum((self.expstresses - calcstresses)**2,0)[i]
            if get_plots:
                data = []
                for ax in axis:
                    cost = np.sum((self.expstresses - calcstresses)**2,0)[ax]
                    data.append((self.expstrains[:,ax], calcstresses[:,ax], cost))
                return data

        return cost

class Convert_experiments():
    '''a class that can deal with biaxial experiments, and writes them 
        into usefull data outputs, for biaxial experiments'''
    def __init__(self):
        pass

    def load_exp(self, filen):
        '''loads in a file'''
        self.data = np.loadtxt(filen, skiprows = 1, dtype =float)
        

    def make_full_strain_stresses(self):
        'calculate the deforamtion fradients from the loaded strains and stresses'
        strains = []
        stresses = []
        dfgrd = np.zeros((3,3))
        for line in self.data:
            dfgrd[0,0] = math.sqrt(line[0])
            dfgrd[1,1] = math.sqrt(line[1])
            dfgrd[2,2] = 1/math.sqrt(line[0]*line[1])
            c = cauchy_green_strain(dfgrd)
            strains.append([c[0,0],c[1,1],c[2,2],c[0,1],c[0,2],c[1,2]])
            stresses.append([ line[2],line[3],0,0,0,0 ])
        return np.array(strains), np.array(stresses)


def plot(x,y,idx):
    '''convenience function to show results of an experiment'''
    xs = x[:,idx]
    ys = y[:,idx]
    plt.xlabel('Cauchy Green Strain')
    plt.ylabel('2Piola Kirchhoff Stress')
    plt.title('Curve Fitting for Reese et al Experiments')
    plt.plot(xs,ys,label='umat Calculation '+str(idx),linewidth=3)


#now make a class whoes sole purpose it is to test paramters for a test
# in a smart fashion

class Kringing():
    '''Brute fore search for most influencial parameters. Was replaced with
       Monte Carlo Counting instead'''
    def __init__(self,axis=False):
        '''An parameter test for material model and certain material'''
        #the axis for cost calculation
        self.axis = axis

    def set_umat(self, umat):
        '''load in the umat'''
        self.u = umat
        self.m = umat.get_model()
    
    def test_params(self, combis = 1):
        '''Investigate the impact of changing single parameters. Start out
           with all of them set to the middle value   '''
        ps = self.m.get_keys()
    
        #store the old parameters as well as set them all to middle val
        backup = {}
        for key in ps:
            default, dtype, mini, maxi, ctr = self.m.getv(key)
            backup[key] = [default, dtype, mini, maxi, ctr]
            if dtype == int:
                val = int((maxi+mini)/2.)
            else:
                val = (maxi + mini) / 2.
            self.m.setv(key, val , dtype, mini, maxi, ctr)
            
 
        combis = combinations(ps, combis)
        #consider to investigate the change due to combination of paramters 
        influence_dic = {}
        for key in ps:
            default, dtype, mini, maxi, ctr = self.m.getv(key)
            if dtype == int:
                stepsize = 1
                steps = maxi - mini
            else:
                steps = 10
                stepsize = (maxi - mini) / float(steps) 
            currentp = mini
            costs = []
            for idx in range(steps):
                currentp += stepsize
                self.m.setv(key, currentp)
                costs.append(self.u.cost_calculations(axis=self.axis))
            influence_dic[key] = np.std(costs)
            #pylab.plot(range(len(costs)),costs,label=key)
            self.m.setv(key, default, dtype, mini, maxi, ctr)
        #pylab.legend(loc='upper right')
        #pylab.show()
        for key in backup:
            default, dtype, mini, maxi, ctr = backup[key]
            self.m.setv(key, default, dtype, mini, maxi, ctr)
        return influence_dic

    def barplot(self,influence_dic):
        '''creates a barplot after parameter testing'''
        keys = influence_dic.keys()
        values = np.array(influence_dic.values())
        idx = np.where(values > 0.5)
        values[idx] = np.log(values[idx])


        N = len(keys)
     
        menMeans = [ val for val in values ]
        menStd =   [ 0 for val in values ]

        ind = np.arange(1,N+1)  # the x locations for the groups
        width = 0.4       # the width of the bars

        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)

        # add some text for labels, title and axes ticks
        ax.set_ylabel('Log of Variance of Costchange')
        ax.set_title('Influence of Individual Parameters over Range of Allowed Values')
        ax.set_xticks(ind+width)
        ax.set_xticklabels( keys )
        #ax.legend( rects1[0], 'Test' )

        def autolabel(rects, m ):
            # attach some text labels
            idx = 0
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%.2f'%round(m[idx],2),
                        ha='center', va='bottom')
                idx += 1
        autolabel(rects1, menMeans)
        plt.setp( ax.xaxis.get_majorticklabels(), rotation=70 )
        plt.show()

class Load_props():
    '''Class that makes it more convenient to load in a properties file and to handle
       the contents, requieres a path to a filename'''
    def __init__(self, fname):
        '''Initializer'''
        self.data = np.loadtxt(fname, str,skiprows=1)
        self.set_prop_mapping()
        self.nprops = len(self.data)
        self.set_start_parameters()    

    def set_prop_mapping(self):
        '''props need to be mapped to entry in props_list for umat call'''
        self.props_map = {}
        for idx in range(len(self.data)):
            self.props_map[idx] = self.data[idx][0]
        
    def get_prop_mapping(self):
        '''gets the mapping between propertie entry and argument index for umat call'''
        return self.props_map


    def set_start_parameters(self):
        '''retrieves the start values for each of the parameters, order as input'''
        self.prop_dic = {}
        for prop,startval,dtype,mini,maxi in self.data:
            if dtype == 'int': #for int
                self.prop_dic[prop] = ( int(float(startval)), int, int(float(mini)), int(float(maxi)), 0)
            else:
                self.prop_dic[prop] = ( float(startval), float, float(mini), float(maxi), 0)
             
            
    def get_start_sortet_parameters(self):
        '''returns the start parameters in the order that the umat want'''
        svals = []
        for idx in range(self.nprops):
            prop = self.props_map[idx]
            svals.append( (prop, self.prop_dic[prop][0]))
        return svals

    def get_start_parameters(self):
        '''returns the dicitonary of properties'''
        return self.prop_dic

class Simulated_Experiments():
    '''Simple class that can perfrom uni or biaxial experiments on a given 
       properties set. The only purpose is to demonstrate how an external 
       function can be used to generate data that can be used for curve fitting.
       This part needs to be replaced with actual calls to computer experiments
       on finite element models'''
    def __init__(self, max_strain = 200):
        '''A class that serves a container for possible tests. Each test
           is associated with a list of strains for parameter
           testing'''
        self.max_strain = max_strain


    def random(self):
        '''generates an totally arbitray strain state, where all but the
           6 component are defined. The last one is calculated based on
           isochoric condition'''
        #!DOE SNOT WORK! 
        end_strain = [ random.randrange(80,self.max_strain) for i in range(6) ]
        #now move along from all one to final state
        start = [ 1 for i in range(6) ]
        strains = []
        steps = [ (end_strain[idx] - start[idx])/(100.-self.max_strain) for idx in range(6) ]
        for idx in range(100, self.max_strain):
           strains.append( [ start[idx] + steps[idx] * idx for idx in range(6) ])
        return strains

    def uniaxial(self, axis = 0):
        '''Strain the Model along the given axis, assume isochoric behaviour
           and isotropy in the other directions as an approximation'''
        if axis == 0:
            strain = lambda x: [ x,1/math.sqrt(x),1/math.sqrt(x),0. ,  0.,   0.  ]
        elif axis == 1: 
            strain = lambda x: [ 1/math.sqrt(x),x,1/math.sqrt(x),0. ,  0.,   0.  ]
        elif axis == 2: 
            strain = lambda x: [ 1/math.sqrt(x),1/math.sqrt(x),x,0. ,  0.,   0.  ]
        return [ strain(x/100.) for x in range(100, self.max_strain) ]


    def biaxial(self,axis1 = 0, axis2 = 1, non_iso = False):
        '''Strain the Model uniform along the given axis, assume 
           isochoric behaviour
           and isotropy in the other direction as an approximation'''
        if non_iso:
            strain = lambda x: [ x,x/1.1,1/float(x*x/1.1),0. ,  0.,   0.  ]
        elif axis1 == 0 and axis2 == 1:
            strain = lambda x: [ x,x,1/float(x**2),0. ,  0.,   0.  ]
        else:
            print 'INVALID TEST'
            return
        return [ strain(x/100.) for x in range(100, self.max_strain) ]

class Optimization():
    '''The meat of the code. A class that can perform simulated annealing of the
       Monte Carlo Method to find local minima of propertie combinations. The cost
       function is defined as the pointwise square function between input stresses and
       calculated stresses'''
    def __init__(self, umat):
        '''An instance of an Obptimization class, this class is able to
           perform changes on a given dictionary of properties and will
           yield the ideal collection of parameters for a given experi-
           mental curve'''
        self.to_optimize = umat.get_props_dic()
        self.mc_steps = 100
        self.temp     = .05 #might need to be adjusted
        self.u     = umat
        self.stepsize = 1/100. #one percent at the beginning?
    
    def set_stepsize(self, stepsize):
        '''define the stepsize for a propertie. (1 if int, 0.01 if float by default)'''
        self.stepsize = stepsize

    def barplot(self):
        '''returns a the names of parameters and the ctr value for a barplot'''
        props = self.u.get_props_dic()
        sort = self.u.m.get_sorted_props()
    
        xs = []
        ys = []
        for key, v in sort:
            xs.append(key)
            ys.append(props[key][-1])
        return xs, ys

    def simulated_annealing(self, howoften = 10, howmuch = 0.9, mc_steps = 1000, exclusion_threshold=False, change_stepsize = 0.9):
        '''use simulated annealing to go lower on the cost surface.
           will repeat monte carlo simulations'''
        self.mc_steps = mc_steps
        for idx in range(howoften):
            self.monte_carlo(exclusion_threshold = exclusion_threshold)
            self.temp *= howmuch
            self.stepsize *= change_stepsize
            print 'Finished Annealing Step', idx, 'of',howoften

    def set_comparison(self, strains, stresses):  
        '''Load in the strain - stress curve that will be used for 
           comparison. The strain will serve as basis for the defor-
           mation gradients that will be used in the evaluation of the
           umat''' 
        self.ref_strains = np.array(strains)
        self.ref_stresses = np.array(stresses)
        self.u.exp_inputs(strains, stresses)
        self.to_optimize = self.u.m.get_keys()        

    def get_props_dic(self):
        '''gets the underlying properties dictionary'''
        return self.u.get_props_dic()

    def set_temperature(self, temp):
        '''Set the starting temperature for simualted annealing'''
        try:
            self.temp = float(temp)
        except:
            print 'INVALID TEMPERATURE!', temp

    def set_duration(self, duration):
        '''Set the number of MC Steps that will be performed before decrease
           in temperature'''
        try:
            self.mc_steps = int(float(duration))
        except:
            print 'INVALID DURATION!', duration

    def set_optimization_set(self, optimize):
        '''define the list of parameters that should be optimized, all otherwise'''
        self.to_optimize = optimize
    
    def set_axis(self, axis):
        '''set the axis you want to investigate'''
        self.axis = axis


    def exclude_props(self, exclude):
        '''takes a list of properties that shall be excluded from oprimization'''
        print 'EXCLUDE FROM SEARCH:', exclude
        for item in exclude:
            self.to_optimize.remove(item)

    def exclude_negative(self, threshold = 50):
        '''search through properties, if one has a ctr > threshold, remove from
           optimization list'''
        exclude = []
        for key in self.get_props_dic():
            val, dtype, mini, maxi, ctr = self.u.getv(key)
            if ctr > threshold:
                exclude.append(key)
        self.exclude_props(exclude)

    def monte_carlo(self, exclusion_threshold = False):
        '''run a monte carlo simulation, it will use simulated annealing
           to find local minima of the property landscape based on cost
           function evaluation for all interesting directions
           Stress can either be in 0, 1 or 0 and 1 direction
           If exclusion is an integer after each MC run parameters are 
           tested for impact. Only necesssary at the beginning'''


        kbtlike= self.temp
        epsilon = 0.1 
        ACCEPTED = 0
        COSTS = []
        old_cost = float('inf')
        for idx in range(self.mc_steps):

            change = random.choice(self.to_optimize)
            val, dtype, mini, maxi, ctr = self.u.getv(change)

            if dtype == int:
                stepsize = 1
            else:
                stepsize = (maxi - mini) * float(self.stepsize) #change in 1% steps
 
            signum = random.choice([-1,1])
            ival = val + signum * stepsize 
           
            if ival < mini or ival > maxi: 
                #went to far, so turn back!
                ival = val + -1* signum * stepsize 
                
    
            self.u.setv(change,ival, dtype, mini, maxi, ctr)            

            #make sure that the correct cost is calculated by specifing
            #the axis
            new_cost = self.u.cost_calculations(get_plots=False, axis = self.axis)
            
            #the monte carlo step
            dE = new_cost - old_cost
            if abs(dE) < epsilon:
                #count the properties that do not change Energy!
                ctr -= 1
            else:
                ctr += 1 
            self.u.setv(change,ival, dtype, mini, maxi, ctr)            
            if dE < 0:
                #accepted step
                old_cost = new_cost
                ACCEPTED += 1
                COSTS.append(old_cost)            

            elif math.exp(- dE /float(kbtlike) ) > random.uniform(0,1):
                #accepted step
                old_cost = new_cost
                ACCEPTED += 1
                COSTS.append(old_cost)            
            else:
                #undo the changes
                self.u.setv(change, val, dtype, mini, maxi, ctr)            

 
            if idx % 50 == 0:
                print 'Step ',idx, 'of', int(self.mc_steps), old_cost
                

        #dict_tester(DDD.get_start_parameters(),self.get_props_dic())
        print 'Acceptence Ratio:', ACCEPTED/float(self.mc_steps)
        if exclusion_threshold:
            self.exclude_negative(exclusion_threshold)

    def plot_refs(self):
        '''returns the x and y components for as many plots as specified in axis'''
        ret = []
        for ax in self.axis:
            ret.append([self.ref_strains[:,ax], self.ref_stresses[:,ax]])
        return ret

def dict_tester(dic_a, dic_b):
    '''sees if mc cost correlates with getting closer to right answers
       dic_a is real value'''
    diff = 0
    for key in dic_a:
            val, dtype, mini, maxi, ctr = dic_a[key]
            val2, dtype, mini, maxi, ctr = dic_b[key]
            diff += ((val - val2))**2
    cost = math.sqrt(diff)
    print 'Current Distance is', cost
    return cost

def get_optimizer(props_file, exp_file):
    '''Return an optmisation object to run command from cmd line on'''
    d = Load_props(props_file)
    u_t = Umat_wrapper(Material_model(d))
    opti = Optimization(u_t)
    data = np.loadtxt(exp_file)
    strains = data[:,0:6]
    stresses = data[:,6:]
    opti.set_comparison(strains, stresses)
    return opti    

font = {'weight' : 'bold',
    'size'   : 18}

plt.rc('font', **font)
        
WIDTH = 1400
HEIGHT = 800
def start_gui():
    '''a Call to this function will start the graphical user interface'''
    root = Tkinter.Tk()
    root.geometry("%dx%d+0+0" % (WIDTH, HEIGHT))
    root.columnconfigure(0,weight=1)
    root.columnconfigure(1,weight=1)
    root.rowconfigure(1,weight=1)
    root.rowconfigure(0,weight=1)
    gui = MC_GUI(root)
    gui.master.title("DCMC Curve Fitting - version .8")
    try:
        #script this for conveniece
        s = Simulated_Experiments(180)
        strains = s.uniaxial(0) #make a uniaxial test - gets strains as if pulling in one direction only
        #strains = s.biaxial()   #make a biaxial test
        d = Load_props('./inputs/faser_input.dat')
        u_t = Umat_wrapper(Material_model(d))
        u_t.test_inputs(strains)
        strains, stresses =  u_t.get_strstres()

        gui.load_props('./inputs/faser_input_corrupt.dat')
        gui.open_exp(strains, stresses)
    except:
        print 'YOU NEED TO COMPILE YOUR UMAT FIRST, USE GUI BUTTON'
    root.mainloop()
if __name__ == '__main__':
    start_gui()

