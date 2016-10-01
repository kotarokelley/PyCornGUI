#!/usr/bin/python

'''
pyCornGUI.py
A GUI editor to extract and edit data from .res files generated
by UNICORN Chromatography software supplied with AKTA Systems.
by Kotaro Kelley 2015 kotarokelley1@gmail.com
'''
import glob
import os
import sys
import traceback
from collections import OrderedDict
    
    
try:
    import wx
    import wx.lib.buttons as buttons
    #from wx.lib.agw import ultimatelistctrl as ULC      # for file browser
except:
    ImportError
    print 'ERROR: wx not found - Abborting Program!'
    sys.exit()
    
try:
    import matplotlib
    matplotlib.use('WXAgg')
    from matplotlib.ticker import AutoMinorLocator
    import matplotlib.pyplot as plt
    plotting = True
    from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
    from matplotlib.backends.backend_wx import NavigationToolbar2Wx
    from matplotlib.figure import Figure
except:
    ImportError
    print("WARNING: Matplotlib not found - Plotting disabled!")
    plotting = False
    
try:
    # import pyCornLibraries
    from pycorn import pc_res3
except:
    ImportError
    print("WARNING: pycorn not found - Plotting disabled!")
    plotting = False

########################################################################
class PlotPanel(wx.Panel):
    """
    Panel for plotting matplotlib graphs. Called by MainPanel and is located in the upper left corner. 
    
    Class Members:
        self.figure(matplotlib.figure.Figure):
            Main matplot figure object. 
        self.axes(subplot):
            Main axes to which other axes are added with the twix method. 
        self.canvas(FigureCanvas):
            This allows our figure to be output and updated as a widget in the PlotPanel. 

    Class Methods: 
        Display:
            Displays a single trace file to screen. 
        mapper:
            Calculate position of given percent relative to min and max. 
        expander:    
            Expand -/+ direction of two values by a percentage of their delta.
        xy_data:
            Takes a data block and returns two lists with x- and y-data      
        smartscale:
            Uses input data and user parameters to determine scaling of x/y-axis.
        Plotter:
            Uses input data and user parameters to plot trace data.
    """
    
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialize the class members. 
        """
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.EXPAND|wx.ALL)
        self.SetSizer(self.sizer)
        self.Fit()

        self.styles = {'UV':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV1_':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV2_':{'color': (0.898,0.086,0.086), 'lw': 1.4, 'ls': "-", 'alpha':1.0},
        'UV3_':{'color': (0.780,0.239,0.902), 'lw': 1.2, 'ls': "-", 'alpha':1.0},
        'Cond':{'color': (1,0.486,0.161), 'lw': 1.4, 'ls': "-", 'alpha':0.75},
        'Inje':{'color': (0.835,0.427,0.616), 'lw': 1.0, 'ls': "-", 'alpha':0.75},
        'UV_Other':{'color': (0.224,0.059,0.137), 'lw': 1.0, 'ls': "-", 'alpha':0.5},
        'UV1_1':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV2_1':{'color': (0.898,0.086,0.086), 'lw': 1.4, 'ls': "-", 'alpha':1.0},
        'UV3_1':{'color': (0.780,0.239,0.902), 'lw': 1.2, 'ls': "-", 'alpha':1.0},
        'Cond1':{'color': (1,0.486,0.161), 'lw': 1.4, 'ls': "-", 'alpha':0.75},
        'UV_Other1':{'color': (0.224,0.059,0.137), 'lw': 1.0, 'ls': "-", 'alpha':0.5},
        'UV1_2':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV2_2':{'color': (0.898,0.086,0.086), 'lw': 1.4, 'ls': "-", 'alpha':1.0},
        'UV3_2':{'color': (0.780,0.239,0.902), 'lw': 1.2, 'ls': "-", 'alpha':1.0},
        'Cond2':{'color': (1,0.486,0.161), 'lw': 1.4, 'ls': "-", 'alpha':0.75},
        'UV_Other2':{'color': (0.224,0.059,0.137), 'lw': 1.0, 'ls': "-", 'alpha':0.5},
        'UV1_3':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV2_3':{'color': (0.898,0.086,0.086), 'lw': 1.4, 'ls': "-", 'alpha':1.0},
        'UV3_3':{'color': (0.780,0.239,0.902), 'lw': 1.2, 'ls': "-", 'alpha':1.0},
        'Cond3':{'color': (1,0.486,0.161), 'lw': 1.4, 'ls': "-", 'alpha':0.75},
        'UV_Other3':{'color': (0.224,0.059,0.137), 'lw': 1.0, 'ls': "-", 'alpha':0.5},
        'UV1_4':{'color': (0.098,0.098,1), 'lw': 1.6, 'ls': "-", 'alpha':1.0},
        'UV2_4':{'color': (0.898,0.086,0.086), 'lw': 1.4, 'ls': "-", 'alpha':1.0},
        'UV3_4':{'color': (0.780,0.239,0.902), 'lw': 1.2, 'ls': "-", 'alpha':1.0},
        'Cond4':{'color': (1,0.486,0.161), 'lw': 1.4, 'ls': "-", 'alpha':0.75},
        'UV_Other4':{'color': (0.224,0.059,0.137), 'lw': 1.0, 'ls': "-", 'alpha':0.5}
        }

        
    #----------------------------------------------------------------------
    def Display(self, fname, traceParms, dispParms, saveParms):         
        """
        Displays a single trace file to screen. 

        Args:
            fname(string):
                The absolute path to file. Must have .res suffix to avoid raising exeption downstream. 
            traceParms:
                Dictionary of parameters used to decide which traces are displayed. 
                This field should come from CommandPanel.GetTraceParms()
            dispParms:
                Dictionary of display parameters. This field should come from CommandPanel.GetDispParms()
            saveParms:
                Dictionary of save parameters. This field for now is determined by the various output buttons in CommandPanel.
                At minimum, it should contain key 'save'. 
        """  
        if not isinstance(fname, list):         # check for single data
            try: 
                fdata = pc_res3(fname)          # load data
                fdata.load()                    # parse data from .res file
                self.figure.clf()                           # clear figure
                self.axes = self.figure.add_subplot(111)    # new axes   
                self.Plotter(fdata, fname, traceParms, dispParms, self.axes)
                if dispParms['title']:
                    self.figure.suptitle(fname.split('/')[-1].split('.')[0])    
                self.canvas.draw()
                self.canvas.Refresh() 
                
                if saveParms['Save']:                       # save figure. 
                    baseName = os.path.basename(fname)      # get just the basename from the absolute path that was passed as fname
                    baseName = baseName[:-4] + '.' + saveParms['format']    # replace .res with .format
                    outDir = self.GetParent().Params['outDir']              # get output directory from MainPanel
                    outPath = os.path.join(outDir, baseName)
                    self.figure.savefig(outPath, bbox_inches='tight', dpi=300)
            except Exception:
                tb = traceback.format_exc()                 # print traceback to diagnose source of error
                print tb
                print "ERROR: could not display data."
        else:                                   # multiple data to be overlayed
            try:
                fdata = [pc_res3(name) for name in fname]   # keep datasets in a list
                for fdat in fdata:
                    fdat.load()
                self.figure.clf()                           # clear figure
                self.axes = self.figure.add_subplot(111)    # new axes   
                self.PlotterOverlay(fdata, fname, traceParms, dispParms, self.axes)
                if dispParms['title']:
                    title = 'Overlay '
                    for name in fname:
                        title += name.split('/')[-1].split('.')[0]
                        title += '\n'
                    self.figure.suptitle(title)    # TODO: make this an option in future through display popup  
                self.canvas.draw()
                self.canvas.Refresh() 

                if saveParms['Save']:                       # save figure. 
                    title = 'Overlay '
                    for name in fname:
                        title += os.path.basename(name)[4]      # get just the basename from the absolute path that was passed as fname
                        title += '_'
                    title += title + '.' + saveParms['format']    # replace .res with .format
                    outDir = self.GetParent().Params['outDir']              # get output directory from MainPanel
                    outPath = os.path.join(outDir, title)
                    self.figure.savefig(outPath, bbox_inches='tight', dpi=300)
                
            except AssertionError: 
                #print "Please select more than one .res file to overlay"
                pass
            except Exception:
                tb = traceback.format_exc()                 # print traceback to diagnose source of error
                print tb
                print "ERROR: could not overlay data"
                
    '''
    The following functions mapper, expander, xy_data, and smartscale were 
    borrowed and modified by Kotaro Kelley on 15/02/10 from:
    PyCORN - script to extract data from .res (results) files generated
    by UNICORN Chromatography software supplied with AKTA Systems
    (c)2014-2015 - Yasar L. Ahmed
    v0.16
    ''' 
    #----------------------------------------------------------------------
    def mapper(self, min_val, max_val, perc):
        
        '''
        Modified from code by Yasar L. Ahmed
        Calculate position of given percent relative to min and max. 
        The order if min and max is not checked. 
        
        Args:
            min_val(float):
                min value
            max_val(foat):
                max value
            perc(float<=1):
                percent
        Returns:
            Position of given percent relative to min and max. 
        '''
        x = abs(max_val - min_val) * perc
        if min_val < 0:
            return (x - abs(min_val))
        else:
            return (x + min_val)
    
    #----------------------------------------------------------------------   
    def expander(self, min_val, max_val, perc):
        '''
        Modified from code by Yasar L. Ahmed
        Expand -/+ direction of two values by a percentage of their delta.
        The order of min and max is not checked. 
                
        Args:
            min_val(float):
                min value
            max_val(foat):
                max value
            perc(float<=1):
                percent
        Returns:
            Tuple of min_value-delta*perc, max_value-delta*perc
                    
        '''
        delta = abs(max_val - min_val)
        x = delta * perc
        return (min_val - x, max_val + x)
    
    #----------------------------------------------------------------------
    def xy_data(self, inp):
        '''
        Modified from code by Yasar L. Ahmed
        Takes a data block and returns two lists with x- and y-data
        
        Args:
            A list of tuples. Only the first two elemets of the tuples will be returned. 
        Returns:
            Two lists containing 0 elements and 1 elements respectively. 
        '''
        x_data = [x[0] for x in inp]
        y_data = [x[1] for x in inp]
        return x_data, y_data
    
    #----------------------------------------------------------------------
    def smartscale(self, inp, args):
        '''
        Modified from code by Yasar L. Ahmed
        Uses input data and user parameters to determine scaling of x/y-axis.
        This function should only take in data from a single run so that the x_values are consistent.
        This will not be checked though. 
        
        Args:
            inp:
                Dictionary of trace data.
            args:
                Dictionary of user input xmin and xmax. 
        Returns:
            Best min/max for x/y
        Raises:
            KeyError:
                When trying to access 'data' keys that don't exist. 
                When trying to access 'Fraction' key that does not exist. 
        '''
        x_data = []         # hold lists of parsed xdata
        y_data = []         # hold lists of parsed ydata    
        try:
            for key in inp.keys():
                xdata, ydata = self.xy_data(inp[key]['data'])   # for each trace, parse data into x, y
                x_data.append(xdata)                            # get all the raw data in selected fdata
                y_data.append(ydata)
        except KeyError:
            print 'ERROR: tried to input invalid data structure to smartscale.'
        try:
            frac_data = inp['Fractions']['data']           # see if there is fraction data available
            frac_x, frac_y = self.xy_data(frac_data)
            frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]   # distance between each data point
            frac_delta.append(frac_delta[-1])                               
        except KeyError:  
            frac_data = None
        
        if 'xmin' in args and args['xmin'] != None:       # requested xmin. Should't get to 2nd statement if args is empty
            plot_x_min = args['xmin']
        elif frac_data:
            plot_x_min = frac_data[:][0]                # no requested xmin, get value from fraction data
        else:
            plot_x_min = min([dat[0] for dat in x_data])                   # if no fraction data, get starting value from the first trace      
        if 'xmax' in args and args['xmax'] != None:                            # requested xmax
            plot_x_max = args['xmax']
        elif frac_data:       
            plot_x_max = frac_data[-1][0] + frac_delta[-1]*2 # no requested xmax, get starting value from first fraction data
        else:
            plot_x_max = max([dat[-1] for dat in x_data])                  # if no fraction data, get starting value from the first trace
                
        if plot_x_min > plot_x_max: # This shlouldn't happen since it is checked when inputing xmin/xmax into args, but just in case.  
            print("Warning: xmin bigger than xmax - adjusting...")      # Also it is safe to assume that traces are in order from min to max. 
            plot_x_min = min([dat[0] for dat in x_data])                   #get starting value from first fraction data
        if plot_x_max < plot_x_min:
            print("Warning: xmax smaller than xmin - adjusting...")
            plot_x_max = max([dat[-1] for dat in x_data])
            
        # optimize y_scaling
        min_y_values = []
        max_y_values = []
        for trace in range(len(x_data)):                   # for each trace
            tmp_x = x_data[trace]
            tmp_y = y_data[trace]
            range_min_lst = [abs(a - plot_x_min) for a in tmp_x]
            range_min_idx = range_min_lst.index(min(range_min_lst))
            range_max_lst = [abs(a - plot_x_max) for a in tmp_x]
            range_max_idx = range_max_lst.index(min(range_max_lst))
            values_in_range = tmp_y[range_min_idx:range_max_idx]
            min_y_values.append(min(values_in_range))
            max_y_values.append(max(values_in_range))
        plot_y_min = min(min_y_values)
        plot_y_max = max(max_y_values)
        #plot_y_min, plot_y_max = self.expander(plot_y_min_tmp, plot_y_max_tmp, 0.085)
        return plot_x_min, plot_x_max, plot_y_min, plot_y_max
    
    #----------------------------------------------------------------------
    def plot(self, inp, ax, parms):
        for key, value in inp.items():
            x_dat, y_dat = self.xy_data(value['data'])
            print("Plotting: " + value['data_name'])
            stl = self.styles[key[:4]]
            p0, = ax.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'],
                            ls=stl['ls'], lw=stl['lw'],alpha=stl['alpha'], linewidth = parms['lineThickness'])
    #----------------------------------------------------------------------
    def Plotter(self,inp,fname, traceParms, dispParms, axes):
        """
        Uses input data and user parameters to plot trace data.
                
        Args:
            inp:
                Dictionary of trace data.
            fname:
                List of file names to be processed. Has to be absolute path. 
            traceParsms:
                Dictionary describing which traces to plot. 
            dispParms:
                Dictionary of display parameters.
            axes:
                Axes object where all outputs  and subplots will be sent to. 
        Raises:
            Exception:
                When no traces are requested. 
        """
        
        if not traceParms:              # raise exception if no traces are requested
            raise Exception('No traces requested') 
                
        inpUV = {key: value for key, value in inp.items() if (key.startswith('UV') and not key.endswith('_0nm'))}   # get all UV traces
        inpUV_280_260 = {}  
        inpUV_220 = {}
        inpUV_Other = {}                        
        
        for key, value in inpUV.items():                            # parse out different traces only if user has requested them
            try:
                try:
                    numkey = float(key.split('_')[-1].split('nm')[0]) # get the numerical value of wavelength used for off peak 220 measurements
                except ValueError:
                    numkey = -1         # this line is for when trace does not have a numerical value. This happens on some single wavelength Akta.
                if '280' in key :
                    if traceParms['280nm']:
                        inpUV_280_260.update({key:value})
                elif '260' in key:
                    if traceParms['260nm']:
                        inpUV_280_260.update({key:value})
                elif 215<=numkey<=225:                              # allow for off peak measuremnt of protein backbone
                    if traceParms['220nm']:
                        inpUV_220.update({key:value})
                elif traceParms['OtherUV']:                         # other wavelengths. 
                    inpUV_Other.update({key:value})
            except Exception:
                tb = traceback.format_exc()
                print tb
                print "ERROR: Something went wront when trying to plot data, %s" , key
        
        inpCond = {key:value for key, value in inp.items() if 'Cond' in key and not key.endswith('%') and traceParms[key[-5:]]}
        
        hNo220, lNo220 = [], []                         # hold values for use in common legend
        h220, l220 = [], []
        hOther, lOther = [], []
        hCond, lCond = [], []
        hOther, lOther = [], []
        
        plot_x_min, plot_x_max = None, None  # if there are non 220 traces then it will be set then. Otherwise, it will be set by 220, and so on. 
        
        if dispParms['show_xlabel']:
            axes.set_xlabel("Elution volume (ml)")
        if dispParms['show_ylabel']:
            axes.set_ylabel("Absorbance (mAu)")
        # show ticks
        axes.tick_params(axis='x', which='both', top='off')
        axes.tick_params(axis='y', which='both', right='off')
        if not dispParms['xticks']:
            axes.tick_params(axis='x', which='both', bottom='off')
        if not dispParms['xticklabel']:
            axes.tick_params(axis='x', which='both', labelbottom='off')
        if not dispParms['yticks']:
            axes.tick_params(axis='y', which='both', left='off', right='off')
        if not dispParms['yticklabel']:
            axes.tick_params(axis='y', which='both', labelleft='off', labelright='off')
            
        if inpUV_280_260:        #skip this plot if no 280 or 260 is not requested
            plot_x_min, plot_x_max, plot_y_min, plot_y_max = self.smartscale(inpUV_280_260, traceParms) # get scale without 220 trace
            plot_y_min , plot_y_max = self.expander(plot_y_min, plot_y_max, 0.085)
            axes.set_xlim(plot_x_min, plot_x_max)
            axes.set_ylim(plot_y_min, plot_y_max)
            
            for key, value in inpUV_280_260.items():
                x_dat, y_dat = self.xy_data(value['data'])
                print("Plotting: " + value['data_name'])
                stl = self.styles[key[:4]]
                p0, = axes.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'],
                                ls=stl['ls'], lw=stl['lw'],alpha=stl['alpha'], linewidth = dispParms['lineThickness'])

            hNo220, lNo220 = axes.get_legend_handles_labels()

        if inpUV_220:               # skip plot of 220 trace if not requested
            par1 = axes.twinx()     # make separate axis for displaying 220 data
            stl = self.styles['UV3_']
            
            if dispParms['show_ylabel']:
                par1.set_ylabel("Absorbance (mAu)", color=stl['color'])
                
            p1_xmin, p1_xmax, p1_ymin, p1_ymax = self.smartscale(inpUV_220, traceParms)
            p1_ymin, p1_ymax = self.expander(p1_ymin, p1_ymax, 0.085)
            if plot_x_min == None or plot_x_max == None:    # x limit has not been set yet 
                plot_x_min, plot_x_max = p1_xmin, p1_xmax   
                par1.set_xlim(plot_x_min, plot_x_max)
                par1.set_ylim(p1_ymin, p1_ymax)
            else:                                           # go with the x limits that have already been set
                par1.set_xlim(plot_x_min, plot_x_max)
                par1.set_ylim(p1_ymin, p1_ymax) 
            for key, value in inpUV_220.items(): 
                x_dat, y_dat = self.xy_data(value['data'])
                print("Plotting: " + value['data_name'])
                p1, = par1.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                                ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'],linewidth = dispParms['lineThickness'])
            # show ticks
            par1.tick_params(axis='x', which='both', top='off')
            par1.tick_params(axis='y', which='both', right='off')
            if not dispParms['xticks']:
                par1.tick_params(axis='x', which='both', bottom='off')
            if not dispParms['xticklabel']:
                par1.tick_params(axis='x', which='both', labelbottom='off')
            if not dispParms['yticks']:
                par1.tick_params(axis='y', which='both', left='off', right='off')
            if not dispParms['yticklabel']:
                par1.tick_params(axis='y', which='both', labelleft='off', labelright='off')

                
            h220, l220 = par1.get_legend_handles_labels()

        if inpCond:                 # check if Cond is requested and if it exists in the data
            par2 = axes.twinx()     # make separate axis for displaying Cond data
            stl = self.styles['Cond']
            
            if dispParms['show_ylabel']:
                par2.set_ylabel('Cond (mS/cm)', color=stl['color'])
            
            p2_xmin, p2_xmax, p2_ymin, p2_ymax = self.smartscale(inpCond, traceParms)
            p2_ymin, p2_ymax = self.expander(p2_ymin, p2_ymax, 0.085)
            if plot_x_min == None or plot_x_max == None:    # x limit has not been set yet
                plot_x_min, plot_x_max = p2_xmin, p2_xmax           
                par2.set_xlim(plot_x_min, plot_x_max)
                par2.set_ylim(p2_ymin, p2_ymax)
            else:                                           # go with the x limits that have already been set
                par2.set_xlim(plot_x_min, plot_x_max)
                par2.set_ylim(p2_ymin, p2_ymax) 
            for key, value in inpCond.items():
                x_dat, y_dat = self.xy_data(value['data'])
                print("Plotting: " + value['data_name'])
                p2, = par2.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                                ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'], linewidth = dispParms['lineThickness'])
            # show ticks
            par2.tick_params(axis='x', which='both', top='off')
            par2.tick_params(axis='y', which='both', right='off')
            if not dispParms['xticks']:
                par2.tick_params(axis='x', which='both', bottom='off')
            if not dispParms['xticklabel']:
                par2.tick_params(axis='x', which='both', labelbottom='off')
            if not dispParms['yticks']:
                par2.tick_params(axis='y', which='both', left='off', right='off')
            if not dispParms['yticklabel']:
                par2.tick_params(axis='y', which='both', labelleft='off', labelright='off')
                
            hCond, lCond = par2.get_legend_handles_labels()
            
        if inpUV_Other:
            # make decision about where to plot this other data
            axisOptions = ['280_260', '220', 'Cond', 'New Axis']       # hard code some output axis
            otherAxes = None                                           # get the user prefered axis
            if dispParms['plototherWith'] == axisOptions[0]:
                otherAxes = axes                                       # 280_260 axes always exists by construction 
            elif dispParms['plototherWith'] == axisOptions[1]:
                try:
                    otherAxes = par1                                   # assign to 220 axes if it exists
                except NameError:
                    otherAxes = axes.twinx()                           # otherwise, make new axes
            elif dispParms['plototherWith'] == axisOptions[2]:
                try:
                    otherAxes = par2                                   # assign to Cond axes if it exists
                except NameError:
                    otherAxes = axes.twinx()                           # otherwise, make new axes
            else:
                otherAxes = axes.twinx()
           
            stl = self.styles['UV_Other']
            
            if dispParms['show_ylabel']:
                otherAxes.set_ylabel('Absorbance (mAu)', color=stl['color'])
                
            other_xmin, other_xmax, other_ymin, other_ymax = self.smartscale(inpUV_Other, traceParms)
            other_ymin, other_ymax = self.expander(other_ymin, other_ymax, 0.085)
            if plot_x_min == None or plot_x_max == None:            # x limit has not been set yet
                plot_x_min, plot_x_max = other_xmin, other_xmax
                otherAxes.set_xlim(plot_x_min, plot_x_max)
                otherAxes.set_ylim(other_ymin, other_ymax)
            else:                                                   # go with the x limits that have already been set
                otherAxes.set_xlim(plot_x_min, plot_x_max)
                otherAxes.set_ylim(other_ymin, other_ymax) 
            for key, value in inpUV_Other.items():
                x_dat, y_dat = self.xy_data(value['data'])
                print("Plotting: " + value['data_name'])
                p1, = otherAxes.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                            ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'],linewidth = dispParms['lineThickness'])
            
            # show ticks
            otherAxes.tick_params(axis='x', which='both', top='off')
            otherAxes.tick_params(axis='y', which='both', right='off')
            if not dispParms['xticks']:
                otherAxes.tick_params(axis='x', which='both', bottom='off')
            if not dispParms['xticklabel']:
                otherAxes.tick_params(axis='x', which='both', labelbottom='off')
            if not dispParms['yticks']:
                otherAxes.tick_params(axis='y', which='both', left='off')
            if not dispParms['yticklabel']:
                otherAxes.tick_params(axis='y', which='both', labelleft='off')
                
            hOther, lOther = otherAxes.get_legend_handles_labels()
            hOther = [i for i in hOther if i not in hNo220]         # avoid redundant labels if plotting on with 280/260 data
            lOther = [i for i in lOther if i not in lNo220]
            
        # plot fraction data
        if dispParms['show_fractions']:
            frac_data = inp['Fractions']['data']
            frac_x, frac_y = self.xy_data(frac_data)
            frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
            frac_delta.append(frac_delta[-1])
            frac_y_pos = self.mapper(axes.get_ylim()[0], axes.get_ylim()[1], 0.015)
            
            if dispParms['eof']:
                frac_data = frac_data[::2]
            
            for ind, i in enumerate(frac_data):
                axes.axvline(x=i[0], ymin=0.065, ymax=0.0, color='r', linewidth=0.85)
                axes.annotate(str(i[1]), xy=(i[0] + frac_delta[ind] * 0.55, frac_y_pos),
                         horizontalalignment='center', verticalalignment='bottom', size=8, rotation=90) 
        # show legend
        if dispParms['show_legend']:
            axes.legend(hNo220+h220+hCond+hOther, lNo220+l220+lCond+lOther, loc='upper right', fontsize=8, fancybox=True )
            axes.xaxis.set_minor_locator(AutoMinorLocator())
            axes.yaxis.set_minor_locator(AutoMinorLocator())
            

            
    #----------------------------------------------------------------------
        
    def PlotterOverlay(self,inp,fname, traceParms, dispParms, axes):
        """
        Uses input data and user parameters to plot trace data.
                
        Args:
            inp:
                Dictionary of trace data.
            fname:
                List of file names to be processed. Has to be absolute path. 
            traceParsms:
                Dictionary describing which traces to plot. 
            dispParms:
                Dictionary of display parameters.
            axes:
                Axes object where all outputs  and subplots will be sent to. 
        Raises:
            Exception:
                When no traces are requested. 
        """
        
        if not traceParms:              # raise exception if no traces are requested
            raise Exception('No traces requested') 
        
        inpUV = [{key: value for key, value in trace.items() if (key.startswith('UV') and not key.endswith('_0nm'))} for trace in inp]   # get all UV traces
        inpUV_280_260 = [{} for i in range(len(inpUV))]
        inpUV_220 = [{} for i in range(len(inpUV))]
        inpUV_Other = [{} for i in range(len(inpUV))]                   
        for i, trace in enumerate(inpUV):
            for key, value in trace.items():                            # parse out different traces only if user has requested them
                try:
                    try:
                        numkey = float(key.split('_')[-1].split('nm')[0]) # get the numerical value of wavelength used for off peak 220 measurements
                    except ValueError:
                        numkey = -1         # this line is for when trace does not have a numerical value. This happens on some single wavelength Akta.
                    if '280' in key :
                        if traceParms['280nm']:
                            inpUV_280_260[i].update({key:value})
                    elif '260' in key:
                        if traceParms['260nm']:
                            inpUV_280_260[i].update({key:value})
                    elif 215<=numkey<=225:                              # allow for off peak measuremnt of protein backbone
                        if traceParms['220nm']:
                            inpUV_220[i].update({key:value})
                    elif traceParms['OtherUV']:                         # other wavelengths. 
                        inpUV_Other[i].update({key:value})
                except Exception:
                    tb = traceback.format_exc()
                    print tb
                    print "ERROR: Something went wront when trying to plot data, %s" , key
                    
        # get rid of emtpy traces. Maybe there is a better way of allocating the list in the lines above? 
        inpUV_280_260 = [trace for trace in inpUV_280_260 if trace]
        inpUV_220 = [trace for trace in inpUV_220 if trace]
        inpUV_Other = [trace for trace in inpUV_Other if trace]
        
        inpCond = [{key:value for key, value in trace.items() if 'Cond' in key and not key.endswith('%') and traceParms[key[-5:]]} for trace in inp]
        # hold values for use in common legend
        hNo220, lNo220, h220, l220, hOther, lOther, hCond, lCond, hOther, lOther = [[] for i in range(10)]                   
        
        plot_x_min, plot_x_max = None, None  # if there are non 220 traces then it will be set then. Otherwise, it will be set by 220, and so on. 
        
        if dispParms['show_xlabel']:
            axes.set_xlabel("Elution volume (ml)")
        if dispParms['show_ylabel']:
            axes.set_ylabel("Absorbance (mAu)")
        # show ticks
        axes.tick_params(axis='x', which='both', top='off')
        axes.tick_params(axis='y', which='both', right='off')
        if not dispParms['xticks']:
            axes.tick_params(axis='x', which='both', bottom='off')
        if not dispParms['xticklabel']:
            axes.tick_params(axis='x', which='both', labelbottom='off')
        if not dispParms['yticks']:
            axes.tick_params(axis='y', which='both', length='off')
        if not dispParms['yticklabel']:
            axes.tick_params(axis='y', which='both', labelleft='off')
        if traceParms['280nm'] or traceParms['260nm']:        #skip this plot if no 280 or 260 is not requested
            plot_x_min, plot_x_max, plot_y_min, plot_y_max = zip(*[self.smartscale(trace, traceParms) for trace in inpUV_280_260]) # get scale without 220 trace
            plot_x_min_idx, plot_x_max_idx = plot_x_min.index(max(plot_x_min)), plot_x_max.index(min(plot_x_max))
            plot_x_min, plot_x_max, plot_y_min, plot_y_max =  plot_x_min[plot_x_min_idx], plot_x_max[plot_x_max_idx], \
                min(plot_y_min), max(plot_y_max)
            axes.set_xlim(plot_x_min, plot_x_max)
            axes.set_ylim(self.expander(plot_y_min, plot_y_max, 0.085))
            counter280, counter260 = 0, 0
            for trace in inpUV_280_260:
                for key, value in trace.items():
                    x_dat, y_dat = self.xy_data(value['data'])
                    if dispParms['linscale']:
                        min_lst = [abs(a - plot_x_min) for a in x_dat]
                        min_idx = min_lst.index(min(min_lst))
                        max_lst = [abs(a - plot_x_max) for a in x_dat]
                        max_idx = max_lst.index(min(max_lst))
                        x_dat = x_dat[min_idx:max_idx]                      # cut the data to region of interest
                        y_dat = y_dat[min_idx:max_idx]
                        y_dat_min = min(y_dat)
                        y_dat_max = max(y_dat)
                        y_dat = [(i-y_dat_min)*(plot_y_max-plot_y_min)/(y_dat_max-y_dat_min)+plot_y_min for i in y_dat]  # linear scaling of data
                    print("Plotting: " + value['data_name'])
                    
                    if '280' in key:
                        if counter280 == 0:
                            stl = self.styles[key[:4]]
                            counter280 += 1
                        elif counter280 < 5:
                            stl = self.styles[key[:4] + str(counter280)]
                            counter280 += 1
                        else:
                            stl = self.styles[key[:4]]   
                            counter280 += 1
   
                    elif '260' in key:
                        if counter260 == 0:
                            stl = self.styles[key[:4]]
                            counter260 += 1                            
                        elif counter260 < 5:
                            stl = self.styles[key[:4] + str(counter280)]
                            counter260 += 1                            
                        else:
                            stl = self.styles[key[:4]]  
                            counter260 += 1     
                                                                
                    #stl = self.styles[key[:4]]
                    p0, = axes.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'],
                                    ls=stl['ls'], lw=stl['lw'],alpha=stl['alpha'], linewidth = dispParms['lineThickness'])

            hNo220, lNo220 = axes.get_legend_handles_labels()
    
        if traceParms['220nm']:               # skip plot of 220 trace if not requested
            par1 = axes.twinx()     # make separate axis for displaying 220 data
            #stl = self.styles['UV3_']

            if dispParms['show_ylabel']:
                par1.set_ylabel("Absorbance (mAu)", color=stl['color'])
                
            par1_x_min, par1_x_max, par1_y_min, par1_y_max = zip(*[self.smartscale(trace, traceParms) for trace in inpUV_220]) # get scale without 220 trace
            par1_x_min_idx, par1_x_max_idx = par1_x_min.index(max(par1_x_min)), par1_x_max.index(min(par1_x_max))
            par1_x_min, par1_x_max, par1_y_min, par1_y_max =  par1_x_min[par1_x_min_idx], par1_x_max[par1_x_max_idx], \
                min(par1_y_min), max(par1_y_max)
            
            if plot_x_min == None or plot_x_max == None:    # x limit has not been set yet 
                plot_x_min, plot_x_max = par1_x_min, par1_x_max   
                par1.set_xlim(plot_x_min, plot_x_max)
                par1.set_ylim(self.expander(par1_y_min, par1_y_max, 0.085))
            else:                                           # go with the x limits that have already been set
                par1.set_xlim(plot_x_min, plot_x_max)
                par1.set_ylim(self.expander(par1_y_min, par1_y_max, 0.085))
            counter220 = 0
            for trace in inpUV_220:
                for key, value in trace.items(): 
                    x_dat, y_dat = self.xy_data(value['data'])
                    if dispParms['linscale']:
                        min_lst = [abs(a - plot_x_min) for a in x_dat]
                        min_idx = min_lst.index(min(min_lst))
                        max_lst = [abs(a - plot_x_max) for a in x_dat]
                        max_idx = max_lst.index(min(max_lst))
                        x_dat = x_dat[min_idx:max_idx]                      # cut the data to region of interest
                        y_dat = y_dat[min_idx:max_idx]
                        y_dat_min = min(y_dat)
                        y_dat_max = max(y_dat)
                        y_dat = [(i-y_dat_min)*(par1_y_max-par1_y_min)/(y_dat_max-y_dat_min)+par1_y_min for i in y_dat]  # linear scaling of data                    
                    print("Plotting: " + value['data_name'])
                    if counter220 == 0:
                        stl = self.styles[key[:4]]
                        counter220 += 1
                    elif counter220 < 5:
                        stl = self.styles[key[:4] + str(counter280)]
                        counter220 += 1
                    else:
                        stl = self.styles[key[:4]]   
                        counter220 += 1
                    p1, = par1.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                                    ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'],linewidth = dispParms['lineThickness'])
            # show ticks
            par1.tick_params(axis='x', which='both', top='off')
            par1.tick_params(axis='y', which='both', right='off')
            if not dispParms['xticks']:
                par1.tick_params(axis='x', which='both', bottom='off')
            if not dispParms['xticklabel']:
                par1.tick_params(axis='x', which='both', labelbottom='off')
            if not dispParms['yticks']:
                par1.tick_params(axis='y', which='both', left='off')
            if not dispParms['yticklabel']:
                par1.tick_params(axis='y', which='both', labelleft='off')
                
            h220, l220 = par1.get_legend_handles_labels()

        if traceParms['Cond']:                 # check if Cond is requested and if it exists in the data
            par2 = axes.twinx()     # make separate axis for displaying Cond data
            #stl = self.styles['Cond']
            
            if dispParms['show_ylabel']:
                par2.set_ylabel('Cond (mS/cm)', color=stl['color'])
            par2_x_min, par2_x_max, par2_y_min, par2_y_max = zip(*[self.smartscale(trace, traceParms) for trace in inpCond]) # get scale without 220 trace
            par2_x_min_idx, par2_x_max_idx = par2_x_min.index(max(par2_x_min)), par2_x_max.index(min(par2_x_max))
            par2_x_min, par2_x_max, par2_y_min, par2_y_max =  par2_x_min[par2_x_min_idx], par2_x_max[par2_x_max_idx], \
                min(par2_y_min), max(par2_y_max)
                
            if plot_x_min == None or plot_x_max == None:    # x limit has not been set yet
                plot_x_min, plot_x_max = par2_x_min, par2_x_max           
                par2.set_xlim(plot_x_min, plot_x_max)
                par2.set_ylim(self.expander(par2_y_min, par2_y_max, 0.085))
            else:                                           # go with the x limits that have already been set
                par2.set_xlim(plot_x_min, plot_x_max)
                par2.set_ylim(self.expander(par2_y_min, par2_y_max, 0.085))
            counterCond = 0
            for trace in inpCond:
                for key, value in trace.items():
                    x_dat, y_dat = self.xy_data(value['data'])
                    if dispParms['linscale']:
                        min_lst = [abs(a - plot_x_min) for a in x_dat]
                        min_idx = min_lst.index(min(min_lst))
                        max_lst = [abs(a - plot_x_max) for a in x_dat]
                        max_idx = max_lst.index(min(max_lst))
                        x_dat = x_dat[min_idx:max_idx]                      # cut the data to region of interest
                        y_dat = y_dat[min_idx:max_idx]
                        y_dat_min = min(y_dat)
                        y_dat_max = max(y_dat)
                        y_dat = [(i-y_dat_min)*(par2_y_max-par2_y_min)/(y_dat_max-y_dat_min)+par2_y_min for i in y_dat]  # linear scaling of data        
                    print("Plotting: " + value['data_name'])
                    if counterCond == 0:
                        stl = self.styles[key[:4]]
                        counterCond += 1
                    elif counterCond < 5:
                        stl = self.styles[key[:4] + str(counter280)]
                        counterCond += 1
                    else:
                        stl = self.styles[key[:4]]   
                        counterCond += 1
                    p2, = par2.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                                    ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'], linewidth = dispParms['lineThickness'])
            hCond, lCond = par2.get_legend_handles_labels()
            
        if traceParms['OtherUV']:
            # make decision about where to plot this other data
            axisOptions = ['280_260', '220', 'Cond', 'New Axis']       # hard code some output axis
            otherAxes = None                                           # get the user prefered axis
            if dispParms['plototherWith'] == axisOptions[0]:
                otherAxes = axes                                       # 280_260 axes always exists by construction 
            elif dispParms['plototherWith'] == axisOptions[1]:
                try:
                    otherAxes = par1                                   # assign to 220 axes if it exists
                except NameError:
                    otherAxes = axes.twinx()                           # otherwise, make new axes
            elif dispParms['plototherWith'] == axisOptions[2]:
                try:
                    otherAxes = par2                                   # assign to Cond axes if it exists
                except NameError:
                    otherAxes = axes.twinx()                           # otherwise, make new axes
            else:
                otherAxes = axes.twinx()
           
            #stl = self.styles['UV_Other']
            
            if dispParms['show_ylabel']:
                otherAxes.set_ylabel('Absorbance (mAu)', color=stl['color'])
                
            other_x_min, other_x_max, other_y_min, other_y_max = zip(*[self.smartscale(trace, traceParms) for trace in inpUV_Other]) # get scale without 220 trace
            other_x_min_idx, other_x_max_idx = other_x_min.index(max(other_x_min)), other_x_max.index(min(other_x_max))
            other_x_min, other_x_max, other_y_min, other_y_max =  other_x_min[other_x_min_idx], other_x_max[other_x_max_idx], \
                min(other_y_min), max(other_y_max)
                
            if plot_x_min == None or plot_x_max == None:            # x limit has not been set yet
                plot_x_min, plot_x_max = other_x_min, other_x_max
                otherAxes.set_xlim(plot_x_min, plot_x_max)
                otherAxes.set_ylim(self.expander(other_y_min, other_y_max, 0.085))
            else:                                                   # go with the x limits that have already been set
                otherAxes.set_xlim(plot_x_min, plot_x_max)
                otherAxes.set_ylim(self.expander(other_y_min, other_y_max, 0.085))
            counterOther = 0
            for trace in inpUV_Other:
                for key, value in trace.items():
                    x_dat, y_dat = self.xy_data(value['data'])
                    if dispParms['linscale']:
                        min_lst = [abs(a - plot_x_min) for a in x_dat]
                        min_idx = min_lst.index(min(min_lst))
                        max_lst = [abs(a - plot_x_max) for a in x_dat]
                        max_idx = max_lst.index(min(max_lst))
                        x_dat = x_dat[min_idx:max_idx]                      # cut the data to region of interest
                        y_dat = y_dat[min_idx:max_idx]
                        y_dat_min = min(y_dat)
                        y_dat_max = max(y_dat)
                        y_dat = [(i-y_dat_min)*(other_y_max-other_y_min)/(y_dat_max-y_dat_min)+other_y_min for i in y_dat]  # linear scaling of data   
                    print("Plotting: " + value['data_name'])
                    if counterOther == 0:
                        stl = self.styles[key[:4]]
                        counterOther += 1
                    elif counterOther < 5:
                        stl = self.styles[key[:4] + str(counter280)]
                        counterOther += 1
                    else:
                        stl = self.styles[key[:4]]   
                        counterOther += 1
                    p1, = otherAxes.plot(x_dat, y_dat, label=value['data_name'], color=stl['color'], 
                                ls=stl['ls'], lw=stl['lw'], alpha=stl['alpha'],linewidth = dispParms['lineThickness'])
            hOther, lOther = otherAxes.get_legend_handles_labels()
            hOther = [i for i in hOther if i not in hNo220]         # avoid redundant labels if plotting on with 280/260 data
            lOther = [i for i in lOther if i not in lNo220]
    
        # plot fraction data
        if dispParms['show_fractions']:             # get fraction data from the first trace
            frac_data = inp[0]['Fractions']['data']
            frac_x, frac_y = self.xy_data(frac_data)
            frac_delta = [abs(a - b) for a, b in zip(frac_x, frac_x[1:])]
            frac_delta.append(frac_delta[-1])
            frac_y_pos = self.mapper(axes.get_ylim()[0], axes.get_ylim()[1], 0.015)
            
            if dispParms['eof']:
                frac_data = frac_data[::2]
            
            for ind, i in enumerate(frac_data):
                axes.axvline(x=i[0], ymin=0.065, ymax=0.0, color='r', linewidth=0.85)
                axes.annotate(str(i[1]), xy=(i[0] + frac_delta[ind] * 0.55, frac_y_pos),
                         horizontalalignment='center', verticalalignment='bottom', size=8, rotation=90) 
            
        # show legend
        if dispParms['show_legend']:
            axes.legend(hNo220+h220+hCond+hOther, lNo220+l220+lCond+lOther, loc='upper right', fontsize=8, fancybox=True )
            axes.xaxis.set_minor_locator(AutoMinorLocator())
            axes.yaxis.set_minor_locator(AutoMinorLocator())

########################################################################
class FilePanel(wx.Panel):
    """
    File navigator and commands to display and output trace files. Instantiated by MainPanel. 
    
    Class Members:
        self.list_ctrl(wx.ListCtrl):
            ListCtrl holding all .res files located in working directory(self.path). 
        self.openBtn(wx.Button):
            Choose a new working directory. 
        self.outputBtn(wx.Button):
            Calls OnOutput to display and output selected files sequentially. 
        self.overlayBtn(wx.Button):
             Calls OnOverlay to overlaty and output selected files simultaneously.
        self.path:
            String representation of the working directory. 
    Class Methods: 
        OnOpenDirectory
            Select new working directory. 
        UpdateDisplay:
            Update the display when a new working directory is chosen. 
        GetSelectedItems:
            Get selected items from self.list_ctrl. 
        OnDoubleClick:
            Display single file. 
        OnOutput:
            Display and output multiple files. 
        OnOverlay
            Overlay two or more files. 
        GetDispParms:
            Helper function to get display parameters. Called by OnDoubleClick, OnOutput, and OnOverlay. 
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Panel.__init__(self, parent)
        self.list_ctrl = wx.ListCtrl(self, style=wx.LB_MULTIPLE | wx.LC_REPORT| wx.BORDER_SUNKEN| wx.LC_SORT_DESCENDING)
        self.list_ctrl.InsertColumn(0, 'Filename', format=wx.LIST_FORMAT_LEFT, width=wx.LIST_AUTOSIZE)
        
        self.openBtn = wx.Button(self, label="Open Folder")             # brows and open files containing .res files. 
        self.openBtn.Bind(wx.EVT_BUTTON, self.OnOpenDirectory)
        self.Bind(wx.EVT_LIST_ITEM_ACTIVATED, self.OnDoubleClick, self.list_ctrl)   # double click on .res file to display trace. 
        
        self.outputBtn = wx.Button(self, label="Output selected \n to file.")
        self.Bind(wx.EVT_BUTTON, self.OnOutput, self.outputBtn)
        
        self.overlayBtn = wx.Button(self, label="Overlay selected \n to file")
        self.Bind(wx.EVT_BUTTON, self.OnOverlay, self.overlayBtn)
        
        vsizer = wx.BoxSizer(wx.VERTICAL)
        vsizer.Add(self.list_ctrl, 1, wx.EXPAND)
        vsizer.Add(self.openBtn, 0, wx.ALL|wx.CENTER, 5)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        hsizer.Add(self.outputBtn, 0, wx.ALL|wx.CENTER)
        hsizer.Add(self.overlayBtn, 0, wx.ALL|wx.CENTER)
        vsizer.Add(hsizer)
        self.SetSizer(vsizer)
        
        self.path = os.getcwd()           # try to load any files in directory program was launched in
        self.UpdateDisplay(self.path)             

    #----------------------------------------------------------------------
    def OnOpenDirectory(self, event):
        """
        Launch pop window to choose a working directory.
        """
        dlg = wx.DirDialog(self, "Choose a directory:")
        if dlg.ShowModal() == wx.ID_OK:
            self.path = dlg.GetPath()
            self.UpdateDisplay(self.path)
        dlg.Destroy()
        
    #----------------------------------------------------------------------
    def UpdateDisplay(self, folder_path):
        """
        Update the listctrl with the .res file names in the working directory. Called by OnOpenDirectory. 
        """
        self.list_ctrl.DeleteAllItems()             # clear the list first
        paths = glob.glob(folder_path + "/*.*")
        counter = 0
        for path in paths:                          # populate list with only .res files in selected directory
            if os.path.basename(path).endswith('.res'):
                self.list_ctrl.InsertStringItem(counter, os.path.basename(path))
                counter += 1
        self.list_ctrl.SetColumnWidth(0, wx.LIST_AUTOSIZE)
            
    #----------------------------------------------------------------------
    def GetSelectedItems(self):
        """    
        Gets the selected items for the list control. Called by OnOutput and OnOverlay. 
        The relative pathnames are returned as a list of strings. 
        """
        selection = []
        index = self.list_ctrl.GetFirstSelected()
        
        if index == -1:        # if nothing is selected return None
            selection= None
        else: 
            selection.append(self.list_ctrl.GetItemText(index))
            while len(selection) != self.list_ctrl.GetSelectedItemCount():
                index = self.list_ctrl.GetNextSelected(index)
                selection.append(self.list_ctrl.GetItemText(index))

        return selection
    
    #----------------------------------------------------------------------
    def OnDoubleClick(self, event):
        """
        Display a single file when user double clicks on any filename in self.list_ctrl. 
        """
        traceParms = self.GetParent().CommandPanel.GetTraceParms()
        dispParms = self.GetDispParms()
        saveParms = {'Save':False}
        
        fname = self.path + '/' + event.GetText()           # send absolute pathname
        self.GetParent().PlotPanel.Display(fname, traceParms, dispParms, saveParms)
        
    #----------------------------------------------------------------------
    def OnOutput(self, event):
        """
        Display and output multiple files sequentially when user clicks on self.outputBtn. 
        
        Raises:
            AssertionError:
                Raised when trying to display and output without selecting files.
            Exception:
                Catches all other possible errors and prints out stack trace. 
                This is for debugging the display panel and will be removed once the display panel is stable. 
        """
        try:
            fnames = self.GetSelectedItems()
            assert not (None in fnames)
            fnames = [self.path + '/' + name for name in fnames]                        # send absolute pathnames

            traceParms = self.GetParent().CommandPanel.GetTraceParms()
            saveParms = {'Save':True, 'format':self.GetParent().Params['outFormat']}    # output to file  
            dispParms = self.GetDispParms()
            for fname in fnames:
                self.GetParent().PlotPanel.Display(fname, traceParms, dispParms, saveParms)
        except AssertionError:
            print "Please load and select a .res file to process."
        except Exception:
            tb = traceback.format_exc()
            print tb

    def OnOverlay(self, event):
        """
        Display and output overlay of multiple files.

        Raises:
            AssertionError:
                Raised when trying to overlay less than 2 files.
            Exception:
                Catches all other possible errors and prints out stack trace. 
                This is for debugging the display panel and will be removed once the display panel is stable.
        """ 
        try:
            fnames = self.GetSelectedItems()
            assert len(fnames)>1 
            fnames = [self.path + '/' + name for name in fnames]        # send absolute pathnames
            traceParms = self.GetParent().CommandPanel.GetTraceParms()
            dispParms = self.GetDispParms()
            saveParms = {'Save':True, 'format':self.GetParent().Params['outFormat']}                 # do output to file  
            self.GetParent().PlotPanel.Display(fnames, traceParms, dispParms, saveParms)
        except AssertionError:
            print "Please select more than one .res file to overlay."
        except Exception:
            tb = traceback.format_exc()
            print tb

    def GetDispParms(self):
        """
        Helper function to gather display parameters from MainPanel.Params. Called by OnDoublClick, OnOutput, and OnOverlay. 
        """
        dispParms = {key:self.GetParent().Params[key] for key in ['lineThickness', 'show_xlabel', 
                                                                 'show_ylabel', 'show_title', 'show_fractions', 
                                                                 'plototherWith', 'show_legend','eof', 'xticks', 
                                                                 'xticklabel', 'yticks', 'yticklabel', 'linscale', 'title']}
        return dispParms


########################################################################
class CommandPanel(wx.Panel):
    """
    Options and commands to process trace files. Instantiated by MainPanel. 
    
    Class Members:
        self.button_280(wx.ToggleButton):
            Output 280nm trace? 
        self.button_260(wx.ToggleButton)
            Output 260nm trace? 
        self.button_220(wx.ToggleButton)
            Output 220nm trace? 
        self.button_OtherUV(wx.ToggleButton)
            Output other wavelength trace? 
        self.button_Cond(wx.ToggleButton)
            Output conductance trace? 
        self.textctrl_Range(wx.TextCtrl)
            Input range in ml to display. 
        self.ToggleList:
            List hodling the UV/Cond output buttons. 
        self.traceParms:
            Dictionary of parameters that will be passed to FilePanel. 
            
    Class Methods:
        GetTraceParms:
            Update selected parameters and output them. Called by FilePanel.
            
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Panel.__init__(self, parent)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        
        self.button_280 = wx.ToggleButton(self,label = "280nm")
        self.button_260 = wx.ToggleButton(self,label = "260nm")
        self.button_220 = wx.ToggleButton(self, label = "220nm")
        self.button_OtherUV = wx.ToggleButton(self, label = "OtherUV")
        self.button_Cond = wx.ToggleButton(self, label = "Cond")
        self.textctrl_Range = wx.TextCtrl(self, name = "Range" )
        
        self.ToggleList =[(self.button_280,), (self.button_260,), (self.button_220,), (self.button_OtherUV,), (self.button_Cond,)]
        sizer.AddMany(self.ToggleList)
        centeredLabel = wx.StaticText(self, -1, 'Range in ml:min-max')
        sizer.Add(centeredLabel)
        sizer.Add(self.textctrl_Range)
        
        self.traceParms = {} # hold all the parementer from the various pressed buttons. 

        self.SetSizer(sizer)
      
    #----------------------------------------------------------------------        
    def Update_TraceParms(self):
        """
        Update selected parameters stored in self.Parms. 
        
        Raises:
            KeyError:
                Raised when incorrect format for range is entered into self.textctrl_Range.
        """
        self.traceParms = {}                     # get the parameters fresh each time 
        for button in self.ToggleList:
            self.traceParms.update({button[0].GetLabel(): button[0].GetValue()})    # get the booleans values describing whether or not button is pushed
        range = None
        if len(self.textctrl_Range.GetValue().split('-')) == 2:                     # The next few lines allows for a few ways to define the range
            range = self.textctrl_Range.GetValue().split('-')
        elif len(self.textctrl_Range.GetValue().split('_')) == 2:
            range = self.textctrl_Range.GetValue().split('_')
        elif len(self.textctrl_Range.GetValue().split(' ')) == 2:
            range = self.textctrl_Range.GetValue().split(' ')      
        try:
            if range and float(range[0]) <= float(range[1]):       # if we have a range designate xmin and xmax only if xmin is smaller than xmax7
                self.traceParms.update({u'xmin':float(range[0]), u'xmax':float(range[1])})
            else:           # no range, no xmin of xmax
                self.traceParms.update({u'xmin':None, u'xmax':None})
        except KeyError:
            print("Warning: Enter correct range. (xmin-xmax)!")
            self.traceParms.update({u'xmin':None, u'xmax':None})

    #----------------------------------------------------------------------
    def GetTraceParms(self):
        """
        Update selected parameters and output them. 

        Returns:
            self.traceParms
        """
        self.Update_TraceParms()
        return self.traceParms

########################################################################
class MenuBar():            # The parent to the menu bar is MainFrame
    """
    Menu bar for the App attached to MainFrame. 
    
    Class Members:
        self.parent:
           Explicitly saved parent since there isn't a parent field in the MenuBar constructor. 
        self.menuBar(wx.MenuBar)
        self.optionsMenu(wx.Menu)              
        self.optionsMenu_SaveOptions(wx.MenuItem) 
        self.optionsMenu_DisplayOptions(wx.MenuItem) 
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Initialze the class members. 
        """
        self.parent = parent
        self.menuBar = wx.MenuBar()                 # the menu bar
       
        self.optionsMenu = wx.Menu()                # options menu tab
        self.optionsMenu_SaveOptions = self.optionsMenu.Append(wx.ID_ANY, 'Save Options')
        self.optionsMenu_DisplayOptions = self.optionsMenu.Append(wx.ID_ANY, 'Display Options')
        self.optionsMenu.Bind(wx.EVT_MENU, self.OnSaveOptions, self.optionsMenu_SaveOptions)
        self.optionsMenu.Bind(wx.EVT_MENU, self.OnDisplayOptions, self.optionsMenu_DisplayOptions)
        
        self.menuBar.Append(self.optionsMenu, '&Options')       # append optionsMenu to menu bar 
        self.parent.SetMenuBar(self.menuBar)        # attach menubar to MainPanel
        
    #----------------------------------------------------------------------
    def OnSaveOptions(self, event):
        """
        Launch pop up to change save parameters.
        """
        popUpSaveOptions = SavePopUp(self.parent)       
        popUpSaveOptions.Show()
        
    #----------------------------------------------------------------------
    def OnDisplayOptions(self, event):    
        """
        Launch pop up to change display parameters.
        """              
        popUpDisplay = DisplayPopUp(self.parent)
        popUpDisplay.Show()

########################################################################
class SavePopUp(wx.Frame):              # The parent of this pop-up is the MainFrame
    """
    A pop up window to view and change save parameters. 
    Should be called by MenuBar so that parameters held in MainPanel can be viewed and changed. 
    
    Class Members:
        self.saveParams:
            Dictionary of save parameters. This is initialized from items in MainPanel.Params.
        self.formatOptions:
            Dictionary of strings describing available output file formats. 
        self.formatCombo
    Class Methods:
        OnOpenDirectory:
            Pop up to select a directory to output files. 
        UpdateOutDir:
            Called by OnOpenDirectory upon directory selection.
        OnFormatCombo:
            Updates output file fomat internal to class and MainPanel.Params by calling SendSaveMenu.
        SendSaveMenu:
            Change save parameters in MainPanel.Params by calling MainPanel.ChangeParams. 

    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Frame.__init__(self, parent, title ='Save Options', size = (600,100),
                           style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER| wx.RESIZE_BOX| wx.MAXIMIZE_BOX))
        panel = wx.Panel(self)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer3 = wx.BoxSizer(wx.HORIZONTAL)
        # internal variables to output to the main panel. 
        self.saveParams = {'outDir':parent.MainPanel.Params['outDir'],
                           'outFormat':parent.MainPanel.Params['outFormat']}
        self.formatOptions = ['pdf','png','jpeg']
        self.formatCombo = wx.ComboBox(panel, -1, self.saveParams['outFormat'], (15,40), 
                                        wx.DefaultSize, self.formatOptions, wx.CB_DROPDOWN)
        sizer2.AddSpacer(10)
        sizer2.Add(wx.StaticText(panel,-1, 'File Format:', (15,15)))
        sizer2.Add(self.formatCombo)
        self.buttonOutDir = wx.Button(panel, label = 'Output Directory:')
        self.outDirTxt = wx.StaticText(panel)
        self.outDirTxt.SetLabel(self.saveParams['outDir'])
        panel.Bind(wx.EVT_BUTTON, self.OnOpenDirectory, self.buttonOutDir)
        panel.Bind(wx.EVT_COMBOBOX, self.OnFormatCombo, self.formatCombo)
           
        sizer3.AddSpacer(10)     
        sizer3.Add(self.buttonOutDir)
        sizer3.Add(self.outDirTxt)
        sizer1.AddSpacer(10)
        sizer1.Add(sizer2,wx.TOP|wx.LEFT, 2)
        sizer1.AddSpacer(10)
        sizer1.Add(sizer3,wx.TOP|wx.LEFT,2)

        panel.SetSizer(sizer1)
        
    #----------------------------------------------------------------------
    def OnOpenDirectory(self, event):
        """
        Pop up to select a directory to output files. 
        """
        dlg = wx.DirDialog(self, "Choose a directory:")
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.UpdateOutDir(path)
        dlg.Destroy()

    #----------------------------------------------------------------------      
    def OnFormatCombo(self, event):
        """
        Update internal output file fomat. 
        Update output file format in MainPanel.Params by calling SendSaveMenu.
        """
        self.saveParams['outFormat'] = self.formatCombo.GetValue()
        self.SendSaveMenu()
        
    #----------------------------------------------------------------------
    def UpdateOutDir(self, path):
        """
        Update internal output directory. 
        Update output directory in MainPanel.Params by calling SendSaveMenu.
        """
        self.saveParams['outDir'] = path
        self.outDirTxt.SetLabel(self.saveParams['outDir'])
        self.SendSaveMenu()


    #----------------------------------------------------------------------      
    def SendSaveMenu(self):
        """
        Change save parameters in MainPanel.Params by calling MainPanel.ChangeParams. 
        """
        self.GetParent().MainPanel.ChangeParams(self.saveParams)
        
########################################################################
class DisplayPopUp(wx.Frame):
    """
    A pop up window to view and change display parameters. 
    Should be called by MenuBar so that parameters held in MainPanel can be viewed and changed. 
    
    Class Members:
        self.DisplayParams:
            Dictionary of display parameters. This is initialized from items in MainPanel.Params.
        self.ThicknessOptions:
            List of string ints describing the thickness of the traces. The string will be converted 
            to ints when calling self.Send_DisplayMenu.
        self.Thickness_Combo 
        self.Show_Xlabel_CkBx 
        self.Show_Ylabel_CkBx 
        self.Show_Fractions_CkBx
        self.Show_Legend_CkBx 
        self.PlotWithOptions:
            List of strings describing which axes to plot other UV traces.
        self.PlotWithOptions_Combo
    Class Methods:
        On_Thickness_Combo
        On_Show_Xlabel_CkBx
        On_Show_Ylabel_CkBx
        On_Show_Fractions_CkBx 
        On_Show_Legend_CkBx
        On_PlotWithOptions_Combo
        Send_DisplayMenu:
            Change display parameters in MainPanel.Params by calling MainPanel.ChangeParams. 
            
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Frame.__init__(self, parent, title = 'Display Options', size = (350,200),
                           style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER| wx.RESIZE_BOX| wx.MAXIMIZE_BOX))
        panel = wx.Panel(self)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        sizer2 = wx.BoxSizer(wx.VERTICAL)
        sizer3 = wx.BoxSizer(wx.HORIZONTAL)
        self.DisplayParams = {key:parent.MainPanel.Params[key] for key in ['lineThickness','show_xlabel', # get current display parameters
                                                        'show_ylabel','show_title', 'show_fractions','plototherWith','show_legend','eof', 
                                                        'xticks', 'xticklabel', 'yticks', 'yticklabel','linscale','title']}
        self.ThicknessOptions = [str(i) for i in range(1,10)]
        self.Thickness_Combo = wx.ComboBox(panel, -1, str(self.DisplayParams['lineThickness']), (15,40), 
                                           wx.DefaultSize, self.ThicknessOptions, wx.CB_DROPDOWN)
        self.Show_Xlabel_CkBx = wx.CheckBox(panel, -1, 'show_xlabel', (15,40), wx.DefaultSize)
        self.Show_Xlabel_CkBx.SetValue(self.DisplayParams['show_xlabel'])
        self.Show_Ylabel_CkBx = wx.CheckBox(panel, -1, 'show_ylabel', (15,40), wx.DefaultSize)
        self.Show_Ylabel_CkBx.SetValue(self.DisplayParams['show_ylabel'])
        self.Show_Xticks_Ckbx = wx.CheckBox(panel, -1, 'show_xticks', (15,40), wx.DefaultSize)
        self.Show_Xticks_Ckbx.SetValue(self.DisplayParams['xticks'])
        self.Show_XtickLabel_Ckbx = wx.CheckBox(panel, -1, 'show_xtick_label', (15,40), wx.DefaultSize)
        self.Show_XtickLabel_Ckbx.SetValue(self.DisplayParams['xticklabel'])
        self.Show_Yticks_Ckbx = wx.CheckBox(panel, -1, 'show_yticks', (15,40), wx.DefaultSize)
        self.Show_Yticks_Ckbx.SetValue(self.DisplayParams['yticks'])
        self.Show_YtickLabel_Ckbx = wx.CheckBox(panel, -1, 'show_ytick_label', (15,40), wx.DefaultSize)
        self.Show_YtickLabel_Ckbx.SetValue(self.DisplayParams['yticklabel'])    
        self.Show_Fractions_CkBx = wx.CheckBox(panel, -1, 'show_fractions', (15,40), wx.DefaultSize)
        self.Show_Fractions_CkBx.SetValue(self.DisplayParams['show_fractions'])
        self.EOF_CkBx = wx.CheckBox(panel, -1, 'every other fraction', (15,40), wx.DefaultSize)
        self.EOF_CkBx.SetValue(self.DisplayParams['eof'])
        self.Show_Legend_CkBx = wx.CheckBox(panel, -1, 'show_legend', (15,40), wx.DefaultSize)
        self.Show_Legend_CkBx.SetValue(self.DisplayParams['show_legend'])
        self.Show_Title_CkBx = wx.CheckBox(panel, -1, 'show_title', (15,40), wx.DefaultSize)
        self.Show_Title_CkBx.SetValue(self.DisplayParams['title'])
        self.LinarScaling_CkBx = wx.CheckBox(panel, -1, 'linear_scaling', (15,40), wx.DefaultSize)
        self.LinarScaling_CkBx.SetValue(self.DisplayParams['linscale'])
        self.PlotWithOptions = ['280_260', '220', 'Cond', 'New Axis']
        self.PlotWithOptions_Combo = wx.ComboBox(panel, -1, self.DisplayParams['plototherWith'], (15,40), 
                                           wx.DefaultSize, self.PlotWithOptions, wx.CB_DROPDOWN)
        self.Color_Btn = wx.Button(panel, -1, 'select_colors', (15,40), wx.DefaultSize)


        panel.Bind(wx.EVT_COMBOBOX, self.On_Thickness_Combo, self.Thickness_Combo)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Xlabel_CkBx, self.Show_Xlabel_CkBx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Ylabel_CkBx, self.Show_Ylabel_CkBx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Xticks_Ckbx, self.Show_Xticks_Ckbx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_XtickLabel_Ckbx, self.Show_XtickLabel_Ckbx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Yticks_Ckbx, self.Show_Yticks_Ckbx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_YtickLabel_Ckbx, self.Show_YtickLabel_Ckbx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Fractions_CkBx, self.Show_Fractions_CkBx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Title_CkBx, self.Show_Title_CkBx)

        panel.Bind(wx.EVT_CHECKBOX, self.On_EOF_CkBx, self.EOF_CkBx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_Show_Legend_CkBx, self.Show_Legend_CkBx)
        panel.Bind(wx.EVT_CHECKBOX, self.On_LinarScaling_CkBx, self.LinarScaling_CkBx)
        panel.Bind(wx.EVT_COMBOBOX, self.On_PlotWithOptions_Combo, self.PlotWithOptions_Combo)
        panel.Bind(wx.EVT_BUTTON, self.On_Color_Btn, self.Color_Btn)

        sizer1.Add(wx.StaticText(panel,-1, 'Line Thickness:', (15,15)))
        sizer1.Add(self.Thickness_Combo)
        sizer1.Add(self.Show_Xlabel_CkBx)
        sizer1.Add(self.Show_Ylabel_CkBx)
        sizer1.Add(self.Show_Xticks_Ckbx)
        sizer1.Add(self.Show_XtickLabel_Ckbx)
        sizer1.Add(self.Show_Yticks_Ckbx)
        sizer1.Add(self.Show_YtickLabel_Ckbx)
        sizer2.Add(self.Show_Fractions_CkBx)
        sizer2.Add(self.Show_Title_CkBx)
        sizer2.Add(self.EOF_CkBx)
        sizer2.Add(self.Show_Legend_CkBx)
        sizer2.Add(self.LinarScaling_CkBx)
        sizer2.Add(wx.StaticText(panel,-1, 'Plot Other UV on axis:', (15,15)))
        sizer2.Add(self.PlotWithOptions_Combo)
        sizer2.Add(self.Color_Btn)
        sizer3.AddSpacer(20)
        sizer3.AddMany([(sizer1, wx.ALL, 5),(sizer2, wx.ALL, 5)])
        panel.SetSizer(sizer3)
    
    #----------------------------------------------------------------------
    def On_Thickness_Combo(self, event):
        """
        Set the thickness of the trace lines.
        """
        self.DisplayParams['lineThickness'] = int(self.Thickness_Combo.GetValue())
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------
    def On_Show_Xlabel_CkBx(self, event):
        """
        Toggle to show xlabel
        """
        self.DisplayParams['show_xlabel'] = self.Show_Xlabel_CkBx.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------
    def On_Show_Ylabel_CkBx(self, event):
        """
        Toggle to show ylabel
        """
        self.DisplayParams['show_ylabel'] = self.Show_Ylabel_CkBx.GetValue()
        self.Send_DisplayMenu()
    
    #----------------------------------------------------------------------
    def On_Show_Xticks_Ckbx(self, event):
        """
        Toggle to show xticks
        """
        self.DisplayParams['xticks'] = self.Show_Xticks_Ckbx.GetValue()
        self.Send_DisplayMenu()

    #----------------------------------------------------------------------
    def On_Show_XtickLabel_Ckbx(self, event):
        """
        Toggle to show xtick labels
        """
        self.DisplayParams['xticklabel'] = self.Show_XtickLabel_Ckbx.GetValue()
        self.Send_DisplayMenu()
  
    #----------------------------------------------------------------------
    def On_Show_Yticks_Ckbx(self, event):
        """
        Toggle to show yticks
        """
        self.DisplayParams['yticks'] = self.Show_Yticks_Ckbx.GetValue()
        self.Send_DisplayMenu()

    #----------------------------------------------------------------------
    def On_Show_YtickLabel_Ckbx(self, event):
        """
        Toggle to show ytick labels
        """
        self.DisplayParams['yticklabel'] = self.Show_YtickLabel_Ckbx.GetValue()
        self.Send_DisplayMenu()

    #----------------------------------------------------------------------  
    def On_Show_Fractions_CkBx(self, event):
        """
        Toggle to show Fractions
        """
        self.DisplayParams['show_fractions'] = self.Show_Fractions_CkBx.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------  
    def On_EOF_CkBx(self, event):
        """
        Toggle to reduce fraction tick marks to half
        """
        self.DisplayParams['eof'] = self.EOF_CkBx.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------  
    def On_Show_Legend_CkBx(self, event):
        """
        Toggle to show Legend
        """
        self.DisplayParams['show_legend'] = self.Show_Legend_CkBx.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------  
    def On_Show_Title_CkBx(self, event):
        """
        Show title at top of plot.
        """
        self.DisplayParams['title'] = self.Show_Title_CkBx.GetValue()
        self.Send_DisplayMenu()

    #----------------------------------------------------------------------  
    def On_PlotWithOptions_Combo(self, event):
        """
        Set which axis to plot other UV tracs.
        """
        self.DisplayParams['plototherWith'] = self.PlotWithOptions_Combo.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------  
    def On_LinarScaling_CkBx(self, event):
        """
        Linearly scale data for overlaying traces.
        """
        self.DisplayParams['linscale'] = self.LinarScaling_CkBx.GetValue()
        self.Send_DisplayMenu()
        
    #----------------------------------------------------------------------
    def Send_DisplayMenu(self):     # calls MainPanel.ChangeParams through MainFrame
        """
        Call MainPanel.ChangeParams through MainFrame.
        """
        self.GetParent().MainPanel.ChangeParams(self.DisplayParams)
        
    #----------------------------------------------------------------------
    def On_Color_Btn(self, event):
        """
        Launch ColorPopUp
        """
        colors = ColorPopUp(self)
        colors.Show()

########################################################################
class ColorPopUp(wx.Frame):
    """
    A pop up window to view and change color parameters. 
    Should be called by DisplayPopUp. 
    
    Class Members:
        
    Class Methods:
    """
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Frame.__init__(self, parent, title = 'Select Colors', size = (500,175),
                           style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER| wx.RESIZE_BOX| wx.MAXIMIZE_BOX))
        panel = wx.Panel(self)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        sizer2 = wx.BoxSizer(wx.VERTICAL)
        sizer3 = wx.BoxSizer(wx.VERTICAL)
        sizer4 = wx.BoxSizer(wx.VERTICAL)
        sizer5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer6 = wx.BoxSizer(wx.VERTICAL)
        

        self.Colors = {key:tuple([int(255*i) for i in parent.GetParent().MainPanel.PlotPanel.styles[key]['color']]) \
                       for key in parent.GetParent().MainPanel.PlotPanel.styles.keys()}     # get the current color settings        
        self.UV_280_Btn = wx.Button(panel, -1, 'UV_280', (15,40), wx.DefaultSize)
        self.UV_280_Btn_1 = wx.Button(panel, -1, 'UV_280_1', (15,40), wx.DefaultSize)
        self.UV_280_Btn_2 = wx.Button(panel, -1, 'UV_280_2', (15,40), wx.DefaultSize)
        self.UV_280_Btn_3 = wx.Button(panel, -1, 'UV_280_3', (15,40), wx.DefaultSize)
        self.UV_280_Btn_4 = wx.Button(panel, -1, 'UV_280_4', (15,40), wx.DefaultSize)

        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_'])
        self.UV_280_Btn_1.SetBackgroundColour(self.Colors['UV1_1'])
        self.UV_280_Btn_2.SetBackgroundColour(self.Colors['UV1_2'])
        self.UV_280_Btn_3.SetBackgroundColour(self.Colors['UV1_3'])
        self.UV_280_Btn_4.SetBackgroundColour(self.Colors['UV1_4'])

        self.UV_260_Btn = wx.Button(panel, -1, 'UV_260', (15,40), wx.DefaultSize)
        self.UV_260_Btn_1 = wx.Button(panel, -1, 'UV_260', (15,40), wx.DefaultSize)
        self.UV_260_Btn_2 = wx.Button(panel, -1, 'UV_260', (15,40), wx.DefaultSize)
        self.UV_260_Btn_3 = wx.Button(panel, -1, 'UV_260', (15,40), wx.DefaultSize)
        self.UV_260_Btn_4 = wx.Button(panel, -1, 'UV_260', (15,40), wx.DefaultSize)

        self.UV_260_Btn.SetBackgroundColour(self.Colors['UV2_'])
        self.UV_260_Btn_1.SetBackgroundColour(self.Colors['UV2_1'])
        self.UV_260_Btn_2.SetBackgroundColour(self.Colors['UV2_2'])
        self.UV_260_Btn_3.SetBackgroundColour(self.Colors['UV2_3'])
        self.UV_260_Btn_4.SetBackgroundColour(self.Colors['UV2_4'])

        self.UV_220_Btn = wx.Button(panel, -1, 'UV_220', (15,40), wx.DefaultSize)
        self.UV_220_Btn_1 = wx.Button(panel, -1, 'UV_220_1', (15,40), wx.DefaultSize)
        self.UV_220_Btn_2 = wx.Button(panel, -1, 'UV_220_2', (15,40), wx.DefaultSize)
        self.UV_220_Btn_3 = wx.Button(panel, -1, 'UV_220_3', (15,40), wx.DefaultSize)
        self.UV_220_Btn_4 = wx.Button(panel, -1, 'UV_220_4', (15,40), wx.DefaultSize)

        self.UV_220_Btn.SetBackgroundColour(self.Colors['UV3_'])
        self.UV_220_Btn_1.SetBackgroundColour(self.Colors['UV3_1'])
        self.UV_220_Btn_2.SetBackgroundColour(self.Colors['UV3_2'])
        self.UV_220_Btn_3.SetBackgroundColour(self.Colors['UV3_3'])
        self.UV_220_Btn_4.SetBackgroundColour(self.Colors['UV3_4'])

        self.UV_Other_Btn = wx.Button(panel, -1, 'UV_Other', (15,40), wx.DefaultSize)
        self.UV_Other_Btn_1 = wx.Button(panel, -1, 'UV_Other1', (15,40), wx.DefaultSize)
        self.UV_Other_Btn_2 = wx.Button(panel, -1, 'UV_Other2', (15,40), wx.DefaultSize)
        self.UV_Other_Btn_3 = wx.Button(panel, -1, 'UV_Other3', (15,40), wx.DefaultSize)
        self.UV_Other_Btn_4 = wx.Button(panel, -1, 'UV_Other4', (15,40), wx.DefaultSize)

        self.UV_Other_Btn.SetBackgroundColour(self.Colors['UV_Other'])
        self.UV_Other_Btn_1.SetBackgroundColour(self.Colors['UV_Other1'])
        self.UV_Other_Btn_2.SetBackgroundColour(self.Colors['UV_Other2'])
        self.UV_Other_Btn_3.SetBackgroundColour(self.Colors['UV_Other3'])
        self.UV_Other_Btn_4.SetBackgroundColour(self.Colors['UV_Other4'])

        self.Cond_Btn = wx.Button(panel, -1, 'Cond', (15,40), wx.DefaultSize)
        self.Cond_Btn_1 = wx.Button(panel, -1, 'Cond1', (15,40), wx.DefaultSize)
        self.Cond_Btn_2 = wx.Button(panel, -1, 'Cond2', (15,40), wx.DefaultSize)
        self.Cond_Btn_3 = wx.Button(panel, -1, 'Cond3', (15,40), wx.DefaultSize)
        self.Cond_Btn_4 = wx.Button(panel, -1, 'Cond4', (15,40), wx.DefaultSize)

        self.Cond_Btn.SetBackgroundColour(self.Colors['Cond'])
        self.Cond_Btn_1.SetBackgroundColour(self.Colors['Cond1'])
        self.Cond_Btn_2.SetBackgroundColour(self.Colors['Cond2'])
        self.Cond_Btn_3.SetBackgroundColour(self.Colors['Cond3'])
        self.Cond_Btn_4.SetBackgroundColour(self.Colors['Cond4'])

        panel.Bind(wx.EVT_BUTTON, self.On_UV_280_Btn, self.UV_280_Btn)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_280_Btn_1, self.UV_280_Btn_1)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_280_Btn_2, self.UV_280_Btn_2)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_280_Btn_3, self.UV_280_Btn_3)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_280_Btn_4, self.UV_280_Btn_4)

        panel.Bind(wx.EVT_BUTTON, self.On_UV_260_Btn, self.UV_260_Btn)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_260_Btn_1, self.UV_260_Btn_1)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_260_Btn_2, self.UV_260_Btn_2)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_260_Btn_3, self.UV_260_Btn_3)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_260_Btn_4, self.UV_260_Btn_4)

        panel.Bind(wx.EVT_BUTTON, self.On_UV_220_Btn, self.UV_220_Btn)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_220_Btn_1, self.UV_220_Btn_1)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_220_Btn_2, self.UV_220_Btn_2)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_220_Btn_3, self.UV_220_Btn_3)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_220_Btn_4, self.UV_220_Btn_4)

        panel.Bind(wx.EVT_BUTTON, self.On_UV_Other_Btn, self.UV_Other_Btn)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_Other_Btn_1, self.UV_Other_Btn_1)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_Other_Btn_2, self.UV_Other_Btn_2)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_Other_Btn_3, self.UV_Other_Btn_3)
        panel.Bind(wx.EVT_BUTTON, self.On_UV_Other_Btn_4, self.UV_Other_Btn_4)

        panel.Bind(wx.EVT_BUTTON, self.On_Cond_Btn, self.Cond_Btn)
        panel.Bind(wx.EVT_BUTTON, self.On_Cond_Btn_1, self.Cond_Btn_1)
        panel.Bind(wx.EVT_BUTTON, self.On_Cond_Btn_2, self.Cond_Btn_2)
        panel.Bind(wx.EVT_BUTTON, self.On_Cond_Btn_3, self.Cond_Btn_3)
        panel.Bind(wx.EVT_BUTTON, self.On_Cond_Btn_4, self.Cond_Btn_4)

        sizer.AddMany([(wx.StaticText(panel,-1,'1st Trace', (15,15))), (wx.StaticText(panel,-1,' ', (10,10))), (self.UV_280_Btn,),(self.UV_260_Btn,),(self.UV_220_Btn,),(self.UV_Other_Btn,),(self.Cond_Btn,)])
        sizer1.AddMany([(wx.StaticText(panel,-1,'2nd Trace', (15,15))), (wx.StaticText(panel,-1,' ', (10,10))), (self.UV_280_Btn_1,),(self.UV_260_Btn_1,),(self.UV_220_Btn_1,),(self.UV_Other_Btn_1,),(self.Cond_Btn_1,)])
        sizer2.AddMany([(wx.StaticText(panel,-1,'3rd Trace', (15,15))), (wx.StaticText(panel,-1,' ', (10,10))), (self.UV_280_Btn_2,),(self.UV_260_Btn_2,),(self.UV_220_Btn_2,),(self.UV_Other_Btn_2,),(self.Cond_Btn_2,)])
        sizer3.AddMany([(wx.StaticText(panel,-1,'4th Trace', (15,15))), (wx.StaticText(panel,-1,' ', (10,10))), (self.UV_280_Btn_3,),(self.UV_260_Btn_3,),(self.UV_220_Btn_3,),(self.UV_Other_Btn_3,),(self.Cond_Btn_3,)])
        sizer4.AddMany([(wx.StaticText(panel,-1,'5th Trace', (15,15))), (wx.StaticText(panel,-1,' ', (10,10))), (self.UV_280_Btn_4,),(self.UV_260_Btn_4,),(self.UV_220_Btn_4,),(self.UV_Other_Btn_4,),(self.Cond_Btn_4,)])
        sizer5.AddSpacer(20)
        sizer5.AddMany([(sizer,), (wx.StaticText(panel,-1,' ', (10,10))), (sizer1,), (wx.StaticText(panel,-1,' ', 
                        (10,10))), (sizer2,), (wx.StaticText(panel,-1,' ', (10,10))), (sizer3,), (wx.StaticText(panel,-1,' ', (10,10))), (sizer4,)])
        panel.SetSizer(sizer5)
        
    #----------------------------------------------------------------------
    def On_UV_280_Btn(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV1_'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV1_'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_280_Btn_1(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV1_1'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV1_1'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_1'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_1'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_280_Btn_2(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV1_2'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV1_2'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_2'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_2'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_280_Btn_3(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV1_3'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV1_3'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_3'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_3'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_280_Btn_4(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV1_4'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV1_4'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_4'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV1_4'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_260_Btn(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV2_'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV2_'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_260_Btn_1(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV2_1'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV2_1'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_1'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_1'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_260_Btn_2(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV2_2'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV2_2'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_2'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_2'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_260_Btn_3(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV2_3'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV2_3'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_3'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_3'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_UV_260_Btn_4(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV2_4'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV2_4'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_4'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV2_4'])
        self.Send_Colors()


    #----------------------------------------------------------------------
    def On_UV_220_Btn(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV3_'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV3_'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_220_Btn_1(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV3_1'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV3_1'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_1'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_1'])
        self.Send_Colors()
    #----------------------------------------------------------------------
    def On_UV_220_Btn_2(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV3_2'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV3_2'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_2'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_2'])
        self.Send_Colors()
    #----------------------------------------------------------------------
    def On_UV_220_Btn_3(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV3_3'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV3_3'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_3'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_3'])
        self.Send_Colors()
    #----------------------------------------------------------------------
    def On_UV_220_Btn_4(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV3_4'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV3_4'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_4'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV3_4'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_Other_Btn(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV_Other'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV_Other'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_Other_Btn_1(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV_Other1'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV_Other1'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other1'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other1'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_Other_Btn_2(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV_Other2'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV_Other2'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other2'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other2'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_Other_Btn_3(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV_Other3'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV_Other3'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other3'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other3'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def On_UV_Other_Btn_4(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['UV_Other4'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['UV_Other4'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other4'])
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['UV_Other4'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_Cond_Btn(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['Cond'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['Cond'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond']) 
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond'])
        self.Send_Colors() 
               
    #----------------------------------------------------------------------
    def On_Cond_Btn_1(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['Cond1'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['Cond1'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond1']) 
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond1'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_Cond_Btn_2(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['Cond2'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['Cond2'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond2']) 
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond2'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_Cond_Btn_3(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['Cond3'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['Cond3'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond3']) 
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond3'])
        self.Send_Colors()
        
    #----------------------------------------------------------------------
    def On_Cond_Btn_4(self, event):
        data = wx.ColourData()
        data.SetChooseFull(True)
        data.SetColour(self.Colors['Cond4'])
        dlg = wx.ColourDialog(self,data)
        if dlg.ShowModal() == wx.ID_OK:
            newColor = dlg.GetColourData()
            self.Colors['Cond4'] = newColor.GetColour().Get()
            self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond4']) 
        dlg.Destroy()
        self.UV_280_Btn.SetBackgroundColour(self.Colors['Cond4'])
        self.Send_Colors()

    #----------------------------------------------------------------------
    def Send_Colors(self):
        for key in self.Colors.keys():
            try:
                self.GetParent().GetParent().MainPanel.PlotPanel.styles[key]['color'] = tuple([col/255. for col in self.Colors[key]])
            except KeyError:                # This exception shouldn't be raised, but just in case. 
                print "Color: %s was not found in MainPanel.Params" % key
        
########################################################################
class MainPanel(wx.Panel):
    """
    The main display. All other displays are contained in pop ups. Instantiated by MainFrame.
    
    Class Members:
        self.PlotPanel(PlotPanel)
        self.FilePanel(FilePanel)
        self.CommandPanel(CommandPanel)
        self.Params:
            Parameters that are used by child panels and pop ups called from MainFrame. 
    Class Methods:
        ChangeParams:
            Called by pop ups to change paramters by user. 
    """
    def __init__(self, parent):
        """
        Use inherited constructor and initialze the class members. 
        """
        wx.Panel.__init__(self, parent)
        
        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        vsizer = wx.BoxSizer(wx.VERTICAL)
        
        self.PlotPanel = PlotPanel(self)
        self.FilePanel= FilePanel(self)
        self.CommandPanel = CommandPanel(self)
             
        vsizer.Add(self.PlotPanel, 1, wx.EXPAND)
        vsizer.Add(self.CommandPanel, 1, wx.ALIGN_CENTER| wx.TOP, 30)
        hsizer.Add(vsizer, 1, wx.EXPAND)
        hsizer.Add(self.FilePanel, 2, wx.EXPAND)
        self.SetSizer(hsizer)
        
        self.Params = {'outDir':os.getcwd(),
                       'outFormat':'pdf', 'lineThickness':2, 'show_xlabel':True, 'show_ylabel':True,
                       'show_title':True, 'show_fractions':True, 'plototherWith':'280_260', 'show_legend':True,
                        'eof':False, 'xticks':True, 'xticklabel':True, 'yticks':True, 'yticklabel':True, 'linscale':False, 'title':True}  
        
        self.Show(True)
    
    def ChangeParams(self, params):
        """
        Args:
            params: 
                A dictionary of parameters. The keys should match self.Params.
        Raises:
            KeyError:    
                If params contains a key that is not in self.Params. 
        """
        for key in params.keys():
            try:
                self.Params[key] = params[key]
            except KeyError:                # This exception shouldn't be raised, but just in case. 
                print "Parameter: %s was not found in MainPanel.Params" % key

########################################################################
class MainFrame(wx.Frame):
    """
    The top level frame/'parent' that contains the rest of the app. 
    
    Class Members:
        self.Menubar(Menubar)
        self.MainPanel(MainPanel)
        
    """

    #----------------------------------------------------------------------
    def __init__(self):
        """
        Use inherited constructor and initialize the main panel and the menu bar. 
        """
        wx.Frame.__init__(self, None, wx.ID_ANY, title="PyCornGUI", size=(900,600), 
                          style=wx.DEFAULT_FRAME_STYLE & ~(wx.RESIZE_BORDER| wx.RESIZE_BOX| wx.MAXIMIZE_BOX)) 
        
        self.Menubar = MenuBar(self)            

        self.MainPanel = MainPanel(self)
        
        self.Show()
        
def run():
    """
    Run GUI    
    """
    app = wx.App(False)
    frame = MainFrame()
    app.MainLoop()  

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    app.MainLoop()
