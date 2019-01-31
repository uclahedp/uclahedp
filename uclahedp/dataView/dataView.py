import sys
import time
import os
import h5py

import numpy as np

from matplotlib.backends.qt_compat import QtCore, QtWidgets, QtGui, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
    
import matplotlib.figure
import matplotlib.cm





class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        


        #self.setStyleSheet(open('stylesheet.css').read())
        
        
        self.dir = os.path.dirname(os.path.realpath(__file__))

        iconfile = os.path.join(self.dir,'program_icon.png' )
        self.setWindowIcon(QtGui.QIcon(iconfile))
        
        self.setWindowTitle("HEDP dataView")
        
        self.filepath = ''
    
        self.axes = []
        self.plottype = 1 #1 or 2 for 1D or 2D plot
        self.cur_axes = [0,0]
        
        self.last_cur_axes = [0,0]
        self.last_axes = []#Used to smoothly transfer axes info when new files loaded  
        
        self.datarange = [None, None]
        
        self.colormap_dict = {'Autumn':'autumn', "Winter":"winter", "Cool":"cool", 
                          "Ocean":"ocean", "Rainbow":"gist_rainbow",
                          "Seismic":"seismic", "RedGrey":"RdGy",
                          "Coolwarm":"coolwarm"}
        

        
        self.layout = QtWidgets.QHBoxLayout(self._main)
        
        #
        #Define Actions
        #
        quitAct = QtWidgets.QAction(" &Quit", self)
        quitAct.triggered.connect(self.close)
        
        loadAct = QtWidgets.QAction(" &Load", self)
        loadAct.triggered.connect(self.fileDialog)
        
        savePlotAct = QtWidgets.QAction(" &Save Plot", self)
        savePlotAct.triggered.connect(self.savePlot)
        
        
        #SETUP MENU
        menubar = self.menuBar()
        #Necessary for OSX, which trys to put menu bar on the top of the screen
        menubar.setNativeMenuBar(False) 
        menubar.addAction(quitAct)
        menubar.addAction(loadAct)
        menubar.addAction(savePlotAct)
        

        self.centerbox = QtWidgets.QVBoxLayout()
        self.layout.addLayout(self.centerbox)
        
        self.rightbox = QtWidgets.QVBoxLayout()
        self.layout.addLayout(self.rightbox)
        
        self.select_ax_box = QtWidgets.QVBoxLayout()
        self.rightbox.addLayout(self.select_ax_box)
        
        #Make divider line
        divFrame = QtWidgets.QFrame()
        divFrame.setFrameShape(QtWidgets.QFrame.HLine)
        self.rightbox.addWidget(divFrame)
        
        self.axesbox = QtWidgets.QVBoxLayout()
        self.rightbox.addLayout(self.axesbox)
        
        
        #Create the plot-type dropdown box
        self.plottype_box = QtWidgets.QHBoxLayout()
        self.centerbox.addLayout(self.plottype_box)
        
        self.plottype_label = QtWidgets.QLabel("Plot Type: ")
        self.plottype_box.addWidget(self.plottype_label)
        
        self.plottype_field = QtWidgets.QComboBox()
        self.plottype_field.currentIndexChanged.connect(self.updatePlotTypeAction)
        self.plottype_field.addItem('1D')
        self.plottype_field.addItem('2D')
        self.plottype_field.show()
        self.plottype_box.addWidget(self.plottype_field)
        
        
        self.colormap_box = QtWidgets.QHBoxLayout()
        self.centerbox.addLayout(self.colormap_box)
        self.colormap_lbl = QtWidgets.QLabel("Colormap: ")
        self.colormap_box.addWidget(self.colormap_lbl)
        self.colormap_field = QtWidgets.QComboBox()
        self.colormap_box.addWidget(self.colormap_field)
        self.colormap_field.currentIndexChanged.connect(self.makePlot)
        for k in self.colormap_dict.keys():
            self.colormap_field.addItem(k)
        
        
        
        
        self.toplabel = QtWidgets.QLabel('')
        self.centerbox.addWidget(self.toplabel)
        
        self.figure = matplotlib.figure.Figure(figsize=(5, 3))
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setMinimumSize(500, 500)
        self.centerbox.addWidget(self.canvas)
        
        
        #Create the datarange box
        
        self.datarange_box = QtWidgets.QHBoxLayout()
        self.centerbox.addLayout(self.datarange_box)
  
        
        self.datarange_auto = QtWidgets.QCheckBox("Autorange?")
        self.datarange_auto.setChecked(True)  
        self.datarange_auto.stateChanged.connect(self.updateDataRangeAction)
        self.datarange_box.addWidget(self.datarange_auto)
        

        
        self.datarange_lbl = QtWidgets.QLabel("Data Range: ")
        self.datarange_box.addWidget(self.datarange_lbl)
        
        self.datarange_a = QtWidgets.QDoubleSpinBox()
        self.datarange_a.editingFinished.connect(self.updateDataRangeAction)
        self.datarange_box.addWidget(self.datarange_a )
        
        self.datarange_b = QtWidgets.QDoubleSpinBox()
        self.datarange_b.editingFinished.connect(self.updateDataRangeAction)
        self.datarange_box.addWidget(self.datarange_b )
        
        self.datarange_unitlbl = QtWidgets.QLabel("")
        self.datarange_box.addWidget(self.datarange_unitlbl)

        
        #Create the first axis dropdown menu
        self.dropdown1_box = QtWidgets.QHBoxLayout()
        self.select_ax_box.addLayout(self.dropdown1_box)
        
        self.dropdown1_label = QtWidgets.QLabel("Axis 1: ")
        self.dropdown1_box.addWidget(self.dropdown1_label)
        
        self.dropdown1 = QtWidgets.QComboBox()
        self.dropdown1.currentIndexChanged.connect(self.updateAxesFieldsAction)
        self.dropdown1_box.addWidget(self.dropdown1)
        
        #Create the second axis dropdown menu
        self.dropdown2_box = QtWidgets.QHBoxLayout()
        self.select_ax_box.addLayout(self.dropdown2_box)
        
        self.dropdown2_label = QtWidgets.QLabel("Axis 2: ")
        self.dropdown2_box.addWidget(self.dropdown2_label)
        
        self.dropdown2 = QtWidgets.QComboBox()
        self.dropdown2.currentIndexChanged.connect(self.updateAxesFieldsAction)
        self.dropdown2_box.addWidget(self.dropdown2)
        
        
        self.updatePlotType()
        self.updateAxesFields()
        
         
        
        
        
        

    def fileDialog(self):
        self.last_axes = self.axes #Copy over any axes to memory
        self.axes = []

        self.last_cur_axes = self.cur_axes
        self.cur_axes = (0,)
        
        print("load")
        opendialog = QtWidgets.QFileDialog()
        opendialog.setNameFilter("HDF5 Files (*.hf5, *.h5)")
        self.filepath = opendialog.getOpenFileName(self, "Select file to open", "", "hdf5 Files (*.hdf5)")[0]
        print(self.filepath)
  
        with h5py.File(self.filepath, 'r') as f:
            temp_axes = ( f['data'].attrs['dimensions'] ) 
            dataunit = f['data'].attrs['unit']
            
            self.datarange_unitlbl.setText(dataunit)
           
            for ind, ax in enumerate(temp_axes):
                ax_dict = {}
                name = ax.decode("utf-8")
                ax_dict['name'] =  name
                ax_dict['ax'] = f[name][:]
                ax_dict['ind'] = ind
                #ax_dict['step'] =  np.mean(np.gradient(ax_dict['ax']))
                ax_dict['valrange'] = ( f[name][0] , f[name][-1] )
                ax_dict['indrange'] = ( 0 ,  len(f[name]) -1 )
                self.axes.append(ax_dict)
                
        self.fillAxesBox()
        self.updatePlotType()
        self.updateAxesFields()
            
    
    def fillAxesBox(self):

        #Remove old widgets
        self.clearLayout(self.axesbox)
        
        #Remove old items from dropdown menus
        self.dropdown1.clear()
        self.dropdown2.clear()

        
        for i, ax in enumerate(self.axes):
            #Take the ax_dict out of the axes array
            ax_dict = self.axes[i]
            
            lbl = (str(ax_dict['name']) + ' : [' + 
                   str(int(ax_dict['indrange'][0])) + ',' + 
                   str(int(ax_dict['indrange'][1])) + ']' )
   
            #Add the axes names to the dropdown menus
            self.dropdown1.addItem(ax_dict['name'])
            self.dropdown2.addItem(ax_dict['name'])

            ax_dict['box'] = QtWidgets.QHBoxLayout()
            self.axesbox.addLayout(ax_dict['box'])
            
            ax_dict['label']  = QtWidgets.QLabel( lbl )  
            ax_dict['box'].addWidget(ax_dict['label'])
            
            #Add the first range field
            ax_dict['field1']  = QtWidgets.QSpinBox()
            ax_dict['field1'].editingFinished.connect(self.updateAxRangeAction)
            ax_dict['field1'].setRange(ax_dict['indrange'][0],ax_dict['indrange'][1])
            ax_dict['field1'].setSingleStep(1)
            ax_dict['field1'].setMinimumWidth(70)
            ax_dict['box'].addWidget(ax_dict['field1'])
        
            
            #Add the second range field (will only show if this dimension is an axis)
            ax_dict['field2']  = QtWidgets.QSpinBox()
            ax_dict['field2'].editingFinished.connect(self.updateAxRangeAction)
            ax_dict['field2'].setRange(ax_dict['indrange'][0],ax_dict['indrange'][1])
            ax_dict['field2'].setSingleStep(1)
            ax_dict['field2'].setMinimumWidth(70)
            ax_dict['box'].addWidget(ax_dict['field2'])

            #Put the ax_dict back into the axes array
            self.axes[i] = ax_dict
        
        #If names of any old axes match those of any new axes
        #attempt to copy over the currently chosen indices
        for i, ax in enumerate(self.axes):
            for j, old_ax in enumerate(self.last_axes):
                if ax['name'] == old_ax['name']:
                    ax['field1'].setValue( old_ax['field1'].value() )
                    ax['field2'].setValue( old_ax['field2'].value() )
                    
        #If new axes match old ones, set the indices so the axes stay
        #the same
        if len(self.last_axes) != 0:
            cur_name = self.last_axes[ self.last_cur_axes[0] ]['name']
            ind = self.dropdown1.findText(cur_name)
            if ind != -1:
                self.dropdown1.setCurrentIndex(ind)
                
            cur_name = self.last_axes[ self.last_cur_axes[1] ]['name']
            ind = self.dropdown2.findText(cur_name)
            if ind != -1:
                self.dropdown2.setCurrentIndex(ind)

        
        
            
            
    def updatePlotType(self):
        ind = self.plottype_field.currentIndex()
        if ind == 0:
            self.plottype = 1
            try:
                self.dropdown2_label.hide()
                self.dropdown2.hide()
            except AttributeError as e:
                pass
            
        elif ind == 1:
            self.plottype = 2
            try:
                self.dropdown2_label.show()
                self.dropdown2.show()
            except AttributeError as e:
                pass

        
    def updatePlotTypeAction(self):
        self.updatePlotType()
        self.updateAxesFields()
        self.makePlot()

            
            
    def updateAxesFields(self):
        try:
            if self.plottype == 1:
                self.cur_axes = (self.dropdown1.currentIndex(), )
            elif self.plottype ==2:
                self.cur_axes = (self.dropdown1.currentIndex(), self.dropdown2.currentIndex())
            
            for i, ax in enumerate(self.axes):
                if 'field2' in ax.keys():
                    if i in self.cur_axes:
                        ax['field2'].show()
                    else:
                        ax['field2'].hide()
        except AttributeError:
            pass
        
    def updateAxesFieldsAction(self):
        self.updateAxesFields()
        self.makePlot()
        
        
    def updateDataRange(self):
        try:
            self.datarange[0] = float(self.datarange_a.text())
            self.datarange[1] = float(self.datarange_b.text())
        except AttributeError:
            pass
        
    def updateDataRangeAction(self):
        self.updateDataRange()
        self.makePlot()
        
        
    def updateAxRange(self):
        pass
        
        
    def updateAxRangeAction(self):
        self.updateAxRange()
        self.makePlot()

        

        
    def validateChoices(self):
        #Try statement catches attribute error for if any of these
        #aren't defined yet
        try:
            
            #Validate file
            if not os.path.isfile(self.filepath):
                self.toplabel.setText("WARNING: Invalid filepath!")
                return False
            elif os.path.splitext(self.filepath)[-1] != '.hdf5':
                self.toplabel.setText("WARNING: Filepath must end in .hdf5!")
                return False
            
            #Validate plot params
            if  self.plottype == 2 and self.dropdown1.currentIndex() == self.dropdown2.currentIndex():
                self.toplabel.setText("WARNING: Axes selected need to be different!")
                return False
            
            if len(self.axes) == 0:
                return False

            for ind, ax_dict in enumerate(self.axes):
                if ind in self.cur_axes:
                    print(ind)
                    if ax_dict['field1'].text()  == ax_dict['field2'].text():
                        self.toplabel.setText("WARNING: Axes range is 0!")
                        return False
                    elif int(ax_dict['field1'].text())  > int(ax_dict['field2'].text()):
                        self.toplabel.setText("WARNING: First range element should be smallest!")
                        return False
                
            self.toplabel.setText("")
        except (AttributeError, KeyError):
            return False
        
        
        return True
    


        
        
    def makePlot(self):
        self.clearCanvas()
        if self.validateChoices():
            if self.plottype == 1:
                self.plot1D()
            elif self.plottype == 2:
                self.plot2D()

            
            
            
    def clearCanvas(self):
        try:
            self.figure.clf()
            self.canvas_ax.clear()
            self.canvas.draw()
        except AttributeError as e:
            pass
        
    def clearLayout(self, layout):
        if layout !=None:
            while layout.count():
                child = layout.takeAt(0)
                if child.widget() is not None:
                    child.widget().deleteLater()
                elif child.layout() is not None:
                    self.clearLayout(child.layout())
        
        

    def plot1D(self):
        print("plot 1d")
        
        #Horizontal axis for this 1D plot - obj
        ax_ind = self.cur_axes[0]

        dslice = []
        
        for i in range(len(self.axes) ):
            d  = self.axes[i]
            if i == ax_ind:
                hname = d['name']
                a = int(d['field1'].text())
                b = int(d['field2'].text())
                dslice.append( slice(a, b, 1) )
                hslice = slice(a,b, 1)
                
            else:
                a = int(d['field1'].text())
                dslice.append( slice(a, a+1, 1) )
                
        with h5py.File(self.filepath, 'r') as f:
            hax = np.squeeze(f[hname][hslice])#already a tuple
            data = np.squeeze(f['data'][tuple(dslice)])
            dataunit  = str(f['data'].attrs['unit'])
            hunit = str(f[hname].attrs['unit'])
            
            
        
            
            

        self.canvas_ax = self.canvas.figure.subplots()
        #t = np.linspace(0, 10, 501)
        self.canvas_ax.plot(hax, data, linestyle='-')
        
        self.canvas_ax.set_xlabel(str(hname) + ' (' + str(hunit) + ')')
        self.canvas_ax.set_ylabel('(' + str(dataunit) + ')')

        
        #Setup axis formats
        self.canvas_ax.ticklabel_format(axis='x', scilimits=(-3,3) )
        self.canvas_ax.ticklabel_format(axis='y', scilimits=(-3,3) )
        
        #Autorange if appropriate
        if not self.datarange_auto.isChecked():
            self.canvas_ax.set_ylim(self.datarange[0], self.datarange[1])
            
        self.canvas.draw()
    
    def plot2D(self):
        
        #Horizontal axis for this 1D plot - obj
        hind = self.cur_axes[0]
        vind = self.cur_axes[1]
  
        dslice = []
        
        for i in range(len(self.axes) ):
            d  = self.axes[i]
            if i == hind:
                hname = d['name']
                a = int(d['field1'].text())
                b = int(d['field2'].text())
                dslice.append( slice(a, b, 1) )
                hslice = slice(a,b, 1)
            elif i == vind:
                vname = d['name']
                a = int(d['field1'].text())
                b = int(d['field2'].text())
                dslice.append( slice(a, b, 1) )
                vslice = slice(a,b, 1)
                
            else:
                a = int(d['field1'].text())
                dslice.append( slice(a, a+1, 1) )
        
        with h5py.File(self.filepath, 'r') as f:
            hax = np.squeeze(f[hname][hslice])#already a tuple
            vax = np.squeeze(f[vname][vslice])#already a tuple
            data = np.squeeze(f['data'][tuple(dslice)])
            dataunit  = str(f['data'].attrs['unit'])
            hunit = f[hname].attrs['unit']
            vunit = f[vname].attrs['unit']
            
        if vind > hind:
            data = data.transpose()
            
        
        
        if not self.datarange_auto.isChecked():
            cmin = self.datarange[0]
            cmax = self.datarange[1]
        else:
            cmin = None
            cmax = None
            
        #levels = np.linspace(cmin, cmax, num=50)
        
        cmap_key = self.colormap_field.currentText()
        cmname = self.colormap_dict[cmap_key]
        colormap = matplotlib.cm.get_cmap(name=cmname)
        
        print("Plotting 2D data")
        self.canvas_ax = self.canvas.figure.subplots()
        
        #self.canvas_ax.contourf(hax, vax, data, levels, cmap=colormap)
        
        cs = self.canvas_ax.imshow(data, origin='upper',
                                   vmin=cmin, vmax=cmax,
                                   aspect='auto', interpolation = 'nearest',
                                   extent=[hax[0], hax[-1], vax[0], vax[-1]],
                                   cmap = colormap)
        
        
        
        
        cbar = self.canvas.figure.colorbar(cs, orientation='horizontal', 
                                           format='%.2f')
        cbar.ax.set_xlabel('(' + str(dataunit) + ')' )
        
        self.canvas_ax.set_xlabel(str(hname) + ' (' + str(hunit) + ')')
        self.canvas_ax.set_ylabel(str(vname) + ' (' + str(vunit) + ')')

        
        #Setup axis formats
        self.canvas_ax.ticklabel_format(axis='x', scilimits=(-3,3) )
        self.canvas_ax.ticklabel_format(axis='y', scilimits=(-3,3) )
        
        
        self.canvas.draw()
        print("Done plotting")
        
        
    def savePlot(self):
        savedialog = QtWidgets.QFileDialog()
        
        
        suggested_name = os.path.splitext(self.filepath)[0] + '.png'

        
        savefile = savedialog.getSaveFileName(self, "Save as: ", suggested_name, "")[0]
        self.figure.savefig(savefile)
        
        
        
        
        

if __name__ == "__main__":
    
    #Check if a QApplicaiton already exists, and don't open a new one if it does
    #This helps avoid kernel crashes on exit.
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
    else:
        print('QApplication instance already exists: %s' % str(app))
    
    w = ApplicationWindow()
    w.show()
    app.exec()