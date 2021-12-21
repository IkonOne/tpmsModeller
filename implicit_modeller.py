import os
import sys
from PySide6.QtWidgets import (
    QApplication,
    QPushButton,
    QComboBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QWidget
)
from PySide6.QtGui import (
    QDoubleValidator,
    QIntValidator,
)

from vtkmodules.vtkFiltersSources import vtkConeSource
from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer
# load implementations for rendering and interaction factory classes
import vtkmodules.vtkRenderingOpenGL2
import vtkmodules.vtkInteractionStyle

from QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

import numpy as np
import vedo
from vedo import *
from vedo.applications import IsosurfaceBrowser

def print_versions():
    import PySide6.QtCore
    print(PySide6.__version__)
    print(PySide6.QtCore.__version__)

def lerp(a, b, t):
    return (1-t)*a + b*t

def invlerp(val, a, b):
    return (val - a) / (b - a)

def map_range(val, min1, max1, min2, max2):
    return lerp(min2, max2, invlerp(val, min1, max1))

def compute_volume(span, resolution, f, g_min=0, g_max=1, offset=0.0):
    x, y, z = span * np.mgrid[g_min:g_max + resolution:resolution,
                              g_min:g_max + resolution:resolution,
                              g_min:g_max + resolution:resolution]
    return f(x,y,z, offset)

def gyroid(x,y,z, offset=0.0):
    def f(x,y,z):
        return np.sin(x)*np.cos(y) + np.sin(y)*np.cos(z) + np.sin(z)*np.cos(x)
    
    def gradient(x,y,z):
        dx = np.cos(x)*np.cos(y) - np.sin(z)*np.sin(x)
        dy = -np.sin(x)*np.sin(y) + np.cos(y)*np.cos(z)
        dz = -np.sin(y)*np.sin(z) + np.cos(z)*np.cos(x)
        grad = np.stack((dx, dy, dz), axis=-1)
        # grad = grad / np.linalg.norm(grad, axis=-1, keepdims=True)
        return grad
    
    if offset == 0.0:
        return f(x,y,z)

    grad = gradient(x,y,z)**2
    grad = np.sqrt(np.sum(grad, axis=-1))
    v = f(x,y,z) - offset*grad

    # Creates a second surface equidistant from the zero level, but in the other direction
    # v = np.maximum(v, np.flip(v))

    # v[[0,v.shape[0]-1],:,:] = 0.01
    # v[:,[0,v.shape[1]-1],:] = 0.01
    # v[:,:,[0,v.shape[2]-1]] = 0.01

    return v

def fks(x,y,z,offset=0):
    def f(x,y,z):
        return np.cos(2*x)*np.sin(y)*np.cos(z) + np.cos(2*y)*np.sin(z)*np.cos(x) + np.cos(2*z)*np.sin(x)*np.cos(y)

    def gradient(x,y,z):
        dx = -2*np.sin(2*x)*np.sin(y)*np.cos(z) + np.cos(2*y)*np.sin(z)*-np.sin(x) + np.cos(2*z)*np.cos(x)*np.cos(y)
        dy = np.cos(2*x)*np.cos(y)*np.cos(z) + -2*np.sin(2*y)*np.sin(z)*np.cos(x) + np.cos(2*z)*np.sin(x)*-np.sin(y)
        dz = np.cos(2*x)*np.sin(y)*-np.sin(z) + np.cos(2*y)*np.cos(z)*np.cos(x) + -2*np.sin(2*z)*np.sin(x)*np.cos(y)
        grad = np.stack((dx, dy, dz), axis=-1)
        # grad = grad / np.linalg.norm(grad, axis=-1, keepdims=True)
        return grad

    if offset == 0.0:
        return f(x,y,z)

    grad = gradient(x,y,z)**2
    grad = np.sqrt(np.sum(grad, axis=-1))
    v = f(x,y,z) - offset*grad

    # Creates a second surface equidistant from the zero level, but in the other direction
    # v = np.maximum(v, np.flip(v))

    # v[[0,v.shape[0]-1],:,:] = 0.01
    # v[:,[0,v.shape[1]-1],:] = 0.01
    # v[:,:,[0,v.shape[2]-1]] = 0.01

    return v

class TPMSRenderer:
    def __init__(self, plotter, fn, periods=1.0, dimensions=50.0, roadWidth=0.5, roadCount=3, resolution=0.05) -> None:
        self._plotter = plotter
        self._isos=[]
        self._fn = fn
        self._resolution = resolution
        self._periods = periods
        self._roadWidth = roadWidth
        self._roadCount = roadCount
        self._dimensions = dimensions
        self._update()
    
    def _update(self):
        if self._fn == None:
            return
        if self._resolution <= 0:
            return
        if self._periods <= 0:
            return
        if self._roadWidth < 0:
            return
        if self._dimensions <= 0:
            return
        
        rad_2_mm = self._dimensions / (self._periods*2*np.pi)
        mm_2_rad = (self._periods*2*np.pi)/self._dimensions

        r = self._resolution
        scale = 2/r
        
        if self._roadCount <= 1:
            halfCount = 0
        elif self._roadCount == 2:
            halfCount = self._roadWidth / 2 * mm_2_rad
        else:
            rw = self._roadWidth * mm_2_rad
            halfCount = self._roadCount*rw/2 - rw
        offsets = np.linspace(-halfCount, halfCount, self._roadCount)
        
        x,y,z = np.mgrid[
            -1:1+r:r,
            -1:1+r:r,
            -1:1+r:r] * np.pi * self._periods
        
        for isos in self._isos:
            self._plotter -= isos

        n = len(offsets)
        if n == 1:
            cols = colorMap(0, 'rainbow', vmin=0, vmax=1)
        elif n == 2:
            cols = colorMap([0,1], 'rainbow', vmin=0, vmax=1)
        else:
            cols = colorMap(range(n), 'rainbow')
        for i in range(n):
            U = self._fn(x,y,z, offsets[i])
            isos = Volume(U, dims=U.shape).isosurface(threshold=0)
            isos.scale(1/scale * self._dimensions)
            isos.c(cols[i])
            self._plotter += isos
            self._isos.append(isos)
        
        self._plotter.window.Render()
             
    def setTPMSfunc(self, fn):
        self._fn = fn
        self._update()
    
    def setDimensions(self, dimensions):
        self._dimensions = float(dimensions)
        self._update()
    
    def setPeriods(self, periods):
        self._periods = float(periods)
        self._update()
    
    def setResolution(self, resolution):
        self._resolution = float(resolution)
        self._update()
    
    def setRoadWidth(self, width):
        self._roadWidth = float(width)
        self._update()

    def setRoadCount(self, roadCount):
        self._roadCount = int(roadCount)
        self._update()

def main(argv):
    # every QT app needs an app
    app = QApplication(['Implicit Surface Model Generator'])

    window = QMainWindow()
    window.setMinimumSize(800, 600)

    frame = QFrame()
    window.setCentralWidget(frame)

    main_layout = QHBoxLayout()
    frame.setLayout(main_layout)

    editor_widget = QWidget()
    editor_layout = QGridLayout()
    editor_layout.setColumnMinimumWidth(0, 50)
    editor_layout.setColumnMinimumWidth(1, 50)
    editor_widget.setLayout(editor_layout)
    main_layout.addWidget(editor_widget, stretch=0)

    row = 0
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("TPMS"))
    combo = QComboBox()
    combo.addItems(['Gyroid', 'Fisher-Koch S'])
    editor_layout.addWidget(combo)

    row = 1
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("Dimensions (mm)"), row, 0)
    dimensions_text = QLineEdit("50.0")
    dimensions_text.setValidator(QDoubleValidator(0.0, 1000.0, 4, dimensions_text))
    editor_layout.addWidget(dimensions_text)

    row = 2
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("Periods"), row, 0)
    period_text = QLineEdit("1.0")
    period_text.setValidator(QDoubleValidator(0.0, 100.0, 4, period_text))
    editor_layout.addWidget(period_text)

    row = 3
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("Road Width (mm)"), row, 0)
    road_width_text = QLineEdit("0.5")
    road_width_text.setValidator(
        QDoubleValidator(0.001, 100.0, 4, road_width_text)
    )
    editor_layout.addWidget(road_width_text)

    row = 4
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("Road Count"), row, 0)
    road_count_text = QLineEdit("3")
    road_count_text.setValidator(
        QIntValidator(0, 50, road_count_text)
    )
    editor_layout.addWidget(road_count_text)

    row = 5
    editor_layout.setRowStretch(row, 0)
    editor_layout.setRowMinimumHeight(row, 20)
    editor_layout.addWidget(QLabel("Resolution"), row, 0)
    resoltuion_text = QLineEdit("0.05")
    resoltuion_text.setValidator(
        QDoubleValidator(0.001, 1.0, 4, resoltuion_text)
    )
    editor_layout.addWidget(resoltuion_text)

    row = 6
    editor_layout.setRowStretch(row,0)
    editor_layout.setRowMinimumHeight(row, 20)
    export_btn = QPushButton('Export')
    def export_tpms(e):
        filename, _ = QFileDialog.getSaveFileName(
            window, 'Export STL', f'{os.getcwd()}/{combo.currentText}_{period_text.text}.stl', 'STL File (*.stl)'
        )
        vedo.io.write(merge(tpms._plotter.actors), filename)
    export_btn.clicked.connect(export_tpms)

    editor_layout.addWidget(export_btn)

    vtk_widget = QVTKRenderWindowInteractor(window)
    main_layout.addWidget(vtk_widget, stretch=1)

    # if you don't want the 'q' key to exit comment this.
    # widget.AddObserver("ExitEvent", lambda o, e, a=app: a.quit())

    vp = Plotter(qtWidget=vtk_widget)

    fns = [gyroid, fks]
    tpms = TPMSRenderer(vp, fns[0], 1.0, 50.0, 0.5)
    combo.currentIndexChanged.connect(lambda idx : tpms.setTPMSfunc(fns[idx]))
    dimensions_text.textChanged.connect(tpms.setDimensions)
    period_text.textChanged.connect(tpms.setPeriods)
    road_width_text.textChanged.connect(tpms.setRoadWidth)
    road_count_text.textChanged.connect(tpms.setRoadCount)
    resoltuion_text.textChanged.connect(tpms.setResolution)

    vp.show(bg='blackboard', axes=8)

    # show the widget
    window.show()

    vtk_widget.Initialize()
    vtk_widget.Start()

    # start event processing
    # Source: https://doc.qt.io/qtforpython/porting_from2.html
    # 'exec_' is deprecated and will be removed in the future.
    # Use 'exec' instead.
    try:
        app.exec()
    except AttributeError:
        app.exec_()

if __name__ == "__main__":
    main(sys.argv)