import os
import sys
from PySide6.QtWidgets import (
    QApplication,
    QPushButton,
    QComboBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QVBoxLayout,
    QWidget
)
from PySide6.QtGui import (
    QDoubleValidator,
    QIntValidator,
)
from PySide6.QtCore import (
    Qt
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
        self._porosity = 0.0
        self.cb_updated = None
        self.update()
    
    def update(self):
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
        
        V = -np.ones_like(x)
        rw = self._roadWidth * mm_2_rad * 0.5
        max_offset = np.max(offsets) + rw
        V = self._fn(x, y, z, max_offset)
        min_offset = np.min(offsets) - rw
        V = np.maximum(V, -self._fn(x, y, z, min_offset))
        empty_voxels = np.count_nonzero(V > 0)
        total_voxels = V.shape[0]*V.shape[1]*V.shape[2]
        self._porosity = empty_voxels / total_voxels
        
        self._plotter.window.Render()

        if self.cb_updated:
            self.cb_updated()
             
    def setTPMSfunc(self, fn):
        self._fn = fn
        self.update()
    
    def setDimensions(self, dimensions):
        self._dimensions = float(dimensions)

    def setPeriods(self, periods):
        self._periods = float(periods)

    def setResolution(self, resolution):
        self._resolution = float(resolution)

    def setRoadWidth(self, width):
        self._roadWidth = float(width)

    def setRoadCount(self, roadCount):
        self._roadCount = int(roadCount)

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
    editor_layout = QVBoxLayout()
    editor_widget.setLayout(editor_layout)
    main_layout.addWidget(editor_widget, stretch=0)

    vtk_widget = QVTKRenderWindowInteractor(window)
    main_layout.addWidget(vtk_widget, stretch=1)

    vp = Plotter(qtWidget=vtk_widget)
    fns = [gyroid, fks]
    tpms = TPMSRenderer(vp, fns[0], 1.0, 50.0, 0.5)

    def export_tpms(e):
        filename, _ = QFileDialog.getSaveFileName(
            window, 'Export STL', f'{os.getcwd()}/{combo.currentText}_{period_text.text}.stl', 'STL File (*.stl)'
        )
        vedo.io.write(merge(tpms._plotter.actors), filename)
    export_btn = QPushButton('Export')
    export_btn.clicked.connect(export_tpms)
    editor_layout.addWidget(export_btn)

    hline = QFrame()
    hline.setFrameShape(QFrame.Shape.HLine)
    hline.setFrameShadow(QFrame.Shadow.Sunken)
    editor_layout.addWidget(hline)

    tpms_layout = QHBoxLayout()
    editor_layout.addLayout(tpms_layout)
    tpms_layout.addWidget(QLabel("TPMS"))
    combo = QComboBox()
    combo.addItems(['Gyroid', 'Fisher-Koch S'])
    tpms_layout.addWidget(combo)

    dimensions_layout = QHBoxLayout()
    editor_layout.addLayout(dimensions_layout)
    dimensions_layout.addWidget(QLabel("Dimensions (mm)"))
    dimensions_text = QLineEdit("50.0")
    dimensions_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    dimensions_text.setValidator(QDoubleValidator(0.0, 1000.0, 4, dimensions_text))
    dimensions_layout.addWidget(dimensions_text)

    period_layout = QHBoxLayout()
    editor_layout.addLayout(period_layout)
    period_layout.addWidget(QLabel("Periods"))
    period_text = QLineEdit("1.0")
    period_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    period_text.setValidator(QDoubleValidator(0.0, 100.0, 4, period_text))
    period_layout.addWidget(period_text)

    road_width_layout = QHBoxLayout()
    editor_layout.addLayout(road_width_layout)
    road_width_layout.addWidget(QLabel("Road Width (mm)"))
    road_width_text = QLineEdit("0.5")
    road_width_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    road_width_text.setValidator(
        QDoubleValidator(0.001, 100.0, 4, road_width_text)
    )
    road_width_layout.addWidget(road_width_text)

    road_count_layout = QHBoxLayout()
    editor_layout.addLayout(road_count_layout)
    road_count_layout.addWidget(QLabel("Road Count"))
    road_count_text = QLineEdit("3")
    road_count_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    road_count_text.setValidator(
        QIntValidator(0, 50, road_count_text)
    )
    road_count_layout.addWidget(road_count_text)

    resolution_layout = QHBoxLayout()
    editor_layout.addLayout(resolution_layout)
    resolution_layout.addWidget(QLabel("Resolution"))
    resolution_text = QLineEdit("0.05")
    resolution_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    resolution_text.setValidator(
        QDoubleValidator(0.001, 1.0, 4, resolution_text)
    )
    resolution_layout.addWidget(resolution_text)

    export_btn = QPushButton('Update')
    export_btn.clicked.connect(tpms.update)
    editor_layout.addWidget(export_btn)

    dimensions_text.returnPressed.connect(export_btn.click)
    period_text.returnPressed.connect(export_btn.click)
    road_width_text.returnPressed.connect(export_btn.click)
    road_count_text.returnPressed.connect(export_btn.click)
    resolution_text.returnPressed.connect(export_btn.click)

    hline = QFrame()
    hline.setFrameShape(QFrame.Shape.HLine)
    hline.setFrameShadow(QFrame.Shadow.Sunken)
    editor_layout.addWidget(hline)

    porosity_layout = QHBoxLayout()
    editor_layout.addLayout(porosity_layout)
    porosity_label = QLabel('Porosity')
    porosity_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
    porosity_layout.addWidget(porosity_label)
    porosity_text = QLabel('0.0')
    porosity_text.setAlignment(Qt.AlignmentFlag.AlignRight)
    porosity_layout.addWidget(porosity_text)

    controls_layout = QHBoxLayout()
    editor_layout.addLayout(controls_layout)
    controls_label = QLabel('''Controls:
        Enter: When editing a textbox, updates model
        R: Center view on model
        W: Wireframe rendering
        S: Solid rendering
        Shift+F: Zoom to mouse on model
        Mouse1: Rotate view around center
        Scroll Wheel: Zoom in/out
        Shift + Mouse1: Move view
    ''')
    controls_label.setMaximumHeight(160)
    controls_layout.addWidget(controls_label)

    # if you don't want the 'q' key to exit comment this.
    # widget.AddObserver("ExitEvent", lambda o, e, a=app: a.quit())

    combo.currentIndexChanged.connect(lambda idx : tpms.setTPMSfunc(fns[idx]))
    dimensions_text.textChanged.connect(tpms.setDimensions)
    period_text.textChanged.connect(tpms.setPeriods)
    road_width_text.textChanged.connect(tpms.setRoadWidth)
    road_count_text.textChanged.connect(tpms.setRoadCount)
    resolution_text.textChanged.connect(tpms.setResolution)

    def handle_tpms_update():
        porosity_text.setText(f'%{tpms._porosity*100:2f}')
    tpms.cb_updated = handle_tpms_update
    tpms.update()

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