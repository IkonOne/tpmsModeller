import QtQuick
import QtQuick.Controls
import QtQuick.Layouts
from QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

Rectangle {
    id: base

    property alias implicitFunctionIndex: implicitFunctionCombo.currentIndex
    property alias periods: periods.text
    property alias surfaceWidth: surfaceWidth.text
    property alias dimensions: dimensions.text
    // property alias resolution: resolution.text

    RowLayout {
        ColumnLayout {
            RowLayout {
                Label { text: "Implicit Surface Function" }
                ComboBox {
                    id : implicitFunctionCombo
                    model : ListModel {
                        ListElement { text: "Gyroid" }
                        ListElement { text: "Fisher Koch S" }
                    }
                }
            }

            RowLayout {
                Label { text: "Periods" }
                TextField {
                    id: periods
                    text: "1.0"
                    validator: DoubleValidator { decimals: 4; bottom: 0.0; top: 100.0 }
                }
            }

            RowLayout {
                Label { text: "Surface Width (mm)" }
                TextField {
                    id: surfaceWidth
                    text: "1.0"
                    validator: DoubleValidator { decimals: 4; bottom: 0.0; top: 10000.0 }
                }
            }


            RowLayout {
                Label { text: "Dimensions (mm)" }
                TextField {
                    id: dimensions
                    text: "100"
                    validator: DoubleValidator { decimals: 4; bottom: 0.0; top: 10000.0 }
                }
            }

            // RowLayout {
            //     Label { text: "Resolution (mm)" }
            //     TextField {
            //         id: resolution
            //         text: "0.01"
            //         validator: DoubleValidator { decimals: 4; bottom: 0.0; top: 10000.0 }
            //     }
            // }

        //     RowLayout {
        //         Button {
        //             text: "Accept"
        //             onClicked: base.accept()
        //         }

        //         Button {
        //             text: "Cancel"
        //             onClicked: base.reject()
        //         }
        //     }
        }

        QVTKRenderWindowInteractor {

        }
    }
}