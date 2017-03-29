#-------------------------------------------------
#
# Project created by QtCreator 2014-06-24T12:13:18
#
#-------------------------------------------------

macx {
  QMAKE_CXXFLAGS+="-g -gdwarf-2 -ftemplate-depth=1024 "
  QMAKE_POST_LINK='/usr/bin/dsymutil RIVET.app/Contents/MacOS/RIVET -o RIVET.app/Contents/MacOS/RIVET.dsym'
}

CONFIG += c++11 debug

QT       +=    core gui \
                widgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = RIVET
TEMPLATE = app

QMAKE_LIBDIR += /usr/local/lib #TODO: figure out how to generalize
LIBS        += -lboost_serialization

SOURCES	+= main.cpp                         \
		visualizationwindow.cpp             \
		dataselectdialog.cpp                \
		dcel/dcel.cpp                       \
		dcel/arrangement.cpp                       \
		interface/control_dot.cpp           \
		#interface/input_manager.cpp         \
		interface/persistence_bar.cpp       \
		interface/persistence_diagram.cpp   \
		interface/persistence_dot.cpp       \
		interface/slice_diagram.cpp         \
		interface/slice_line.cpp            \
	    math/bool_array.cpp                 \
		math/index_matrix.cpp               \
		math/map_matrix.cpp                 \
		#math/multi_betti.cpp                \
		#math/simplex_tree.cpp               \
		#math/st_node.cpp                    \
		dcel/barcode.cpp               \
		dcel/barcode_template.cpp           \
		dcel/anchor.cpp                     \
		dcel/arrangement_message.cpp               \
		dcel/grades.cpp                     \
		#math/persistence_updater.cpp        \
		math/template_points_matrix.cpp          \
		math/template_point.cpp                   \
		interface/progressdialog.cpp        \
		computationthread.cpp               \
		interface/aboutmessagebox.cpp       \
		interface/configuredialog.cpp       \
		interface/config_parameters.cpp     \
		interface/file_input_reader.cpp \
    #driver.cpp \
    interface/file_writer.cpp \
    debug.cpp \
    timer.cpp \
    interface/console_interaction.cpp \
    numerics.cpp \
    dendrogram_data.cpp \
    dcel/arrangement_builder.cpp \
    math/persistence_updater.cpp \
    math/simplex_tree.cpp \
    math/st_node.cpp \
    computation.cpp \
    interface/input_manager.cpp \
    math/multi_betti.cpp \
    # front_end/qneblock.cpp \
    # front_end/qneconnection.cpp \
    # front_end/qnemainwindow.cpp \
    # front_end/qneport.cpp \
    # front_end/qnodeseditor.cpp
    qnodeseditor/qneblock.cpp \
    qnodeseditor/qneconnection.cpp \
    qnodeseditor/qnemainwindow.cpp \
    qnodeseditor/qneport.cpp \
    qnodeseditor/qnodeseditor.cpp \
    math/bifiltration_data.cpp


HEADERS  += visualizationwindow.h			\
		dataselectdialog.h					\
		dcel/dcel.h							\
		dcel/arrangement.h							\
		interface/control_dot.h				\
		interface/input_manager.h			\
		interface/persistence_bar.h			\
		interface/persistence_diagram.h		\
		interface/persistence_dot.h			\
		interface/slice_diagram.h			\
		interface/slice_line.h				\
		math/bool_array.h                 \
		math/index_matrix.h					\
		math/map_matrix.h					\
		math/multi_betti.h					\
		math/simplex_tree.h					\
		math/st_node.h						\
		dcel/barcode.h	    				\
		dcel/barcode_template.h				\
		dcel/anchor.h						\
		dcel/grades.h                       \
		math/persistence_updater.h			\
		math/template_points_matrix.h			\
		math/template_point.h \
    interface/progressdialog.h \
    computationthread.h \
    interface/input_parameters.h \
    interface/aboutmessagebox.h \
    interface/configuredialog.h \
    interface/config_parameters.h \
    interface/file_input_reader.h \
    #driver.h \
    interface/file_writer.h \
    cutgraph.h \
    dcel/serialization.h \
    interface/console_interaction.h \
    numerics.h \
    dendrogram_data.h \
    time_root.h \
    computing_s.h \
    Graph.h \
    subset.h \
    dcel/arrangement_builder.h \
    dcel/arrangement_message.h \
    computation.h \
    # front_end/qneblock.h \
    # front_end/qneconnection.h \
    # front_end/qneport.h \
    # front_end/qnodeseditor.h \
    # front_end/qnemainwindow.h
    dendrogram_viz.h \
    qnodeseditor/IntervalTree.h \
    qnodeseditor/qnemainwindow.h \
    qnodeseditor/qnodeseditor.h \
    qnodeseditor/qneconnection.h \
    qnodeseditor/qneport.h \
    qnodeseditor/qneblock.h \
    math/bifiltration_data.h

FORMS   += visualizationwindow.ui			\
		dataselectdialog.ui \
    interface/progressdialog.ui \
    interface/aboutmessagebox.ui \
    interface/configuredialog.ui

# QMAKE_MAC_SDK = macosx10.12

CONFIG += c++11
QMAKE_CFLAGS += -std=c++11 -stdlib=libc++ -mmacosx-version-min=10.8
QMAKE_CXXFLAGS += -std=c++11 -stdlib=libc++ -mmacosx-version-min=10.8

QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9


LIBS += -L"/usr/local/Cellar/boost/1.60.0_2/lib"
INCLUDEPATH += "/usr/local/Cellar/boost/1.60.0_2/include"
LIBS += -L"/usr/local/Cellar/boost/1.60.0_2/lib" -lboost_random

INCLUDEPATH += "../qt4/QtCore"
