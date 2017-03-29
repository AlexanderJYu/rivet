/**********************************************************************
Copyright 2014-2016 The RIVET Developers. See the COPYRIGHT file at
the top-level directory of this distribution.

This file is part of RIVET.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
**********************************************************************/

#include "visualizationwindow.h"
#include "ui_visualizationwindow.h"

#include "dcel/arrangement_message.h"
#include "dcel/barcode.h"
#include "dcel/barcode_template.h"
#include "interface/config_parameters.h"
#include "interface/file_writer.h"
#include "numerics.h"
#include "qnodeseditor/qnodeseditor.h"

#include <QDateTime>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include <QMessageBox>
#include <QTime>

#include <algorithm>
#include <fstream>
#include <sstream>

const QString VisualizationWindow::DEFAULT_SAVE_DIR_KEY("default_save_dir");

VisualizationWindow::VisualizationWindow(InputParameters& params)
    : QMainWindow()
    , ui(new Ui::VisualizationWindow)
    , verbosity(params.verbosity)
    , data_selected(false)
    , unsaved_data(false)
    , input_params(params)
    , config_params()
    , ds_dialog(input_params, this)
    , grades()
    , angle_precise(0)
    , offset_precise(0)
    , template_points()
    , cthread(input_params)
    , prog_dialog(this)
    , line_selection_ready(false)
    , slice_diagram(&config_params, grades.x, grades.y, this)
    , slice_update_lock(false)
    , p_diagram(&config_params, this)
    , persistence_diagram_drawn(false)
    , dendrogram_diagram(NULL,&config_params,this)
    , dendrogram_arrangement_received(false)
{
    ui->setupUi(this);

    //set up the slice diagram
    ui->sliceView->setScene(&slice_diagram);
    //    ui->sliceView->setDragMode(QGraphicsView::ScrollHandDrag);
    ui->sliceView->scale(1, -1);
    ui->sliceView->setRenderHint(QPainter::Antialiasing);

    //set up the persistence diagram scene
    ui->pdView->setScene(&p_diagram);
    ui->pdView->scale(1, -1);
    ui->pdView->setRenderHint(QPainter::Antialiasing);

    //set up the dendrogram diagram scene
    std::cout << "setting up the dendrogram diagram" << std::endl;
    ui->dendrogramView->setScene(&dendrogram_diagram);
    ui->dendrogramView->scale(1,-1);
    ui->dendrogramView->setRenderHint(QPainter::Antialiasing);
    ui->xposBox->setMinimum(-100.00); // MAGIC NUMBER
    std::cout << "done setting up the dendrogram diagram" << std::endl;

    ui->clusterText->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Fixed);
    // Style just to make it clear that the widget is
    // being resized to fit the parent, it doesn't "overflow"
    ui->clusterText->setFrameShape(QFrame::Box);
    ui->clusterText->setFrameShadow(QFrame::Raised);

    //connect signal from DataSelectDialog to start the computation
    QObject::connect(&ds_dialog, &DataSelectDialog::dataSelected, this, &VisualizationWindow::start_computation);

    //connect signals from ComputationThread to slots in VisualizationWindow
    QObject::connect(&cthread, &ComputationThread::advanceProgressStage, &prog_dialog, &ProgressDialog::advanceToNextStage);
    QObject::connect(&cthread, &ComputationThread::setProgressMaximum, &prog_dialog, &ProgressDialog::setStageMaximum);
    QObject::connect(&cthread, &ComputationThread::setCurrentProgress, &prog_dialog, &ProgressDialog::updateProgress);
    QObject::connect(&cthread, &ComputationThread::templatePointsReady, this, &VisualizationWindow::paint_template_points);
    QObject::connect(&cthread, &ComputationThread::arrangementReady, this, &VisualizationWindow::augmented_arrangement_ready);
    QObject::connect(&cthread, &ComputationThread::dendrogramArrangementReady, this, &VisualizationWindow::dendrogram_arrangement_ready);
    QObject::connect(&cthread, &ComputationThread::finished, &prog_dialog, &ProgressDialog::setComputationFinished);

    //connect signals and slots for the diagrams
    QObject::connect(&slice_diagram, &SliceDiagram::set_line_control_elements, this, &VisualizationWindow::set_line_parameters);
    QObject::connect(&slice_diagram, &SliceDiagram::persistence_bar_selected, &p_diagram, &PersistenceDiagram::receive_dot_selection);
    QObject::connect(&slice_diagram, &SliceDiagram::persistence_bar_deselected, &p_diagram, &PersistenceDiagram::receive_dot_deselection);
    QObject::connect(&p_diagram, &PersistenceDiagram::persistence_dot_selected, &slice_diagram, &SliceDiagram::receive_bar_selection);
    QObject::connect(&p_diagram, &PersistenceDiagram::persistence_dot_secondary_selection, &slice_diagram, &SliceDiagram::receive_bar_secondary_selection);
    QObject::connect(&p_diagram, &PersistenceDiagram::persistence_dot_deselected, &slice_diagram, &SliceDiagram::receive_bar_deselection);

    std::cout << "setting up the dendrogram connections" << std::endl;
    QObject::connect(&dendrogram_diagram, &QNEMainWindow::set_cluster_label, this, &VisualizationWindow::set_cluster_label);
    QObject::connect(&dendrogram_diagram, &QNEMainWindow::set_xpos_box, this, &VisualizationWindow::set_xpos_box);
    std::cout << "done setting up the dendrogram connections" << std::endl;
    //connect other signals and slots
    QObject::connect(&prog_dialog, &ProgressDialog::stopComputation, &cthread, &ComputationThread::terminate); ///TODO: don't use QThread::terminate()! modify ComputationThread so that it can stop gracefully and clean up after itself
}

VisualizationWindow::~VisualizationWindow()
{
    delete ui;
}

//slot that starts the persistent homology computation in a new thread
void VisualizationWindow::start_computation()
{
    data_selected = true;

    //show the progress box
    prog_dialog.show();
    prog_dialog.activateWindow();
    prog_dialog.raise();

    //start the computation in a new thread
    cthread.compute();

} //end start_computation()

//this slot is signaled when the xi support points are ready to be drawn
void VisualizationWindow::paint_template_points(std::shared_ptr<TemplatePointsMessage> points)
{
    qDebug() << "VisualizationWindow: Received template points";

    template_points = points;

    //first load our local copies of the data
    grades = Grades(template_points->x_exact, template_points->y_exact);

    //send xi support points to the SliceDiagram
    slice_diagram.clear_points();
    for (auto point : template_points->template_points)
        slice_diagram.add_point(grades.x[point.x], grades.y[point.y], point.zero, point.one, point.two);

    //create the SliceDiagram
    if (!slice_diagram.is_created()) {
        slice_diagram.create_diagram(
            QString::fromStdString(template_points->x_label),
            QString::fromStdString(template_points->y_label),
            grades.x.front(), grades.x.back(),
            grades.y.front(), grades.y.back(),
            ui->normCoordCheckBox->isChecked(), template_points->homology_dimensions);
    }

    // create slice_line in dendrogramView
    qDebug() << "create slice_line in dendrogramView";
    std::cout << "create slice_line in dendrogramView" << std::endl;
    slice_line = slice_diagram.get_slice_line();
    QObject::connect(slice_line, &SliceLine::slice_line_released, this, &VisualizationWindow::slice_line_released);

    //create right control dot
    dot_right = slice_diagram.get_dot_right();
    QObject::connect(dot_right, &ControlDot::right_dot_released, this, &VisualizationWindow::right_dot_released);

    //create left control dot
    dot_left = slice_diagram.get_dot_left();
    QObject::connect(dot_left, &ControlDot::left_dot_released, this, &VisualizationWindow::left_dot_released);
    std::cout << "finished slice_line stuff" << std::endl;

    //enable control items
    ui->BettiLabel->setEnabled(true);
    ui->xi0CheckBox->setEnabled(true);
    ui->xi1CheckBox->setEnabled(true);
    ui->xi2CheckBox->setEnabled(true);
    ui->normCoordCheckBox->setEnabled(true);

    //update offset extents
    ///TODO: maybe these extents should be updated dynamically, based on the slope of the slice line
    ui->offsetSpinBox->setMinimum(grades.min_offset());
    ui->offsetSpinBox->setMaximum(grades.max_offset());

    //update status
    line_selection_ready = true;
    ui->statusBar->showMessage("bigraded Betti number visualization ready");
}

//this slot is signaled when the dendrogram arrangement is ready
void VisualizationWindow::dendrogram_arrangement_ready(std::shared_ptr<ArrangementMessage> dendrogram_arrangement)
{
    //receive the arrangement
    qDebug() << "dendrogram_arrangement_ready";
    dendrogram_arrangement_received = true;
    this->dendrogram_arrangement = dendrogram_arrangement;
}

//this slot is signaled when the augmented arrangement is ready
void VisualizationWindow::augmented_arrangement_ready(std::shared_ptr<ArrangementMessage> arrangement)
{
    qDebug() << "augmented_arrangement_ready";
    //receive the arrangement
    this->arrangement = arrangement;

    //TESTING: print arrangement info and verify consistency
    //    arrangement->print_stats();
    //    arrangement->test_consistency();

    //inialize persistence diagram
    p_diagram.create_diagram(QString::fromStdString(input_params.shortName), input_params.dim);

    //get the barcode
    BarcodeTemplate dbc = arrangement->get_barcode_template(angle_precise, offset_precise);
    barcode = dbc.rescale(angle_precise, offset_precise, template_points->template_points, grades);

    //TESTING
    barcode->print();

    //draw the barcode
    double zero_coord = rivet::numeric::project_zero(angle_precise, offset_precise, grades.x[0], grades.y[0]);
    p_diagram.set_barcode(zero_coord, *barcode);
    p_diagram.resize_diagram(slice_diagram.get_slice_length(), slice_diagram.get_pd_scale());

    slice_diagram.draw_barcode(*barcode, zero_coord, ui->barcodeCheckBox->isChecked());

    //enable slice diagram control items
    slice_diagram.enable_slice_line();
    ui->angleLabel->setEnabled(true);
    ui->angleDoubleSpinBox->setEnabled(true);
    ui->offsetLabel->setEnabled(true);
    ui->offsetSpinBox->setEnabled(true);
    ui->barcodeCheckBox->setEnabled(true);

    //update status
    if (verbosity >= 2) {
        qDebug() << "COMPUTATION FINISHED; READY FOR INTERACTIVITY.";
    }
    persistence_diagram_drawn = true;
    ui->statusBar->showMessage("ready for interactive barcode exploration");

    //Enable save menu item
    ui->actionSave->setEnabled(true);
    //if an output file has been specified, then save the arrangement
    if (!input_params.outputFile.empty())
        save_arrangement(QString::fromStdString(input_params.outputFile));
    //TODO: we don't have file reading tools here anymore, so we don't know what kind of file it was
    //Have to rely on console to either a) always save (to tmp file if needed), or b) tell us filetype in the output.
    //    else if(input_params.raw_data)
    //        unsaved_data = true;

} //end augmented_arrangement_ready()

void VisualizationWindow::on_angleDoubleSpinBox_valueChanged(double angle)
{
    if (line_selection_ready && !slice_update_lock) {
        angle_precise = angle;
        slice_diagram.update_line(angle_precise, ui->offsetSpinBox->value());
    }

    update_persistence_diagram();
    qDebug() << "update_dendrogram_diagram 1";
    if (dendrogram_arrangement_received)
    {
        qDebug() << "dendrogram_arrangement_received";
        update_dendrogram_diagram();
    }
}

void VisualizationWindow::on_offsetSpinBox_valueChanged(double offset)
{
    if (line_selection_ready && !slice_update_lock) {
        offset_precise = offset;
        slice_diagram.update_line(ui->angleDoubleSpinBox->value(), offset_precise);
    }

    update_persistence_diagram();
    qDebug() << "update_dendrogram_diagram 2";
    if (dendrogram_arrangement_received)
    {
        qDebug() << "dendrogram_arrangement_received";
        update_dendrogram_diagram();
    }
}

void VisualizationWindow::on_normCoordCheckBox_clicked(bool checked)
{
    if (line_selection_ready) {
        slice_diagram.set_normalized_coords(checked);
        slice_diagram.resize_diagram();
        if (persistence_diagram_drawn)
            p_diagram.resize_diagram(slice_diagram.get_slice_length(), slice_diagram.get_pd_scale());
    }
}

void VisualizationWindow::on_barcodeCheckBox_clicked(bool checked)
{
    if (line_selection_ready)
        slice_diagram.toggle_barcode(checked);
}

void VisualizationWindow::on_xi0CheckBox_toggled(bool checked)
{
    if (line_selection_ready)
        slice_diagram.toggle_xi0_points(checked);
}

void VisualizationWindow::on_xi1CheckBox_toggled(bool checked)
{
    if (line_selection_ready)
        slice_diagram.toggle_xi1_points(checked);
}

void VisualizationWindow::on_xi2CheckBox_toggled(bool checked)
{
    if (line_selection_ready)
        slice_diagram.toggle_xi2_points(checked);
}

enum {is_forest_connector = -2}; // MAGIC NUMBER

void VisualizationWindow::project_template_dendrogram(Time_root& tr)
{
    int x,y;
    if (tr.r == is_forest_connector)
    {
        x = 0;
        y = 0;
    }
    else
    {
        x = tr.bigrade.first;
        y = tr.bigrade.second;
    }
    TemplatePoint pt(x,y,0,0,0);
    tr.t = project(pt, angle_precise, offset_precise);
    std::vector<std::shared_ptr<Time_root>> new_children;
    for (auto child : tr.children)
    {
        Time_root new_child = *child;
        std::shared_ptr<Time_root> new_child_ptr = std::make_shared<Time_root>(new_child);
        project_template_dendrogram(*new_child_ptr);
        new_children.push_back(new_child_ptr);
    }
    tr.children = new_children;
    return;
}

void VisualizationWindow::update_dendrogram_diagram()
{
    qDebug() << "******** UPDATE_DENDROGRAM_DIAGRAM *********";
    //projected_grades_of_appearance(angle_precise, offset_precise);
    qDebug() << "angle_precise = " << angle_precise;
    qDebug() << "offset_precise = " << offset_precise;
    const Time_root& tr = dendrogram_arrangement->get_dendrogram_template(angle_precise, offset_precise);
    qDebug() << "recursively_print_tr(tr) in update_dendrogram_diagram";
    Time_root::recursively_print_tr(tr);
    //std::unordered_map<std::pair<double,int>, Time_root> new_time_root_to_tr;
    //Time_root new_last_upper_tr;
    Time_root projected_tr = tr; //shallow copy hopefully
    project_template_dendrogram(projected_tr);
    qDebug() << "recursively_print_tr(projected_tr) in update_dendrogram_diagram";
    Time_root::recursively_print_tr(projected_tr, true);

    dendrogram_diagram.removeEventFilter(nodesEditor);
    qDebug() << "removeEventFilter success";
    dendrogram_diagram.update_dendrogram(projected_tr);
    qDebug() << "update_dendrogram success";




    // attempt at continuous updates
    /*QNEMainWindow* d_diagram = new QNEMainWindow(projected_tr,
                                                 &config_params,
                                                 this);*/

    /*QNodesEditor* nodesEditor = new QNodesEditor(this, dendrogram_diagram.get_block_x_scale(),
                                                 dendrogram_diagram.get_block_y_scale(),
                                                 &config_params);*/
    nodesEditor = new QNodesEditor(this, dendrogram_diagram.get_block_x_scale(),
                                   dendrogram_diagram.get_block_y_scale(),
                                   &config_params);

    nodesEditor->install(&dendrogram_diagram);

    /*
    ui->xposBox->setMinimum(-100.00); // MAGIC NUMBER
    dendrogram_diagram = new QNEMainWindow(projected_tr,
                                           &config_params,
                                           this);

    //dendrogram_diagram = new QNEMainWindow(graph, vertex_appearance, new_time_root_to_tr, new_last_upper_tr, &config_params, this);
    QObject::connect(dendrogram_diagram, &QNEMainWindow::set_cluster_label, this, &VisualizationWindow::set_cluster_label);
    QObject::connect(dendrogram_diagram, &QNEMainWindow::set_xpos_box, this, &VisualizationWindow::set_xpos_box);


    ui->dendrogramView->setScene(dendrogram_diagram);
    ui->dendrogramView->setRenderHint(QPainter::Antialiasing);

    QNodesEditor* nodesEditor = new QNodesEditor(this, dendrogram_diagram->get_block_x_scale(),
                                                 dendrogram_diagram->get_block_y_scale(),
                                                 &config_params);
    nodesEditor->install(dendrogram_diagram);*/
}

double VisualizationWindow::project(TemplatePoint& pt, double angle, double offset)
{
    if(angle == 0)  //then line is horizontal
    {
        if( grades.y[pt.y] <= offset)   //then point is below the line, so projection exists
            return grades.x[pt.x];
        else    //then no projection
            return INFTY;
    }
    else if(angle == 90)    //then line is vertical
    {
        if( grades.x[pt.x] <= -1*offset)   //then point is left of the line, so projection exists
            return grades.y[pt.y];
        else    //then no projection
            return INFTY;
    }
    //if we get here, then line is neither horizontal nor vertical
    double radians = angle*PI/180;
    double x = grades.x[pt.x];
    double y = grades.y[pt.y];

    if( y > x*tan(radians) + offset/cos(radians) )	//then point is above line
        return y/sin(radians) - offset/tan(radians); //project right

    return x/cos(radians) + offset*tan(radians); //project up
}//end project()

//updates the persistence diagram and barcode after a change in the slice line
void VisualizationWindow::update_persistence_diagram()
{
    if (persistence_diagram_drawn) {
        //get the barcode
        if (verbosity >= 4) {
            qDebug() << "  QUERY: angle =" << angle_precise << ", offset =" << offset_precise;
        }
        BarcodeTemplate dbc = arrangement->get_barcode_template(angle_precise, offset_precise);
        barcode = dbc.rescale(angle_precise, offset_precise, template_points->template_points, grades);

        //TESTING
        //qDebug() << "  XI SUPPORT VECTOR:";
        //for (unsigned i = 0; i < template_points->template_points.size(); i++) {
        //    TemplatePoint p = template_points->template_points[i];
        //    qDebug().nospace() << "    [" << i << "]: (" << p.x << "," << p.y << ") --> (" << grades.x[p.x] << "," << grades.y[p.y] << ")";
        //}
        if (verbosity >= 4) {
            dbc.print();
            barcode->print();
        }

        double zero_coord = rivet::numeric::project_zero(angle_precise, offset_precise, grades.x[0], grades.y[0]);

        //draw the barcode
        p_diagram.update_diagram(slice_diagram.get_slice_length(), slice_diagram.get_pd_scale(), zero_coord, *barcode);
        slice_diagram.update_barcode(*barcode, zero_coord, ui->barcodeCheckBox->isChecked());
    }
}

void VisualizationWindow::right_dot_released()
{
    //std::cout << "UPDATED DENDROGRAM_DIAGRAM" << std::endl;
    //update_dendrogram_diagram();
}

void VisualizationWindow::left_dot_released()
{
    //std::cout << "UPDATED DENDROGRAM_DIAGRAM" << std::endl;
    //update_dendrogram_diagram();
}

void VisualizationWindow::slice_line_released()
{
    //std::cout << "UPDATED DENDROGRAM_DIAGRAM" << std::endl;
    //update_dendrogram_diagram();
}

// sets xpos label for dendrogram
void VisualizationWindow::set_xpos_box(double xpos)
{
    //std::cout << " ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ in set_xpos_box ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
    //std::cout << "xpos = " << xpos << std::endl;

    slice_update_lock = true;

    ui->xposBox->setValue(xpos);

    slice_update_lock = false;
}

void VisualizationWindow::set_cluster_label(std::string cluster)
{
    slice_update_lock = true;

    ui->clusterText->setText(QString::fromStdString(cluster));

    slice_update_lock = false;
}

void VisualizationWindow::set_line_parameters(double angle, double offset)
{
    slice_update_lock = true;

    //correct for slight numerical errors that the interface might introduce
    if (angle < 0 && angle > -45)
        angle = 0;
    if (angle > 90 || angle < -40)
        angle = 90;

    //store values internally
    angle_precise = angle;
    offset_precise = offset;

    //update UI elements (values will be truncated)
    ui->angleDoubleSpinBox->setValue(angle);
    ui->offsetSpinBox->setValue(offset);

    slice_update_lock = false;

    update_persistence_diagram();
    qDebug() << "update_dendrogram_diagram 3";
    if (dendrogram_arrangement_received)
    {
        qDebug() << "dendrogram_arrangement_received";
        update_dendrogram_diagram();
    }
}

void VisualizationWindow::showEvent(QShowEvent* event)
{
    QMainWindow::showEvent(event);
    if (!data_selected) {
        ds_dialog.exec(); //show();
    }
}

void VisualizationWindow::resizeEvent(QResizeEvent* /*unused*/)
{
    if (line_selection_ready) {
        slice_diagram.resize_diagram();

        if (persistence_diagram_drawn)
            p_diagram.resize_diagram(slice_diagram.get_slice_length(), slice_diagram.get_pd_scale());
    }
}

void VisualizationWindow::closeEvent(QCloseEvent* event)
{
    if (unsaved_data) {
        QMessageBox::StandardButton reallyExit;
        reallyExit = QMessageBox::question(this, "Exit?", "Are you sure you want to exit? Your augmented arrangement has not been saved and will be lost!", QMessageBox::Yes | QMessageBox::No);

        if (reallyExit == QMessageBox::Yes) {
            event->accept();
            //        qDebug() << "User has closed RIVET.";
        } else
            event->ignore();
    } else {
        event->accept();
    }
}

void VisualizationWindow::on_actionExit_triggered()
{
    close();
}

void VisualizationWindow::on_actionAbout_triggered()
{
    aboutBox.show();
}

void VisualizationWindow::on_actionConfigure_triggered()
{
    configBox = new ConfigureDialog(config_params, input_params, this);
    configBox->exec();

    if (line_selection_ready) {
        slice_diagram.receive_parameter_change(QString::fromStdString(template_points->x_label),
            QString::fromStdString(template_points->y_label));

        if (persistence_diagram_drawn)
            p_diagram.receive_parameter_change();
    }

    delete configBox;
}

void VisualizationWindow::on_actionSave_persistence_diagram_as_image_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this, "Export persistence diagram as image",
        suggestedName("persist_offset_" + QString::number(offset_precise)
            + "_angle_" + QString::number(angle_precise) + ".png"),
        "PNG Image (*.png)");
    if (!fileName.isNull()) {
        QSettings settings;
        settings.setValue(DEFAULT_SAVE_DIR_KEY, QFileInfo(fileName).absolutePath());
        QPixmap pixMap = ui->pdView->grab();
        pixMap.save(fileName, "PNG");
    }
    ///TODO: error handling?
}

void VisualizationWindow::on_actionSave_line_selection_window_as_image_triggered()
{

    QString fileName = QFileDialog::getSaveFileName(this, "Export line selection window as image",
        suggestedName("line_offset_" + QString::number(offset_precise)
            + "_angle_" + QString::number(angle_precise) + ".png"),
        "PNG Image (*.png)");
    if (!fileName.isNull()) {
        QSettings settings;
        settings.setValue(DEFAULT_SAVE_DIR_KEY, QFileInfo(fileName).absolutePath());
        QPixmap pixMap = ui->sliceView->grab();
        pixMap.save(fileName, "PNG");
    }
    ///TODO: error handling?
}

void VisualizationWindow::on_actionSave_triggered()
{

    QString fileName = QFileDialog::getSaveFileName(this, "Save computed data", suggestedName("rivet"));
    if (!fileName.isNull()) {
        QSettings settings;
        settings.setValue(DEFAULT_SAVE_DIR_KEY, QFileInfo(fileName).absolutePath());
        save_arrangement(fileName);
    }
    ///TODO: error handling?
} //end on_actionSave_triggered()

void VisualizationWindow::save_arrangement(const QString& filename)
{
    try {
        write_boost_file(filename, input_params, *template_points, *arrangement);
    } catch (std::exception& e) {
        QMessageBox errorBox(QMessageBox::Warning, "Error",
            QString("Unable to write file: ").append(filename).append(": ").append(e.what()));
        errorBox.exec();
    }
} //end save_arrangement()

void VisualizationWindow::on_actionOpen_triggered()
{
    ///TODO: get user confirmation and clear the existing data structures

    QMessageBox msgBox;
    msgBox.setText("This feature is not implemented yet.");
    msgBox.exec();

    ///TODO: open the data select dialog box and load new data
} //end on_actionOpen_triggered()

QString VisualizationWindow::suggestedName(QString extension)
{
    QSettings settings;
    auto name = QString::fromStdString(input_params.fileName + ".H"
                    + std::to_string(input_params.dim)
                    + "_" + std::to_string(input_params.x_bins)
                    + "_" + std::to_string(input_params.y_bins)
                    + ".")
        + extension;
    auto suggested = QDir(settings.value(DEFAULT_SAVE_DIR_KEY).toString()).filePath(name);
    return suggested;
}
