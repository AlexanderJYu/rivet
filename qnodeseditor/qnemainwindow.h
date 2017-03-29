/* Copyright (c) 2012, STANISLAW ADASZEWSKI
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of STANISLAW ADASZEWSKI nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL STANISLAW ADASZEWSKI BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */

#ifndef QNEMAINWINDOW_H
#define QNEMAINWINDOW_H

#include <QMainWindow>
#include <QtWidgets>
#include "Graph.h"
#include "time_root.h"
#include "qnodeseditor/qneblock.h"
#include "interface/control_dot.h"
#include "interface/slice_line.h"
#include <set>
#include <boost/unordered_map.hpp>
#include "IntervalTree.h"

class QNodesEditor;
/*
class QNEMainWindow : public QMainWindow
{
    Q_OBJECT*/

class QNEMainWindow : public QGraphicsScene
{
    Q_OBJECT

public:
    /*explicit QNEMainWindow(DenseGRAPH<Edge>* graph,
                           boost::unordered::unordered_map<unsigned int, double> vertex_appearance,
                           std::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                           Time_root& last_upper_tr,
                           ConfigParameters* params,
                           QObject *parent = 0);*/

    explicit QNEMainWindow(const Time_root& tr,
                           ConfigParameters* params,
                           QObject *parent = 0);
	~QNEMainWindow();

    void update_line(double angle, double offset);  //updates the line, in response to a change in the controls in the VisualizationWindow
    void update_window_controls();   //computes new angle and offset in response to a change in the line, emits signal for the VisualizationWindow

    void resize_diagram(/*double xmin, double xmax, double ymin, double ymax*/);

    void draw_dendrogram(QGraphicsScene* scene,
                    /*DenseGRAPH<Edge>* graph,*/
                    /*std::unordered_map<std::pair<int,int>, Time_root>& time_root_to_tr,*/
                    const Time_root& current_tr,
                    std::pair<int,int> low_high,
                    QNEPort* current_in,
                    double block_x_scale,
                    double block_y_scale,
                    std::vector<Interval<std::string, double>>& intervals);

    //int populate_num_leaves(Time_root& tr);
    std::string populate_time_root_to_complabel(const Time_root& tr);

    // update_dendrogram
    void update_dendrogram(const Time_root& tr);


    // compute data_xmin and data_xmax
    void compute_data_parameters(const Time_root& tr, std::shared_ptr<int> n);

    double get_block_x_scale() { return block_x_scale; }
    double get_block_y_scale() { return block_y_scale ; }

    //double get_slice_length();  //gets the length of the slice, for scaling the persistence diagram
    //double get_pd_scale();      //gets the number of pixels per unit, for the persistence diagram
    //double get_zero();          //gets the coordinate on the slice line which we consider "zero" for the persistence diagram

private slots:
	void saveFile();
	void loadFile();
	void addBlock();

signals:
    //void set_line_control_elements(double angle, double offset);    //sends updates to, e.g., the VisualizationWindow
    void set_xpos_box(double xpos); // sends updates to Visualization Window to display x position
    void set_cluster_label(std::string cluster); //sends updates to Visualization Window to display cluster
private:

    std::set<double> edge_wts;

	QNodesEditor *nodesEditor;
    QMenu *fileMenu;
    QGraphicsView *view;
    QGraphicsScene *scene;

    //parameters
    ConfigParameters* config_params;

    QGraphicsRectItem* control_rect;            //control dots live on this rectangle

    ControlDot* dot_left;
    ControlDot* dot_right;

    SliceLine* slice_line;
    QGraphicsLineItem* highlight_line;

    const double PI;   //used in get_pd_scale() when the slice line is vertical

    double scale_x, scale_y;            //x- and y-scales for drawing data points
    double block_x_scale, block_y_scale;   //x- and y- scales for drawing QNEblocks
    const int padding;  //distance between xi support point area and control rectangle (on the top and right sides)
    double data_xmin, data_xmax, data_ymin, data_ymax;  //min and max coordinates of the data
    int diagram_width, diagram_height;  //pixel size of the diagram
    double line_slope;   //slope of the slice line in data units
    bool line_vert;      //true if the line is vertical, false otherwise
    double line_pos;     //relative position of left endpoint of line: 0 is lower left corner, positive values (up to 1) are along left side, negative values (to -1) are along bottom edge of box
    //std::unordered_map<double, std::string> xpos_to_cluster;
    boost::unordered::unordered_map<std::pair<double,int>, std::string> time_root_to_complabel;
    std::vector<Interval<std::string, double>> intervals;
    IntervalTree<std::string, double> interval_tree;
    //QNEBlock *b;
    //QNEPort* current_in;
    //QNEPort* current_out;
};

#endif // QNEMAINWINDOW_H
