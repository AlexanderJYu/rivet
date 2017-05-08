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

#include "qnemainwindow.h"
#include "ui_qnemainwindow.h"

#include "qneblock.h"
#include "qnodeseditor.h"
#include "qneconnection.h"
//#include "qnodeseditor/dendrogram_data.h"
//#include "IntervalTree.h"
//#include "interface/control_dot.h"
//#include "interface/slice_line.h"


#include <QGraphicsScene>
#include <QGraphicsView>
#include <QFileDialog>

#include "qneport.h"
#include <sstream>
#include <string>
#include <cmath>    //c++ version of math.h; includes overloaded absolute value functions

#ifdef DEBUG_BUILD
    #define PRINTDEBUG(debugString) std::cout << (debugString) << std::endl;
#else
    #define PRINTDEBUG(debugString)
#endif

enum {is_forest_connector = -2,
      dendrogramViewWidth = 500,
      dendrogramViewHeight = 550}; // MAGIC NUMBER

void QNEMainWindow::draw_dendrogram(QGraphicsScene* scene,
                     /*DenseGRAPH<Edge>* graph *//* dont need this paramter probs ,*/
                     const Time_root& current_tr,
                     std::pair<int,int> low_high,
                     QNEPort* current_in,
                     double block_x_scale,
                     double block_y_scale,
                     std::vector<Interval<std::string, double>>& intervals)
{
    if (current_tr.children.size() == 0)
        return;

    QNEConnection *conn;
    //std::vector<Time_root> children = current_tr.children;
    //int top_of_partition = low_high.first;
    int bottom_of_partition = low_high.first;

    for (std::shared_ptr<Time_root> child_ptr : current_tr.children)
    {
        Time_root& child = *child_ptr;

        //int num_children = child.num_children;
        //double block_x_scale = 5.0 / (graph->V());
        //double block_y_scale = (-5.0) / (graph->V());
        int num_leaves = child.num_leaves;
        int xpos = child.t * 100 * block_x_scale;
        //int ypos = top_of_partition + (num_leaves * 100)/2;
        int ypos = bottom_of_partition + (num_leaves * 100 * (-1 * block_y_scale))/2;

        //std::cout << "NNNNNNNNNNNNNNNNNNNN" << std::endl;
        //std::cout << "low = " << low_high.first << std::endl;
        //std::cout << "high = " << low_high.second << std::endl;
        //std::cout << "bottom_of_partition = " << bottom_of_partition << std::endl;
        //std::cout << "num_leaves = " << num_leaves << std::endl;
        //std::cout << "(" << child.t << "," << child.r << ")" << std::endl;
        //std::cout << "xpos = " << xpos << std::endl;
        //std::cout << "ypos = " << ypos << std::endl;


        //std::pair<int,int> next_low_high(top_of_partition, top_of_partition + (num_leaves * 100));
        std::pair<int,int> next_low_high(bottom_of_partition, bottom_of_partition + (num_leaves * 100)* (-1 * block_y_scale));
        //top_of_partition += num_leaves * 100;
        bottom_of_partition += num_leaves * 100 * (-1 * block_y_scale);

        if (child.t != std::numeric_limits<double>::infinity())
        {
            QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
            scene->addItem(b);

            std::stringstream ss;
            ss << "(" << child.t << "," << child.r << ")";
            QString in_str =  QString::fromStdString(ss.str());

            std::stringstream ss_label;
            ss_label << "[" << time_root_to_complabel[std::pair<double,int>(child.t, child.r)];
            std::string complabel = ss_label.str();
            complabel.erase(complabel.size() - 1);
            complabel.append("]");

            //xpos_to_cluster[child.t].append(complabel);
            intervals.push_back(Interval<std::string,double>(child.t, current_tr.t, complabel));


            //std::cout << "child.comp_vertices.size() = " << child.comp_vertices.size() << std::endl;
            QString out_str = QString::fromStdString(child.birth_label_str);
            QNEPort* in = (b->addPort(in_str,false/*,true*/));
            QNEPort* out = (b->addPort(out_str,true));
            out->setLabelVisibility(false);
            in->setLabelVisibility(false);
            b->setPos(xpos,ypos);
            if (current_tr.r != is_forest_connector && current_tr.t != std::numeric_limits<double>::infinity())
            {
                conn = new QNEConnection(0);
                scene->addItem(conn);

                //QVector<QNEPort*> ports = vertex_blocks[v].ports();
                //QGraphicsItem *item_v = itemAt(vertex_pos[v].first,vertex_pos[v].second);
                conn->setPort1(out);
                conn->setPort2(current_in);
                conn->setPos1(out->scenePos());
                conn->setPos2(current_in->scenePos());
                conn->updatePath();
            }
            draw_dendrogram(scene, child, next_low_high, in, block_x_scale, block_y_scale, intervals);

        }
        else
        {
            QNEPort* in; // dummy variable to fulfill draw_dendrogram argument
            draw_dendrogram(scene, child, next_low_high, in, block_x_scale, block_y_scale, intervals);
        }
    }
}

//resizes diagram to fill the QGraphicsView
void QNEMainWindow::resize_diagram(/*double xmin, double xmax, double ymin, double ymax*/)
{
    /*
    //parameters
    int scene_padding = 30; //pixels
    int text_padding = 5;   //pixels*/

    //get dimensions of the QGraphicsView
    QList<QGraphicsView*> view_list = views();
    //std::cout << "view_list.size() = " << view_list.size() << std::endl;
    //int view_width = view_list[1]->width();
    //int view_height = view_list[1]->height();
    int view_width = dendrogramViewWidth; // MAGIC NUMBER
    int view_height = dendrogramViewHeight; // MAGIC NUMBER



    //determine scale
    //double left_text_width = std::max(data_ymin_text->boundingRect().width(), data_ymax_text->boundingRect().width());
    double diagram_max_width = view_width - padding; //- 2*scene_padding - text_padding - left_text_width;
    //double lower_text_height = std::max(data_xmin_text->boundingRect().height(), data_xmax_text->boundingRect().height());
    double diagram_max_height = view_height - padding; //- 2*scene_padding - text_padding - lower_text_height;
    /*
    data_xmin = xmin;
    data_xmax = xmax;
    data_ymin = ymin;
    data_ymax = ymax;*/

    if(data_xmax > data_xmin)
         scale_x = diagram_max_width/(data_xmax - data_xmin);
    else    //then there is only one x-grade
        scale_x = 1;                 ///IS THIS WHAT WE WANT???

    if(data_ymax > data_ymin)
        scale_y = diagram_max_height/(data_ymax - data_ymin);
    else    //then there is only one x-grade
        scale_y = 1;                 ///IS THIS WHAT WE WANT???
    /*
    if(!normalized_coords)  //then we want scale_x and scale_y to be the same (choose the smaller of the two)
    {
        if(scale_y < scale_x)
            scale_x = scale_y;
        else
            scale_y = scale_x;
    }*/

    //determine diagram size

    // FIGURE OUT EXACTLY WHAT FORMULAS I NEED HERE
    diagram_width = scale_x*(data_xmax - data_xmin);  //units: pixels
    diagram_height = scale_y*(data_ymax - data_ymin); //units: pixels

    //reposition reference objects
    control_rect->setRect(data_xmin, 0, diagram_width - 50, diagram_height - 50);
    /*gray_line_vertical->setLine(diagram_width, 0, diagram_width, diagram_height);
    gray_line_horizontal->setLine(0, diagram_height, diagram_width, diagram_height);

    data_xmin_text->setPos(data_xmin_text->boundingRect().width()/(-2), -1*text_padding);
    data_xmax_text->setPos(diagram_width - data_xmax_text->boundingRect().width()/2, -1*text_padding);
    data_ymin_text->setPos(-1*text_padding - data_ymin_text->boundingRect().width(), data_ymin_text->boundingRect().height()/2);
    data_ymax_text->setPos(-1*text_padding - data_ymax_text->boundingRect().width(), diagram_height + data_ymax_text->boundingRect().height()/2);

    x_label->setPos((diagram_width - x_label->boundingRect().width())/2, -1*text_padding);
    y_label->setPos(-1*text_padding - y_label->boundingRect().height(), (diagram_height - y_label->boundingRect().width())/2);*/


    //reposition dimension rectangles
    //redraw_dim_rects();

    //reposition xi points
    //redraw_dots();

    //reposition slice line
    slice_line->update_bounds(diagram_width, diagram_height - 50, 0/*padding*/, data_xmin);

    double x = 0, y = 0;
    x = data_xmin + (-1*line_pos*diagram_width);
    /*
    if(line_pos < 0)    //then left-bottom endpoint is along bottom edge of box
        x = -1*line_pos*diagram_width;
    else                //then left-bottom endpoint is along left edge of box
        y = line_pos*diagram_height;*/
    slice_line->update_position(x, y, line_vert, line_slope*scale_y/scale_x);

    //clear selection (because resizing window might combine or split dots in the upper strip of the persistence diagram)
    /*clear_selection();
    highlight_line->hide();

    //reposition highlighting
    if(primary_selected.size() > 0)
        update_highlight();*/

    // I MIGHT NEED THIS ACTUALLY
    //set scene rectangle (necessary to prevent auto-scrolling)
    double scene_rect_x = data_xmin; //-left_text_width - text_padding;
    double scene_rect_y = -25;//-lower_text_height - text_padding;
    double scene_rect_w = diagram_width; //+ padding; // + text_padding + left_text_width;
    double scene_rect_h = diagram_height; //+ padding; // + text_padding + lower_text_height;
    setSceneRect(scene_rect_x, scene_rect_y, scene_rect_w, scene_rect_h);
}//end resize_diagram()

//updates controls in the VisualizationWindow in response to a change in the line
void QNEMainWindow::update_window_controls()
{
    //refresh the scene to avoid artifacts from old lines, which otherwise can occur when the user moves the line quickly
    update(sceneRect());    //NOTE: this updates more items than necessary, but that is fine as long as it is fast

    //update SliceDiagram data values
    line_vert = slice_line->is_vertical();
    line_slope = slice_line->get_slope()*scale_x/scale_y;   //convert pixel units to data units
    line_pos = -1*(slice_line->pos().x() - data_xmin)/diagram_width;
    /* UNCOMMENT POSSIBLY
    if(slice_line->pos().x() > 0)
    {
        line_pos = (slice_line->pos().x() - data_xmin)/diagram_width;
    }
    else
    {
        line_pos = (slice_line->pos().x() - data_xmin)/diagram_width;
    }*/
    /* UNCOMMENT POSSIBLY
    if(slice_line->pos().x() > 0)
    {
        if(diagram_width > 0)
            line_pos = -1*slice_line->pos().x()/diagram_width;
        else    //can this ever happen?
            line_pos = 0;
    }
    else
    {
        if(diagram_height > 0)
            line_pos = slice_line->pos().y()/diagram_height;
        else
            line_pos = 0;
    }*/

    //update VisualizatoinWindow control objects
    //defaults for vertical line
    /*double angle = 90;
    double offset = -1*(slice_line->pos().x()/scale_x + data_xmin); //data units

    //handle non-vertical line
    if(!line_vert)
    {
        angle = atan(line_slope);   //radians

        double y_intercept = (slice_line->pos().y()/scale_y + data_ymin - line_slope * (slice_line->pos().x()/scale_x + data_xmin) );   //data units
        offset = cos(angle) * y_intercept;

        angle = angle*180/PI;   //convert to degrees
    }

    //send updates
    emit set_line_control_elements(angle, offset);*/
    double xpos = slice_line->pos().x()/(100 * block_x_scale);

    /*
    auto floor_xpos_iter = edge_wts.lower_bound(xpos);
    double floor_xpos;
    if (floor_xpos_iter != edge_wts.begin())
         floor_xpos = *(std::prev(edge_wts.lower_bound(xpos)));
    else
        floor_xpos = *(edge_wts.lower_bound(xpos));*/
    //double floor_xpos = *(std::prev(edge_wts.lower_bound(xpos)));

    emit set_xpos_box(xpos);

    //std::cout << "floor_xpos = " << floor_xpos << std::endl;
    //std::cout << "xpos_to_cluster[floor_xpos] = " << xpos_to_cluster[floor_xpos] << std::endl;
    //std::cout << "all elements of edge_wts:";
    //for (double x : edge_wts)
    //    std::cout << x << ",";

    std::vector<Interval<std::string, double>> intervals_containing_xpos;
    interval_tree.findOverlapping(xpos, xpos, intervals_containing_xpos);
    //std::cout << "interval_tree.intervals.size() = " << interval_tree.intervals.size() << std::endl;
    std::string cluster = "";
    //std::cout << "======= ALL INTERVALS IN TREE =======";
    /*for (Interval<std::string, double> interval : interval_tree.intervals)
    {
        std::cout << "interval.start = " << interval.start << std::endl;
        std::cout << "interval.stop = " << interval.stop << std::endl;
        std::cout << "interval.value = " << interval.value << std::endl;
        std::cout << "----------------------" << std::endl;
    }*/
    for (Interval<std::string, double> interval : intervals_containing_xpos)
    {
        //std::cout << "interval.start = " << interval.start << std::endl;
        //std::cout << "interval.stop = " << interval.stop << std::endl;
        //std::cout << "interval.value = " << interval.value << std::endl;
        //std::cout << "----------------------" << std::endl;
        cluster.append(interval.value);
    }

    //std::cout << "intervals_containing_xpos.size() = " << intervals_containing_xpos.size() << std::endl;
    //std::cout << "cluster = " << cluster << std::endl;
    //std::cout << "xpos = " << xpos << std::endl;
    emit set_cluster_label(cluster);

    //emit set_cluster_label(xpos_to_cluster[floor_xpos]);

    //since the line has changed, the highlighting is no longer valid
    //highlight_line->hide();
}//end update_window_controls()

std::string QNEMainWindow::populate_time_root_to_complabel(const Time_root& tr)
{
    std::stringstream ss;
    ss << tr.birth_label_str;
    if (tr.children.size() == 0)
    {
        time_root_to_complabel[std::pair<double,int>(tr.t,tr.r)] = ss.str();
        return ss.str();
    }
    for (auto child : tr.children)
    {
        ss << populate_time_root_to_complabel(*child);
    }
    time_root_to_complabel[std::pair<double,int>(tr.t,tr.r)] = ss.str();
    return ss.str();
}

// set data_xmax and data_xmin in the units of the data (to be converted later)
void QNEMainWindow::compute_data_parameters(const Time_root& tr,
                                            std::shared_ptr<int> n)
{
    *n = *n + 1;
    if (tr.t > finite_data_xmax && tr.t < std::numeric_limits<double>::infinity())
        finite_data_xmax = tr.t;
    if (tr.t > data_xmax)
        data_xmax = tr.t;
    if (tr.t < data_xmin)
        data_xmin = tr.t;
    for (auto child : tr.children)
        compute_data_parameters(*child, n);
    return;
}

void QNEMainWindow::update_dendrogram(const Time_root& tr)
{
    // reset everything and clear diagram
    this->clear();
    time_root_to_complabel.clear();
    populate_time_root_to_complabel(tr);
    std::shared_ptr<int> n = std::make_shared<int>(0);
    data_xmin = 0;
    data_xmax = 0;
    finite_data_xmax = 0;
    compute_data_parameters(tr, n);
    intervals.clear();

    //std::cout << "*********************** finished compute_data_paramters ***********************" << std::endl;
    //std::cout << "data_xmin = " << data_xmin << std::endl;
    //std::cout << "data_xmax = " << data_xmax << std::endl;
    //std::cout << "n = " << *n << std::endl;

    bool has_infinity_nodes = false;
    if (data_xmax == std::numeric_limits<double>::infinity())
        has_infinity_nodes = true;

    if (has_infinity_nodes)
        data_xmax = finite_data_xmax;

    block_x_scale = 3.0 / (data_xmax - data_xmin);
    block_y_scale = (-5.0) / ((double) *n);

    data_xmin = data_xmin * 100 * block_x_scale;
    data_xmax = data_xmax * 100 * block_x_scale;

    data_ymax = (*n) * 100;


    //pens and brushes
    QPen blackPen(Qt::black);
    blackPen.setWidth(2);

    control_rect = addRect(QRectF(), blackPen);  //0,0,diagram_width + padding,diagram_height + padding, blackPen);
    //add control objects
    line_vert = true;
    line_slope = (data_ymax - data_ymin)/(data_xmax - data_xmin);    //slope in data units
    line_pos = 0;
    //line_pos = data_xmin;   //start the line at the lower left corner of the box

    slice_line = new SliceLine(this, config_params);
    addItem(slice_line);

    dot_left = new ControlDot(slice_line, true, config_params, false, true);
    addItem(dot_left);

    dot_right = new ControlDot(slice_line, false, config_params, true);
    addItem(dot_right);

    slice_line->setDots(dot_left, dot_right);

    //data_ymax = 100 * graph->V();
    //data_xmax = 100 * mst.back()->wt;

    //data_ymax = 800;
    //data_xmax = 800;

    //block_x_scale = 3.0 / (mst.back()->wt);
    //block_y_scale = (-4.0) / (graph->V());

    // adjust dimensions of diagram
    resize_diagram();
    update_window_controls();

    if (tr.r != is_forest_connector && tr.t != std::numeric_limits<double>::infinity())
    {
        QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
        this->addItem(b);

        std::stringstream ss;
        ss << "(" << tr.t << "," << tr.r << ")";
        QString in_str =  QString::fromStdString(ss.str());
        /*ss.str("");
        ss << "[";
        for (int v : last_upper_tr.comp_vertices)
        {
            ss << v << ",";
        }
        std::string ss_str = ss.str();
        ss_str.erase(ss_str.size() - 1);
        ss_str.append("]");
        QString out_str = QString::fromStdString(ss_str);*/
        QString out_str = QString::fromStdString(tr.birth_label_str);
        //std::cout << "out_str = " << ss.str() << std::endl;

        // UNCOMMENT LATER xpos_to_cluster[last_upper_tr.t].append(ss_str);

        QNEPort* current_in = (b->addPort(in_str,false));
        QNEPort* current_out = (b->addPort(out_str,true));
        current_in->setLabelVisibility(false);
        current_out->setLabelVisibility(false);
        //current_in = (b->addPort(in_str,false));
        //current_out = (b->addPort(out_str,true));
        /*
        int total_leaves = 0;
        //total_leaves = time_root_to_tr[0].size();
        //std::cout << "******************* COMPUTE TOTAL_LEAVES ***********************" << std::endl;
        for (auto kv : time_root_to_tr)
        {
            Time_root tr = kv.second;
            //std::cout << "tr = (" << tr.t << "," << tr.r << ")" << std::endl;
            //std::cout << "tr.num_leaves = " << tr.num_leaves << std::endl;
            //std::cout << "tr.children.size() = " << tr.children.size() << std::endl;
            //std::cout << "************************" << std::endl;
            if (tr.children.size() == 0)
                total_leaves += 1;
        }

        //std::pair<int,int> low_high(0, graph->V() * 100);
        std::pair<int,int> low_high(0, total_leaves * 100 * -1 * block_y_scale);*/
        std::pair<int,int> low_high(0, tr.num_leaves * 100 * -1 * block_y_scale);
        int xpos = tr.t * 100 * block_x_scale;
        int ypos = (low_high.first + low_high.second)/2;
        b->setPos(xpos,ypos);

        //std::cout << "*^*^*^*^*^**^**^*^*^" << std::endl;
        //std::cout << "total_leaves = " << total_leaves << std::endl;
        //std::cout << "xpos_to_cluster[0] = "
                  //<< xpos_to_cluster[0] << std::endl;
        //std::cout << "xpos_to_cluster[1] = "
                  //<< xpos_to_cluster[1] << std::endl;
        //std::cout << "xpos_to_cluster[2] = "
                  //<< xpos_to_cluster[2] << std::endl;

        //std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM" << std::endl;
        //std::cout << "block_x_scale = " << block_x_scale << std::endl;
        //std::cout << "block_y_scale = " << block_y_scale << std::endl;
        draw_dendrogram(this, /*graph,*/ tr, low_high, current_in, block_x_scale, block_y_scale, intervals);
    }
    else
    {
        std::pair<int,int> low_high(0, tr.num_leaves * 100 * -1 * block_y_scale);
        QNEPort* current_in; // dummy variable to fulfill draw_dendrogram argument
        draw_dendrogram(this, /*graph,*/ tr, low_high, current_in, block_x_scale, block_y_scale, intervals);
    }
    interval_tree = IntervalTree<std::string, double>(intervals);
}

QNEMainWindow::QNEMainWindow(/*const Time_root& tr,*/
                             ConfigParameters* params,
                             QObject *parent) :
    QGraphicsScene(parent), dot_left(), dot_right(), slice_line(), config_params(params), PI(3.14159265358979323846),
    padding(30), data_xmin(0), data_xmax(0), data_ymin(0), data_ymax(1000), time_root_to_complabel()
{
    // Add slice line to diagram
    /*Time_root nonconst_tr = tr;
    std::cout << "QNEMainWindow constructor" << std::endl;
    std::cout << ">>>>>>>>>>>>>>> tr before populate_num_leaves: " << std::endl;
    Time_root::recursively_print_tr(tr);
    populate_num_leaves(tr);
    std::cout << ">>>>>>>>>>>>>>> tr after populate_num_leaves: " << std::endl;
    Time_root::recursively_print_tr(tr);*/

    /*
    populate_time_root_to_complabel(tr);
    std::shared_ptr<int> n = std::make_shared<int>(0);
    compute_data_parameters(tr, n);

    std::cout << "*********************** finished compute_data_paramters ***********************" << std::endl;
    std::cout << "data_xmin = " << data_xmin << std::endl;
    std::cout << "data_xmax = " << data_xmax << std::endl;
    std::cout << "n = " << *n << std::endl;

    block_x_scale = 3.0 / (data_xmax - data_xmin);
    block_y_scale = (-5.0) / ((double) *n);

    data_xmin = data_xmin * 100 * block_x_scale;
    data_xmax = data_xmax * 100 * block_x_scale;

    data_ymax = (*n) * 100;


    //pens and brushes
    QPen blackPen(Qt::black);
    blackPen.setWidth(2);

    control_rect = addRect(QRectF(), blackPen);  //0,0,diagram_width + padding,diagram_height + padding, blackPen);
    //add control objects
    line_vert = true;
    line_slope = (data_ymax - data_ymin)/(data_xmax - data_xmin);    //slope in data units
    line_pos = 0;
    //line_pos = data_xmin;   //start the line at the lower left corner of the box

    slice_line = new SliceLine(this, config_params);
    addItem(slice_line);

    dot_left = new ControlDot(slice_line, true, config_params, false, true);
    addItem(dot_left);

    dot_right = new ControlDot(slice_line, false, config_params, true);
    addItem(dot_right);

    slice_line->setDots(dot_left, dot_right);

    //data_ymax = 100 * graph->V();
    //data_xmax = 100 * mst.back()->wt;

    //data_ymax = 800;
    //data_xmax = 800;

    //block_x_scale = 3.0 / (mst.back()->wt);
    //block_y_scale = (-4.0) / (graph->V());

    // adjust dimensions of diagram
    resize_diagram();
    update_window_controls();

    if (tr.r != is_forest_connector)
    {
        QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
        this->addItem(b);

        std::stringstream ss;
        ss << "(" << tr.t << "," << tr.r << ")";
        QString in_str =  QString::fromStdString(ss.str());

        QString out_str = QString::fromStdString(tr.birth_label_str);
        //std::cout << "out_str = " << ss.str() << std::endl;

        // UNCOMMENT LATER xpos_to_cluster[last_upper_tr.t].append(ss_str);

        QNEPort* current_in = (b->addPort(in_str,false));
        QNEPort* current_out = (b->addPort(out_str,true));
        current_in->setLabelVisibility(false);
        current_out->setLabelVisibility(false);
        //current_in = (b->addPort(in_str,false));
        //current_out = (b->addPort(out_str,true));

        std::pair<int,int> low_high(0, tr.num_leaves * 100 * -1 * block_y_scale);
        int xpos = tr.t * 100 * block_x_scale;
        int ypos = (low_high.first + low_high.second)/2;
        b->setPos(xpos,ypos);

        //std::cout << "*^*^*^*^*^**^**^*^*^" << std::endl;
        //std::cout << "total_leaves = " << total_leaves << std::endl;
        //std::cout << "xpos_to_cluster[0] = "
                  //<< xpos_to_cluster[0] << std::endl;
        //std::cout << "xpos_to_cluster[1] = "
                  //<< xpos_to_cluster[1] << std::endl;
        //std::cout << "xpos_to_cluster[2] = "
                  //<< xpos_to_cluster[2] << std::endl;

        //std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM" << std::endl;
        //std::cout << "block_x_scale = " << block_x_scale << std::endl;
        //std::cout << "block_y_scale = " << block_y_scale << std::endl;
        draw_dendrogram(this, tr, low_high, current_in, block_x_scale, block_y_scale, intervals);
    }
    else
    {
        std::pair<int,int> low_high(0, tr.num_leaves * 100 * -1 * block_y_scale);
        QNEPort* current_in; // dummy variable to fulfill draw_dendrogram argument
        draw_dendrogram(this, tr, low_high, current_in, block_x_scale, block_y_scale, intervals);
    }
    interval_tree = IntervalTree<std::string, double>(intervals);*/
}
/*
QNEMainWindow::QNEMainWindow(DenseGRAPH<Edge>* graph,
                             boost::unordered::unordered_map<unsigned int, double> vertex_appearance,
                             std::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                             Time_root& last_upper_tr,
                             ConfigParameters* params,
                             QObject *parent) :
    //QMainWindow(parent)
    QGraphicsScene(parent), dot_left(), dot_right(), slice_line(), padding(30), PI(3.14159265358979323846),
    config_params(params), data_xmin(0), data_xmax(1000), data_ymin(0), data_ymax(1000), xpos_to_cluster()

{
    //scene = new QGraphicsScene();


    if (graph->V() == 0)
    {
        graph = new DenseGRAPH<Edge>( 9, true );
        //std::vector<int> vertex_appearance(9);
        boost::unordered::unordered_map<unsigned int, double> vertex_appearance;
        vertex_appearance[0] = 0;
        vertex_appearance[1] = 0;
        vertex_appearance[2] = 0;
        vertex_appearance[3] = 0;
        vertex_appearance[4] = 1;
        vertex_appearance[5] = 1;

        vertex_appearance[6] = 0;
        vertex_appearance[7] = 1;
        vertex_appearance[8] = 0;
        //vertex_appearance[6] = 0;
        //vertex_appearance[7] = 1;
        //vertex_appearance[8] = 1;
        //vertex_appearance[9] = 0;

        graph->insert( EdgePtr( new Edge( 0, 1, 1.0 ) ) );
        graph->insert( EdgePtr( new Edge( 2, 3, 0.0 ) ) );
        graph->insert( EdgePtr( new Edge( 4, 5, 2.0 ) ) );
        graph->insert( EdgePtr( new Edge( 0, 2, 2.0 ) ) );
        graph->insert( EdgePtr( new Edge( 0, 4, 3.0 ) ) );

        graph->insert( EdgePtr( new Edge( 8, 6, 1.0 ) ) );
        graph->insert( EdgePtr( new Edge( 7, 6, 1.0 ) ) );

        graph->insert( EdgePtr( new Edge( 6, 0, 2.0 ) ) );

    }
    // Add slice line to diagram

    //pens and brushes
    QPen blackPen(Qt::black);
    blackPen.setWidth(2);

    control_rect = addRect(QRectF(), blackPen);  //0,0,diagram_width + padding,diagram_height + padding, blackPen);
    //add control objects
    line_vert = true;  //IS IT POSSIBLE THAT THE INITIAL LINE COULD BE VERTICAL???????????????????????????????????
    line_slope = (data_ymax - data_ymin)/(data_xmax - data_xmin);    //slope in data units
    line_pos = 0;   //start the line at the lower left corner of the box

    slice_line = new SliceLine(this, config_params);
    addItem(slice_line);

    dot_left = new ControlDot(slice_line, true, config_params, false, true);
    addItem(dot_left);

    dot_right = new ControlDot(slice_line, false, config_params, true);
    addItem(dot_right);

    slice_line->setDots(dot_left, dot_right);

    // Add blocks to diagram

    std::vector<EdgePtr> mst;
    //std::unordered_map<std::pair<int,int>, Time_root> time_root_to_tr;
    //Time_root last_upper_tr;

    Dendrogram_data::KruskalMST( graph, mst, time_root_to_tr, last_upper_tr, vertex_appearance);

    // add vertex_appearance times to edge_wts, which suggests its a misnomer

    for (auto kv : vertex_appearance)
        edge_wts.insert(kv.second);


    //std::cout << "GGGGGGGGGGGGGGGGGGGGGGGGGGGG" << std::endl;

    for (auto e : mst)
    {
        edge_wts.insert(e->wt);
        //std::cout << "e->wt = " << e->wt << std::endl;
    }

    data_xmax = 100 * graph->V();
    data_ymax = 100 * mst.back()->wt;

    block_x_scale = 3.0 / (mst.back()->wt);
    block_y_scale = (-4.0) / (graph->V());

    // adjust dimensions of diagram
    resize_diagram();
    update_window_controls();

    QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
    this->addItem(b);

    std::stringstream ss;
    ss << "(" << last_upper_tr.t << "," << last_upper_tr.r << ")";
    QString in_str =  QString::fromStdString(ss.str());
    ss.str("");
    ss << "[";
    for (int v : last_upper_tr.comp_vertices)
    {
        ss << v << ",";
    }
    std::string ss_str = ss.str();
    ss_str.erase(ss_str.size() - 1);
    ss_str.append("]");
    QString out_str = QString::fromStdString(ss_str);
    //std::cout << "out_str = " << ss.str() << std::endl;

    xpos_to_cluster[last_upper_tr.t].append(ss_str);

    QNEPort* current_in = (b->addPort(in_str,false));
    QNEPort* current_out = (b->addPort(out_str,true));
    current_in->setLabelVisibility(false);
    current_out->setLabelVisibility(false);
    //current_in = (b->addPort(in_str,false));
    //current_out = (b->addPort(out_str,true));
    int total_leaves = 0;
    //total_leaves = time_root_to_tr[0].size();
    //std::cout << "******************* COMPUTE TOTAL_LEAVES ***********************" << std::endl;
    for (auto kv : time_root_to_tr)
    {
        Time_root tr = kv.second;
        //std::cout << "tr = (" << tr.t << "," << tr.r << ")" << std::endl;
        //std::cout << "tr.num_leaves = " << tr.num_leaves << std::endl;
        //std::cout << "tr.children.size() = " << tr.children.size() << std::endl;
        //std::cout << "************************" << std::endl;
        if (tr.children.size() == 0)
            total_leaves += 1;
    }

    //std::pair<int,int> low_high(0, graph->V() * 100);
    std::pair<int,int> low_high(0, total_leaves * 100 * -1 * block_y_scale);
    int xpos = last_upper_tr.t * 100 * block_x_scale;
    int ypos = (low_high.first + low_high.second)/2;
    b->setPos(xpos,ypos);

    //std::cout << "*^*^*^*^*^**^**^*^*^" << std::endl;
    //std::cout << "total_leaves = " << total_leaves << std::endl;
    //std::cout << "xpos_to_cluster[0] = "
              //<< xpos_to_cluster[0] << std::endl;
    //std::cout << "xpos_to_cluster[1] = "
              //<< xpos_to_cluster[1] << std::endl;
    //std::cout << "xpos_to_cluster[2] = "
              //<< xpos_to_cluster[2] << std::endl;

    //std::cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM" << std::endl;
    //std::cout << "block_x_scale = " << block_x_scale << std::endl;
    //std::cout << "block_y_scale = " << block_y_scale << std::endl;

    draw_dendrogram(this, graph, last_upper_tr, low_high, current_in, block_x_scale, block_y_scale, intervals);
    interval_tree = IntervalTree<std::string, double>(intervals);

}*/

QNEMainWindow::~QNEMainWindow()
{

}

void QNEMainWindow::saveFile()
{
	QString fname = QFileDialog::getSaveFileName();
	if (fname.isEmpty())
		return;

	QFile f(fname);
	f.open(QFile::WriteOnly);
	QDataStream ds(&f);
	nodesEditor->save(ds);
}

void QNEMainWindow::loadFile()
{
	QString fname = QFileDialog::getOpenFileName();
	if (fname.isEmpty())
		return;

	QFile f(fname);
	f.open(QFile::ReadOnly);
	QDataStream ds(&f);
	nodesEditor->load(ds);
}

void QNEMainWindow::addBlock()
{
    QNEBlock *b = new QNEBlock(0);

    scene->addItem(b);
	static const char* names[] = {"Vin", "Voutsadfasdf", "Imin", "Imax", "mul", "add", "sub", "div", "Conv", "FFT"};
	for (int i = 0; i < 4 + rand() % 3; i++)
	{
		b->addPort(names[rand() % 10], rand() % 2, 0, 0);
        b->setPos(view->sceneRect().center().toPoint());
	}
}
