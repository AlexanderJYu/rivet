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

#include <QGraphicsScene>
#include <QFileDialog>

#include "qneport.h"
#include <sstream>
#include <string>

void draw_dendrogram(QGraphicsScene* scene,
                     DenseGRAPH<Edge>* graph,
                     /*std::unordered_map<std::pair<int,int>, Time_root>& time_root_to_tr,*/
                     boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                     Time_root& current_tr,
                     std::pair<int,int> low_high,
                     QNEPort* current_in)
{
    if (current_tr.t == 0)
        return;

    QNEConnection *conn;
    std::vector<Time_root> children = current_tr.children;
    int top_of_partition = low_high.first;

    for (Time_root child : children)
    {
        int num_children = child.num_children;
        int xpos = child.t * 100;
        int ypos = top_of_partition + (num_children * 100)/2;
        std::pair<int,int> next_low_high(top_of_partition, top_of_partition + (num_children * 100));
        top_of_partition += num_children * 100;
        QNEBlock *b = new QNEBlock(0);
        scene->addItem(b);

        std::stringstream ss;
        ss << "(" << child.t << "," << child.r << ")";
        QString in_str =  QString::fromStdString(ss.str());

        ss.str("");
        //for (int v : child.comp_vertices)
        for (int v : child.birth_label)
        {
            ss << v << ",";
        }
        QString out_str = QString::fromStdString(ss.str());
        //std::cout << "out_str = " << ss.str() << std::endl;
        //std::cout << "child.t = " << child.t << " , " << "child.r = " << child.r << std::endl;
        //std::cout << "child.comp_vertices.size() = " << child.comp_vertices.size() << std::endl;

        QNEPort* in = (b->addPort(in_str,false));
        QNEPort* out = (b->addPort(out_str,true));
        b->setPos(xpos,ypos);

        //int v = it->get()->v;
        //int w = it->get()->w;
        conn = new QNEConnection(0);
        scene->addItem(conn);

        //QVector<QNEPort*> ports = vertex_blocks[v].ports();
        //QGraphicsItem *item_v = itemAt(vertex_pos[v].first,vertex_pos[v].second);
        conn->setPort1(out);
        conn->setPort2(current_in);
        conn->setPos1(out->scenePos());
        conn->setPos2(current_in->scenePos());
        conn->updatePath();

        draw_dendrogram(scene, graph, time_root_to_tr, child, next_low_high, in);
    }
}

QNEMainWindow::QNEMainWindow(DenseGRAPH<Edge>* graph,
                             /*std::unordered_map<std::pair<int,int>, Time_root>& time_root_to_tr,*/
                             boost::unordered::unordered_map<std::pair<double,int>, Time_root>& time_root_to_tr,
                             Time_root& last_upper_tr,
                             QWidget *parent) :
    QMainWindow(parent)
{
    scene = new QGraphicsScene();

    /*QAction *quitAct = new QAction(tr("&Quit"), this);
    quitAct->setShortcuts(QKeySequence::Quit);
    quitAct->setStatusTip(tr("Quit the application"));
    connect(quitAct, SIGNAL(triggered()), qApp, SLOT(quit()));

    QAction *loadAct = new QAction(tr("&Load"), this);
    loadAct->setShortcuts(QKeySequence::Open);
    loadAct->setStatusTip(tr("Open a file"));
    connect(loadAct, SIGNAL(triggered()), this, SLOT(loadFile()));

    QAction *saveAct = new QAction(tr("&Save"), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save a file"));
    connect(saveAct, SIGNAL(triggered()), this, SLOT(saveFile()));

    QAction *addAct = new QAction(tr("&Add"), this);
    addAct->setStatusTip(tr("Add a block"));
    connect(addAct, SIGNAL(triggered()), this, SLOT(addBlock()));

    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(addAct);
    fileMenu->addAction(loadAct);
    fileMenu->addAction(saveAct);
    fileMenu->addSeparator();
    fileMenu->addAction(quitAct);*/

    setWindowTitle(tr("Node Editor"));


    QDockWidget *dock = new QDockWidget(tr("Nodes"), this);
    dock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
    view = new QGraphicsView(dock);
    view->setScene(scene);

    view->setRenderHint(QPainter::Antialiasing, true);

    dock->setWidget(view);
    addDockWidget(Qt::LeftDockWidgetArea, dock);

    nodesEditor = new QNodesEditor(this);
    nodesEditor->install(scene);

    QNEBlock *b = new QNEBlock(0);
    scene->addItem(b);

    std::stringstream ss;
    ss << "(" << last_upper_tr.t << "," << last_upper_tr.r << ")";
    QString in_str =  QString::fromStdString(ss.str());

    ss.str("");
    for (int v : last_upper_tr.birth_label)
    {
        ss << v << ",";
    }
    QString out_str = QString::fromStdString(ss.str());
    std::cout << "out_str = " << ss.str() << std::endl;

    QNEPort* current_in = (b->addPort(in_str,false));
    QNEPort* current_out = (b->addPort(out_str,true));
    std::pair<int,int> low_high(0, graph->V() * 100);
    int xpos = last_upper_tr.t * 100;
    int ypos = (low_high.first + low_high.second)/2;
    b->setPos(xpos,ypos);

    draw_dendrogram(scene, graph, time_root_to_tr, last_upper_tr, low_high, current_in);
    /*
    int vertex_offset = 0;
    std::vector<std::pair<int,int>> vertex_pos (graph->V());
    std::vector<QNEBlock> vertex_blocks;
    std::vector<std::pair<QNEPort*,QNEPort*>> vertex_ports (graph->V()); //first is in, second is out

    for(int u = 0; u < graph->V(); ++u)
    {
        QNEBlock *b = new QNEBlock(0);
        scene->addItem(b);
        //b->addPort((QString) u,0,QNEPort::NamePort);
        QNEPort* in = (b->addPort("in",false));
        QNEPort* out = (b->addPort("out",true));
        b->setPos(vertex_offset,vertex_offset);
        vertex_pos[u] = std::pair<int,int>(vertex_offset,vertex_offset);
        vertex_offset += 100;
        //vertex_blocks[u] = b;
        std::pair<QNEPort*,QNEPort*> in_out (in,out);
        vertex_ports[u] = in_out;

    }*/

    //QNEConnection *conn;


    /*
    for ( std::vector<EdgePtr>::const_iterator it = mst.begin(); it != mst.end(); ++it )
    {
        int v = it->get()->v;
        int w = it->get()->w;
        conn = new QNEConnection(0);
        scene->addItem(conn);

        //QVector<QNEPort*> ports = vertex_blocks[v].ports();
        //QGraphicsItem *item_v = itemAt(vertex_pos[v].first,vertex_pos[v].second);
        conn->setPort1(vertex_ports[v].first);
        conn->setPort2(vertex_ports[w].second);
        conn->setPos1(vertex_ports[v].first->scenePos());
        conn->setPos2(vertex_ports[w].second->scenePos());
        conn->updatePath();

        //conn->setPort1((QNEPort*) item);
        //conn->setPos1(item->scenePos());
        //conn->setPos2(me->scenePos());
        //conn->updatePath();
    }

    for ( std::vector<EdgePtr>::const_iterator it = mst.begin(); it != mst.end(); ++it )
    {
        std::cout << it->get()->v << " "
                  << it->get()->w << " "
                  << it->get()->wt << std::endl;
    }*/

    /*
    QNEBlock *b = new QNEBlock(0);
    scene->addItem(b);
    b->addPort("test", 0, QNEPort::NamePort);
    b->addPort("TestBlock", 0, QNEPort::TypePort);
    b->addInputPort("in1");
    b->addInputPort("in2");
    b->addInputPort("in3");
    b->addOutputPort("out1");
    b->addOutputPort("out2");
    b->addOutputPort("out3");

    b = b->clone();
    b->setPos(150, 0);

    b = b->clone();
    b->setPos(150, 150);*/
}

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
