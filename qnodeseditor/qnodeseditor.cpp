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

#include "qnodeseditor.h"

#include <QGraphicsScene>
#include <QEvent>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsItem>
#include <QFontMetrics>
#include <QPainter>
//#include <QGraphicsPathItem>

#include "qneport.h"
#include "qneconnection.h"
#include "qneblock.h"
#include <iostream>

QNodesEditor::QNodesEditor(QObject *parent, double block_x_scale, double block_y_scale, ConfigParameters *params) :
    QObject(parent)
{
	conn = 0;
    this->block_x_scale = block_x_scale;
    this->block_y_scale = block_y_scale;
    this->config_params = params;
}

void QNodesEditor::install(QGraphicsScene *s)
{
	s->installEventFilter(this);
	scene = s;
}

QGraphicsItem* QNodesEditor::itemAt(const QPointF &pos)
{
	QList<QGraphicsItem*> items = scene->items(QRectF(pos - QPointF(1,1), QSize(3,3)));

	foreach(QGraphicsItem *item, items)
		if (item->type() > QGraphicsItem::UserType)
			return item;

	return 0;
}

bool QNodesEditor::eventFilter(QObject *o, QEvent *e)
{
	QGraphicsSceneMouseEvent *me = (QGraphicsSceneMouseEvent*) e;

	switch ((int) e->type())
	{
	case QEvent::GraphicsSceneMousePress:
	{

		switch ((int) me->button())
		{
		case Qt::LeftButton:
		{
			QGraphicsItem *item = itemAt(me->scenePos());
            if (item && item->type() == QNEPort::Type)
			{
                /*
                conn = new QNEConnection(0);
                scene->addItem(conn);
				conn->setPort1((QNEPort*) item);
				conn->setPos1(item->scenePos());
				conn->setPos2(me->scenePos());
                conn->updatePath();

                return true;*/
                item = item->parentItem();
            }
            if (item && item->type() == QNEBlock::Type)
			{
                QNEBlock* old_b = qgraphicsitem_cast<QNEBlock *>(item);
                QVector<QNEPort*> ports = old_b->ports();
                //QRectF* rect = qgraphicsitem_cast<QRectF *> b->path().elementAt(0);
                QString in_str, out_str;
                QVector<QNEConnection*> in_connections;
                QVector<QNEConnection*> out_connections;
                QFontMetrics fm(scene->font());
                int w;
                bool outportIsVisible;
                QNEPort* outport;
                QNEPort* inport;
                for (auto port: ports)
                {
                    if (port->isOutput())
                    {
                        out_str = port->portName();
                        out_connections = port->connections();
                        w = fm.width(port->portName());
                        outportIsVisible = port->getLabel()->isVisible();
                        outport = port;
                    }
                    else
                    {
                        in_str = port->portName();
                        in_connections = port->connections();
                        inport = port;
                    }
                }

                int xpos = old_b->scenePos().x();
                int ypos = old_b->scenePos().y();

                //QNEPort* in, out;

                if (/*item->isSelected()*/outportIsVisible)
                    //item->setVisible(true);
                {
                    QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
                    scene->addItem(b);
                    b->setPos(xpos,ypos);

                    QNEPort* in = (b->addPort(in_str,false));
                    QNEPort* out = (b->addPort(out_str,true));

                    out->setLabelVisibility(false);
                    in->setLabelVisibility(false);

                    // redraw edges
                    for (auto in_conn : in_connections)
                    {
                        QNEConnection* new_conn = new QNEConnection(0);
                        scene->addItem(new_conn);
                        new_conn->setPort1(in_conn->port1());
                        new_conn->setPort2(in);
                        new_conn->setPos1(in_conn->port1()->scenePos());
                        new_conn->setPos2(in->scenePos());

                        new_conn->updatePath();
                    }
                    //std::cout << "WIDTH RESIZE 5" << std::endl;
                    for (auto out_conn : out_connections)
                    {
                        //std::cout << "WIDTH RESIZE 6" << std::endl;
                        QNEConnection* new_conn = new QNEConnection(0);
                        //conn = new QNEConnection(0);
                        scene->addItem(new_conn);

                        //std::cout << "WIDTH RESIZE 7" << std::endl;

                        new_conn->setPort1(out);
                        //std::cout << "WIDTH RESIZE 7.25" << std::endl;
                        //std::cout << "stuff = " << out_conn->port2()->isOutput() << std::endl;
                        new_conn->setPort2(out_conn->port2());
                        //std::cout << "WIDTH RESIZE 7.5" << std::endl;
                        new_conn->setPos1(out->scenePos());
                        new_conn->setPos2(out_conn->port2()->scenePos());

                        //std::cout << "WIDTH RESIZE 8" << std::endl;

                        new_conn->updatePath();

                        //std::cout << "WIDTH RESIZE 9" << std::endl;
                    }
                    //std::cout << "WIDTH RESIZE 10" << std::endl;
                    delete item;
                }
                else
                {
                    //std::cout << "ELSE CASE" << std::endl;
                    outport->setLabelVisibility(true);
                    inport->setLabelVisibility(true);
                    /*
                    if (w > old_b->getWidth() - old_b->getHorzMargin())
                    {
                        QNEBlock *b = new QNEBlock(0, block_x_scale, block_y_scale, config_params);
                        scene->addItem(b);
                        b->setPos(xpos,ypos);
                        std::cout << "WIDTH RESIZE" << std::endl;

                        QNEPort* in = (b->addPort(in_str,false, true));
                        QNEPort* out = (b->addPort(out_str,true, true));
                        //out->setLabelVisibility(true);
                        //in->setLabelVisibility(true);

                        // redraw edges
                        for (auto in_conn : in_connections)
                        {
                            QNEConnection* new_conn = new QNEConnection(0);
                            scene->addItem(new_conn);
                            new_conn->setPort1(in_conn->port1());
                            new_conn->setPort2(in);
                            new_conn->setPos1(in_conn->port1()->scenePos());
                            new_conn->setPos2(in->scenePos());

                            new_conn->updatePath();
                        }
                        std::cout << "WIDTH RESIZE 5" << std::endl;
                        for (auto out_conn : out_connections)
                        {
                            std::cout << "WIDTH RESIZE 6" << std::endl;
                            QNEConnection* new_conn = new QNEConnection(0);
                            //conn = new QNEConnection(0);
                            scene->addItem(new_conn);

                            std::cout << "WIDTH RESIZE 7" << std::endl;

                            new_conn->setPort1(out);
                            std::cout << "WIDTH RESIZE 7.25" << std::endl;
                            std::cout << "stuff = " << out_conn->port2()->isOutput() << std::endl;
                            new_conn->setPort2(out_conn->port2());
                            std::cout << "WIDTH RESIZE 7.5" << std::endl;
                            new_conn->setPos1(out->scenePos());
                            new_conn->setPos2(out_conn->port2()->scenePos());

                            std::cout << "WIDTH RESIZE 8" << std::endl;

                            new_conn->updatePath();

                            std::cout << "WIDTH RESIZE 9" << std::endl;
                        }
                        std::cout << "WIDTH RESIZE 10" << std::endl;
                        delete item;
                    }
                    else
                    {
                        outport->setLabelVisibility(true);
                        inport->setLabelVisibility(true);
                    }*/

                }
			}
			break;
		}
        /*
		case Qt::RightButton:
		{
			QGraphicsItem *item = itemAt(me->scenePos());
			if (item && (item->type() == QNEConnection::Type || item->type() == QNEBlock::Type))
				delete item;
			// if (selBlock == (QNEBlock*) item)
				// selBlock = 0;
			break;
        }*/
		}
	}
	case QEvent::GraphicsSceneMouseMove:
	{
		if (conn)
		{
			conn->setPos2(me->scenePos());
			conn->updatePath();
			return true;
		}
		break;
	}
	case QEvent::GraphicsSceneMouseRelease:
	{
		if (conn && me->button() == Qt::LeftButton)
		{
			QGraphicsItem *item = itemAt(me->scenePos());
			if (item && item->type() == QNEPort::Type)
			{
				QNEPort *port1 = conn->port1();
				QNEPort *port2 = (QNEPort*) item;

				if (port1->block() != port2->block() && port1->isOutput() != port2->isOutput() && !port1->isConnected(port2))
				{
					conn->setPos2(port2->scenePos());
					conn->setPort2(port2);
					conn->updatePath();
					conn = 0;
					return true;
				}
			}

			delete conn;
			conn = 0;
			return true;
		}
		break;
	}
	}
	return QObject::eventFilter(o, e);
}

void QNodesEditor::save(QDataStream &ds)
{
	foreach(QGraphicsItem *item, scene->items())
		if (item->type() == QNEBlock::Type)
		{
			ds << item->type();
			((QNEBlock*) item)->save(ds);
		}

	foreach(QGraphicsItem *item, scene->items())
		if (item->type() == QNEConnection::Type)
		{
			ds << item->type();
			((QNEConnection*) item)->save(ds);
		}
}

void QNodesEditor::load(QDataStream &ds)
{
	scene->clear();

	QMap<quint64, QNEPort*> portMap;

	while (!ds.atEnd())
	{
		int type;
		ds >> type;
		if (type == QNEBlock::Type)
		{
            QNEBlock *block = new QNEBlock(0);
            scene->addItem(block);
			block->load(ds, portMap);
		} else if (type == QNEConnection::Type)
		{
            QNEConnection *conn = new QNEConnection(0);
            scene->addItem(conn);
			conn->load(ds, portMap);
		}
	}
}
