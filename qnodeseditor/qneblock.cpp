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

#include "qneblock.h"

#include <QPen>
#include <QGraphicsScene>
#include <QFontMetrics>
#include <QPainter>

#include "qneport.h"
#include <iostream>

QNEBlock::QNEBlock(QGraphicsItem *parent, double block_x_scale, double block_y_scale, ConfigParameters *params) : QGraphicsPathItem(parent)
{
	QPainterPath p;

    //p.addRoundedRect(-50, -15, 0, 100, 5, 5);
	setPath(p);
    //setPen(QPen(Qt::darkGreen));
    setPen(params->xi0color);
    //setPen(params->xi0color);
    //setBrush(Qt::green);
    setBrush(params->xi0color);
    //setFlag(QGraphicsItem::ItemIsMovable);
	setFlag(QGraphicsItem::ItemIsSelectable);
    horzMargin = 10;
    vertMargin = 10;
    width = horzMargin / block_x_scale;
    height = vertMargin / block_y_scale;
    this->block_x_scale = block_x_scale;
    this->block_y_scale = block_y_scale;
    config_params = params;

    QTransform transform;                    // scale
    transform.scale(block_x_scale, block_y_scale);
    //std::cout << "========= SCALES===========" << std::endl;
    //std::cout << "block_x_scale = " << block_x_scale << std::endl;
    //std::cout << "block_y_scale = " << block_y_scale << std::endl;
    this->setTransform(transform);
    this->update();


}

QNEPort* QNEBlock::addPort(const QString &name, bool isOutput, bool shouldResize, int flags, int ptr)
{
    QNEPort *port = new QNEPort(this, block_x_scale, block_y_scale, config_params);
	port->setName(name);
	port->setIsOutput(isOutput);
	port->setNEBlock(this);
	port->setPortFlags(flags);
	port->setPtr(ptr);


	QFontMetrics fm(scene()->font());
	int w = fm.width(name);
	int h = fm.height();
	// port->setPos(0, height + h/2);
    /*if (shouldResize && w > width - horzMargin)
        width = w + horzMargin;*/
    //height += h;

	QPainterPath p;
    //std::cout << "======= DIMENSIONS =========" << std::endl;
    //std::cout << "width = " << width << std::endl;
    //std::cout << "height = " << height << std::endl;
    p.addRoundedRect(-width/2, -(height)/2, width, height, 5, 5);
	setPath(p);

    int y = -(height) / 2 /*+ vertMargin + port->radius()*/;
    foreach(QGraphicsItem *port_, childItems()) {
		if (port_->type() != QNEPort::Type)
			continue;

		QNEPort *port = (QNEPort*) port_;
		if (port->isOutput())
            port->setPos(/*width/2 + port->radius()*/0, y + height/2 /*+ 2*port->radius()*/);
		else
            port->setPos(/*-width/2 - port->radius()*/0, y + height/2 /*+ 2*port->radius()*/);
        //y += h;
	}

	return port;
}

/*void QNEBlock::updatePath()
{
    QPainterPath p;

    //QPointF pos1(m_port1->scenePos());
    //QPointF pos2(m_port2->scenePos());

    p.moveTo(pos1);

    qreal dx = pos2.x() - pos1.x();
    qreal dy = pos2.y() - pos1.y();

    QPointF ctr1(pos1.x() + dx * 0.25, pos1.y() + dy * 0.1);
    QPointF ctr2(pos1.x() + dx * 0.75, pos1.y() + dy * 0.9);

    p.cubicTo(ctr1, ctr2, pos2);

    setPath(p);
}*/

void QNEBlock::setWidth(qreal w)
{
    width = w;
}

void QNEBlock::setHeight(qreal h)
{
    height = h;
}

void QNEBlock::addInputPort(const QString &name)
{
	addPort(name, false);
}

void QNEBlock::addOutputPort(const QString &name)
{
	addPort(name, true);
}

void QNEBlock::addInputPorts(const QStringList &names)
{
	foreach(QString n, names)
		addInputPort(n);
}

void QNEBlock::addOutputPorts(const QStringList &names)
{
	foreach(QString n, names)
		addOutputPort(n);
}

void QNEBlock::save(QDataStream &ds)
{
	ds << pos();

	int count(0);

    foreach(QGraphicsItem *port_, childItems())
	{
		if (port_->type() != QNEPort::Type)
			continue;

		count++;
	}

	ds << count;

    foreach(QGraphicsItem *port_, childItems())
	{
		if (port_->type() != QNEPort::Type)
			continue;

		QNEPort *port = (QNEPort*) port_;
		ds << (quint64) port;
		ds << port->portName();
		ds << port->isOutput();
		ds << port->portFlags();
	}
}

void QNEBlock::load(QDataStream &ds, QMap<quint64, QNEPort*> &portMap)
{
	QPointF p;
	ds >> p;
	setPos(p);
	int count;
	ds >> count;
	for (int i = 0; i < count; i++)
	{
		QString name;
		bool output;
		int flags;
		quint64 ptr;

		ds >> ptr;
		ds >> name;
		ds >> output;
		ds >> flags;
		portMap[ptr] = addPort(name, output, flags, ptr);
	}
}

#include <QStyleOptionGraphicsItem>

void QNEBlock::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
	Q_UNUSED(option)
	Q_UNUSED(widget)

	if (isSelected()) {
        painter->setPen(QPen(Qt::black));
        painter->setBrush(config_params->xi0color);
	} else {
        //painter->setPen(QPen(Qt::darkGreen));
        //painter->setBrush(Qt::green);
        painter->setPen(config_params->xi0color);
        painter->setBrush(config_params->xi0color);
	}

	painter->drawPath(path());
}

QNEBlock* QNEBlock::clone()
{
    QNEBlock *b = new QNEBlock(0);
    this->scene()->addItem(b);

	foreach(QGraphicsItem *port_, childItems())
	{
		if (port_->type() == QNEPort::Type)
		{
			QNEPort *port = (QNEPort*) port_;
			b->addPort(port->portName(), port->isOutput(), port->portFlags(), port->ptr());
		}
	}

	return b;
}

QVector<QNEPort*> QNEBlock::ports()
{
	QVector<QNEPort*> res;
	foreach(QGraphicsItem *port_, childItems())
	{
		if (port_->type() == QNEPort::Type)
			res.append((QNEPort*) port_);
	}
	return res;
}

QVariant QNEBlock::itemChange(GraphicsItemChange change, const QVariant &value)
{

    Q_UNUSED(change);

	return value;
}

