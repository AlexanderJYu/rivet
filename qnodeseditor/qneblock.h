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

#ifndef QNEBLOCK_H
#define QNEBLOCK_H

#include <QGraphicsPathItem>
#include "interface/config_parameters.h"

class QNEPort;

class QNEBlock : public QGraphicsPathItem
{
public:
	enum { Type = QGraphicsItem::UserType + 3 };

    QNEBlock(QGraphicsItem *parent = 0, double block_x_scale = 1, double block_y_scale = -1, ConfigParameters *params = 0);

    QNEPort* addPort(const QString &name, bool isOutput, bool shouldResize = false, int flags = 0, int ptr = 0);
	void addInputPort(const QString &name);
	void addOutputPort(const QString &name);
	void addInputPorts(const QStringList &names);
	void addOutputPorts(const QStringList &names);
    void setWidth(qreal w);
    void setHeight(qreal h);
    int getWidth() {return width; }
    int getHeight() {return height; }
    int getHorzMargin() {return horzMargin; }
	void save(QDataStream&);
	void load(QDataStream&, QMap<quint64, QNEPort*> &portMap);
	void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

	QNEBlock* clone();
	QVector<QNEPort*> ports();

	int type() const { return Type; }

protected:
	QVariant itemChange(GraphicsItemChange change, const QVariant &value);

private:
    double horzMargin;
    double vertMargin;
    double width;
    double height;
    double block_x_scale;
    double block_y_scale;
    ConfigParameters* config_params;
};

#endif // QNEBLOCK_H
