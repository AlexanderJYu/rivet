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

//
// Created by Matthew L. Wright, 2014
//

#ifndef CONTROL_DOT_H
#define CONTROL_DOT_H

struct ConfigParameters;
class SliceLine;

#include <QGraphicsItem>
#include <QPainter>
#include <QRectF>
#include <QVariant>

class ControlDot : public QObject, public QGraphicsItem {
    Q_OBJECT
public:
    //ControlDot(SliceLine* line, bool left_bottom, ConfigParameters* params);
    ControlDot(SliceLine* line, bool left_bottom, ConfigParameters* params, bool isRightDendrogramDot = false,
               bool isLeftDendrogramDot = false);
    QRectF boundingRect() const;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*);
    QVariant itemChange(GraphicsItemChange change, const QVariant& value);

    void set_position(const QPointF& newpos);

signals:
    void right_dot_released();
    void left_dot_released();

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

private:
    SliceLine* slice_line;
    ConfigParameters* config_params;

    bool pressed;
    bool left_bottom; //TRUE if this is a left-bottom control dot, FALSE if this is a right-top control dot
    bool update_lock; //TRUE when the dot is being moved as result of external input; to avoid update loops
    bool isRightDendrogramDot; // TRUE if this is the right-top dot in the dendrogram window
    bool isLeftDendrogramDot; // TRUE if this is the left-bottom dot in the dendrogram window
};

#endif // CONTROL_DOT_H
