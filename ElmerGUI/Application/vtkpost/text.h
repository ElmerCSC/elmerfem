/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland    *
 *                                                                           *
 *  This program is free software; you can redistribute it and/or            *
 *  modify it under the terms of the GNU General Public License              *
 *  as published by the Free Software Foundation; either version 2           *
 *  of the License, or (at your option) any later version.                   *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful,          *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU General Public License for more details.                             *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License        *
 *  along with this program (in file fem/GPL-2); if not, write to the        *
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,         *
 *  Boston, MA 02110-1301, USA.                                              *
 *                                                                           *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 *  ElmerGUI text                                                            *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly, Juha Ruokolainen and Peter Råback                   *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                 *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 15 Mar 2008                                               *
 *                                                                           *
 *****************************************************************************/

#ifndef TEXT_H
#define TEXT_H

#include <QWidget>
#include "ui_text.h"

class VtkPost;

class Text : public QDialog
{
  Q_OBJECT

public:
  Text(QWidget *parent = 0);
  ~Text();

  Ui::textDialog ui;

  void draw(VtkPost*);

public slots:
  void SetMessage(QString);
  void SetPosX(int);
  void SetPosY(int);
  void SetLeft();
  void SetCentered();
  void SetRight();
  void SetSize(int);
  void SetBold(bool);
  void SetItalic(bool);
  void SetShadow(bool);
  void SetRed(double);
  void SetGreen(double);
  void SetBlue(double);
  void SetRGB(double, double, double);

signals:
  void drawTextSignal();
  void hideTextSignal();

private slots:
  void applyButtonClicked();
  void cancelButtonClicked();
  void okButtonClicked();

private:
  double red;
  double green;
  double blue;

};

#endif // TEXT_H
