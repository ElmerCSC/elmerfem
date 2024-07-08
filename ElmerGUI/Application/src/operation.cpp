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
 *  ElmerGUI operation_t                                                     *
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

#include <iostream>
#include "operation.h"

using namespace std;

operation_t::operation_t()
{
  next = 0;
  type = 0;
  angle = 0.0;
  selected = 0;
  select_set = 0;
}

operation_t::~operation_t()
{
}

void operation_t::appendToProject(QDomDocument *projectDoc, QDomElement *ops)
{
  operation_t *p = this->next;
  for(int index = 0; p; p = p->next, index++) {
    QDomElement op = projectDoc->createElement("operation");
    op.setAttribute("index", QString::number(index));
    ops->appendChild(op);
    
    QDomElement type = projectDoc->createElement("type");
    QDomText typeValue = projectDoc->createTextNode(QString::number(p->type));
    type.appendChild(typeValue);
    op.appendChild(type);

    QDomElement angle = projectDoc->createElement("angle");
    QDomText angleValue = projectDoc->createTextNode(QString::number(p->angle));
    angle.appendChild(angleValue);
    op.appendChild(angle);

    QDomElement selected = projectDoc->createElement("selected");
    selected.setAttribute("lists", QString::number(p->selected));
    op.appendChild(selected);

    for(int list = 0; list < p->selected; list++) {
      QDomElement selection = projectDoc->createElement("list");
      QDomText selectionValue = projectDoc->createTextNode(QString::number(p->select_set[list]));
      selection.appendChild(selectionValue);
      selected.appendChild(selection);
    }    
  }
}


int operation_t::readFromProject(QDomDocument *projectDoc, QDomElement *ops)
{
  operation_t *p = this->next;
  operation_t *q = NULL;

  while(p != NULL) {
    if(p->select_set != NULL)
      delete [] p->select_set;
    q = p->next;
    if(p != NULL)
      delete p;
    p = q;
  }

  int operations = 0;
  q = this;
  q->next = NULL;

  QDomElement op = ops->firstChildElement("operation");

  for( ; !op.isNull(); op = op.nextSiblingElement()) {
    operations++;

    p = new operation_t;
    p->next = NULL;
    q->next = p;
    q = p;

    p->type = op.firstChildElement("type").text().toInt();
    p->angle = op.firstChildElement("angle").text().toDouble();

    QDomElement selected = op.firstChildElement("selected");
    p->selected = selected.attribute("lists").toInt();

    p->select_set = new int[p->selected];

    QDomElement selection = selected.firstChildElement("list");

    for(int list = 0; !selection.isNull(); selection = selection.nextSiblingElement(), list++) {
      if(list >= p->selected) {
	cout << "Project loader: load operations: index out of bounds" << endl;
	return 0;
      }
      p->select_set[list] = selection.text().toInt();
    }
  }

  return operations;
}
