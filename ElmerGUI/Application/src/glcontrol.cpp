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
 *  ElmerGUI glcontrol                                                       *
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

#include <QtGui>
#include <iostream>
#include "glcontrol.h"

using namespace std;

GLcontrol::GLcontrol(QWidget *parent)
  : QDialog(parent)
{
  ui.setupUi(this);

  connect(ui.okButton, SIGNAL(clicked()), 
	  this, SLOT(okButtonClicked()));

  connect(ui.closeButton, SIGNAL(clicked()), 
	  this, SLOT(close()));

  connect(ui.ambientDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.diffuseDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.specularDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.posxSpinBox, SIGNAL(editingFinished()),
	  this, SLOT(okButtonClicked()));

  connect(ui.posySpinBox, SIGNAL(editingFinished()),
	  this, SLOT(okButtonClicked()));

  connect(ui.poszSpinBox, SIGNAL(editingFinished()),
	  this, SLOT(okButtonClicked()));

  connect(ui.matAmbientDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.matDiffuseDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.matSpecularDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  connect(ui.matShininessDial, SIGNAL(sliderReleased()),
	  this, SLOT(okButtonClicked()));

  setWindowIcon(QIcon(":/icons/Mesh3D.png"));
}

GLcontrol::~GLcontrol()
{
}

void GLcontrol::okButtonClicked()
{
  GLfloat light_ambient[4];
  GLfloat light_diffuse[4];
  GLfloat light_specular[4];
  GLfloat light_position[4];
  GLfloat mat_ambient[4];
  GLfloat mat_diffuse[4];
  GLfloat mat_specular[4];
  GLfloat high_shininess[1];

  int ambientMaximum = ui.ambientDial->maximum();
  int diffuseMaximum = ui.diffuseDial->maximum();
  int specularMaximum = ui.specularDial->maximum();

  int ambientValue = ui.ambientDial->value();
  int diffuseValue = ui.diffuseDial->value();
  int specularValue = ui.specularDial->value();

  int posxValue = ui.posxSpinBox->value();
  int posyValue = ui.posySpinBox->value();
  int poszValue = ui.poszSpinBox->value();

  int matAmbientMaximum = ui.matAmbientDial->maximum();
  int matDiffuseMaximum = ui.matDiffuseDial->maximum();
  int matSpecularMaximum = ui.matSpecularDial->maximum();
  int matShininessMaximum = ui.matShininessDial->maximum();

  int matAmbientValue = ui.matAmbientDial->value();
  int matDiffuseValue = ui.matDiffuseDial->value();
  int matSpecularValue = ui.matSpecularDial->value();
  int matShininessValue = ui.matShininessDial->value();

  ambient = (GLfloat)ambientValue / (GLfloat)ambientMaximum;
  diffuse = (GLfloat)diffuseValue / (GLfloat)diffuseMaximum;
  specular = (GLfloat)specularValue / (GLfloat)specularMaximum;

  posx = (GLfloat)posxValue;
  posy = (GLfloat)posyValue;
  posz = (GLfloat)poszValue;

  matAmbient = (GLfloat)matAmbientValue / (GLfloat)matAmbientMaximum;
  matDiffuse = (GLfloat)matDiffuseValue / (GLfloat)matDiffuseMaximum;
  matSpecular = (GLfloat)matSpecularValue / (GLfloat)matSpecularMaximum;
  matShininess = (GLfloat)matShininessValue / (GLfloat)matShininessMaximum;
  matShininess *= 100.0;

  light_ambient[0] = ambient;
  light_ambient[1] = ambient;
  light_ambient[2] = ambient;
  light_ambient[3] = 1.0;

  light_diffuse[0] = diffuse;
  light_diffuse[1] = diffuse;
  light_diffuse[2] = diffuse;
  light_diffuse[3] = 1.0;

  light_specular[0] = specular;
  light_specular[1] = specular;
  light_specular[2] = specular;
  light_specular[3] = 1.0;

  light_position[0] = posx;
  light_position[1] = posy;
  light_position[2] = posz;
  light_position[3] = 0.0;

  mat_ambient[0] = matAmbient;
  mat_ambient[1] = matAmbient;
  mat_ambient[2] = matAmbient;
  mat_ambient[3] = 1.0;

  mat_diffuse[0] = matDiffuse;
  mat_diffuse[1] = matDiffuse;
  mat_diffuse[2] = matDiffuse;
  mat_diffuse[3] = 1.0;

  mat_specular[0] = matSpecular;
  mat_specular[1] = matSpecular;
  mat_specular[2] = matSpecular;
  mat_specular[3] = 1.0;

  high_shininess[0] = matShininess;

  glPushMatrix();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHT0);  

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glPopMatrix();

  glWidget->updateGL();
}
