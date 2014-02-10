/*****************************************************************************
 *                                                                           *
 *  Elmer, A Finite Element Software for Multiphysical Problems              *
 *                                                                           *
 *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland   *
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
 *  ElmerClips                                                               *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Mikko Lyly                                                      *
 *  Email:   Juha.Ruokolainen@csc.fi                                         *
 *  Web:     http://www.csc.fi/elmer                                         *
 *  Address: CSC - IT Center for Science Ltd.                                *
 *           Keilaranta 14                                                   *
 *           02101 Espoo, Finland                                            *
 *                                                                           *
 *  Original Date: 14 Nov 2010                                               *
 *                                                                           *
 *****************************************************************************/
#ifndef ENCODER_H
#define ENCODER_H

#include <QtGui>

#ifndef INT64_C
#define INT64_C(c) (c ## LL)
#define UINT64_C(c) (c ## ULL)
#endif

extern "C" {
#include <libavcodec/avcodec.h>
#include <libswscale/swscale.h>
}

class Encoder : public QThread
{
 Q_OBJECT

 public:
  Encoder(QObject *parent = 0);
  void setUrls(const QList<QUrl> &urls);
  void setResolutions(const QList<int> &resolutions);
  void setQuality(int quality);
  void run();

 signals:
  void information(const QString &fileName);
  void progress(int value);

 private:
  void findImages(const QStringList &list);
  bool isImage(const QFileInfo &info) const;
  void sortImages();
  void compressImages(int targetWidth);
  bool convertToYUV(const QImage &image, int widthYUV, int heightYUV);

  QList<QUrl> urls;
  QList<int> resolutions;
  int quality;
  QStringList imageFileList;
  AVFrame *frameYUV;
  AVFrame *frameRGB;
  QByteArray bufferYUV;
  QByteArray bufferMPG;
  int totalFrames;
  int currentFrame;
};

#endif // ENCODER_H
