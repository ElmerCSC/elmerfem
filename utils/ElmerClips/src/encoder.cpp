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
#include "encoder.h"

Encoder::Encoder(QObject *parent) : QThread(parent), quality(2)
{
  frameRGB = avcodec_alloc_frame();
  frameYUV = avcodec_alloc_frame();
  bufferMPG.resize(1000000);
}

void Encoder::setUrls(const QList<QUrl> &urls)
{
  this->urls = urls;
}

void Encoder::setResolutions(const QList<int> &resolutions)
{
  this->resolutions = resolutions;
}

void Encoder::setQuality(int quality)
{
  this->quality = quality;
}

void Encoder::run()
{
  QStringList fileNameList;
  
  foreach(const QUrl &url, urls)
    fileNameList << url.toLocalFile();
  
  findImages(fileNameList);
  
  totalFrames = imageFileList.count() * resolutions.count();

  currentFrame = 0;

  emit progress(0);

  foreach(int resolution, resolutions)
    compressImages(resolution);

  emit information("DONE");

  emit progress(0);
}

void Encoder::findImages(const QStringList &fileNameList)
{
  imageFileList.clear();
  
  foreach(const QString &fileName, fileNameList) {
    QFileInfo fileInfo(fileName);
    
    if(fileInfo.isFile() && isImage(fileInfo))
      imageFileList << fileName;
    
    if(fileInfo.isDir()) {
      QDirIterator iterator(fileName, QDirIterator::Subdirectories);
      
      while(iterator.hasNext()) {
	QString current = iterator.next();
	
	if(isImage(iterator.fileInfo()))
	  imageFileList << current;
      }
    }
  }

  if(imageFileList.isEmpty()) {
    emit information("ERROR: No image files found");
    return;
  }

  sortImages();
}

bool Encoder::isImage(const QFileInfo &info) const
{
  QString suffix = info.suffix().toLower();

  if((suffix == "png") || 
     (suffix == "jpg") || 
     (suffix == "jpeg") ||
     (suffix == "tiff") ||
     (suffix == "gif"))
    return true;
  
  return false;
}

void Encoder::sortImages()
{
  QRegExp rx("[0-9]{1,10}");

  QMap<int, QString> tmpMap;

  foreach(const QString &imageFile, imageFileList) {
    QFileInfo info(imageFile);

    QString fileName = info.fileName();

    int index = rx.indexIn(fileName);

    if(index < 0) {
      emit information("ERROR: Unable to sort images");
      imageFileList.clear();
      return;
    }

    int order = rx.cap().toInt();

    tmpMap[order] = imageFile;
  }

  imageFileList = tmpMap.values();
}

void Encoder::compressImages(int targetWidth)
{
  if(imageFileList.isEmpty())
    return;

  if((!frameRGB) || (!frameYUV)) {
    emit information("ERROR: Out of resources");
    return;
  }

  // Load first image to determine the frame size:
  //-----------------------------------------------
  int widthYUV = targetWidth;

  widthYUV -= widthYUV % 2;

  QImage tmpImage(imageFileList[0]);

  tmpImage = tmpImage.scaledToWidth(widthYUV);

  int heightYUV = tmpImage.height();

  heightYUV -= heightYUV % 2;

  int pixels = widthYUV * heightYUV;

  if(pixels < 1) {
    emit information("ERROR: Illegal image size");
    return;
  }

  // Initialize avcodec:
  //---------------------
  avcodec_init();

  avcodec_register_all();

  // Init encoder:
  //---------------
  CodecID codec_id = CODEC_ID_MPEG1VIDEO;

  AVCodec *codec = avcodec_find_encoder(codec_id);

  if(!codec) {
    emit information("ERROR: Unable to initialize encoder");
    return;
  }

  // Init context:
  //---------------
  AVCodecContext *context = avcodec_alloc_context();

  if(!context) {
    emit information("ERROR: Unable to initialize encoder");
    return;
  }

  context->codec_id = codec_id;
  context->codec_type = CODEC_TYPE_VIDEO;
  context->qmin = 2;
  context->qmax = qMax(2, qMin(31, quality));
  context->width = widthYUV;
  context->height = heightYUV;
  context->time_base.num = 1;
  context->time_base.den = 25;
  context->gop_size = 15;
  context->max_b_frames = 2;
  context->pix_fmt = PIX_FMT_YUV420P;

  if(avcodec_open(context, codec) < 0) {
    emit information("ERROR: Unable to initialize encoder");
    avcodec_close(context);
    return;
  }

  // Open the output file:
  //-----------------------
  QFileInfo info(imageFileList[0]);

  QString fileName = info.absolutePath()
    + "/" + info.dir().dirName()
    + "_" + QString::number(widthYUV)
    + "x" + QString::number(heightYUV)
    + ".mpg";

  emit information("FILE: " + fileName);

  QFile file(fileName);

  if(!file.open(QFile::WriteOnly)) {
    emit information("ERROR: Unable to open output file");
    avcodec_close(context);
    return;
  }

  // Compress frames:
  //------------------
  int bytes = 0;

  foreach(const QString &imageFile, imageFileList) {
    ++currentFrame;

    emit progress(100 * double(currentFrame) / totalFrames);

    emit information(imageFile);

    QImage image(imageFile);

    if(!convertToYUV(image, widthYUV, heightYUV)) {
      emit information("ERROR: Unable to scale image");
      continue;
    }
    
    bytes = avcodec_encode_video(context, (uint8_t *)bufferMPG.data(),
				 bufferMPG.size(), frameYUV);

    if(file.write(bufferMPG, bytes) != bytes) {
      emit information("ERROR: Unable to write to output file");
      file.close();
      avcodec_close(context);
      return;
    }
  }

  // Get the delayed frames:
  //-------------------------
  for( ; bytes; ) {
    bytes = avcodec_encode_video(context, (uint8_t *)bufferMPG.data(),
				 bufferMPG.size(), NULL);

    if(file.write(bufferMPG, bytes) != bytes) {
      emit information("ERROR: Unable to write to output file");
      file.close();
      avcodec_close(context);
      return;
    }
  }

  // MPEG1 end code:
  //-----------------
  bufferMPG[0] = (char)0x00;
  bufferMPG[1] = (char)0x00;
  bufferMPG[2] = (char)0x01;
  bufferMPG[3] = (char)0xb7;

  if(file.write(bufferMPG, 4) != 4)
      emit information("ERROR: Unable to write to output file");

  // Done:
  //-------
  file.close();
  avcodec_close(context);
}

bool Encoder::convertToYUV(const QImage &image,
			   int widthYUV, int heightYUV)
{
  // Allocate scaler context:
  //--------------------------
  int widthRGB = image.width();
  int heightRGB = image.height();

  SwsContext *context = sws_getContext(widthRGB, heightRGB, PIX_FMT_RGB24,
				       widthYUV, heightYUV, PIX_FMT_YUV420P,
				       SWS_LANCZOS, NULL, NULL, NULL);
  if(!context) {
    emit information("ERROR: Out of resources");
    return false;
  }

  // Prepare the RGB frame:
  //------------------------
  QImage imageRGB = image.convertToFormat(QImage::Format_RGB888);

  frameRGB->data[0] = imageRGB.bits();
  frameRGB->data[1] = imageRGB.bits();
  frameRGB->data[2] = imageRGB.bits();

  frameRGB->linesize[0] = imageRGB.bytesPerLine();
  frameRGB->linesize[1] = imageRGB.bytesPerLine();
  frameRGB->linesize[2] = imageRGB.bytesPerLine();

  // Prepare the YUV frame:
  //------------------------
  int pixelsYUV = widthYUV * heightYUV;
  int bytesYUV = 3 * pixelsYUV / 2;

  if(bufferYUV.size() < bytesYUV)
    bufferYUV.resize(bytesYUV);

  frameYUV->data[0] = (uint8_t *)bufferYUV.data();
  frameYUV->data[1] = (uint8_t *)bufferYUV.data() + pixelsYUV;
  frameYUV->data[2] = (uint8_t *)bufferYUV.data() + pixelsYUV + pixelsYUV / 4;

  frameYUV->linesize[0] = widthYUV;
  frameYUV->linesize[1] = widthYUV / 2;
  frameYUV->linesize[2] = widthYUV / 2;
  
  // Scale && transform color space:
  //---------------------------------
  sws_scale(context, frameRGB->data, frameRGB->linesize,
	    0, heightRGB, frameYUV->data, frameYUV->linesize);

  av_free(context);

  return true;
}
