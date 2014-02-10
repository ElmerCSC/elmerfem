
#include <QBitmap>
#include <QPainterPath>
#include <QVariant>
#include <QObject>
/****************************************************************************
** Meta object code from reading C++ file 'qbitarray.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qbitarray.h"
class PythonQtQBitArrayWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QBitArray() { return QBitArray(); }
QVariant new_QBitArray(int arg0,bool arg1) { return QBitArray(arg0,arg1); }
QVariant new_QBitArray(int arg0) { return QBitArray(arg0); }
QVariant new_QBitArray(const QBitArray & arg0) { return QBitArray(arg0); }
int size(QBitArray* obj)  const  {return obj->size(); }
int count(QBitArray* obj)  const  {return obj->count(); }
int count(QBitArray* obj,bool arg0)  const  {return obj->count(arg0); }
bool isEmpty(QBitArray* obj)  const  {return obj->isEmpty(); }
bool isNull(QBitArray* obj)  const  {return obj->isNull(); }
void resize(QBitArray* obj,int arg0)  {obj->resize(arg0); }
void detach(QBitArray* obj)  {obj->detach(); }
bool isDetached(QBitArray* obj)  const  {return obj->isDetached(); }
void clear(QBitArray* obj)  {obj->clear(); }
bool testBit(QBitArray* obj,int arg0)  const  {return obj->testBit(arg0); }
void setBit(QBitArray* obj,int arg0)  {obj->setBit(arg0); }
void setBit(QBitArray* obj,int arg0,bool arg1)  {obj->setBit(arg0,arg1); }
void clearBit(QBitArray* obj,int arg0)  {obj->clearBit(arg0); }
bool toggleBit(QBitArray* obj,int arg0)  {return obj->toggleBit(arg0); }
bool at(QBitArray* obj,int arg0)  const  {return obj->at(arg0); }
bool fill(QBitArray* obj,bool arg0,int arg1)  {return obj->fill(arg0,arg1); }
bool fill(QBitArray* obj,bool arg0)  {return obj->fill(arg0); }
void fill(QBitArray* obj,bool arg0,int arg1,int arg2)  {obj->fill(arg0,arg1,arg2); }
void truncate(QBitArray* obj,int arg0)  {obj->truncate(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qdatetime.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qdatetime.h"
class PythonQtQDateWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QDate() { return QDate(); }
QVariant new_QDate(int arg0,int arg1,int arg2) { return QDate(arg0,arg1,arg2); }
bool isNull(QDate* obj)  const  {return obj->isNull(); }
bool isValid(QDate* obj)  const  {return obj->isValid(); }
int year(QDate* obj)  const  {return obj->year(); }
int month(QDate* obj)  const  {return obj->month(); }
int day(QDate* obj)  const  {return obj->day(); }
int dayOfWeek(QDate* obj)  const  {return obj->dayOfWeek(); }
int dayOfYear(QDate* obj)  const  {return obj->dayOfYear(); }
int daysInMonth(QDate* obj)  const  {return obj->daysInMonth(); }
int daysInYear(QDate* obj)  const  {return obj->daysInYear(); }
int weekNumber(QDate* obj,int * arg0)  const  {return obj->weekNumber(arg0); }
int weekNumber(QDate* obj)  const  {return obj->weekNumber(); }
QString static_QDate_shortMonthName(int arg0)  {return QDate::shortMonthName(arg0); }
QString static_QDate_shortDayName(int arg0)  {return QDate::shortDayName(arg0); }
QString static_QDate_longMonthName(int arg0)  {return QDate::longMonthName(arg0); }
QString static_QDate_longDayName(int arg0)  {return QDate::longDayName(arg0); }
QString toString(QDate* obj,Qt::DateFormat arg0)  const  {return obj->toString(arg0); }
QString toString(QDate* obj)  const  {return obj->toString(); }
QString toString(QDate* obj,const QString & arg0)  const  {return obj->toString(arg0); }
bool setYMD(QDate* obj,int arg0,int arg1,int arg2)  {return obj->setYMD(arg0,arg1,arg2); }
bool setDate(QDate* obj,int arg0,int arg1,int arg2)  {return obj->setDate(arg0,arg1,arg2); }
QDate addDays(QDate* obj,int arg0)  const  {return obj->addDays(arg0); }
QDate addMonths(QDate* obj,int arg0)  const  {return obj->addMonths(arg0); }
QDate addYears(QDate* obj,int arg0)  const  {return obj->addYears(arg0); }
int daysTo(QDate* obj,const QDate & arg0)  const  {return obj->daysTo(arg0); }
QDate static_QDate_currentDate()  {return QDate::currentDate(); }
QDate static_QDate_fromString(const QString & arg0,Qt::DateFormat arg1)  {return QDate::fromString(arg0,arg1); }
QDate static_QDate_fromString(const QString & arg0)  {return QDate::fromString(arg0); }
QDate static_QDate_fromString(const QString & arg0,const QString & arg1)  {return QDate::fromString(arg0,arg1); }
bool static_QDate_isValid(int arg0,int arg1,int arg2)  {return QDate::isValid(arg0,arg1,arg2); }
bool static_QDate_isLeapYear(int arg0)  {return QDate::isLeapYear(arg0); }
uint static_QDate_gregorianToJulian(int arg0,int arg1,int arg2)  {return QDate::gregorianToJulian(arg0,arg1,arg2); }
void static_QDate_julianToGregorian(uint arg0,int & arg1,int & arg2,int & arg3)  {QDate::julianToGregorian(arg0,arg1,arg2,arg3); }
QDate fromJulianDay(QDate* obj,int arg0)  {return obj->fromJulianDay(arg0); }
int toJulianDay(QDate* obj)  const  {return obj->toJulianDay(); }

};

class PythonQtQTimeWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QTime() { return QTime(); }
bool isNull(QTime* obj)  const  {return obj->isNull(); }
bool isValid(QTime* obj)  const  {return obj->isValid(); }
int hour(QTime* obj)  const  {return obj->hour(); }
int minute(QTime* obj)  const  {return obj->minute(); }
int second(QTime* obj)  const  {return obj->second(); }
int msec(QTime* obj)  const  {return obj->msec(); }
QString toString(QTime* obj,Qt::DateFormat arg0)  const  {return obj->toString(arg0); }
QString toString(QTime* obj)  const  {return obj->toString(); }
QString toString(QTime* obj,const QString & arg0)  const  {return obj->toString(arg0); }
bool setHMS(QTime* obj,int arg0,int arg1,int arg2,int arg3)  {return obj->setHMS(arg0,arg1,arg2,arg3); }
bool setHMS(QTime* obj,int arg0,int arg1,int arg2)  {return obj->setHMS(arg0,arg1,arg2); }
QTime addSecs(QTime* obj,int arg0)  const  {return obj->addSecs(arg0); }
int secsTo(QTime* obj,const QTime & arg0)  const  {return obj->secsTo(arg0); }
QTime addMSecs(QTime* obj,int arg0)  const  {return obj->addMSecs(arg0); }
int msecsTo(QTime* obj,const QTime & arg0)  const  {return obj->msecsTo(arg0); }
QTime static_QTime_currentTime()  {return QTime::currentTime(); }
QTime static_QTime_fromString(const QString & arg0,Qt::DateFormat arg1)  {return QTime::fromString(arg0,arg1); }
QTime static_QTime_fromString(const QString & arg0)  {return QTime::fromString(arg0); }
QTime static_QTime_fromString(const QString & arg0,const QString & arg1)  {return QTime::fromString(arg0,arg1); }
bool static_QTime_isValid(int arg0,int arg1,int arg2,int arg3)  {return QTime::isValid(arg0,arg1,arg2,arg3); }
bool static_QTime_isValid(int arg0,int arg1,int arg2)  {return QTime::isValid(arg0,arg1,arg2); }
void start(QTime* obj)  {obj->start(); }
int restart(QTime* obj)  {return obj->restart(); }
int elapsed(QTime* obj)  const  {return obj->elapsed(); }

};

class PythonQtQDateTimeWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QDateTime() { return QDateTime(); }
QVariant new_QDateTime(const QDate & arg0) { return QDateTime(arg0); }
QVariant new_QDateTime(const QDate & arg0,const QTime & arg1,Qt::TimeSpec arg2) { return QDateTime(arg0,arg1,arg2); }
QVariant new_QDateTime(const QDate & arg0,const QTime & arg1) { return QDateTime(arg0,arg1); }
QVariant new_QDateTime(const QDateTime & arg0) { return QDateTime(arg0); }
bool isNull(QDateTime* obj)  const  {return obj->isNull(); }
bool isValid(QDateTime* obj)  const  {return obj->isValid(); }
QDate date(QDateTime* obj)  const  {return obj->date(); }
QTime time(QDateTime* obj)  const  {return obj->time(); }
Qt::TimeSpec timeSpec(QDateTime* obj)  const  {return obj->timeSpec(); }
uint toTime_t(QDateTime* obj)  const  {return obj->toTime_t(); }
void setDate(QDateTime* obj,const QDate & arg0)  {obj->setDate(arg0); }
void setTime(QDateTime* obj,const QTime & arg0)  {obj->setTime(arg0); }
void setTimeSpec(QDateTime* obj,Qt::TimeSpec arg0)  {obj->setTimeSpec(arg0); }
void setTime_t(QDateTime* obj,uint arg0)  {obj->setTime_t(arg0); }
QString toString(QDateTime* obj,Qt::DateFormat arg0)  const  {return obj->toString(arg0); }
QString toString(QDateTime* obj)  const  {return obj->toString(); }
QString toString(QDateTime* obj,const QString & arg0)  const  {return obj->toString(arg0); }
QDateTime addDays(QDateTime* obj,int arg0)  const  {return obj->addDays(arg0); }
QDateTime addMonths(QDateTime* obj,int arg0)  const  {return obj->addMonths(arg0); }
QDateTime addYears(QDateTime* obj,int arg0)  const  {return obj->addYears(arg0); }
QDateTime addSecs(QDateTime* obj,int arg0)  const  {return obj->addSecs(arg0); }
QDateTime addMSecs(QDateTime* obj,qint64 arg0)  const  {return obj->addMSecs(arg0); }
QDateTime toTimeSpec(QDateTime* obj,Qt::TimeSpec arg0)  const  {return obj->toTimeSpec(arg0); }
QDateTime toLocalTime(QDateTime* obj)  const  {return obj->toLocalTime(); }
QDateTime toUTC(QDateTime* obj)  const  {return obj->toUTC(); }
int daysTo(QDateTime* obj,const QDateTime & arg0)  const  {return obj->daysTo(arg0); }
int secsTo(QDateTime* obj,const QDateTime & arg0)  const  {return obj->secsTo(arg0); }
QDateTime static_QDateTime_currentDateTime()  {return QDateTime::currentDateTime(); }
QDateTime static_QDateTime_fromString(const QString & arg0,Qt::DateFormat arg1)  {return QDateTime::fromString(arg0,arg1); }
QDateTime static_QDateTime_fromString(const QString & arg0)  {return QDateTime::fromString(arg0); }
QDateTime static_QDateTime_fromString(const QString & arg0,const QString & arg1)  {return QDateTime::fromString(arg0,arg1); }
QDateTime static_QDateTime_fromTime_t(uint arg0)  {return QDateTime::fromTime_t(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qurl.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qurl.h"
class PythonQtQUrlWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(ParsingMode FormattingOption )
enum ParsingMode {TolerantMode = QUrl::TolerantMode, 
StrictMode = QUrl::StrictMode }; 
enum FormattingOption {None = QUrl::None, 
RemoveScheme = QUrl::RemoveScheme, 
RemovePassword = QUrl::RemovePassword, 
RemoveUserInfo = QUrl::RemoveUserInfo, 
RemovePort = QUrl::RemovePort, 
RemoveAuthority = QUrl::RemoveAuthority, 
RemovePath = QUrl::RemovePath, 
RemoveQuery = QUrl::RemoveQuery, 
RemoveFragment = QUrl::RemoveFragment, 
StripTrailingSlash = QUrl::StripTrailingSlash }; 
Q_DECLARE_FLAGS(FormattingOptions, FormattingOption)
public slots:
QVariant new_QUrl(const QString & arg0) { return QUrl(arg0); }
QVariant new_QUrl(const QString & arg0,ParsingMode arg1) { return QUrl(arg0,(QUrl::ParsingMode)arg1); }
QVariant new_QUrl(const QUrl & arg0) { return QUrl(arg0); }
void setUrl(QUrl* obj,const QString & arg0)  {obj->setUrl(arg0); }
void setUrl(QUrl* obj,const QString & arg0,ParsingMode arg1)  {obj->setUrl(arg0,(QUrl::ParsingMode)arg1); }
void setEncodedUrl(QUrl* obj,const QByteArray & arg0)  {obj->setEncodedUrl(arg0); }
void setEncodedUrl(QUrl* obj,const QByteArray & arg0,ParsingMode arg1)  {obj->setEncodedUrl(arg0,(QUrl::ParsingMode)arg1); }
bool isValid(QUrl* obj)  const  {return obj->isValid(); }
bool isEmpty(QUrl* obj)  const  {return obj->isEmpty(); }
void clear(QUrl* obj)  {obj->clear(); }
void setScheme(QUrl* obj,const QString & arg0)  {obj->setScheme(arg0); }
QString scheme(QUrl* obj)  const  {return obj->scheme(); }
void setAuthority(QUrl* obj,const QString & arg0)  {obj->setAuthority(arg0); }
QString authority(QUrl* obj)  const  {return obj->authority(); }
void setUserInfo(QUrl* obj,const QString & arg0)  {obj->setUserInfo(arg0); }
QString userInfo(QUrl* obj)  const  {return obj->userInfo(); }
void setUserName(QUrl* obj,const QString & arg0)  {obj->setUserName(arg0); }
QString userName(QUrl* obj)  const  {return obj->userName(); }
void setPassword(QUrl* obj,const QString & arg0)  {obj->setPassword(arg0); }
QString password(QUrl* obj)  const  {return obj->password(); }
void setHost(QUrl* obj,const QString & arg0)  {obj->setHost(arg0); }
QString host(QUrl* obj)  const  {return obj->host(); }
void setPort(QUrl* obj,int arg0)  {obj->setPort(arg0); }
int port(QUrl* obj)  const  {return obj->port(); }
int port(QUrl* obj,int arg0)  const  {return obj->port(arg0); }
void setPath(QUrl* obj,const QString & arg0)  {obj->setPath(arg0); }
QString path(QUrl* obj)  const  {return obj->path(); }
bool hasQuery(QUrl* obj)  const  {return obj->hasQuery(); }
void setEncodedQuery(QUrl* obj,const QByteArray & arg0)  {obj->setEncodedQuery(arg0); }
QByteArray encodedQuery(QUrl* obj)  const  {return obj->encodedQuery(); }
void setQueryDelimiters(QUrl* obj,char arg0,char arg1)  {obj->setQueryDelimiters(arg0,arg1); }
char queryValueDelimiter(QUrl* obj)  const  {return obj->queryValueDelimiter(); }
char queryPairDelimiter(QUrl* obj)  const  {return obj->queryPairDelimiter(); }
void setQueryItems(QUrl* obj,const QList<QPair<QString,QString> > & arg0)  {obj->setQueryItems(arg0); }
void addQueryItem(QUrl* obj,const QString & arg0,const QString & arg1)  {obj->addQueryItem(arg0,arg1); }
QList<QPair<QString,QString> > queryItems(QUrl* obj)  const  {return obj->queryItems(); }
bool hasQueryItem(QUrl* obj,const QString & arg0)  const  {return obj->hasQueryItem(arg0); }
QString queryItemValue(QUrl* obj,const QString & arg0)  const  {return obj->queryItemValue(arg0); }
QStringList allQueryItemValues(QUrl* obj,const QString & arg0)  const  {return obj->allQueryItemValues(arg0); }
void removeQueryItem(QUrl* obj,const QString & arg0)  {obj->removeQueryItem(arg0); }
void removeAllQueryItems(QUrl* obj,const QString & arg0)  {obj->removeAllQueryItems(arg0); }
void setFragment(QUrl* obj,const QString & arg0)  {obj->setFragment(arg0); }
QString fragment(QUrl* obj)  const  {return obj->fragment(); }
bool hasFragment(QUrl* obj)  const  {return obj->hasFragment(); }
QUrl resolved(QUrl* obj,const QUrl & arg0)  const  {return obj->resolved(arg0); }
bool isRelative(QUrl* obj)  const  {return obj->isRelative(); }
bool isParentOf(QUrl* obj,const QUrl & arg0)  const  {return obj->isParentOf(arg0); }
QUrl static_QUrl_fromLocalFile(const QString & arg0)  {return QUrl::fromLocalFile(arg0); }
QString toLocalFile(QUrl* obj)  const  {return obj->toLocalFile(); }
QString toString(QUrl* obj,FormattingOptions arg0)  const  {return obj->toString((QUrl::FormattingOptions)QFlag(arg0)); }
QString toString(QUrl* obj)  const  {return obj->toString(); }
QByteArray toEncoded(QUrl* obj,FormattingOptions arg0)  const  {return obj->toEncoded((QUrl::FormattingOptions)QFlag(arg0)); }
QByteArray toEncoded(QUrl* obj)  const  {return obj->toEncoded(); }
QUrl static_QUrl_fromEncoded(const QByteArray & arg0)  {return QUrl::fromEncoded(arg0); }
QUrl static_QUrl_fromEncoded(const QByteArray & arg0,ParsingMode arg1)  {return QUrl::fromEncoded(arg0,(QUrl::ParsingMode)arg1); }
void detach(QUrl* obj)  {obj->detach(); }
bool isDetached(QUrl* obj)  const  {return obj->isDetached(); }
QString static_QUrl_fromPercentEncoding(const QByteArray & arg0)  {return QUrl::fromPercentEncoding(arg0); }
QByteArray static_QUrl_toPercentEncoding(const QString & arg0,const QByteArray & arg1,const QByteArray & arg2)  {return QUrl::toPercentEncoding(arg0,arg1,arg2); }
QByteArray static_QUrl_toPercentEncoding(const QString & arg0,const QByteArray & arg1)  {return QUrl::toPercentEncoding(arg0,arg1); }
QByteArray static_QUrl_toPercentEncoding(const QString & arg0)  {return QUrl::toPercentEncoding(arg0); }
QString static_QUrl_fromPunycode(const QByteArray & arg0)  {return QUrl::fromPunycode(arg0); }
QByteArray static_QUrl_toPunycode(const QString & arg0)  {return QUrl::toPunycode(arg0); }
QString static_QUrl_fromAce(const QByteArray & arg0)  {return QUrl::fromAce(arg0); }
QByteArray static_QUrl_toAce(const QString & arg0)  {return QUrl::toAce(arg0); }
QStringList static_QUrl_idnWhitelist()  {return QUrl::idnWhitelist(); }
void static_QUrl_setIdnWhitelist(const QStringList & arg0)  {QUrl::setIdnWhitelist(arg0); }
QString errorString(QUrl* obj)  const  {return obj->errorString(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qlocale.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qlocale.h"
class PythonQtQLocaleWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(Language Country FormatType NumberOption )
enum Language {C = QLocale::C, 
Abkhazian = QLocale::Abkhazian, 
Afan = QLocale::Afan, 
Afar = QLocale::Afar, 
Afrikaans = QLocale::Afrikaans, 
Albanian = QLocale::Albanian, 
Amharic = QLocale::Amharic, 
Arabic = QLocale::Arabic, 
Armenian = QLocale::Armenian, 
Assamese = QLocale::Assamese, 
Aymara = QLocale::Aymara, 
Azerbaijani = QLocale::Azerbaijani, 
Bashkir = QLocale::Bashkir, 
Basque = QLocale::Basque, 
Bengali = QLocale::Bengali, 
Bhutani = QLocale::Bhutani, 
Bihari = QLocale::Bihari, 
Bislama = QLocale::Bislama, 
Breton = QLocale::Breton, 
Bulgarian = QLocale::Bulgarian, 
Burmese = QLocale::Burmese, 
Byelorussian = QLocale::Byelorussian, 
Cambodian = QLocale::Cambodian, 
Catalan = QLocale::Catalan, 
Chinese = QLocale::Chinese, 
Corsican = QLocale::Corsican, 
Croatian = QLocale::Croatian, 
Czech = QLocale::Czech, 
Danish = QLocale::Danish, 
Dutch = QLocale::Dutch, 
English = QLocale::English, 
Esperanto = QLocale::Esperanto, 
Estonian = QLocale::Estonian, 
Faroese = QLocale::Faroese, 
FijiLanguage = QLocale::FijiLanguage, 
Finnish = QLocale::Finnish, 
French = QLocale::French, 
Frisian = QLocale::Frisian, 
Gaelic = QLocale::Gaelic, 
Galician = QLocale::Galician, 
Georgian = QLocale::Georgian, 
German = QLocale::German, 
Greek = QLocale::Greek, 
Greenlandic = QLocale::Greenlandic, 
Guarani = QLocale::Guarani, 
Gujarati = QLocale::Gujarati, 
Hausa = QLocale::Hausa, 
Hebrew = QLocale::Hebrew, 
Hindi = QLocale::Hindi, 
Hungarian = QLocale::Hungarian, 
Icelandic = QLocale::Icelandic, 
Indonesian = QLocale::Indonesian, 
Interlingua = QLocale::Interlingua, 
Interlingue = QLocale::Interlingue, 
Inuktitut = QLocale::Inuktitut, 
Inupiak = QLocale::Inupiak, 
Irish = QLocale::Irish, 
Italian = QLocale::Italian, 
Japanese = QLocale::Japanese, 
Javanese = QLocale::Javanese, 
Kannada = QLocale::Kannada, 
Kashmiri = QLocale::Kashmiri, 
Kazakh = QLocale::Kazakh, 
Kinyarwanda = QLocale::Kinyarwanda, 
Kirghiz = QLocale::Kirghiz, 
Korean = QLocale::Korean, 
Kurdish = QLocale::Kurdish, 
Kurundi = QLocale::Kurundi, 
Laothian = QLocale::Laothian, 
Latin = QLocale::Latin, 
Latvian = QLocale::Latvian, 
Lingala = QLocale::Lingala, 
Lithuanian = QLocale::Lithuanian, 
Macedonian = QLocale::Macedonian, 
Malagasy = QLocale::Malagasy, 
Malay = QLocale::Malay, 
Malayalam = QLocale::Malayalam, 
Maltese = QLocale::Maltese, 
Maori = QLocale::Maori, 
Marathi = QLocale::Marathi, 
Moldavian = QLocale::Moldavian, 
Mongolian = QLocale::Mongolian, 
NauruLanguage = QLocale::NauruLanguage, 
Nepali = QLocale::Nepali, 
Norwegian = QLocale::Norwegian, 
Occitan = QLocale::Occitan, 
Oriya = QLocale::Oriya, 
Pashto = QLocale::Pashto, 
Persian = QLocale::Persian, 
Polish = QLocale::Polish, 
Portuguese = QLocale::Portuguese, 
Punjabi = QLocale::Punjabi, 
Quechua = QLocale::Quechua, 
RhaetoRomance = QLocale::RhaetoRomance, 
Romanian = QLocale::Romanian, 
Russian = QLocale::Russian, 
Samoan = QLocale::Samoan, 
Sangho = QLocale::Sangho, 
Sanskrit = QLocale::Sanskrit, 
Serbian = QLocale::Serbian, 
SerboCroatian = QLocale::SerboCroatian, 
Sesotho = QLocale::Sesotho, 
Setswana = QLocale::Setswana, 
Shona = QLocale::Shona, 
Sindhi = QLocale::Sindhi, 
Singhalese = QLocale::Singhalese, 
Siswati = QLocale::Siswati, 
Slovak = QLocale::Slovak, 
Slovenian = QLocale::Slovenian, 
Somali = QLocale::Somali, 
Spanish = QLocale::Spanish, 
Sundanese = QLocale::Sundanese, 
Swahili = QLocale::Swahili, 
Swedish = QLocale::Swedish, 
Tagalog = QLocale::Tagalog, 
Tajik = QLocale::Tajik, 
Tamil = QLocale::Tamil, 
Tatar = QLocale::Tatar, 
Telugu = QLocale::Telugu, 
Thai = QLocale::Thai, 
Tibetan = QLocale::Tibetan, 
Tigrinya = QLocale::Tigrinya, 
TongaLanguage = QLocale::TongaLanguage, 
Tsonga = QLocale::Tsonga, 
Turkish = QLocale::Turkish, 
Turkmen = QLocale::Turkmen, 
Twi = QLocale::Twi, 
Uigur = QLocale::Uigur, 
Ukrainian = QLocale::Ukrainian, 
Urdu = QLocale::Urdu, 
Uzbek = QLocale::Uzbek, 
Vietnamese = QLocale::Vietnamese, 
Volapuk = QLocale::Volapuk, 
Welsh = QLocale::Welsh, 
Wolof = QLocale::Wolof, 
Xhosa = QLocale::Xhosa, 
Yiddish = QLocale::Yiddish, 
Yoruba = QLocale::Yoruba, 
Zhuang = QLocale::Zhuang, 
Zulu = QLocale::Zulu, 
Nynorsk = QLocale::Nynorsk, 
Bosnian = QLocale::Bosnian, 
Divehi = QLocale::Divehi, 
Manx = QLocale::Manx, 
Cornish = QLocale::Cornish, 
LastLanguage = QLocale::LastLanguage }; 
enum Country {AnyCountry = QLocale::AnyCountry, 
Afghanistan = QLocale::Afghanistan, 
Albania = QLocale::Albania, 
Algeria = QLocale::Algeria, 
AmericanSamoa = QLocale::AmericanSamoa, 
Andorra = QLocale::Andorra, 
Angola = QLocale::Angola, 
Anguilla = QLocale::Anguilla, 
Antarctica = QLocale::Antarctica, 
AntiguaAndBarbuda = QLocale::AntiguaAndBarbuda, 
Argentina = QLocale::Argentina, 
Armenia = QLocale::Armenia, 
Aruba = QLocale::Aruba, 
Australia = QLocale::Australia, 
Austria = QLocale::Austria, 
Azerbaijan = QLocale::Azerbaijan, 
Bahamas = QLocale::Bahamas, 
Bahrain = QLocale::Bahrain, 
Bangladesh = QLocale::Bangladesh, 
Barbados = QLocale::Barbados, 
Belarus = QLocale::Belarus, 
Belgium = QLocale::Belgium, 
Belize = QLocale::Belize, 
Benin = QLocale::Benin, 
Bermuda = QLocale::Bermuda, 
Bhutan = QLocale::Bhutan, 
Bolivia = QLocale::Bolivia, 
BosniaAndHerzegowina = QLocale::BosniaAndHerzegowina, 
Botswana = QLocale::Botswana, 
BouvetIsland = QLocale::BouvetIsland, 
Brazil = QLocale::Brazil, 
BritishIndianOceanTerritory = QLocale::BritishIndianOceanTerritory, 
BruneiDarussalam = QLocale::BruneiDarussalam, 
Bulgaria = QLocale::Bulgaria, 
BurkinaFaso = QLocale::BurkinaFaso, 
Burundi = QLocale::Burundi, 
Cambodia = QLocale::Cambodia, 
Cameroon = QLocale::Cameroon, 
Canada = QLocale::Canada, 
CapeVerde = QLocale::CapeVerde, 
CaymanIslands = QLocale::CaymanIslands, 
CentralAfricanRepublic = QLocale::CentralAfricanRepublic, 
Chad = QLocale::Chad, 
Chile = QLocale::Chile, 
China = QLocale::China, 
ChristmasIsland = QLocale::ChristmasIsland, 
CocosIslands = QLocale::CocosIslands, 
Colombia = QLocale::Colombia, 
Comoros = QLocale::Comoros, 
DemocraticRepublicOfCongo = QLocale::DemocraticRepublicOfCongo, 
PeoplesRepublicOfCongo = QLocale::PeoplesRepublicOfCongo, 
CookIslands = QLocale::CookIslands, 
CostaRica = QLocale::CostaRica, 
IvoryCoast = QLocale::IvoryCoast, 
Croatia = QLocale::Croatia, 
Cuba = QLocale::Cuba, 
Cyprus = QLocale::Cyprus, 
CzechRepublic = QLocale::CzechRepublic, 
Denmark = QLocale::Denmark, 
Djibouti = QLocale::Djibouti, 
Dominica = QLocale::Dominica, 
DominicanRepublic = QLocale::DominicanRepublic, 
EastTimor = QLocale::EastTimor, 
Ecuador = QLocale::Ecuador, 
Egypt = QLocale::Egypt, 
ElSalvador = QLocale::ElSalvador, 
EquatorialGuinea = QLocale::EquatorialGuinea, 
Eritrea = QLocale::Eritrea, 
Estonia = QLocale::Estonia, 
Ethiopia = QLocale::Ethiopia, 
FalklandIslands = QLocale::FalklandIslands, 
FaroeIslands = QLocale::FaroeIslands, 
FijiCountry = QLocale::FijiCountry, 
Finland = QLocale::Finland, 
France = QLocale::France, 
MetropolitanFrance = QLocale::MetropolitanFrance, 
FrenchGuiana = QLocale::FrenchGuiana, 
FrenchPolynesia = QLocale::FrenchPolynesia, 
FrenchSouthernTerritories = QLocale::FrenchSouthernTerritories, 
Gabon = QLocale::Gabon, 
Gambia = QLocale::Gambia, 
Georgia = QLocale::Georgia, 
Germany = QLocale::Germany, 
Ghana = QLocale::Ghana, 
Gibraltar = QLocale::Gibraltar, 
Greece = QLocale::Greece, 
Greenland = QLocale::Greenland, 
Grenada = QLocale::Grenada, 
Guadeloupe = QLocale::Guadeloupe, 
Guam = QLocale::Guam, 
Guatemala = QLocale::Guatemala, 
Guinea = QLocale::Guinea, 
GuineaBissau = QLocale::GuineaBissau, 
Guyana = QLocale::Guyana, 
Haiti = QLocale::Haiti, 
HeardAndMcDonaldIslands = QLocale::HeardAndMcDonaldIslands, 
Honduras = QLocale::Honduras, 
HongKong = QLocale::HongKong, 
Hungary = QLocale::Hungary, 
Iceland = QLocale::Iceland, 
India = QLocale::India, 
Indonesia = QLocale::Indonesia, 
Iran = QLocale::Iran, 
Iraq = QLocale::Iraq, 
Ireland = QLocale::Ireland, 
Israel = QLocale::Israel, 
Italy = QLocale::Italy, 
Jamaica = QLocale::Jamaica, 
Japan = QLocale::Japan, 
Jordan = QLocale::Jordan, 
Kazakhstan = QLocale::Kazakhstan, 
Kenya = QLocale::Kenya, 
Kiribati = QLocale::Kiribati, 
DemocraticRepublicOfKorea = QLocale::DemocraticRepublicOfKorea, 
RepublicOfKorea = QLocale::RepublicOfKorea, 
Kuwait = QLocale::Kuwait, 
Kyrgyzstan = QLocale::Kyrgyzstan, 
Lao = QLocale::Lao, 
Latvia = QLocale::Latvia, 
Lebanon = QLocale::Lebanon, 
Lesotho = QLocale::Lesotho, 
Liberia = QLocale::Liberia, 
LibyanArabJamahiriya = QLocale::LibyanArabJamahiriya, 
Liechtenstein = QLocale::Liechtenstein, 
Lithuania = QLocale::Lithuania, 
Luxembourg = QLocale::Luxembourg, 
Macau = QLocale::Macau, 
Macedonia = QLocale::Macedonia, 
Madagascar = QLocale::Madagascar, 
Malawi = QLocale::Malawi, 
Malaysia = QLocale::Malaysia, 
Maldives = QLocale::Maldives, 
Mali = QLocale::Mali, 
Malta = QLocale::Malta, 
MarshallIslands = QLocale::MarshallIslands, 
Martinique = QLocale::Martinique, 
Mauritania = QLocale::Mauritania, 
Mauritius = QLocale::Mauritius, 
Mayotte = QLocale::Mayotte, 
Mexico = QLocale::Mexico, 
Micronesia = QLocale::Micronesia, 
Moldova = QLocale::Moldova, 
Monaco = QLocale::Monaco, 
Mongolia = QLocale::Mongolia, 
Montserrat = QLocale::Montserrat, 
Morocco = QLocale::Morocco, 
Mozambique = QLocale::Mozambique, 
Myanmar = QLocale::Myanmar, 
Namibia = QLocale::Namibia, 
NauruCountry = QLocale::NauruCountry, 
Nepal = QLocale::Nepal, 
Netherlands = QLocale::Netherlands, 
NetherlandsAntilles = QLocale::NetherlandsAntilles, 
NewCaledonia = QLocale::NewCaledonia, 
NewZealand = QLocale::NewZealand, 
Nicaragua = QLocale::Nicaragua, 
Niger = QLocale::Niger, 
Nigeria = QLocale::Nigeria, 
Niue = QLocale::Niue, 
NorfolkIsland = QLocale::NorfolkIsland, 
NorthernMarianaIslands = QLocale::NorthernMarianaIslands, 
Norway = QLocale::Norway, 
Oman = QLocale::Oman, 
Pakistan = QLocale::Pakistan, 
Palau = QLocale::Palau, 
PalestinianTerritory = QLocale::PalestinianTerritory, 
Panama = QLocale::Panama, 
PapuaNewGuinea = QLocale::PapuaNewGuinea, 
Paraguay = QLocale::Paraguay, 
Peru = QLocale::Peru, 
Philippines = QLocale::Philippines, 
Pitcairn = QLocale::Pitcairn, 
Poland = QLocale::Poland, 
Portugal = QLocale::Portugal, 
PuertoRico = QLocale::PuertoRico, 
Qatar = QLocale::Qatar, 
Reunion = QLocale::Reunion, 
Romania = QLocale::Romania, 
RussianFederation = QLocale::RussianFederation, 
Rwanda = QLocale::Rwanda, 
SaintKittsAndNevis = QLocale::SaintKittsAndNevis, 
StLucia = QLocale::StLucia, 
StVincentAndTheGrenadines = QLocale::StVincentAndTheGrenadines, 
Samoa = QLocale::Samoa, 
SanMarino = QLocale::SanMarino, 
SaoTomeAndPrincipe = QLocale::SaoTomeAndPrincipe, 
SaudiArabia = QLocale::SaudiArabia, 
Senegal = QLocale::Senegal, 
Seychelles = QLocale::Seychelles, 
SierraLeone = QLocale::SierraLeone, 
Singapore = QLocale::Singapore, 
Slovakia = QLocale::Slovakia, 
Slovenia = QLocale::Slovenia, 
SolomonIslands = QLocale::SolomonIslands, 
Somalia = QLocale::Somalia, 
SouthAfrica = QLocale::SouthAfrica, 
SouthGeorgiaAndTheSouthSandwichIslands = QLocale::SouthGeorgiaAndTheSouthSandwichIslands, 
Spain = QLocale::Spain, 
SriLanka = QLocale::SriLanka, 
StHelena = QLocale::StHelena, 
StPierreAndMiquelon = QLocale::StPierreAndMiquelon, 
Sudan = QLocale::Sudan, 
Suriname = QLocale::Suriname, 
SvalbardAndJanMayenIslands = QLocale::SvalbardAndJanMayenIslands, 
Swaziland = QLocale::Swaziland, 
Sweden = QLocale::Sweden, 
Switzerland = QLocale::Switzerland, 
SyrianArabRepublic = QLocale::SyrianArabRepublic, 
Taiwan = QLocale::Taiwan, 
Tajikistan = QLocale::Tajikistan, 
Tanzania = QLocale::Tanzania, 
Thailand = QLocale::Thailand, 
Togo = QLocale::Togo, 
Tokelau = QLocale::Tokelau, 
TongaCountry = QLocale::TongaCountry, 
TrinidadAndTobago = QLocale::TrinidadAndTobago, 
Tunisia = QLocale::Tunisia, 
Turkey = QLocale::Turkey, 
Turkmenistan = QLocale::Turkmenistan, 
TurksAndCaicosIslands = QLocale::TurksAndCaicosIslands, 
Tuvalu = QLocale::Tuvalu, 
Uganda = QLocale::Uganda, 
Ukraine = QLocale::Ukraine, 
UnitedArabEmirates = QLocale::UnitedArabEmirates, 
UnitedKingdom = QLocale::UnitedKingdom, 
UnitedStates = QLocale::UnitedStates, 
UnitedStatesMinorOutlyingIslands = QLocale::UnitedStatesMinorOutlyingIslands, 
Uruguay = QLocale::Uruguay, 
Uzbekistan = QLocale::Uzbekistan, 
Vanuatu = QLocale::Vanuatu, 
VaticanCityState = QLocale::VaticanCityState, 
Venezuela = QLocale::Venezuela, 
VietNam = QLocale::VietNam, 
BritishVirginIslands = QLocale::BritishVirginIslands, 
USVirginIslands = QLocale::USVirginIslands, 
WallisAndFutunaIslands = QLocale::WallisAndFutunaIslands, 
WesternSahara = QLocale::WesternSahara, 
Yemen = QLocale::Yemen, 
Yugoslavia = QLocale::Yugoslavia, 
Zambia = QLocale::Zambia, 
Zimbabwe = QLocale::Zimbabwe, 
SerbiaAndMontenegro = QLocale::SerbiaAndMontenegro, 
LastCountry = QLocale::LastCountry }; 
enum FormatType {LongFormat = QLocale::LongFormat, 
ShortFormat = QLocale::ShortFormat }; 
enum NumberOption {OmitGroupSeparator = QLocale::OmitGroupSeparator, 
RejectGroupSeparator = QLocale::RejectGroupSeparator }; 
Q_DECLARE_FLAGS(NumberOptions, NumberOption)
public slots:
QVariant new_QLocale(const QString & arg0) { return QLocale(arg0); }
QVariant new_QLocale(Language arg0,Country arg1) { return QLocale((QLocale::Language)arg0,(QLocale::Country)arg1); }
QVariant new_QLocale(Language arg0) { return QLocale((QLocale::Language)arg0); }
QVariant new_QLocale(const QLocale & arg0) { return QLocale(arg0); }
Language language(QLocale* obj)  const  {return (PythonQtQLocaleWrapper::Language)obj->language(); }
Country country(QLocale* obj)  const  {return (PythonQtQLocaleWrapper::Country)obj->country(); }
QString name(QLocale* obj)  const  {return obj->name(); }
short toShort(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toShort(arg0,arg1,arg2); }
short toShort(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toShort(arg0,arg1); }
short toShort(QLocale* obj,const QString & arg0)  const  {return obj->toShort(arg0); }
ushort toUShort(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toUShort(arg0,arg1,arg2); }
ushort toUShort(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toUShort(arg0,arg1); }
ushort toUShort(QLocale* obj,const QString & arg0)  const  {return obj->toUShort(arg0); }
int toInt(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toInt(arg0,arg1,arg2); }
int toInt(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toInt(arg0,arg1); }
int toInt(QLocale* obj,const QString & arg0)  const  {return obj->toInt(arg0); }
uint toUInt(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toUInt(arg0,arg1,arg2); }
uint toUInt(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toUInt(arg0,arg1); }
uint toUInt(QLocale* obj,const QString & arg0)  const  {return obj->toUInt(arg0); }
qlonglong toLongLong(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toLongLong(arg0,arg1,arg2); }
qlonglong toLongLong(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toLongLong(arg0,arg1); }
qlonglong toLongLong(QLocale* obj,const QString & arg0)  const  {return obj->toLongLong(arg0); }
qlonglong toULongLong(QLocale* obj,const QString & arg0,bool * arg1,int arg2)  const  {return obj->toULongLong(arg0,arg1,arg2); }
qlonglong toULongLong(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toULongLong(arg0,arg1); }
qlonglong toULongLong(QLocale* obj,const QString & arg0)  const  {return obj->toULongLong(arg0); }
float toFloat(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toFloat(arg0,arg1); }
float toFloat(QLocale* obj,const QString & arg0)  const  {return obj->toFloat(arg0); }
double toDouble(QLocale* obj,const QString & arg0,bool * arg1)  const  {return obj->toDouble(arg0,arg1); }
double toDouble(QLocale* obj,const QString & arg0)  const  {return obj->toDouble(arg0); }
QString toString(QLocale* obj,qlonglong arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,qulonglong arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,short arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,ushort arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,int arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,uint arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,double arg0,char arg1,int arg2)  const  {return obj->toString(arg0,arg1,arg2); }
QString toString(QLocale* obj,double arg0,char arg1)  const  {return obj->toString(arg0,arg1); }
QString toString(QLocale* obj,double arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,float arg0,char arg1,int arg2)  const  {return obj->toString(arg0,arg1,arg2); }
QString toString(QLocale* obj,float arg0,char arg1)  const  {return obj->toString(arg0,arg1); }
QString toString(QLocale* obj,float arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,const QDate & arg0,const QString & arg1)  const  {return obj->toString(arg0,arg1); }
QString toString(QLocale* obj,const QDate & arg0,FormatType arg1)  const  {return obj->toString(arg0,(QLocale::FormatType)arg1); }
QString toString(QLocale* obj,const QDate & arg0)  const  {return obj->toString(arg0); }
QString toString(QLocale* obj,const QTime & arg0,const QString & arg1)  const  {return obj->toString(arg0,arg1); }
QString toString(QLocale* obj,const QTime & arg0,FormatType arg1)  const  {return obj->toString(arg0,(QLocale::FormatType)arg1); }
QString toString(QLocale* obj,const QTime & arg0)  const  {return obj->toString(arg0); }
QString dateFormat(QLocale* obj,FormatType arg0)  const  {return obj->dateFormat((QLocale::FormatType)arg0); }
QString dateFormat(QLocale* obj)  const  {return obj->dateFormat(); }
QString timeFormat(QLocale* obj,FormatType arg0)  const  {return obj->timeFormat((QLocale::FormatType)arg0); }
QString timeFormat(QLocale* obj)  const  {return obj->timeFormat(); }
QChar decimalPoint(QLocale* obj)  const  {return obj->decimalPoint(); }
QChar groupSeparator(QLocale* obj)  const  {return obj->groupSeparator(); }
QChar percent(QLocale* obj)  const  {return obj->percent(); }
QChar zeroDigit(QLocale* obj)  const  {return obj->zeroDigit(); }
QChar negativeSign(QLocale* obj)  const  {return obj->negativeSign(); }
QChar exponential(QLocale* obj)  const  {return obj->exponential(); }
QString monthName(QLocale* obj,int arg0,FormatType arg1)  const  {return obj->monthName(arg0,(QLocale::FormatType)arg1); }
QString monthName(QLocale* obj,int arg0)  const  {return obj->monthName(arg0); }
QString dayName(QLocale* obj,int arg0,FormatType arg1)  const  {return obj->dayName(arg0,(QLocale::FormatType)arg1); }
QString dayName(QLocale* obj,int arg0)  const  {return obj->dayName(arg0); }
QString static_QLocale_languageToString(Language arg0)  {return QLocale::languageToString((QLocale::Language)arg0); }
QString static_QLocale_countryToString(Country arg0)  {return QLocale::countryToString((QLocale::Country)arg0); }
void static_QLocale_setDefault(const QLocale & arg0)  {QLocale::setDefault(arg0); }
QLocale static_QLocale_c()  {return QLocale::c(); }
QLocale static_QLocale_system()  {return QLocale::system(); }
void setNumberOptions(QLocale* obj,NumberOptions arg0)  {obj->setNumberOptions((QLocale::NumberOptions)QFlag(arg0)); }
NumberOptions numberOptions(QLocale* obj)  const  {return (PythonQtQLocaleWrapper::NumberOptions)QFlag(obj->numberOptions()); }

};

class PythonQtQSystemLocaleWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(QueryType )
enum QueryType {LanguageId = QSystemLocale::LanguageId, 
CountryId = QSystemLocale::CountryId, 
DecimalPoint = QSystemLocale::DecimalPoint, 
GroupSeparator = QSystemLocale::GroupSeparator, 
ZeroDigit = QSystemLocale::ZeroDigit, 
NegativeSign = QSystemLocale::NegativeSign, 
DateFormatLong = QSystemLocale::DateFormatLong, 
DateFormatShort = QSystemLocale::DateFormatShort, 
TimeFormatLong = QSystemLocale::TimeFormatLong, 
TimeFormatShort = QSystemLocale::TimeFormatShort, 
DayNameLong = QSystemLocale::DayNameLong, 
DayNameShort = QSystemLocale::DayNameShort, 
MonthNameLong = QSystemLocale::MonthNameLong, 
MonthNameShort = QSystemLocale::MonthNameShort, 
DateToStringLong = QSystemLocale::DateToStringLong, 
DateToStringShort = QSystemLocale::DateToStringShort, 
TimeToStringLong = QSystemLocale::TimeToStringLong, 
TimeToStringShort = QSystemLocale::TimeToStringShort }; 
public slots:
void delete_QSystemLocale(QSystemLocale* obj) { delete obj; }
QSystemLocale* new_QSystemLocale() { return new QSystemLocale(); }
QVariant query(QSystemLocale* obj,QueryType arg0,QVariant arg1)  const  {return obj->query((QSystemLocale::QueryType)arg0,arg1); }
QLocale fallbackLocale(QSystemLocale* obj)  const  {return obj->fallbackLocale(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qrect.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qrect.h"
class PythonQtQRectWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QRect() { return QRect(); }
QVariant new_QRect(const QPoint & arg0,const QPoint & arg1) { return QRect(arg0,arg1); }
QVariant new_QRect(const QPoint & arg0,const QSize & arg1) { return QRect(arg0,arg1); }
QVariant new_QRect(int arg0,int arg1,int arg2,int arg3) { return QRect(arg0,arg1,arg2,arg3); }
bool isNull(QRect* obj)  const  {return obj->isNull(); }
bool isEmpty(QRect* obj)  const  {return obj->isEmpty(); }
bool isValid(QRect* obj)  const  {return obj->isValid(); }
int left(QRect* obj)  const  {return obj->left(); }
int top(QRect* obj)  const  {return obj->top(); }
int right(QRect* obj)  const  {return obj->right(); }
int bottom(QRect* obj)  const  {return obj->bottom(); }
QRect normalized(QRect* obj)  const  {return obj->normalized(); }
int x(QRect* obj)  const  {return obj->x(); }
int y(QRect* obj)  const  {return obj->y(); }
void setLeft(QRect* obj,int arg0)  {obj->setLeft(arg0); }
void setTop(QRect* obj,int arg0)  {obj->setTop(arg0); }
void setRight(QRect* obj,int arg0)  {obj->setRight(arg0); }
void setBottom(QRect* obj,int arg0)  {obj->setBottom(arg0); }
void setX(QRect* obj,int arg0)  {obj->setX(arg0); }
void setY(QRect* obj,int arg0)  {obj->setY(arg0); }
void setTopLeft(QRect* obj,const QPoint & arg0)  {obj->setTopLeft(arg0); }
void setBottomRight(QRect* obj,const QPoint & arg0)  {obj->setBottomRight(arg0); }
void setTopRight(QRect* obj,const QPoint & arg0)  {obj->setTopRight(arg0); }
void setBottomLeft(QRect* obj,const QPoint & arg0)  {obj->setBottomLeft(arg0); }
QPoint topLeft(QRect* obj)  const  {return obj->topLeft(); }
QPoint bottomRight(QRect* obj)  const  {return obj->bottomRight(); }
QPoint topRight(QRect* obj)  const  {return obj->topRight(); }
QPoint bottomLeft(QRect* obj)  const  {return obj->bottomLeft(); }
QPoint center(QRect* obj)  const  {return obj->center(); }
void moveLeft(QRect* obj,int arg0)  {obj->moveLeft(arg0); }
void moveTop(QRect* obj,int arg0)  {obj->moveTop(arg0); }
void moveRight(QRect* obj,int arg0)  {obj->moveRight(arg0); }
void moveBottom(QRect* obj,int arg0)  {obj->moveBottom(arg0); }
void moveTopLeft(QRect* obj,const QPoint & arg0)  {obj->moveTopLeft(arg0); }
void moveBottomRight(QRect* obj,const QPoint & arg0)  {obj->moveBottomRight(arg0); }
void moveTopRight(QRect* obj,const QPoint & arg0)  {obj->moveTopRight(arg0); }
void moveBottomLeft(QRect* obj,const QPoint & arg0)  {obj->moveBottomLeft(arg0); }
void moveCenter(QRect* obj,const QPoint & arg0)  {obj->moveCenter(arg0); }
void translate(QRect* obj,int arg0,int arg1)  {obj->translate(arg0,arg1); }
void translate(QRect* obj,const QPoint & arg0)  {obj->translate(arg0); }
QRect translated(QRect* obj,int arg0,int arg1)  const  {return obj->translated(arg0,arg1); }
QRect translated(QRect* obj,const QPoint & arg0)  const  {return obj->translated(arg0); }
void moveTo(QRect* obj,int arg0,int arg1)  {obj->moveTo(arg0,arg1); }
void moveTo(QRect* obj,const QPoint & arg0)  {obj->moveTo(arg0); }
void setRect(QRect* obj,int arg0,int arg1,int arg2,int arg3)  {obj->setRect(arg0,arg1,arg2,arg3); }
void getRect(QRect* obj,int * arg0,int * arg1,int * arg2,int * arg3)  const  {obj->getRect(arg0,arg1,arg2,arg3); }
void setCoords(QRect* obj,int arg0,int arg1,int arg2,int arg3)  {obj->setCoords(arg0,arg1,arg2,arg3); }
void getCoords(QRect* obj,int * arg0,int * arg1,int * arg2,int * arg3)  const  {obj->getCoords(arg0,arg1,arg2,arg3); }
void adjust(QRect* obj,int arg0,int arg1,int arg2,int arg3)  {obj->adjust(arg0,arg1,arg2,arg3); }
QRect adjusted(QRect* obj,int arg0,int arg1,int arg2,int arg3)  const  {return obj->adjusted(arg0,arg1,arg2,arg3); }
QSize size(QRect* obj)  const  {return obj->size(); }
int width(QRect* obj)  const  {return obj->width(); }
int height(QRect* obj)  const  {return obj->height(); }
void setWidth(QRect* obj,int arg0)  {obj->setWidth(arg0); }
void setHeight(QRect* obj,int arg0)  {obj->setHeight(arg0); }
void setSize(QRect* obj,const QSize & arg0)  {obj->setSize(arg0); }
bool contains(QRect* obj,const QPoint & arg0,bool arg1)  const  {return obj->contains(arg0,arg1); }
bool contains(QRect* obj,const QPoint & arg0)  const  {return obj->contains(arg0); }
bool contains(QRect* obj,int arg0,int arg1)  const  {return obj->contains(arg0,arg1); }
bool contains(QRect* obj,int arg0,int arg1,bool arg2)  const  {return obj->contains(arg0,arg1,arg2); }
bool contains(QRect* obj,const QRect & arg0,bool arg1)  const  {return obj->contains(arg0,arg1); }
bool contains(QRect* obj,const QRect & arg0)  const  {return obj->contains(arg0); }
QRect unite(QRect* obj,const QRect & arg0)  const  {return obj->unite(arg0); }
QRect united(QRect* obj,const QRect & arg0)  const  {return obj->united(arg0); }
QRect intersect(QRect* obj,const QRect & arg0)  const  {return obj->intersect(arg0); }
QRect intersected(QRect* obj,const QRect & arg0)  const  {return obj->intersected(arg0); }
bool intersects(QRect* obj,const QRect & arg0)  const  {return obj->intersects(arg0); }

};

class PythonQtQRectFWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QRectF() { return QRectF(); }
QVariant new_QRectF(const QPointF & arg0,const QSizeF & arg1) { return QRectF(arg0,arg1); }
QVariant new_QRectF(qreal arg0,qreal arg1,qreal arg2,qreal arg3) { return QRectF(arg0,arg1,arg2,arg3); }
QVariant new_QRectF(const QRect & arg0) { return QRectF(arg0); }
bool isNull(QRectF* obj)  const  {return obj->isNull(); }
bool isEmpty(QRectF* obj)  const  {return obj->isEmpty(); }
bool isValid(QRectF* obj)  const  {return obj->isValid(); }
QRectF normalized(QRectF* obj)  const  {return obj->normalized(); }
qreal left(QRectF* obj)  const  {return obj->left(); }
qreal top(QRectF* obj)  const  {return obj->top(); }
qreal right(QRectF* obj)  const  {return obj->right(); }
qreal bottom(QRectF* obj)  const  {return obj->bottom(); }
qreal x(QRectF* obj)  const  {return obj->x(); }
qreal y(QRectF* obj)  const  {return obj->y(); }
void setLeft(QRectF* obj,qreal arg0)  {obj->setLeft(arg0); }
void setTop(QRectF* obj,qreal arg0)  {obj->setTop(arg0); }
void setRight(QRectF* obj,qreal arg0)  {obj->setRight(arg0); }
void setBottom(QRectF* obj,qreal arg0)  {obj->setBottom(arg0); }
void setX(QRectF* obj,qreal arg0)  {obj->setX(arg0); }
void setY(QRectF* obj,qreal arg0)  {obj->setY(arg0); }
QPointF topLeft(QRectF* obj)  const  {return obj->topLeft(); }
QPointF bottomRight(QRectF* obj)  const  {return obj->bottomRight(); }
QPointF topRight(QRectF* obj)  const  {return obj->topRight(); }
QPointF bottomLeft(QRectF* obj)  const  {return obj->bottomLeft(); }
QPointF center(QRectF* obj)  const  {return obj->center(); }
void setTopLeft(QRectF* obj,const QPointF & arg0)  {obj->setTopLeft(arg0); }
void setBottomRight(QRectF* obj,const QPointF & arg0)  {obj->setBottomRight(arg0); }
void setTopRight(QRectF* obj,const QPointF & arg0)  {obj->setTopRight(arg0); }
void setBottomLeft(QRectF* obj,const QPointF & arg0)  {obj->setBottomLeft(arg0); }
void moveLeft(QRectF* obj,qreal arg0)  {obj->moveLeft(arg0); }
void moveTop(QRectF* obj,qreal arg0)  {obj->moveTop(arg0); }
void moveRight(QRectF* obj,qreal arg0)  {obj->moveRight(arg0); }
void moveBottom(QRectF* obj,qreal arg0)  {obj->moveBottom(arg0); }
void moveTopLeft(QRectF* obj,const QPointF & arg0)  {obj->moveTopLeft(arg0); }
void moveBottomRight(QRectF* obj,const QPointF & arg0)  {obj->moveBottomRight(arg0); }
void moveTopRight(QRectF* obj,const QPointF & arg0)  {obj->moveTopRight(arg0); }
void moveBottomLeft(QRectF* obj,const QPointF & arg0)  {obj->moveBottomLeft(arg0); }
void moveCenter(QRectF* obj,const QPointF & arg0)  {obj->moveCenter(arg0); }
void translate(QRectF* obj,qreal arg0,qreal arg1)  {obj->translate(arg0,arg1); }
void translate(QRectF* obj,const QPointF & arg0)  {obj->translate(arg0); }
QRectF translated(QRectF* obj,qreal arg0,qreal arg1)  const  {return obj->translated(arg0,arg1); }
QRectF translated(QRectF* obj,const QPointF & arg0)  const  {return obj->translated(arg0); }
void moveTo(QRectF* obj,qreal arg0,qreal arg1)  {obj->moveTo(arg0,arg1); }
void moveTo(QRectF* obj,const QPointF & arg0)  {obj->moveTo(arg0); }
void setRect(QRectF* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->setRect(arg0,arg1,arg2,arg3); }
void getRect(QRectF* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3)  const  {obj->getRect(arg0,arg1,arg2,arg3); }
void setCoords(QRectF* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->setCoords(arg0,arg1,arg2,arg3); }
void getCoords(QRectF* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3)  const  {obj->getCoords(arg0,arg1,arg2,arg3); }
void adjust(QRectF* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->adjust(arg0,arg1,arg2,arg3); }
QRectF adjusted(QRectF* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  const  {return obj->adjusted(arg0,arg1,arg2,arg3); }
QSizeF size(QRectF* obj)  const  {return obj->size(); }
qreal width(QRectF* obj)  const  {return obj->width(); }
qreal height(QRectF* obj)  const  {return obj->height(); }
void setWidth(QRectF* obj,qreal arg0)  {obj->setWidth(arg0); }
void setHeight(QRectF* obj,qreal arg0)  {obj->setHeight(arg0); }
void setSize(QRectF* obj,const QSizeF & arg0)  {obj->setSize(arg0); }
bool contains(QRectF* obj,const QPointF & arg0)  const  {return obj->contains(arg0); }
bool contains(QRectF* obj,qreal arg0,qreal arg1)  const  {return obj->contains(arg0,arg1); }
bool contains(QRectF* obj,const QRectF & arg0)  const  {return obj->contains(arg0); }
QRectF unite(QRectF* obj,const QRectF & arg0)  const  {return obj->unite(arg0); }
QRectF united(QRectF* obj,const QRectF & arg0)  const  {return obj->united(arg0); }
QRectF intersect(QRectF* obj,const QRectF & arg0)  const  {return obj->intersect(arg0); }
QRectF intersected(QRectF* obj,const QRectF & arg0)  const  {return obj->intersected(arg0); }
bool intersects(QRectF* obj,const QRectF & arg0)  const  {return obj->intersects(arg0); }
QRect toRect(QRectF* obj)  const  {return obj->toRect(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qsize.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qsize.h"
class PythonQtQSizeWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QSize() { return QSize(); }
QVariant new_QSize(int arg0,int arg1) { return QSize(arg0,arg1); }
bool isNull(QSize* obj)  const  {return obj->isNull(); }
bool isEmpty(QSize* obj)  const  {return obj->isEmpty(); }
bool isValid(QSize* obj)  const  {return obj->isValid(); }
int width(QSize* obj)  const  {return obj->width(); }
int height(QSize* obj)  const  {return obj->height(); }
void setWidth(QSize* obj,int arg0)  {obj->setWidth(arg0); }
void setHeight(QSize* obj,int arg0)  {obj->setHeight(arg0); }
void transpose(QSize* obj)  {obj->transpose(); }
void scale(QSize* obj,int arg0,int arg1,Qt::AspectRatioMode arg2)  {obj->scale(arg0,arg1,arg2); }
void scale(QSize* obj,const QSize & arg0,Qt::AspectRatioMode arg1)  {obj->scale(arg0,arg1); }
QSize expandedTo(QSize* obj,const QSize & arg0)  const  {return obj->expandedTo(arg0); }
QSize boundedTo(QSize* obj,const QSize & arg0)  const  {return obj->boundedTo(arg0); }
void rwidth(QSize* obj)  {obj->rwidth(); }
void rheight(QSize* obj)  {obj->rheight(); }

};

class PythonQtQSizeFWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QSizeF() { return QSizeF(); }
QVariant new_QSizeF(const QSize & arg0) { return QSizeF(arg0); }
QVariant new_QSizeF(qreal arg0,qreal arg1) { return QSizeF(arg0,arg1); }
bool isNull(QSizeF* obj)  const  {return obj->isNull(); }
bool isEmpty(QSizeF* obj)  const  {return obj->isEmpty(); }
bool isValid(QSizeF* obj)  const  {return obj->isValid(); }
qreal width(QSizeF* obj)  const  {return obj->width(); }
qreal height(QSizeF* obj)  const  {return obj->height(); }
void setWidth(QSizeF* obj,qreal arg0)  {obj->setWidth(arg0); }
void setHeight(QSizeF* obj,qreal arg0)  {obj->setHeight(arg0); }
void transpose(QSizeF* obj)  {obj->transpose(); }
void scale(QSizeF* obj,qreal arg0,qreal arg1,Qt::AspectRatioMode arg2)  {obj->scale(arg0,arg1,arg2); }
void scale(QSizeF* obj,const QSizeF & arg0,Qt::AspectRatioMode arg1)  {obj->scale(arg0,arg1); }
QSizeF expandedTo(QSizeF* obj,const QSizeF & arg0)  const  {return obj->expandedTo(arg0); }
QSizeF boundedTo(QSizeF* obj,const QSizeF & arg0)  const  {return obj->boundedTo(arg0); }
void rwidth(QSizeF* obj)  {obj->rwidth(); }
void rheight(QSizeF* obj)  {obj->rheight(); }
QSize toSize(QSizeF* obj)  const  {return obj->toSize(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qline.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qline.h"
class PythonQtQLineWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QLine() { return QLine(); }
QVariant new_QLine(const QPoint & arg0,const QPoint & arg1) { return QLine(arg0,arg1); }
QVariant new_QLine(int arg0,int arg1,int arg2,int arg3) { return QLine(arg0,arg1,arg2,arg3); }
bool isNull(QLine* obj)  const  {return obj->isNull(); }
QPoint p1(QLine* obj)  const  {return obj->p1(); }
QPoint p2(QLine* obj)  const  {return obj->p2(); }
int x1(QLine* obj)  const  {return obj->x1(); }
int y1(QLine* obj)  const  {return obj->y1(); }
int x2(QLine* obj)  const  {return obj->x2(); }
int y2(QLine* obj)  const  {return obj->y2(); }
int dx(QLine* obj)  const  {return obj->dx(); }
int dy(QLine* obj)  const  {return obj->dy(); }
void translate(QLine* obj,const QPoint & arg0)  {obj->translate(arg0); }
void translate(QLine* obj,int arg0,int arg1)  {obj->translate(arg0,arg1); }

};

class PythonQtQLineFWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(IntersectType )
enum IntersectType {NoIntersection = QLineF::NoIntersection, 
BoundedIntersection = QLineF::BoundedIntersection, 
UnboundedIntersection = QLineF::UnboundedIntersection }; 
public slots:
QVariant new_QLineF() { return QLineF(); }
QVariant new_QLineF(const QPointF & arg0,const QPointF & arg1) { return QLineF(arg0,arg1); }
QVariant new_QLineF(qreal arg0,qreal arg1,qreal arg2,qreal arg3) { return QLineF(arg0,arg1,arg2,arg3); }
QVariant new_QLineF(const QLine & arg0) { return QLineF(arg0); }
int isNull(QLineF* obj)  const  {return obj->isNull(); }
QPointF p1(QLineF* obj)  const  {return obj->p1(); }
QPointF p2(QLineF* obj)  const  {return obj->p2(); }
qreal x1(QLineF* obj)  const  {return obj->x1(); }
qreal y1(QLineF* obj)  const  {return obj->y1(); }
qreal x2(QLineF* obj)  const  {return obj->x2(); }
qreal y2(QLineF* obj)  const  {return obj->y2(); }
qreal dx(QLineF* obj)  const  {return obj->dx(); }
qreal dy(QLineF* obj)  const  {return obj->dy(); }
qreal length(QLineF* obj)  const  {return obj->length(); }
void setLength(QLineF* obj,qreal arg0)  {obj->setLength(arg0); }
QLineF unitVector(QLineF* obj)  const  {return obj->unitVector(); }
QLineF normalVector(QLineF* obj)  const  {return obj->normalVector(); }
IntersectType intersect(QLineF* obj,const QLineF & arg0,QPointF * arg1)  const  {return (PythonQtQLineFWrapper::IntersectType)obj->intersect(arg0,arg1); }
qreal angle(QLineF* obj,const QLineF & arg0)  const  {return obj->angle(arg0); }
QPointF pointAt(QLineF* obj,qreal arg0)  const  {return obj->pointAt(arg0); }
void translate(QLineF* obj,const QPointF & arg0)  {obj->translate(arg0); }
void translate(QLineF* obj,qreal arg0,qreal arg1)  {obj->translate(arg0,arg1); }
QLine toLine(QLineF* obj)  const  {return obj->toLine(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qpoint.h'
**
** Created: Thu 12. Apr 14:07:28 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qpoint.h"
class PythonQtQPointWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QPoint() { return QPoint(); }
QVariant new_QPoint(int arg0,int arg1) { return QPoint(arg0,arg1); }
bool isNull(QPoint* obj)  const  {return obj->isNull(); }
int x(QPoint* obj)  const  {return obj->x(); }
int y(QPoint* obj)  const  {return obj->y(); }
void setX(QPoint* obj,int arg0)  {obj->setX(arg0); }
void setY(QPoint* obj,int arg0)  {obj->setY(arg0); }
int manhattanLength(QPoint* obj)  const  {return obj->manhattanLength(); }
void rx(QPoint* obj)  {obj->rx(); }
void ry(QPoint* obj)  {obj->ry(); }

};

class PythonQtQPointFWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QPointF() { return QPointF(); }
QVariant new_QPointF(const QPoint & arg0) { return QPointF(arg0); }
QVariant new_QPointF(qreal arg0,qreal arg1) { return QPointF(arg0,arg1); }
bool isNull(QPointF* obj)  const  {return obj->isNull(); }
qreal x(QPointF* obj)  const  {return obj->x(); }
qreal y(QPointF* obj)  const  {return obj->y(); }
void setX(QPointF* obj,qreal arg0)  {obj->setX(arg0); }
void setY(QPointF* obj,qreal arg0)  {obj->setY(arg0); }
void rx(QPointF* obj)  {obj->rx(); }
void ry(QPointF* obj)  {obj->ry(); }
QPoint toPoint(QPointF* obj)  const  {return obj->toPoint(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qregexp.h'
**
** Created: Thu 12. Apr 14:07:29 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qregexp.h"
class PythonQtQRegExpWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(PatternSyntax CaretMode )
enum PatternSyntax {RegExp = QRegExp::RegExp, 
Wildcard = QRegExp::Wildcard, 
FixedString = QRegExp::FixedString, 
RegExp2 = QRegExp::RegExp2 }; 
enum CaretMode {CaretAtZero = QRegExp::CaretAtZero, 
CaretAtOffset = QRegExp::CaretAtOffset, 
CaretWontMatch = QRegExp::CaretWontMatch }; 
public slots:
QVariant new_QRegExp() { return QRegExp(); }
QVariant new_QRegExp(const QString & arg0,Qt::CaseSensitivity arg1,PatternSyntax arg2) { return QRegExp(arg0,arg1,(QRegExp::PatternSyntax)arg2); }
QVariant new_QRegExp(const QString & arg0,Qt::CaseSensitivity arg1) { return QRegExp(arg0,arg1); }
QVariant new_QRegExp(const QString & arg0) { return QRegExp(arg0); }
QVariant new_QRegExp(const QRegExp & arg0) { return QRegExp(arg0); }
bool isEmpty(QRegExp* obj)  const  {return obj->isEmpty(); }
bool isValid(QRegExp* obj)  const  {return obj->isValid(); }
QString pattern(QRegExp* obj)  const  {return obj->pattern(); }
void setPattern(QRegExp* obj,const QString & arg0)  {obj->setPattern(arg0); }
Qt::CaseSensitivity caseSensitivity(QRegExp* obj)  const  {return obj->caseSensitivity(); }
void setCaseSensitivity(QRegExp* obj,Qt::CaseSensitivity arg0)  {obj->setCaseSensitivity(arg0); }
PatternSyntax patternSyntax(QRegExp* obj)  const  {return (PythonQtQRegExpWrapper::PatternSyntax)obj->patternSyntax(); }
void setPatternSyntax(QRegExp* obj,PatternSyntax arg0)  {obj->setPatternSyntax((QRegExp::PatternSyntax)arg0); }
bool isMinimal(QRegExp* obj)  const  {return obj->isMinimal(); }
void setMinimal(QRegExp* obj,bool arg0)  {obj->setMinimal(arg0); }
bool exactMatch(QRegExp* obj,const QString & arg0)  const  {return obj->exactMatch(arg0); }
int indexIn(QRegExp* obj,const QString & arg0,int arg1,CaretMode arg2)  const  {return obj->indexIn(arg0,arg1,(QRegExp::CaretMode)arg2); }
int indexIn(QRegExp* obj,const QString & arg0,int arg1)  const  {return obj->indexIn(arg0,arg1); }
int indexIn(QRegExp* obj,const QString & arg0)  const  {return obj->indexIn(arg0); }
int lastIndexIn(QRegExp* obj,const QString & arg0,int arg1,CaretMode arg2)  const  {return obj->lastIndexIn(arg0,arg1,(QRegExp::CaretMode)arg2); }
int lastIndexIn(QRegExp* obj,const QString & arg0,int arg1)  const  {return obj->lastIndexIn(arg0,arg1); }
int lastIndexIn(QRegExp* obj,const QString & arg0)  const  {return obj->lastIndexIn(arg0); }
int matchedLength(QRegExp* obj)  const  {return obj->matchedLength(); }
int numCaptures(QRegExp* obj)  const  {return obj->numCaptures(); }
QStringList capturedTexts(QRegExp* obj)  {return obj->capturedTexts(); }
QString cap(QRegExp* obj,int arg0)  {return obj->cap(arg0); }
QString cap(QRegExp* obj)  {return obj->cap(); }
int pos(QRegExp* obj,int arg0)  {return obj->pos(arg0); }
int pos(QRegExp* obj)  {return obj->pos(); }
QString errorString(QRegExp* obj)  {return obj->errorString(); }
QString static_QRegExp_escape(const QString & arg0)  {return QRegExp::escape(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qfont.h'
**
** Created: Thu 12. Apr 14:07:29 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qfont.h"
class PythonQtQFontWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(StyleHint StyleStrategy Weight Style Stretch )
enum StyleHint {Helvetica = QFont::Helvetica, 
SansSerif = QFont::SansSerif, 
Times = QFont::Times, 
Serif = QFont::Serif, 
Courier = QFont::Courier, 
TypeWriter = QFont::TypeWriter, 
OldEnglish = QFont::OldEnglish, 
Decorative = QFont::Decorative, 
System = QFont::System, 
AnyStyle = QFont::AnyStyle }; 
enum StyleStrategy {PreferDefault = QFont::PreferDefault, 
PreferBitmap = QFont::PreferBitmap, 
PreferDevice = QFont::PreferDevice, 
PreferOutline = QFont::PreferOutline, 
ForceOutline = QFont::ForceOutline, 
PreferMatch = QFont::PreferMatch, 
PreferQuality = QFont::PreferQuality, 
PreferAntialias = QFont::PreferAntialias, 
NoAntialias = QFont::NoAntialias, 
OpenGLCompatible = QFont::OpenGLCompatible, 
NoFontMerging = QFont::NoFontMerging }; 
enum Weight {Light = QFont::Light, 
Normal = QFont::Normal, 
DemiBold = QFont::DemiBold, 
Bold = QFont::Bold, 
Black = QFont::Black }; 
enum Style {StyleNormal = QFont::StyleNormal, 
StyleItalic = QFont::StyleItalic, 
StyleOblique = QFont::StyleOblique }; 
enum Stretch {UltraCondensed = QFont::UltraCondensed, 
ExtraCondensed = QFont::ExtraCondensed, 
Condensed = QFont::Condensed, 
SemiCondensed = QFont::SemiCondensed, 
Unstretched = QFont::Unstretched, 
SemiExpanded = QFont::SemiExpanded, 
Expanded = QFont::Expanded, 
ExtraExpanded = QFont::ExtraExpanded, 
UltraExpanded = QFont::UltraExpanded }; 
public slots:
QVariant new_QFont() { return QFont(); }
QVariant new_QFont(const QString & arg0,int arg1,int arg2,bool arg3) { return QFont(arg0,arg1,arg2,arg3); }
QVariant new_QFont(const QString & arg0,int arg1,int arg2) { return QFont(arg0,arg1,arg2); }
QVariant new_QFont(const QString & arg0,int arg1) { return QFont(arg0,arg1); }
QVariant new_QFont(const QString & arg0) { return QFont(arg0); }
QVariant new_QFont(const QFont & arg0,QPaintDevice * arg1) { return QFont(arg0,arg1); }
QVariant new_QFont(const QFont & arg0) { return QFont(arg0); }
QString family(QFont* obj)  const  {return obj->family(); }
void setFamily(QFont* obj,const QString & arg0)  {obj->setFamily(arg0); }
int pointSize(QFont* obj)  const  {return obj->pointSize(); }
void setPointSize(QFont* obj,int arg0)  {obj->setPointSize(arg0); }
qreal pointSizeF(QFont* obj)  const  {return obj->pointSizeF(); }
void setPointSizeF(QFont* obj,qreal arg0)  {obj->setPointSizeF(arg0); }
int pixelSize(QFont* obj)  const  {return obj->pixelSize(); }
void setPixelSize(QFont* obj,int arg0)  {obj->setPixelSize(arg0); }
int weight(QFont* obj)  const  {return obj->weight(); }
void setWeight(QFont* obj,int arg0)  {obj->setWeight(arg0); }
bool bold(QFont* obj)  const  {return obj->bold(); }
void setBold(QFont* obj,bool arg0)  {obj->setBold(arg0); }
void setStyle(QFont* obj,Style arg0)  {obj->setStyle((QFont::Style)arg0); }
Style style(QFont* obj)  const  {return (PythonQtQFontWrapper::Style)obj->style(); }
bool italic(QFont* obj)  const  {return obj->italic(); }
void setItalic(QFont* obj,bool arg0)  {obj->setItalic(arg0); }
bool underline(QFont* obj)  const  {return obj->underline(); }
void setUnderline(QFont* obj,bool arg0)  {obj->setUnderline(arg0); }
bool overline(QFont* obj)  const  {return obj->overline(); }
void setOverline(QFont* obj,bool arg0)  {obj->setOverline(arg0); }
bool strikeOut(QFont* obj)  const  {return obj->strikeOut(); }
void setStrikeOut(QFont* obj,bool arg0)  {obj->setStrikeOut(arg0); }
bool fixedPitch(QFont* obj)  const  {return obj->fixedPitch(); }
void setFixedPitch(QFont* obj,bool arg0)  {obj->setFixedPitch(arg0); }
bool kerning(QFont* obj)  const  {return obj->kerning(); }
void setKerning(QFont* obj,bool arg0)  {obj->setKerning(arg0); }
StyleHint styleHint(QFont* obj)  const  {return (PythonQtQFontWrapper::StyleHint)obj->styleHint(); }
StyleStrategy styleStrategy(QFont* obj)  const  {return (PythonQtQFontWrapper::StyleStrategy)obj->styleStrategy(); }
void setStyleHint(QFont* obj,StyleHint arg0,StyleStrategy arg1)  {obj->setStyleHint((QFont::StyleHint)arg0,(QFont::StyleStrategy)arg1); }
void setStyleHint(QFont* obj,StyleHint arg0)  {obj->setStyleHint((QFont::StyleHint)arg0); }
void setStyleStrategy(QFont* obj,StyleStrategy arg0)  {obj->setStyleStrategy((QFont::StyleStrategy)arg0); }
int stretch(QFont* obj)  const  {return obj->stretch(); }
void setStretch(QFont* obj,int arg0)  {obj->setStretch(arg0); }
bool rawMode(QFont* obj)  const  {return obj->rawMode(); }
void setRawMode(QFont* obj,bool arg0)  {obj->setRawMode(arg0); }
bool exactMatch(QFont* obj)  const  {return obj->exactMatch(); }
bool isCopyOf(QFont* obj,const QFont & arg0)  const  {return obj->isCopyOf(arg0); }
Qt::HANDLE handle(QFont* obj)  const  {return obj->handle(); }
void setRawName(QFont* obj,const QString & arg0)  {obj->setRawName(arg0); }
QString rawName(QFont* obj)  const  {return obj->rawName(); }
QString key(QFont* obj)  const  {return obj->key(); }
QString toString(QFont* obj)  const  {return obj->toString(); }
bool fromString(QFont* obj,const QString & arg0)  {return obj->fromString(arg0); }
QString static_QFont_substitute(const QString & arg0)  {return QFont::substitute(arg0); }
QStringList static_QFont_substitutes(const QString & arg0)  {return QFont::substitutes(arg0); }
QStringList static_QFont_substitutions()  {return QFont::substitutions(); }
void static_QFont_insertSubstitution(const QString & arg0,const QString & arg1)  {QFont::insertSubstitution(arg0,arg1); }
void static_QFont_insertSubstitutions(const QString & arg0,const QStringList & arg1)  {QFont::insertSubstitutions(arg0,arg1); }
void static_QFont_removeSubstitution(const QString & arg0)  {QFont::removeSubstitution(arg0); }
void static_QFont_initialize()  {QFont::initialize(); }
void static_QFont_cleanup()  {QFont::cleanup(); }
void static_QFont_cacheStatistics()  {QFont::cacheStatistics(); }
QString defaultFamily(QFont* obj)  const  {return obj->defaultFamily(); }
QString lastResortFamily(QFont* obj)  const  {return obj->lastResortFamily(); }
QString lastResortFont(QFont* obj)  const  {return obj->lastResortFont(); }
QFont resolve(QFont* obj,const QFont & arg0)  const  {return obj->resolve(arg0); }
uint resolve(QFont* obj)  const  {return obj->resolve(); }
void resolve(QFont* obj,uint arg0)  {obj->resolve(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qpixmap.h'
**
** Created: Thu 12. Apr 14:07:29 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qpixmap.h"
class PythonQtQPixmapWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QPixmap() { return QPixmap(); }
QVariant new_QPixmap(int arg0,int arg1) { return QPixmap(arg0,arg1); }
QVariant new_QPixmap(const QSize & arg0) { return QPixmap(arg0); }
QVariant new_QPixmap(const QString & arg0,const char * arg1,Qt::ImageConversionFlags arg2) { return QPixmap(arg0,arg1,arg2); }
QVariant new_QPixmap(const QString & arg0,const char * arg1) { return QPixmap(arg0,arg1); }
QVariant new_QPixmap(const QString & arg0) { return QPixmap(arg0); }
QVariant new_QPixmap(const QPixmap & arg0) { return QPixmap(arg0); }
bool isNull(QPixmap* obj)  const  {return obj->isNull(); }
int devType(QPixmap* obj)  const  {return obj->devType(); }
int width(QPixmap* obj)  const  {return obj->width(); }
int height(QPixmap* obj)  const  {return obj->height(); }
QSize size(QPixmap* obj)  const  {return obj->size(); }
QRect rect(QPixmap* obj)  const  {return obj->rect(); }
int depth(QPixmap* obj)  const  {return obj->depth(); }
int static_QPixmap_defaultDepth()  {return QPixmap::defaultDepth(); }
void fill(QPixmap* obj,const QColor & arg0)  {obj->fill(arg0); }
void fill(QPixmap* obj)  {obj->fill(); }
void fill(QPixmap* obj,const QWidget * arg0,const QPoint & arg1)  {obj->fill(arg0,arg1); }
void fill(QPixmap* obj,const QWidget * arg0,int arg1,int arg2)  {obj->fill(arg0,arg1,arg2); }
QBitmap mask(QPixmap* obj)  const  {return obj->mask(); }
void setMask(QPixmap* obj,const QBitmap & arg0)  {obj->setMask(arg0); }
QPixmap alphaChannel(QPixmap* obj)  const  {return obj->alphaChannel(); }
void setAlphaChannel(QPixmap* obj,const QPixmap & arg0)  {obj->setAlphaChannel(arg0); }
bool hasAlpha(QPixmap* obj)  const  {return obj->hasAlpha(); }
bool hasAlphaChannel(QPixmap* obj)  const  {return obj->hasAlphaChannel(); }
QBitmap createHeuristicMask(QPixmap* obj,bool arg0)  const  {return obj->createHeuristicMask(arg0); }
QBitmap createHeuristicMask(QPixmap* obj)  const  {return obj->createHeuristicMask(); }
QBitmap createMaskFromColor(QPixmap* obj,const QColor & arg0)  const  {return obj->createMaskFromColor(arg0); }
QPixmap static_QPixmap_grabWindow(WId arg0,int arg1,int arg2,int arg3,int arg4)  {return QPixmap::grabWindow(arg0,arg1,arg2,arg3,arg4); }
QPixmap static_QPixmap_grabWindow(WId arg0,int arg1,int arg2,int arg3)  {return QPixmap::grabWindow(arg0,arg1,arg2,arg3); }
QPixmap static_QPixmap_grabWindow(WId arg0,int arg1,int arg2)  {return QPixmap::grabWindow(arg0,arg1,arg2); }
QPixmap static_QPixmap_grabWindow(WId arg0,int arg1)  {return QPixmap::grabWindow(arg0,arg1); }
QPixmap static_QPixmap_grabWindow(WId arg0)  {return QPixmap::grabWindow(arg0); }
QPixmap static_QPixmap_grabWidget(QWidget * arg0,const QRect & arg1)  {return QPixmap::grabWidget(arg0,arg1); }
QPixmap grabWidget(QPixmap* obj,QWidget * arg0,int arg1,int arg2,int arg3,int arg4)  {return obj->grabWidget(arg0,arg1,arg2,arg3,arg4); }
QPixmap grabWidget(QPixmap* obj,QWidget * arg0,int arg1,int arg2,int arg3)  {return obj->grabWidget(arg0,arg1,arg2,arg3); }
QPixmap grabWidget(QPixmap* obj,QWidget * arg0,int arg1,int arg2)  {return obj->grabWidget(arg0,arg1,arg2); }
QPixmap grabWidget(QPixmap* obj,QWidget * arg0,int arg1)  {return obj->grabWidget(arg0,arg1); }
QPixmap grabWidget(QPixmap* obj,QWidget * arg0)  {return obj->grabWidget(arg0); }
QPixmap scaled(QPixmap* obj,int arg0,int arg1,Qt::AspectRatioMode arg2,Qt::TransformationMode arg3)  const  {return obj->scaled(arg0,arg1,arg2,arg3); }
QPixmap scaled(QPixmap* obj,int arg0,int arg1,Qt::AspectRatioMode arg2)  const  {return obj->scaled(arg0,arg1,arg2); }
QPixmap scaled(QPixmap* obj,int arg0,int arg1)  const  {return obj->scaled(arg0,arg1); }
QPixmap scaled(QPixmap* obj,const QSize & arg0,Qt::AspectRatioMode arg1,Qt::TransformationMode arg2)  const  {return obj->scaled(arg0,arg1,arg2); }
QPixmap scaled(QPixmap* obj,const QSize & arg0,Qt::AspectRatioMode arg1)  const  {return obj->scaled(arg0,arg1); }
QPixmap scaled(QPixmap* obj,const QSize & arg0)  const  {return obj->scaled(arg0); }
QPixmap scaledToWidth(QPixmap* obj,int arg0,Qt::TransformationMode arg1)  const  {return obj->scaledToWidth(arg0,arg1); }
QPixmap scaledToWidth(QPixmap* obj,int arg0)  const  {return obj->scaledToWidth(arg0); }
QPixmap scaledToHeight(QPixmap* obj,int arg0,Qt::TransformationMode arg1)  const  {return obj->scaledToHeight(arg0,arg1); }
QPixmap scaledToHeight(QPixmap* obj,int arg0)  const  {return obj->scaledToHeight(arg0); }
QPixmap transformed(QPixmap* obj,const QMatrix & arg0,Qt::TransformationMode arg1)  const  {return obj->transformed(arg0,arg1); }
QPixmap transformed(QPixmap* obj,const QMatrix & arg0)  const  {return obj->transformed(arg0); }
QMatrix static_QPixmap_trueMatrix(const QMatrix & arg0,int arg1,int arg2)  {return QPixmap::trueMatrix(arg0,arg1,arg2); }
QImage toImage(QPixmap* obj)  const  {return obj->toImage(); }
QPixmap static_QPixmap_fromImage(const QImage & arg0,Qt::ImageConversionFlags arg1)  {return QPixmap::fromImage(arg0,arg1); }
QPixmap static_QPixmap_fromImage(const QImage & arg0)  {return QPixmap::fromImage(arg0); }
bool load(QPixmap* obj,const QString & arg0,const char * arg1,Qt::ImageConversionFlags arg2)  {return obj->load(arg0,arg1,arg2); }
bool load(QPixmap* obj,const QString & arg0,const char * arg1)  {return obj->load(arg0,arg1); }
bool load(QPixmap* obj,const QString & arg0)  {return obj->load(arg0); }
bool loadFromData(QPixmap* obj,const uchar * arg0,uint arg1,const char * arg2,Qt::ImageConversionFlags arg3)  {return obj->loadFromData(arg0,arg1,arg2,arg3); }
bool loadFromData(QPixmap* obj,const uchar * arg0,uint arg1,const char * arg2)  {return obj->loadFromData(arg0,arg1,arg2); }
bool loadFromData(QPixmap* obj,const uchar * arg0,uint arg1)  {return obj->loadFromData(arg0,arg1); }
bool loadFromData(QPixmap* obj,const QByteArray & arg0,const char * arg1,Qt::ImageConversionFlags arg2)  {return obj->loadFromData(arg0,arg1,arg2); }
bool loadFromData(QPixmap* obj,const QByteArray & arg0,const char * arg1)  {return obj->loadFromData(arg0,arg1); }
bool loadFromData(QPixmap* obj,const QByteArray & arg0)  {return obj->loadFromData(arg0); }
bool save(QPixmap* obj,const QString & arg0,const char * arg1,int arg2)  const  {return obj->save(arg0,arg1,arg2); }
bool save(QPixmap* obj,const QString & arg0,const char * arg1)  const  {return obj->save(arg0,arg1); }
bool save(QPixmap* obj,const QString & arg0)  const  {return obj->save(arg0); }
bool save(QPixmap* obj,QIODevice * arg0,const char * arg1,int arg2)  const  {return obj->save(arg0,arg1,arg2); }
bool save(QPixmap* obj,QIODevice * arg0,const char * arg1)  const  {return obj->save(arg0,arg1); }
bool save(QPixmap* obj,QIODevice * arg0)  const  {return obj->save(arg0); }
QPixmap copy(QPixmap* obj,int arg0,int arg1,int arg2,int arg3)  const  {return obj->copy(arg0,arg1,arg2,arg3); }
QPixmap copy(QPixmap* obj,const QRect & arg0)  const  {return obj->copy(arg0); }
QPixmap copy(QPixmap* obj)  const  {return obj->copy(); }
int serialNumber(QPixmap* obj)  const  {return obj->serialNumber(); }
bool isDetached(QPixmap* obj)  const  {return obj->isDetached(); }
void detach(QPixmap* obj)  {obj->detach(); }
bool isQBitmap(QPixmap* obj)  const  {return obj->isQBitmap(); }
QPaintEngine* paintEngine(QPixmap* obj)  const  {return obj->paintEngine(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qbrush.h'
**
** Created: Thu 12. Apr 14:07:29 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qbrush.h"
class PythonQtQBrushWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QBrush() { return QBrush(); }
QVariant new_QBrush(Qt::BrushStyle arg0) { return QBrush(arg0); }
QVariant new_QBrush(const QColor & arg0,Qt::BrushStyle arg1) { return QBrush(arg0,arg1); }
QVariant new_QBrush(const QColor & arg0) { return QBrush(arg0); }
QVariant new_QBrush(Qt::GlobalColor arg0,Qt::BrushStyle arg1) { return QBrush(arg0,arg1); }
QVariant new_QBrush(Qt::GlobalColor arg0) { return QBrush(arg0); }
QVariant new_QBrush(const QColor & arg0,const QPixmap & arg1) { return QBrush(arg0,arg1); }
QVariant new_QBrush(Qt::GlobalColor arg0,const QPixmap & arg1) { return QBrush(arg0,arg1); }
QVariant new_QBrush(const QPixmap & arg0) { return QBrush(arg0); }
QVariant new_QBrush(const QImage & arg0) { return QBrush(arg0); }
QVariant new_QBrush(const QBrush & arg0) { return QBrush(arg0); }
QVariant new_QBrush(const QGradient & arg0) { return QBrush(arg0); }
Qt::BrushStyle style(QBrush* obj)  const  {return obj->style(); }
void setStyle(QBrush* obj,Qt::BrushStyle arg0)  {obj->setStyle(arg0); }
void matrix(QBrush* obj)  const  {obj->matrix(); }
void setMatrix(QBrush* obj,const QMatrix & arg0)  {obj->setMatrix(arg0); }
QPixmap texture(QBrush* obj)  const  {return obj->texture(); }
void setTexture(QBrush* obj,const QPixmap & arg0)  {obj->setTexture(arg0); }
QImage textureImage(QBrush* obj)  const  {return obj->textureImage(); }
void setTextureImage(QBrush* obj,const QImage & arg0)  {obj->setTextureImage(arg0); }
void color(QBrush* obj)  const  {obj->color(); }
void setColor(QBrush* obj,const QColor & arg0)  {obj->setColor(arg0); }
void setColor(QBrush* obj,Qt::GlobalColor arg0)  {obj->setColor(arg0); }
const QGradient* gradient(QBrush* obj)  const  {return obj->gradient(); }
bool isOpaque(QBrush* obj)  const  {return obj->isOpaque(); }

};

class PythonQtQBrushDataWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QBrushData(QBrushData* obj) { delete obj; }

};

class PythonQtQGradientWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(Type Spread CoordinateMode )
enum Type {LinearGradient = QGradient::LinearGradient, 
RadialGradient = QGradient::RadialGradient, 
ConicalGradient = QGradient::ConicalGradient, 
NoGradient = QGradient::NoGradient }; 
enum Spread {PadSpread = QGradient::PadSpread, 
ReflectSpread = QGradient::ReflectSpread, 
RepeatSpread = QGradient::RepeatSpread }; 
enum CoordinateMode {LogicalMode = QGradient::LogicalMode, 
StretchToDeviceMode = QGradient::StretchToDeviceMode }; 
public slots:
void delete_QGradient(QGradient* obj) { delete obj; }
QGradient* new_QGradient() { return new QGradient(); }
Type type(QGradient* obj)  const  {return (PythonQtQGradientWrapper::Type)obj->type(); }
void setSpread(QGradient* obj,Spread arg0)  {obj->setSpread((QGradient::Spread)arg0); }
Spread spread(QGradient* obj)  const  {return (PythonQtQGradientWrapper::Spread)obj->spread(); }
void setColorAt(QGradient* obj,qreal arg0,const QColor & arg1)  {obj->setColorAt(arg0,arg1); }
void setStops(QGradient* obj,const QGradientStops & arg0)  {obj->setStops(arg0); }
QGradientStops stops(QGradient* obj)  const  {return obj->stops(); }
CoordinateMode coordinateMode(QGradient* obj)  const  {return (PythonQtQGradientWrapper::CoordinateMode)obj->coordinateMode(); }
void setCoordinateMode(QGradient* obj,CoordinateMode arg0)  {obj->setCoordinateMode((QGradient::CoordinateMode)arg0); }

};

class PythonQtQLinearGradientWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QLinearGradient(QLinearGradient* obj) { delete obj; }
QLinearGradient* new_QLinearGradient() { return new QLinearGradient(); }
QLinearGradient* new_QLinearGradient(const QPointF & arg0,const QPointF & arg1) { return new QLinearGradient(arg0,arg1); }
QLinearGradient* new_QLinearGradient(qreal arg0,qreal arg1,qreal arg2,qreal arg3) { return new QLinearGradient(arg0,arg1,arg2,arg3); }
QPointF start(QLinearGradient* obj)  const  {return obj->start(); }
void setStart(QLinearGradient* obj,const QPointF & arg0)  {obj->setStart(arg0); }
void setStart(QLinearGradient* obj,qreal arg0,qreal arg1)  {obj->setStart(arg0,arg1); }
QPointF finalStop(QLinearGradient* obj)  const  {return obj->finalStop(); }
void setFinalStop(QLinearGradient* obj,const QPointF & arg0)  {obj->setFinalStop(arg0); }
void setFinalStop(QLinearGradient* obj,qreal arg0,qreal arg1)  {obj->setFinalStop(arg0,arg1); }

};

class PythonQtQRadialGradientWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QRadialGradient(QRadialGradient* obj) { delete obj; }
QRadialGradient* new_QRadialGradient() { return new QRadialGradient(); }
QRadialGradient* new_QRadialGradient(const QPointF & arg0,qreal arg1,const QPointF & arg2) { return new QRadialGradient(arg0,arg1,arg2); }
QRadialGradient* new_QRadialGradient(qreal arg0,qreal arg1,qreal arg2,qreal arg3,qreal arg4) { return new QRadialGradient(arg0,arg1,arg2,arg3,arg4); }
QRadialGradient* new_QRadialGradient(const QPointF & arg0,qreal arg1) { return new QRadialGradient(arg0,arg1); }
QRadialGradient* new_QRadialGradient(qreal arg0,qreal arg1,qreal arg2) { return new QRadialGradient(arg0,arg1,arg2); }
QPointF center(QRadialGradient* obj)  const  {return obj->center(); }
void setCenter(QRadialGradient* obj,const QPointF & arg0)  {obj->setCenter(arg0); }
void setCenter(QRadialGradient* obj,qreal arg0,qreal arg1)  {obj->setCenter(arg0,arg1); }
QPointF focalPoint(QRadialGradient* obj)  const  {return obj->focalPoint(); }
void setFocalPoint(QRadialGradient* obj,const QPointF & arg0)  {obj->setFocalPoint(arg0); }
void setFocalPoint(QRadialGradient* obj,qreal arg0,qreal arg1)  {obj->setFocalPoint(arg0,arg1); }
qreal radius(QRadialGradient* obj)  const  {return obj->radius(); }
void setRadius(QRadialGradient* obj,qreal arg0)  {obj->setRadius(arg0); }

};

class PythonQtQConicalGradientWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QConicalGradient(QConicalGradient* obj) { delete obj; }
QConicalGradient* new_QConicalGradient() { return new QConicalGradient(); }
QConicalGradient* new_QConicalGradient(const QPointF & arg0,qreal arg1) { return new QConicalGradient(arg0,arg1); }
QConicalGradient* new_QConicalGradient(qreal arg0,qreal arg1,qreal arg2) { return new QConicalGradient(arg0,arg1,arg2); }
QPointF center(QConicalGradient* obj)  const  {return obj->center(); }
void setCenter(QConicalGradient* obj,const QPointF & arg0)  {obj->setCenter(arg0); }
void setCenter(QConicalGradient* obj,qreal arg0,qreal arg1)  {obj->setCenter(arg0,arg1); }
qreal angle(QConicalGradient* obj)  const  {return obj->angle(); }
void setAngle(QConicalGradient* obj,qreal arg0)  {obj->setAngle(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qcolor.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qcolor.h"
class PythonQtQColorWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(Spec )
enum Spec {Invalid = QColor::Invalid, 
Rgb = QColor::Rgb, 
Hsv = QColor::Hsv, 
Cmyk = QColor::Cmyk }; 
public slots:
QVariant new_QColor() { return QColor(); }
QVariant new_QColor(Qt::GlobalColor arg0) { return QColor(arg0); }
QVariant new_QColor(int arg0,int arg1,int arg2,int arg3) { return QColor(arg0,arg1,arg2,arg3); }
QVariant new_QColor(int arg0,int arg1,int arg2) { return QColor(arg0,arg1,arg2); }
QVariant new_QColor(QRgb arg0) { return QColor(arg0); }
QVariant new_QColor(const QString & arg0) { return QColor(arg0); }
QVariant new_QColor(const char * arg0) { return QColor(arg0); }
QVariant new_QColor(const QColor & arg0) { return QColor(arg0); }
QVariant new_QColor(Spec arg0) { return QColor((QColor::Spec)arg0); }
bool isValid(QColor* obj)  const  {return obj->isValid(); }
QString name(QColor* obj)  const  {return obj->name(); }
void setNamedColor(QColor* obj,const QString & arg0)  {obj->setNamedColor(arg0); }
QStringList static_QColor_colorNames()  {return QColor::colorNames(); }
Spec spec(QColor* obj)  const  {return (PythonQtQColorWrapper::Spec)obj->spec(); }
int alpha(QColor* obj)  const  {return obj->alpha(); }
void setAlpha(QColor* obj,int arg0)  {obj->setAlpha(arg0); }
qreal alphaF(QColor* obj)  const  {return obj->alphaF(); }
void setAlphaF(QColor* obj,qreal arg0)  {obj->setAlphaF(arg0); }
int red(QColor* obj)  const  {return obj->red(); }
int green(QColor* obj)  const  {return obj->green(); }
int blue(QColor* obj)  const  {return obj->blue(); }
void setRed(QColor* obj,int arg0)  {obj->setRed(arg0); }
void setGreen(QColor* obj,int arg0)  {obj->setGreen(arg0); }
void setBlue(QColor* obj,int arg0)  {obj->setBlue(arg0); }
qreal redF(QColor* obj)  const  {return obj->redF(); }
qreal greenF(QColor* obj)  const  {return obj->greenF(); }
qreal blueF(QColor* obj)  const  {return obj->blueF(); }
void setRedF(QColor* obj,qreal arg0)  {obj->setRedF(arg0); }
void setGreenF(QColor* obj,qreal arg0)  {obj->setGreenF(arg0); }
void setBlueF(QColor* obj,qreal arg0)  {obj->setBlueF(arg0); }
void getRgb(QColor* obj,int * arg0,int * arg1,int * arg2,int * arg3)  const  {obj->getRgb(arg0,arg1,arg2,arg3); }
void getRgb(QColor* obj,int * arg0,int * arg1,int * arg2)  const  {obj->getRgb(arg0,arg1,arg2); }
void setRgb(QColor* obj,int arg0,int arg1,int arg2,int arg3)  {obj->setRgb(arg0,arg1,arg2,arg3); }
void setRgb(QColor* obj,int arg0,int arg1,int arg2)  {obj->setRgb(arg0,arg1,arg2); }
void getRgbF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3)  const  {obj->getRgbF(arg0,arg1,arg2,arg3); }
void getRgbF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2)  const  {obj->getRgbF(arg0,arg1,arg2); }
void setRgbF(QColor* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->setRgbF(arg0,arg1,arg2,arg3); }
void setRgbF(QColor* obj,qreal arg0,qreal arg1,qreal arg2)  {obj->setRgbF(arg0,arg1,arg2); }
QRgb rgba(QColor* obj)  const  {return obj->rgba(); }
void setRgba(QColor* obj,QRgb arg0)  {obj->setRgba(arg0); }
QRgb rgb(QColor* obj)  const  {return obj->rgb(); }
void setRgb(QColor* obj,QRgb arg0)  {obj->setRgb(arg0); }
int hue(QColor* obj)  const  {return obj->hue(); }
int saturation(QColor* obj)  const  {return obj->saturation(); }
int value(QColor* obj)  const  {return obj->value(); }
qreal hueF(QColor* obj)  const  {return obj->hueF(); }
qreal saturationF(QColor* obj)  const  {return obj->saturationF(); }
qreal valueF(QColor* obj)  const  {return obj->valueF(); }
void getHsv(QColor* obj,int * arg0,int * arg1,int * arg2,int * arg3)  const  {obj->getHsv(arg0,arg1,arg2,arg3); }
void getHsv(QColor* obj,int * arg0,int * arg1,int * arg2)  const  {obj->getHsv(arg0,arg1,arg2); }
void setHsv(QColor* obj,int arg0,int arg1,int arg2,int arg3)  {obj->setHsv(arg0,arg1,arg2,arg3); }
void setHsv(QColor* obj,int arg0,int arg1,int arg2)  {obj->setHsv(arg0,arg1,arg2); }
void getHsvF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3)  const  {obj->getHsvF(arg0,arg1,arg2,arg3); }
void getHsvF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2)  const  {obj->getHsvF(arg0,arg1,arg2); }
void setHsvF(QColor* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->setHsvF(arg0,arg1,arg2,arg3); }
void setHsvF(QColor* obj,qreal arg0,qreal arg1,qreal arg2)  {obj->setHsvF(arg0,arg1,arg2); }
int cyan(QColor* obj)  const  {return obj->cyan(); }
int magenta(QColor* obj)  const  {return obj->magenta(); }
int yellow(QColor* obj)  const  {return obj->yellow(); }
int black(QColor* obj)  const  {return obj->black(); }
qreal cyanF(QColor* obj)  const  {return obj->cyanF(); }
qreal magentaF(QColor* obj)  const  {return obj->magentaF(); }
qreal yellowF(QColor* obj)  const  {return obj->yellowF(); }
qreal blackF(QColor* obj)  const  {return obj->blackF(); }
void getCmyk(QColor* obj,int * arg0,int * arg1,int * arg2,int * arg3,int * arg4)  {obj->getCmyk(arg0,arg1,arg2,arg3,arg4); }
void getCmyk(QColor* obj,int * arg0,int * arg1,int * arg2,int * arg3)  {obj->getCmyk(arg0,arg1,arg2,arg3); }
void setCmyk(QColor* obj,int arg0,int arg1,int arg2,int arg3,int arg4)  {obj->setCmyk(arg0,arg1,arg2,arg3,arg4); }
void setCmyk(QColor* obj,int arg0,int arg1,int arg2,int arg3)  {obj->setCmyk(arg0,arg1,arg2,arg3); }
void getCmykF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3,qreal * arg4)  {obj->getCmykF(arg0,arg1,arg2,arg3,arg4); }
void getCmykF(QColor* obj,qreal * arg0,qreal * arg1,qreal * arg2,qreal * arg3)  {obj->getCmykF(arg0,arg1,arg2,arg3); }
void setCmykF(QColor* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3,qreal arg4)  {obj->setCmykF(arg0,arg1,arg2,arg3,arg4); }
void setCmykF(QColor* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {obj->setCmykF(arg0,arg1,arg2,arg3); }
QColor toRgb(QColor* obj)  const  {return obj->toRgb(); }
QColor toHsv(QColor* obj)  const  {return obj->toHsv(); }
QColor toCmyk(QColor* obj)  const  {return obj->toCmyk(); }
QColor convertTo(QColor* obj,Spec arg0)  const  {return obj->convertTo((QColor::Spec)arg0); }
QColor static_QColor_fromRgb(QRgb arg0)  {return QColor::fromRgb(arg0); }
QColor static_QColor_fromRgba(QRgb arg0)  {return QColor::fromRgba(arg0); }
QColor static_QColor_fromRgb(int arg0,int arg1,int arg2,int arg3)  {return QColor::fromRgb(arg0,arg1,arg2,arg3); }
QColor static_QColor_fromRgb(int arg0,int arg1,int arg2)  {return QColor::fromRgb(arg0,arg1,arg2); }
QColor static_QColor_fromRgbF(qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {return QColor::fromRgbF(arg0,arg1,arg2,arg3); }
QColor static_QColor_fromRgbF(qreal arg0,qreal arg1,qreal arg2)  {return QColor::fromRgbF(arg0,arg1,arg2); }
QColor static_QColor_fromHsv(int arg0,int arg1,int arg2,int arg3)  {return QColor::fromHsv(arg0,arg1,arg2,arg3); }
QColor static_QColor_fromHsv(int arg0,int arg1,int arg2)  {return QColor::fromHsv(arg0,arg1,arg2); }
QColor static_QColor_fromHsvF(qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {return QColor::fromHsvF(arg0,arg1,arg2,arg3); }
QColor static_QColor_fromHsvF(qreal arg0,qreal arg1,qreal arg2)  {return QColor::fromHsvF(arg0,arg1,arg2); }
QColor static_QColor_fromCmyk(int arg0,int arg1,int arg2,int arg3,int arg4)  {return QColor::fromCmyk(arg0,arg1,arg2,arg3,arg4); }
QColor static_QColor_fromCmyk(int arg0,int arg1,int arg2,int arg3)  {return QColor::fromCmyk(arg0,arg1,arg2,arg3); }
QColor static_QColor_fromCmykF(qreal arg0,qreal arg1,qreal arg2,qreal arg3,qreal arg4)  {return QColor::fromCmykF(arg0,arg1,arg2,arg3,arg4); }
QColor static_QColor_fromCmykF(qreal arg0,qreal arg1,qreal arg2,qreal arg3)  {return QColor::fromCmykF(arg0,arg1,arg2,arg3); }
QColor light(QColor* obj,int arg0)  const  {return obj->light(arg0); }
QColor light(QColor* obj)  const  {return obj->light(); }
QColor dark(QColor* obj,int arg0)  const  {return obj->dark(arg0); }
QColor dark(QColor* obj)  const  {return obj->dark(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qpalette.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qpalette.h"
class PythonQtQPaletteWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(ColorGroup ColorRole )
enum ColorGroup {Active = QPalette::Active, 
Disabled = QPalette::Disabled, 
Inactive = QPalette::Inactive, 
NColorGroups = QPalette::NColorGroups, 
Current = QPalette::Current, 
All = QPalette::All, 
Normal = QPalette::Normal }; 
enum ColorRole {WindowText = QPalette::WindowText, 
Button = QPalette::Button, 
Light = QPalette::Light, 
Midlight = QPalette::Midlight, 
Dark = QPalette::Dark, 
Mid = QPalette::Mid, 
Text = QPalette::Text, 
BrightText = QPalette::BrightText, 
ButtonText = QPalette::ButtonText, 
Base = QPalette::Base, 
Window = QPalette::Window, 
Shadow = QPalette::Shadow, 
Highlight = QPalette::Highlight, 
HighlightedText = QPalette::HighlightedText, 
Link = QPalette::Link, 
LinkVisited = QPalette::LinkVisited, 
AlternateBase = QPalette::AlternateBase, 
NoRole = QPalette::NoRole, 
NColorRoles = QPalette::NColorRoles, 
Foreground = QPalette::Foreground, 
Background = QPalette::Background }; 
public slots:
QVariant new_QPalette() { return QPalette(); }
QVariant new_QPalette(const QColor & arg0) { return QPalette(arg0); }
QVariant new_QPalette(Qt::GlobalColor arg0) { return QPalette(arg0); }
QVariant new_QPalette(const QColor & arg0,const QColor & arg1) { return QPalette(arg0,arg1); }
QVariant new_QPalette(const QBrush & arg0,const QBrush & arg1,const QBrush & arg2,const QBrush & arg3,const QBrush & arg4,const QBrush & arg5,const QBrush & arg6,const QBrush & arg7,const QBrush & arg8) { return QPalette(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8); }
QVariant new_QPalette(const QColor & arg0,const QColor & arg1,const QColor & arg2,const QColor & arg3,const QColor & arg4,const QColor & arg5,const QColor & arg6) { return QPalette(arg0,arg1,arg2,arg3,arg4,arg5,arg6); }
QVariant new_QPalette(const QPalette & arg0) { return QPalette(arg0); }
ColorGroup currentColorGroup(QPalette* obj)  const  {return (PythonQtQPaletteWrapper::ColorGroup)obj->currentColorGroup(); }
void setCurrentColorGroup(QPalette* obj,ColorGroup arg0)  {obj->setCurrentColorGroup((QPalette::ColorGroup)arg0); }
void color(QPalette* obj,ColorGroup arg0,ColorRole arg1)  const  {obj->color((QPalette::ColorGroup)arg0,(QPalette::ColorRole)arg1); }
void brush(QPalette* obj,ColorGroup arg0,ColorRole arg1)  const  {obj->brush((QPalette::ColorGroup)arg0,(QPalette::ColorRole)arg1); }
void setColor(QPalette* obj,ColorGroup arg0,ColorRole arg1,const QColor & arg2)  {obj->setColor((QPalette::ColorGroup)arg0,(QPalette::ColorRole)arg1,arg2); }
void setColor(QPalette* obj,ColorRole arg0,const QColor & arg1)  {obj->setColor((QPalette::ColorRole)arg0,arg1); }
void setBrush(QPalette* obj,ColorRole arg0,const QBrush & arg1)  {obj->setBrush((QPalette::ColorRole)arg0,arg1); }
bool isBrushSet(QPalette* obj,ColorGroup arg0,ColorRole arg1)  const  {return obj->isBrushSet((QPalette::ColorGroup)arg0,(QPalette::ColorRole)arg1); }
void setBrush(QPalette* obj,ColorGroup arg0,ColorRole arg1,const QBrush & arg2)  {obj->setBrush((QPalette::ColorGroup)arg0,(QPalette::ColorRole)arg1,arg2); }
void setColorGroup(QPalette* obj,ColorGroup arg0,const QBrush & arg1,const QBrush & arg2,const QBrush & arg3,const QBrush & arg4,const QBrush & arg5,const QBrush & arg6,const QBrush & arg7,const QBrush & arg8,const QBrush & arg9)  {obj->setColorGroup((QPalette::ColorGroup)arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9); }
bool isEqual(QPalette* obj,ColorGroup arg0,ColorGroup arg1)  const  {return obj->isEqual((QPalette::ColorGroup)arg0,(QPalette::ColorGroup)arg1); }
void color(QPalette* obj,ColorRole arg0)  const  {obj->color((QPalette::ColorRole)arg0); }
void brush(QPalette* obj,ColorRole arg0)  const  {obj->brush((QPalette::ColorRole)arg0); }
void foreground(QPalette* obj)  const  {obj->foreground(); }
void windowText(QPalette* obj)  const  {obj->windowText(); }
void button(QPalette* obj)  const  {obj->button(); }
void light(QPalette* obj)  const  {obj->light(); }
void dark(QPalette* obj)  const  {obj->dark(); }
void mid(QPalette* obj)  const  {obj->mid(); }
void text(QPalette* obj)  const  {obj->text(); }
void base(QPalette* obj)  const  {obj->base(); }
void alternateBase(QPalette* obj)  const  {obj->alternateBase(); }
void background(QPalette* obj)  const  {obj->background(); }
void window(QPalette* obj)  const  {obj->window(); }
void midlight(QPalette* obj)  const  {obj->midlight(); }
void brightText(QPalette* obj)  const  {obj->brightText(); }
void buttonText(QPalette* obj)  const  {obj->buttonText(); }
void shadow(QPalette* obj)  const  {obj->shadow(); }
void highlight(QPalette* obj)  const  {obj->highlight(); }
void highlightedText(QPalette* obj)  const  {obj->highlightedText(); }
void link(QPalette* obj)  const  {obj->link(); }
void linkVisited(QPalette* obj)  const  {obj->linkVisited(); }
bool isCopyOf(QPalette* obj,const QPalette & arg0)  const  {return obj->isCopyOf(arg0); }
int serialNumber(QPalette* obj)  const  {return obj->serialNumber(); }
QPalette resolve(QPalette* obj,const QPalette & arg0)  const  {return obj->resolve(arg0); }
uint resolve(QPalette* obj)  const  {return obj->resolve(); }
void resolve(QPalette* obj,uint arg0)  {obj->resolve(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qicon.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qicon.h"
class PythonQtQIconWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(Mode State )
enum Mode {Normal = QIcon::Normal, 
Disabled = QIcon::Disabled, 
Active = QIcon::Active, 
Selected = QIcon::Selected }; 
enum State {On = QIcon::On, 
Off = QIcon::Off }; 
public slots:
QVariant new_QIcon() { return QIcon(); }
QVariant new_QIcon(const QPixmap & arg0) { return QIcon(arg0); }
QVariant new_QIcon(const QIcon & arg0) { return QIcon(arg0); }
QVariant new_QIcon(const QString & arg0) { return QIcon(arg0); }
QVariant new_QIcon(QIconEngine * arg0) { return QIcon(arg0); }
QPixmap pixmap(QIcon* obj,const QSize & arg0,Mode arg1,State arg2)  const  {return obj->pixmap(arg0,(QIcon::Mode)arg1,(QIcon::State)arg2); }
QPixmap pixmap(QIcon* obj,const QSize & arg0,Mode arg1)  const  {return obj->pixmap(arg0,(QIcon::Mode)arg1); }
QPixmap pixmap(QIcon* obj,const QSize & arg0)  const  {return obj->pixmap(arg0); }
QPixmap pixmap(QIcon* obj,int arg0,int arg1,Mode arg2,State arg3)  const  {return obj->pixmap(arg0,arg1,(QIcon::Mode)arg2,(QIcon::State)arg3); }
QPixmap pixmap(QIcon* obj,int arg0,int arg1,Mode arg2)  const  {return obj->pixmap(arg0,arg1,(QIcon::Mode)arg2); }
QPixmap pixmap(QIcon* obj,int arg0,int arg1)  const  {return obj->pixmap(arg0,arg1); }
QPixmap pixmap(QIcon* obj,int arg0,Mode arg1,State arg2)  const  {return obj->pixmap(arg0,(QIcon::Mode)arg1,(QIcon::State)arg2); }
QPixmap pixmap(QIcon* obj,int arg0,Mode arg1)  const  {return obj->pixmap(arg0,(QIcon::Mode)arg1); }
QPixmap pixmap(QIcon* obj,int arg0)  const  {return obj->pixmap(arg0); }
QSize actualSize(QIcon* obj,const QSize & arg0,Mode arg1,State arg2)  const  {return obj->actualSize(arg0,(QIcon::Mode)arg1,(QIcon::State)arg2); }
QSize actualSize(QIcon* obj,const QSize & arg0,Mode arg1)  const  {return obj->actualSize(arg0,(QIcon::Mode)arg1); }
QSize actualSize(QIcon* obj,const QSize & arg0)  const  {return obj->actualSize(arg0); }
void paint(QIcon* obj,QPainter * arg0,const QRect & arg1,Qt::Alignment arg2,Mode arg3,State arg4)  const  {obj->paint(arg0,arg1,arg2,(QIcon::Mode)arg3,(QIcon::State)arg4); }
void paint(QIcon* obj,QPainter * arg0,const QRect & arg1,Qt::Alignment arg2,Mode arg3)  const  {obj->paint(arg0,arg1,arg2,(QIcon::Mode)arg3); }
void paint(QIcon* obj,QPainter * arg0,const QRect & arg1,Qt::Alignment arg2)  const  {obj->paint(arg0,arg1,arg2); }
void paint(QIcon* obj,QPainter * arg0,const QRect & arg1)  const  {obj->paint(arg0,arg1); }
void paint(QIcon* obj,QPainter * arg0,int arg1,int arg2,int arg3,int arg4,Qt::Alignment arg5,Mode arg6,State arg7)  const  {obj->paint(arg0,arg1,arg2,arg3,arg4,arg5,(QIcon::Mode)arg6,(QIcon::State)arg7); }
void paint(QIcon* obj,QPainter * arg0,int arg1,int arg2,int arg3,int arg4,Qt::Alignment arg5,Mode arg6)  const  {obj->paint(arg0,arg1,arg2,arg3,arg4,arg5,(QIcon::Mode)arg6); }
void paint(QIcon* obj,QPainter * arg0,int arg1,int arg2,int arg3,int arg4,Qt::Alignment arg5)  const  {obj->paint(arg0,arg1,arg2,arg3,arg4,arg5); }
void paint(QIcon* obj,QPainter * arg0,int arg1,int arg2,int arg3,int arg4)  const  {obj->paint(arg0,arg1,arg2,arg3,arg4); }
bool isNull(QIcon* obj)  const  {return obj->isNull(); }
bool isDetached(QIcon* obj)  const  {return obj->isDetached(); }
int serialNumber(QIcon* obj)  const  {return obj->serialNumber(); }
void addPixmap(QIcon* obj,const QPixmap & arg0,Mode arg1,State arg2)  {obj->addPixmap(arg0,(QIcon::Mode)arg1,(QIcon::State)arg2); }
void addPixmap(QIcon* obj,const QPixmap & arg0,Mode arg1)  {obj->addPixmap(arg0,(QIcon::Mode)arg1); }
void addPixmap(QIcon* obj,const QPixmap & arg0)  {obj->addPixmap(arg0); }
void addFile(QIcon* obj,const QString & arg0,const QSize & arg1,Mode arg2,State arg3)  {obj->addFile(arg0,arg1,(QIcon::Mode)arg2,(QIcon::State)arg3); }
void addFile(QIcon* obj,const QString & arg0,const QSize & arg1,Mode arg2)  {obj->addFile(arg0,arg1,(QIcon::Mode)arg2); }
void addFile(QIcon* obj,const QString & arg0,const QSize & arg1)  {obj->addFile(arg0,arg1); }
void addFile(QIcon* obj,const QString & arg0)  {obj->addFile(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qimage.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qimage.h"
class PythonQtQImageTextKeyLangWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QImageTextKeyLang(QImageTextKeyLang* obj) { delete obj; }
QImageTextKeyLang* new_QImageTextKeyLang(const char * arg0,const char * arg1) { return new QImageTextKeyLang(arg0,arg1); }

};

class PythonQtQImageWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(InvertMode Format )
enum InvertMode {InvertRgb = QImage::InvertRgb, 
InvertRgba = QImage::InvertRgba }; 
enum Format {Format_Invalid = QImage::Format_Invalid, 
Format_Mono = QImage::Format_Mono, 
Format_MonoLSB = QImage::Format_MonoLSB, 
Format_Indexed8 = QImage::Format_Indexed8, 
Format_RGB32 = QImage::Format_RGB32, 
Format_ARGB32 = QImage::Format_ARGB32, 
Format_ARGB32_Premultiplied = QImage::Format_ARGB32_Premultiplied, 
Format_RGB16 = QImage::Format_RGB16, 
NImageFormats = QImage::NImageFormats }; 
public slots:
QVariant new_QImage() { return QImage(); }
QVariant new_QImage(const QSize & arg0,Format arg1) { return QImage(arg0,(QImage::Format)arg1); }
QVariant new_QImage(int arg0,int arg1,Format arg2) { return QImage(arg0,arg1,(QImage::Format)arg2); }
QVariant new_QImage(uchar * arg0,int arg1,int arg2,Format arg3) { return QImage(arg0,arg1,arg2,(QImage::Format)arg3); }
QVariant new_QImage(const uchar * arg0,int arg1,int arg2,Format arg3) { return QImage(arg0,arg1,arg2,(QImage::Format)arg3); }
QVariant new_QImage(const QString & arg0,const char * arg1) { return QImage(arg0,arg1); }
QVariant new_QImage(const QString & arg0) { return QImage(arg0); }
QVariant new_QImage(const char * arg0,const char * arg1) { return QImage(arg0,arg1); }
QVariant new_QImage(const char * arg0) { return QImage(arg0); }
QVariant new_QImage(const QImage & arg0) { return QImage(arg0); }
bool isNull(QImage* obj)  const  {return obj->isNull(); }
int devType(QImage* obj)  const  {return obj->devType(); }
void detach(QImage* obj)  {obj->detach(); }
bool isDetached(QImage* obj)  const  {return obj->isDetached(); }
QImage copy(QImage* obj,const QRect & arg0)  const  {return obj->copy(arg0); }
QImage copy(QImage* obj)  const  {return obj->copy(); }
QImage copy(QImage* obj,int arg0,int arg1,int arg2,int arg3)  const  {return obj->copy(arg0,arg1,arg2,arg3); }
Format format(QImage* obj)  const  {return (PythonQtQImageWrapper::Format)obj->format(); }
QImage convertToFormat(QImage* obj,Format arg0,Qt::ImageConversionFlags arg1)  const  {return obj->convertToFormat((QImage::Format)arg0,arg1); }
QImage convertToFormat(QImage* obj,Format arg0)  const  {return obj->convertToFormat((QImage::Format)arg0); }
QImage convertToFormat(QImage* obj,Format arg0,const QVector<QRgb> & arg1,Qt::ImageConversionFlags arg2)  const  {return obj->convertToFormat((QImage::Format)arg0,arg1,arg2); }
QImage convertToFormat(QImage* obj,Format arg0,const QVector<QRgb> & arg1)  const  {return obj->convertToFormat((QImage::Format)arg0,arg1); }
int width(QImage* obj)  const  {return obj->width(); }
int height(QImage* obj)  const  {return obj->height(); }
QSize size(QImage* obj)  const  {return obj->size(); }
QRect rect(QImage* obj)  const  {return obj->rect(); }
int depth(QImage* obj)  const  {return obj->depth(); }
int numColors(QImage* obj)  const  {return obj->numColors(); }
QRgb color(QImage* obj,int arg0)  const  {return obj->color(arg0); }
void setColor(QImage* obj,int arg0,QRgb arg1)  {obj->setColor(arg0,arg1); }
void setNumColors(QImage* obj,int arg0)  {obj->setNumColors(arg0); }
bool allGray(QImage* obj)  const  {return obj->allGray(); }
bool isGrayscale(QImage* obj)  const  {return obj->isGrayscale(); }
uchar* bits(QImage* obj)  {return obj->bits(); }
const uchar* bits(QImage* obj)  const  {return obj->bits(); }
int numBytes(QImage* obj)  const  {return obj->numBytes(); }
uchar* scanLine(QImage* obj,int arg0)  {return obj->scanLine(arg0); }
const uchar* scanLine(QImage* obj,int arg0)  const  {return obj->scanLine(arg0); }
int bytesPerLine(QImage* obj)  const  {return obj->bytesPerLine(); }
bool valid(QImage* obj,int arg0,int arg1)  const  {return obj->valid(arg0,arg1); }
bool valid(QImage* obj,const QPoint & arg0)  const  {return obj->valid(arg0); }
int pixelIndex(QImage* obj,int arg0,int arg1)  const  {return obj->pixelIndex(arg0,arg1); }
int pixelIndex(QImage* obj,const QPoint & arg0)  const  {return obj->pixelIndex(arg0); }
QRgb pixel(QImage* obj,int arg0,int arg1)  const  {return obj->pixel(arg0,arg1); }
QRgb pixel(QImage* obj,const QPoint & arg0)  const  {return obj->pixel(arg0); }
void setPixel(QImage* obj,int arg0,int arg1,uint arg2)  {obj->setPixel(arg0,arg1,arg2); }
void setPixel(QImage* obj,const QPoint & arg0,uint arg1)  {obj->setPixel(arg0,arg1); }
QVector<QRgb> colorTable(QImage* obj)  const  {return obj->colorTable(); }
void setColorTable(QImage* obj,const QVector<QRgb> arg0)  {obj->setColorTable(arg0); }
void fill(QImage* obj,uint arg0)  {obj->fill(arg0); }
bool hasAlphaChannel(QImage* obj)  const  {return obj->hasAlphaChannel(); }
void setAlphaChannel(QImage* obj,const QImage & arg0)  {obj->setAlphaChannel(arg0); }
QImage alphaChannel(QImage* obj)  const  {return obj->alphaChannel(); }
QImage createAlphaMask(QImage* obj,Qt::ImageConversionFlags arg0)  const  {return obj->createAlphaMask(arg0); }
QImage createAlphaMask(QImage* obj)  const  {return obj->createAlphaMask(); }
QImage createHeuristicMask(QImage* obj,bool arg0)  const  {return obj->createHeuristicMask(arg0); }
QImage createHeuristicMask(QImage* obj)  const  {return obj->createHeuristicMask(); }
QImage scaled(QImage* obj,int arg0,int arg1,Qt::AspectRatioMode arg2,Qt::TransformationMode arg3)  const  {return obj->scaled(arg0,arg1,arg2,arg3); }
QImage scaled(QImage* obj,int arg0,int arg1,Qt::AspectRatioMode arg2)  const  {return obj->scaled(arg0,arg1,arg2); }
QImage scaled(QImage* obj,int arg0,int arg1)  const  {return obj->scaled(arg0,arg1); }
QImage scaled(QImage* obj,const QSize & arg0,Qt::AspectRatioMode arg1,Qt::TransformationMode arg2)  const  {return obj->scaled(arg0,arg1,arg2); }
QImage scaled(QImage* obj,const QSize & arg0,Qt::AspectRatioMode arg1)  const  {return obj->scaled(arg0,arg1); }
QImage scaled(QImage* obj,const QSize & arg0)  const  {return obj->scaled(arg0); }
QImage scaledToWidth(QImage* obj,int arg0,Qt::TransformationMode arg1)  const  {return obj->scaledToWidth(arg0,arg1); }
QImage scaledToWidth(QImage* obj,int arg0)  const  {return obj->scaledToWidth(arg0); }
QImage scaledToHeight(QImage* obj,int arg0,Qt::TransformationMode arg1)  const  {return obj->scaledToHeight(arg0,arg1); }
QImage scaledToHeight(QImage* obj,int arg0)  const  {return obj->scaledToHeight(arg0); }
QImage transformed(QImage* obj,const QMatrix & arg0,Qt::TransformationMode arg1)  const  {return obj->transformed(arg0,arg1); }
QImage transformed(QImage* obj,const QMatrix & arg0)  const  {return obj->transformed(arg0); }
QMatrix static_QImage_trueMatrix(const QMatrix & arg0,int arg1,int arg2)  {return QImage::trueMatrix(arg0,arg1,arg2); }
QImage mirrored(QImage* obj,bool arg0,bool arg1)  const  {return obj->mirrored(arg0,arg1); }
QImage mirrored(QImage* obj,bool arg0)  const  {return obj->mirrored(arg0); }
QImage mirrored(QImage* obj)  const  {return obj->mirrored(); }
QImage rgbSwapped(QImage* obj)  const  {return obj->rgbSwapped(); }
void invertPixels(QImage* obj,InvertMode arg0)  {obj->invertPixels((QImage::InvertMode)arg0); }
void invertPixels(QImage* obj)  {obj->invertPixels(); }
bool load(QImage* obj,QIODevice * arg0,const char * arg1)  {return obj->load(arg0,arg1); }
bool load(QImage* obj,const QString & arg0,const char * arg1)  {return obj->load(arg0,arg1); }
bool load(QImage* obj,const QString & arg0)  {return obj->load(arg0); }
bool loadFromData(QImage* obj,const uchar * arg0,int arg1,const char * arg2)  {return obj->loadFromData(arg0,arg1,arg2); }
bool loadFromData(QImage* obj,const uchar * arg0,int arg1)  {return obj->loadFromData(arg0,arg1); }
bool loadFromData(QImage* obj,const QByteArray & arg0,const char * arg1)  {return obj->loadFromData(arg0,arg1); }
bool loadFromData(QImage* obj,const QByteArray & arg0)  {return obj->loadFromData(arg0); }
bool save(QImage* obj,const QString & arg0,const char * arg1,int arg2)  const  {return obj->save(arg0,arg1,arg2); }
bool save(QImage* obj,const QString & arg0,const char * arg1)  const  {return obj->save(arg0,arg1); }
bool save(QImage* obj,const QString & arg0)  const  {return obj->save(arg0); }
bool save(QImage* obj,QIODevice * arg0,const char * arg1,int arg2)  const  {return obj->save(arg0,arg1,arg2); }
bool save(QImage* obj,QIODevice * arg0,const char * arg1)  const  {return obj->save(arg0,arg1); }
bool save(QImage* obj,QIODevice * arg0)  const  {return obj->save(arg0); }
QImage static_QImage_fromData(const uchar * arg0,int arg1,const char * arg2)  {return QImage::fromData(arg0,arg1,arg2); }
QImage static_QImage_fromData(const uchar * arg0,int arg1)  {return QImage::fromData(arg0,arg1); }
QImage static_QImage_fromData(const QByteArray & arg0,const char * arg1)  {return QImage::fromData(arg0,arg1); }
QImage static_QImage_fromData(const QByteArray & arg0)  {return QImage::fromData(arg0); }
int serialNumber(QImage* obj)  const  {return obj->serialNumber(); }
QPaintEngine* paintEngine(QImage* obj)  const  {return obj->paintEngine(); }
int dotsPerMeterX(QImage* obj)  const  {return obj->dotsPerMeterX(); }
int dotsPerMeterY(QImage* obj)  const  {return obj->dotsPerMeterY(); }
void setDotsPerMeterX(QImage* obj,int arg0)  {obj->setDotsPerMeterX(arg0); }
void setDotsPerMeterY(QImage* obj,int arg0)  {obj->setDotsPerMeterY(arg0); }
QPoint offset(QImage* obj)  const  {return obj->offset(); }
void setOffset(QImage* obj,const QPoint & arg0)  {obj->setOffset(arg0); }
QStringList textKeys(QImage* obj)  const  {return obj->textKeys(); }
QString text(QImage* obj,const QString & arg0)  const  {return obj->text(arg0); }
QString text(QImage* obj)  const  {return obj->text(); }
void setText(QImage* obj,const QString & arg0,const QString & arg1)  {obj->setText(arg0,arg1); }
QString text(QImage* obj,const char * arg0,const char * arg1)  const  {return obj->text(arg0,arg1); }
QString text(QImage* obj,const char * arg0)  const  {return obj->text(arg0); }
QList<QImageTextKeyLang> textList(QImage* obj)  const  {return obj->textList(); }
QStringList textLanguages(QImage* obj)  const  {return obj->textLanguages(); }
QString text(QImage* obj,const QImageTextKeyLang & arg0)  const  {return obj->text(arg0); }
void setText(QImage* obj,const char * arg0,const char * arg1,const QString & arg2)  {obj->setText(arg0,arg1,arg2); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qpolygon.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qpolygon.h"
class PythonQtQPolygonWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QPolygon() { return QPolygon(); }
QVariant new_QPolygon(int arg0) { return QPolygon(arg0); }
QVariant new_QPolygon(const QPolygon & arg0) { return QPolygon(arg0); }
QVariant new_QPolygon(const QVector<QPoint> & arg0) { return QPolygon(arg0); }
QVariant new_QPolygon(int arg0,const int * arg1) { return QPolygon(arg0,arg1); }
void translate(QPolygon* obj,int arg0,int arg1)  {obj->translate(arg0,arg1); }
void translate(QPolygon* obj,const QPoint & arg0)  {obj->translate(arg0); }
QRect boundingRect(QPolygon* obj)  const  {return obj->boundingRect(); }
void point(QPolygon* obj,int arg0,int * arg1,int * arg2)  const  {obj->point(arg0,arg1,arg2); }
QPoint point(QPolygon* obj,int arg0)  const  {return obj->point(arg0); }
void setPoint(QPolygon* obj,int arg0,int arg1,int arg2)  {obj->setPoint(arg0,arg1,arg2); }
void setPoint(QPolygon* obj,int arg0,const QPoint & arg1)  {obj->setPoint(arg0,arg1); }

};

class PythonQtQPolygonFWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
void delete_QPolygonF(QPolygonF* obj) { delete obj; }
QPolygonF* new_QPolygonF() { return new QPolygonF(); }
QPolygonF* new_QPolygonF(int arg0) { return new QPolygonF(arg0); }
QPolygonF* new_QPolygonF(const QPolygonF & arg0) { return new QPolygonF(arg0); }
QPolygonF* new_QPolygonF(const QVector<QPointF> & arg0) { return new QPolygonF(arg0); }
QPolygonF* new_QPolygonF(const QPolygon & arg0) { return new QPolygonF(arg0); }
void translate(QPolygonF* obj,qreal arg0,qreal arg1)  {obj->translate(arg0,arg1); }
void translate(QPolygonF* obj,const QPointF & arg0)  {obj->translate(arg0); }
QPolygon toPolygon(QPolygonF* obj)  const  {return obj->toPolygon(); }
bool isClosed(QPolygonF* obj)  const  {return obj->isClosed(); }
QRectF boundingRect(QPolygonF* obj)  const  {return obj->boundingRect(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qregion.h'
**
** Created: Thu 12. Apr 14:07:30 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qregion.h"
class PythonQtQRegionWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(RegionType )
enum RegionType {Rectangle = QRegion::Rectangle, 
Ellipse = QRegion::Ellipse }; 
public slots:
QVariant new_QRegion() { return QRegion(); }
QVariant new_QRegion(int arg0,int arg1,int arg2,int arg3,RegionType arg4) { return QRegion(arg0,arg1,arg2,arg3,(QRegion::RegionType)arg4); }
QVariant new_QRegion(int arg0,int arg1,int arg2,int arg3) { return QRegion(arg0,arg1,arg2,arg3); }
QVariant new_QRegion(const QRect & arg0,RegionType arg1) { return QRegion(arg0,(QRegion::RegionType)arg1); }
QVariant new_QRegion(const QRect & arg0) { return QRegion(arg0); }
QVariant new_QRegion(const QPolygon & arg0,Qt::FillRule arg1) { return QRegion(arg0,arg1); }
QVariant new_QRegion(const QPolygon & arg0) { return QRegion(arg0); }
QVariant new_QRegion(const QRegion & arg0) { return QRegion(arg0); }
QVariant new_QRegion(const QBitmap & arg0) { return QRegion(arg0); }
bool isEmpty(QRegion* obj)  const  {return obj->isEmpty(); }
bool contains(QRegion* obj,const QPoint & arg0)  const  {return obj->contains(arg0); }
bool contains(QRegion* obj,const QRect & arg0)  const  {return obj->contains(arg0); }
void translate(QRegion* obj,int arg0,int arg1)  {obj->translate(arg0,arg1); }
void translate(QRegion* obj,const QPoint & arg0)  {obj->translate(arg0); }
QRegion translated(QRegion* obj,int arg0,int arg1)  const  {return obj->translated(arg0,arg1); }
QRegion translated(QRegion* obj,const QPoint & arg0)  const  {return obj->translated(arg0); }
QRegion unite(QRegion* obj,const QRegion & arg0)  const  {return obj->unite(arg0); }
QRegion intersect(QRegion* obj,const QRegion & arg0)  const  {return obj->intersect(arg0); }
QRegion subtract(QRegion* obj,const QRegion & arg0)  const  {return obj->subtract(arg0); }
QRegion eor(QRegion* obj,const QRegion & arg0)  const  {return obj->eor(arg0); }
QRegion united(QRegion* obj,const QRegion & arg0)  const  {return obj->united(arg0); }
QRegion intersected(QRegion* obj,const QRegion & arg0)  const  {return obj->intersected(arg0); }
QRegion subtracted(QRegion* obj,const QRegion & arg0)  const  {return obj->subtracted(arg0); }
QRegion xored(QRegion* obj,const QRegion & arg0)  const  {return obj->xored(arg0); }
bool intersects(QRegion* obj,const QRegion & arg0)  const  {return obj->intersects(arg0); }
bool intersects(QRegion* obj,const QRect & arg0)  const  {return obj->intersects(arg0); }
QRect boundingRect(QRegion* obj)  const  {return obj->boundingRect(); }
QVector<QRect> rects(QRegion* obj)  const  {return obj->rects(); }
void setRects(QRegion* obj,const QRect * arg0,int arg1)  {obj->setRects(arg0,arg1); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qbitmap.h'
**
** Created: Thu 12. Apr 14:07:31 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qbitmap.h"
class PythonQtQBitmapWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QBitmap() { return QBitmap(); }
QVariant new_QBitmap(const QPixmap & arg0) { return QBitmap(arg0); }
QVariant new_QBitmap(int arg0,int arg1) { return QBitmap(arg0,arg1); }
QVariant new_QBitmap(const QSize & arg0) { return QBitmap(arg0); }
QVariant new_QBitmap(const QString & arg0,const char * arg1) { return QBitmap(arg0,arg1); }
QVariant new_QBitmap(const QString & arg0) { return QBitmap(arg0); }
void clear(QBitmap* obj)  {obj->clear(); }
QBitmap static_QBitmap_fromImage(const QImage & arg0,Qt::ImageConversionFlags arg1)  {return QBitmap::fromImage(arg0,arg1); }
QBitmap static_QBitmap_fromImage(const QImage & arg0)  {return QBitmap::fromImage(arg0); }
QBitmap static_QBitmap_fromData(const QSize & arg0,const uchar * arg1,QImage::Format arg2)  {return QBitmap::fromData(arg0,arg1,arg2); }
QBitmap static_QBitmap_fromData(const QSize & arg0,const uchar * arg1)  {return QBitmap::fromData(arg0,arg1); }
QBitmap transformed(QBitmap* obj,const QMatrix & arg0)  const  {return obj->transformed(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qcursor.h'
**
** Created: Thu 12. Apr 14:07:31 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qcursor.h"
class PythonQtQCursorWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QCursor() { return QCursor(); }
QVariant new_QCursor(Qt::CursorShape arg0) { return QCursor(arg0); }
QVariant new_QCursor(const QBitmap & arg0,const QBitmap & arg1,int arg2,int arg3) { return QCursor(arg0,arg1,arg2,arg3); }
QVariant new_QCursor(const QBitmap & arg0,const QBitmap & arg1,int arg2) { return QCursor(arg0,arg1,arg2); }
QVariant new_QCursor(const QBitmap & arg0,const QBitmap & arg1) { return QCursor(arg0,arg1); }
QVariant new_QCursor(const QPixmap & arg0,int arg1,int arg2) { return QCursor(arg0,arg1,arg2); }
QVariant new_QCursor(const QPixmap & arg0,int arg1) { return QCursor(arg0,arg1); }
QVariant new_QCursor(const QPixmap & arg0) { return QCursor(arg0); }
QVariant new_QCursor(const QCursor & arg0) { return QCursor(arg0); }
Qt::CursorShape shape(QCursor* obj)  const  {return obj->shape(); }
void setShape(QCursor* obj,Qt::CursorShape arg0)  {obj->setShape(arg0); }
const QBitmap* bitmap(QCursor* obj)  const  {return obj->bitmap(); }
const QBitmap* mask(QCursor* obj)  const  {return obj->mask(); }
QPixmap pixmap(QCursor* obj)  const  {return obj->pixmap(); }
QPoint hotSpot(QCursor* obj)  const  {return obj->hotSpot(); }
QPoint static_QCursor_pos()  {return QCursor::pos(); }
void static_QCursor_setPos(int arg0,int arg1)  {QCursor::setPos(arg0,arg1); }
void static_QCursor_setPos(const QPoint & arg0)  {QCursor::setPos(arg0); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qsizepolicy.h'
**
** Created: Thu 12. Apr 14:07:33 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qsizepolicy.h"
class PythonQtQSizePolicyWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(PolicyFlag Policy )
enum PolicyFlag {GrowFlag = QSizePolicy::GrowFlag, 
ExpandFlag = QSizePolicy::ExpandFlag, 
ShrinkFlag = QSizePolicy::ShrinkFlag, 
IgnoreFlag = QSizePolicy::IgnoreFlag }; 
enum Policy {Fixed = QSizePolicy::Fixed, 
Minimum = QSizePolicy::Minimum, 
Maximum = QSizePolicy::Maximum, 
Preferred = QSizePolicy::Preferred, 
MinimumExpanding = QSizePolicy::MinimumExpanding, 
Expanding = QSizePolicy::Expanding, 
Ignored = QSizePolicy::Ignored }; 
public slots:
QVariant new_QSizePolicy() { return QSizePolicy(); }
Policy horizontalPolicy(QSizePolicy* obj)  const  {return (PythonQtQSizePolicyWrapper::Policy)obj->horizontalPolicy(); }
Policy verticalPolicy(QSizePolicy* obj)  const  {return (PythonQtQSizePolicyWrapper::Policy)obj->verticalPolicy(); }
void setHorizontalPolicy(QSizePolicy* obj,Policy arg0)  {obj->setHorizontalPolicy((QSizePolicy::Policy)arg0); }
void setVerticalPolicy(QSizePolicy* obj,Policy arg0)  {obj->setVerticalPolicy((QSizePolicy::Policy)arg0); }
Qt::Orientations expandingDirections(QSizePolicy* obj)  const  {return obj->expandingDirections(); }
void setHeightForWidth(QSizePolicy* obj,bool arg0)  {obj->setHeightForWidth(arg0); }
bool hasHeightForWidth(QSizePolicy* obj)  const  {return obj->hasHeightForWidth(); }
int horizontalStretch(QSizePolicy* obj)  const  {return obj->horizontalStretch(); }
int verticalStretch(QSizePolicy* obj)  const  {return obj->verticalStretch(); }
void setHorizontalStretch(QSizePolicy* obj,uchar arg0)  {obj->setHorizontalStretch(arg0); }
void setVerticalStretch(QSizePolicy* obj,uchar arg0)  {obj->setVerticalStretch(arg0); }
void transpose(QSizePolicy* obj)  {obj->transpose(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qkeysequence.h'
**
** Created: Thu 12. Apr 14:07:34 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qkeysequence.h"
class PythonQtQKeySequenceWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(StandardKey SequenceMatch SequenceFormat )
enum StandardKey {UnknownKey = QKeySequence::UnknownKey, 
HelpContents = QKeySequence::HelpContents, 
WhatsThis = QKeySequence::WhatsThis, 
Open = QKeySequence::Open, 
Close = QKeySequence::Close, 
Save = QKeySequence::Save, 
New = QKeySequence::New, 
Delete = QKeySequence::Delete, 
Cut = QKeySequence::Cut, 
Copy = QKeySequence::Copy, 
Paste = QKeySequence::Paste, 
Undo = QKeySequence::Undo, 
Redo = QKeySequence::Redo, 
Back = QKeySequence::Back, 
Forward = QKeySequence::Forward, 
Refresh = QKeySequence::Refresh, 
ZoomIn = QKeySequence::ZoomIn, 
ZoomOut = QKeySequence::ZoomOut, 
Print = QKeySequence::Print, 
AddTab = QKeySequence::AddTab, 
NextChild = QKeySequence::NextChild, 
PreviousChild = QKeySequence::PreviousChild, 
Find = QKeySequence::Find, 
FindNext = QKeySequence::FindNext, 
FindPrevious = QKeySequence::FindPrevious, 
Replace = QKeySequence::Replace, 
SelectAll = QKeySequence::SelectAll, 
Bold = QKeySequence::Bold, 
Italic = QKeySequence::Italic, 
Underline = QKeySequence::Underline, 
MoveToNextChar = QKeySequence::MoveToNextChar, 
MoveToPreviousChar = QKeySequence::MoveToPreviousChar, 
MoveToNextWord = QKeySequence::MoveToNextWord, 
MoveToPreviousWord = QKeySequence::MoveToPreviousWord, 
MoveToNextLine = QKeySequence::MoveToNextLine, 
MoveToPreviousLine = QKeySequence::MoveToPreviousLine, 
MoveToNextPage = QKeySequence::MoveToNextPage, 
MoveToPreviousPage = QKeySequence::MoveToPreviousPage, 
MoveToStartOfLine = QKeySequence::MoveToStartOfLine, 
MoveToEndOfLine = QKeySequence::MoveToEndOfLine, 
MoveToStartOfBlock = QKeySequence::MoveToStartOfBlock, 
MoveToEndOfBlock = QKeySequence::MoveToEndOfBlock, 
MoveToStartOfDocument = QKeySequence::MoveToStartOfDocument, 
MoveToEndOfDocument = QKeySequence::MoveToEndOfDocument, 
SelectNextChar = QKeySequence::SelectNextChar, 
SelectPreviousChar = QKeySequence::SelectPreviousChar, 
SelectNextWord = QKeySequence::SelectNextWord, 
SelectPreviousWord = QKeySequence::SelectPreviousWord, 
SelectNextLine = QKeySequence::SelectNextLine, 
SelectPreviousLine = QKeySequence::SelectPreviousLine, 
SelectNextPage = QKeySequence::SelectNextPage, 
SelectPreviousPage = QKeySequence::SelectPreviousPage, 
SelectStartOfLine = QKeySequence::SelectStartOfLine, 
SelectEndOfLine = QKeySequence::SelectEndOfLine, 
SelectStartOfBlock = QKeySequence::SelectStartOfBlock, 
SelectEndOfBlock = QKeySequence::SelectEndOfBlock, 
SelectStartOfDocument = QKeySequence::SelectStartOfDocument, 
SelectEndOfDocument = QKeySequence::SelectEndOfDocument, 
DeleteStartOfWord = QKeySequence::DeleteStartOfWord, 
DeleteEndOfWord = QKeySequence::DeleteEndOfWord, 
DeleteEndOfLine = QKeySequence::DeleteEndOfLine }; 
enum SequenceMatch {NoMatch = QKeySequence::NoMatch, 
PartialMatch = QKeySequence::PartialMatch, 
ExactMatch = QKeySequence::ExactMatch }; 
enum SequenceFormat {NativeText = QKeySequence::NativeText, 
PortableText = QKeySequence::PortableText }; 
public slots:
QVariant new_QKeySequence() { return QKeySequence(); }
QVariant new_QKeySequence(const QString & arg0) { return QKeySequence(arg0); }
QVariant new_QKeySequence(int arg0,int arg1,int arg2,int arg3) { return QKeySequence(arg0,arg1,arg2,arg3); }
QVariant new_QKeySequence(int arg0,int arg1,int arg2) { return QKeySequence(arg0,arg1,arg2); }
QVariant new_QKeySequence(int arg0,int arg1) { return QKeySequence(arg0,arg1); }
QVariant new_QKeySequence(int arg0) { return QKeySequence(arg0); }
QVariant new_QKeySequence(const QKeySequence & arg0) { return QKeySequence(arg0); }
QVariant new_QKeySequence(StandardKey arg0) { return QKeySequence((QKeySequence::StandardKey)arg0); }
uint count(QKeySequence* obj)  const  {return obj->count(); }
bool isEmpty(QKeySequence* obj)  const  {return obj->isEmpty(); }
QString toString(QKeySequence* obj,SequenceFormat arg0)  const  {return obj->toString((QKeySequence::SequenceFormat)arg0); }
QString toString(QKeySequence* obj)  const  {return obj->toString(); }
QKeySequence static_QKeySequence_fromString(const QString & arg0,SequenceFormat arg1)  {return QKeySequence::fromString(arg0,(QKeySequence::SequenceFormat)arg1); }
QKeySequence static_QKeySequence_fromString(const QString & arg0)  {return QKeySequence::fromString(arg0); }
SequenceMatch matches(QKeySequence* obj,const QKeySequence & arg0)  const  {return (PythonQtQKeySequenceWrapper::SequenceMatch)obj->matches(arg0); }
QKeySequence static_QKeySequence_mnemonic(const QString & arg0)  {return QKeySequence::mnemonic(arg0); }
QList<QKeySequence> static_QKeySequence_keyBindings(StandardKey arg0)  {return QKeySequence::keyBindings((QKeySequence::StandardKey)arg0); }
bool isDetached(QKeySequence* obj)  const  {return obj->isDetached(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qpen.h'
**
** Created: Thu 12. Apr 14:07:35 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qpen.h"
class PythonQtQPenWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QPen() { return QPen(); }
QVariant new_QPen(Qt::PenStyle arg0) { return QPen(arg0); }
QVariant new_QPen(const QColor & arg0) { return QPen(arg0); }
QVariant new_QPen(const QBrush & arg0,qreal arg1,Qt::PenStyle arg2,Qt::PenCapStyle arg3,Qt::PenJoinStyle arg4) { return QPen(arg0,arg1,arg2,arg3,arg4); }
QVariant new_QPen(const QBrush & arg0,qreal arg1,Qt::PenStyle arg2,Qt::PenCapStyle arg3) { return QPen(arg0,arg1,arg2,arg3); }
QVariant new_QPen(const QBrush & arg0,qreal arg1,Qt::PenStyle arg2) { return QPen(arg0,arg1,arg2); }
QVariant new_QPen(const QBrush & arg0,qreal arg1) { return QPen(arg0,arg1); }
QVariant new_QPen(const QPen & arg0) { return QPen(arg0); }
Qt::PenStyle style(QPen* obj)  const  {return obj->style(); }
void setStyle(QPen* obj,Qt::PenStyle arg0)  {obj->setStyle(arg0); }
QVector<qreal> dashPattern(QPen* obj)  const  {return obj->dashPattern(); }
void setDashPattern(QPen* obj,const QVector<qreal> & arg0)  {obj->setDashPattern(arg0); }
qreal miterLimit(QPen* obj)  const  {return obj->miterLimit(); }
void setMiterLimit(QPen* obj,qreal arg0)  {obj->setMiterLimit(arg0); }
qreal widthF(QPen* obj)  const  {return obj->widthF(); }
void setWidthF(QPen* obj,qreal arg0)  {obj->setWidthF(arg0); }
int width(QPen* obj)  const  {return obj->width(); }
void setWidth(QPen* obj,int arg0)  {obj->setWidth(arg0); }
QColor color(QPen* obj)  const  {return obj->color(); }
void setColor(QPen* obj,const QColor & arg0)  {obj->setColor(arg0); }
QBrush brush(QPen* obj)  const  {return obj->brush(); }
void setBrush(QPen* obj,const QBrush & arg0)  {obj->setBrush(arg0); }
bool isSolid(QPen* obj)  const  {return obj->isSolid(); }
Qt::PenCapStyle capStyle(QPen* obj)  const  {return obj->capStyle(); }
void setCapStyle(QPen* obj,Qt::PenCapStyle arg0)  {obj->setCapStyle(arg0); }
Qt::PenJoinStyle joinStyle(QPen* obj)  const  {return obj->joinStyle(); }
void setJoinStyle(QPen* obj,Qt::PenJoinStyle arg0)  {obj->setJoinStyle(arg0); }
bool isDetached(QPen* obj)  {return obj->isDetached(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qtextformat.h'
**
** Created: Thu 12. Apr 14:07:35 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qtextformat.h"
class PythonQtQTextLengthWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(Type )
enum Type {VariableLength = QTextLength::VariableLength, 
FixedLength = QTextLength::FixedLength, 
PercentageLength = QTextLength::PercentageLength }; 
public slots:
QVariant new_QTextLength() { return QTextLength(); }
QVariant new_QTextLength(Type arg0,qreal arg1) { return QTextLength((QTextLength::Type)arg0,arg1); }
Type type(QTextLength* obj)  const  {return (PythonQtQTextLengthWrapper::Type)obj->type(); }
qreal value(QTextLength* obj,qreal arg0)  const  {return obj->value(arg0); }
qreal rawValue(QTextLength* obj)  const  {return obj->rawValue(); }

};

class PythonQtQTextFormatWrapper : public QObject  {
  Q_OBJECT

public:
Q_ENUMS(FormatType Property ObjectTypes PageBreakFlag )
enum FormatType {InvalidFormat = QTextFormat::InvalidFormat, 
BlockFormat = QTextFormat::BlockFormat, 
CharFormat = QTextFormat::CharFormat, 
ListFormat = QTextFormat::ListFormat, 
TableFormat = QTextFormat::TableFormat, 
FrameFormat = QTextFormat::FrameFormat, 
UserFormat = QTextFormat::UserFormat }; 
enum Property {ObjectIndex = QTextFormat::ObjectIndex, 
CssFloat = QTextFormat::CssFloat, 
LayoutDirection = QTextFormat::LayoutDirection, 
OutlinePen = QTextFormat::OutlinePen, 
BackgroundBrush = QTextFormat::BackgroundBrush, 
ForegroundBrush = QTextFormat::ForegroundBrush, 
BlockAlignment = QTextFormat::BlockAlignment, 
BlockTopMargin = QTextFormat::BlockTopMargin, 
BlockBottomMargin = QTextFormat::BlockBottomMargin, 
BlockLeftMargin = QTextFormat::BlockLeftMargin, 
BlockRightMargin = QTextFormat::BlockRightMargin, 
TextIndent = QTextFormat::TextIndent, 
BlockIndent = QTextFormat::BlockIndent, 
BlockNonBreakableLines = QTextFormat::BlockNonBreakableLines, 
BlockTrailingHorizontalRulerWidth = QTextFormat::BlockTrailingHorizontalRulerWidth, 
FontFamily = QTextFormat::FontFamily, 
FontPointSize = QTextFormat::FontPointSize, 
FontSizeAdjustment = QTextFormat::FontSizeAdjustment, 
FontSizeIncrement = QTextFormat::FontSizeIncrement, 
FontWeight = QTextFormat::FontWeight, 
FontItalic = QTextFormat::FontItalic, 
FontUnderline = QTextFormat::FontUnderline, 
FontOverline = QTextFormat::FontOverline, 
FontStrikeOut = QTextFormat::FontStrikeOut, 
FontFixedPitch = QTextFormat::FontFixedPitch, 
FontPixelSize = QTextFormat::FontPixelSize, 
TextUnderlineColor = QTextFormat::TextUnderlineColor, 
TextVerticalAlignment = QTextFormat::TextVerticalAlignment, 
TextOutline = QTextFormat::TextOutline, 
TextUnderlineStyle = QTextFormat::TextUnderlineStyle, 
IsAnchor = QTextFormat::IsAnchor, 
AnchorHref = QTextFormat::AnchorHref, 
AnchorName = QTextFormat::AnchorName, 
ObjectType = QTextFormat::ObjectType, 
ListStyle = QTextFormat::ListStyle, 
ListIndent = QTextFormat::ListIndent, 
FrameBorder = QTextFormat::FrameBorder, 
FrameMargin = QTextFormat::FrameMargin, 
FramePadding = QTextFormat::FramePadding, 
FrameWidth = QTextFormat::FrameWidth, 
FrameHeight = QTextFormat::FrameHeight, 
TableColumns = QTextFormat::TableColumns, 
TableColumnWidthConstraints = QTextFormat::TableColumnWidthConstraints, 
TableCellSpacing = QTextFormat::TableCellSpacing, 
TableCellPadding = QTextFormat::TableCellPadding, 
TableHeaderRowCount = QTextFormat::TableHeaderRowCount, 
TableCellRowSpan = QTextFormat::TableCellRowSpan, 
TableCellColumnSpan = QTextFormat::TableCellColumnSpan, 
ImageName = QTextFormat::ImageName, 
ImageWidth = QTextFormat::ImageWidth, 
ImageHeight = QTextFormat::ImageHeight, 
FullWidthSelection = QTextFormat::FullWidthSelection, 
PageBreakPolicy = QTextFormat::PageBreakPolicy, 
UserProperty = QTextFormat::UserProperty }; 
enum ObjectTypes {NoObject = QTextFormat::NoObject, 
ImageObject = QTextFormat::ImageObject, 
TableObject = QTextFormat::TableObject, 
UserObject = QTextFormat::UserObject }; 
enum PageBreakFlag {PageBreak_Auto = QTextFormat::PageBreak_Auto, 
PageBreak_AlwaysBefore = QTextFormat::PageBreak_AlwaysBefore, 
PageBreak_AlwaysAfter = QTextFormat::PageBreak_AlwaysAfter }; 
Q_DECLARE_FLAGS(PageBreakFlags, PageBreakFlag)
public slots:
QVariant new_QTextFormat(int arg0) { return QTextFormat(arg0); }
QVariant new_QTextFormat(const QTextFormat & arg0) { return QTextFormat(arg0); }
void merge(QTextFormat* obj,const QTextFormat & arg0)  {obj->merge(arg0); }
bool isValid(QTextFormat* obj)  const  {return obj->isValid(); }
int type(QTextFormat* obj)  const  {return obj->type(); }
int objectIndex(QTextFormat* obj)  const  {return obj->objectIndex(); }
void setObjectIndex(QTextFormat* obj,int arg0)  {obj->setObjectIndex(arg0); }
QVariant property(QTextFormat* obj,int arg0)  const  {return obj->property(arg0); }
void setProperty(QTextFormat* obj,int arg0,const QVariant & arg1)  {obj->setProperty(arg0,arg1); }
void clearProperty(QTextFormat* obj,int arg0)  {obj->clearProperty(arg0); }
bool hasProperty(QTextFormat* obj,int arg0)  const  {return obj->hasProperty(arg0); }
bool boolProperty(QTextFormat* obj,int arg0)  const  {return obj->boolProperty(arg0); }
int intProperty(QTextFormat* obj,int arg0)  const  {return obj->intProperty(arg0); }
qreal doubleProperty(QTextFormat* obj,int arg0)  const  {return obj->doubleProperty(arg0); }
QString stringProperty(QTextFormat* obj,int arg0)  const  {return obj->stringProperty(arg0); }
QColor colorProperty(QTextFormat* obj,int arg0)  const  {return obj->colorProperty(arg0); }
QPen penProperty(QTextFormat* obj,int arg0)  const  {return obj->penProperty(arg0); }
QBrush brushProperty(QTextFormat* obj,int arg0)  const  {return obj->brushProperty(arg0); }
QTextLength lengthProperty(QTextFormat* obj,int arg0)  const  {return obj->lengthProperty(arg0); }
QVector<QTextLength> lengthVectorProperty(QTextFormat* obj,int arg0)  const  {return obj->lengthVectorProperty(arg0); }
void setProperty(QTextFormat* obj,int arg0,const QVector<QTextLength> & arg1)  {obj->setProperty(arg0,arg1); }
QMap<int,QVariant> properties(QTextFormat* obj)  const  {return obj->properties(); }
void setObjectType(QTextFormat* obj,int arg0)  {obj->setObjectType(arg0); }
int objectType(QTextFormat* obj)  const  {return obj->objectType(); }
bool isCharFormat(QTextFormat* obj)  const  {return obj->isCharFormat(); }
bool isBlockFormat(QTextFormat* obj)  const  {return obj->isBlockFormat(); }
bool isListFormat(QTextFormat* obj)  const  {return obj->isListFormat(); }
bool isFrameFormat(QTextFormat* obj)  const  {return obj->isFrameFormat(); }
bool isImageFormat(QTextFormat* obj)  const  {return obj->isImageFormat(); }
bool isTableFormat(QTextFormat* obj)  const  {return obj->isTableFormat(); }
QTextBlockFormat toBlockFormat(QTextFormat* obj)  const  {return obj->toBlockFormat(); }
QTextCharFormat toCharFormat(QTextFormat* obj)  const  {return obj->toCharFormat(); }
QTextListFormat toListFormat(QTextFormat* obj)  const  {return obj->toListFormat(); }
QTextTableFormat toTableFormat(QTextFormat* obj)  const  {return obj->toTableFormat(); }
QTextFrameFormat toFrameFormat(QTextFormat* obj)  const  {return obj->toFrameFormat(); }
QTextImageFormat toImageFormat(QTextFormat* obj)  const  {return obj->toImageFormat(); }
void setLayoutDirection(QTextFormat* obj,Qt::LayoutDirection arg0)  {obj->setLayoutDirection(arg0); }
Qt::LayoutDirection layoutDirection(QTextFormat* obj)  const  {return obj->layoutDirection(); }
void setBackground(QTextFormat* obj,const QBrush & arg0)  {obj->setBackground(arg0); }
QBrush background(QTextFormat* obj)  const  {return obj->background(); }
void clearBackground(QTextFormat* obj)  {obj->clearBackground(); }
void setForeground(QTextFormat* obj,const QBrush & arg0)  {obj->setForeground(arg0); }
QBrush foreground(QTextFormat* obj)  const  {return obj->foreground(); }
void clearForeground(QTextFormat* obj)  {obj->clearForeground(); }

};

/****************************************************************************
** Meta object code from reading C++ file 'qmatrix.h'
**
** Created: Thu 12. Apr 14:07:35 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.2.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "qmatrix.h"
class PythonQtQMatrixWrapper : public QObject  {
  Q_OBJECT

public:
public slots:
QVariant new_QMatrix() { return QMatrix(); }
QVariant new_QMatrix(qreal arg0,qreal arg1,qreal arg2,qreal arg3,qreal arg4,qreal arg5) { return QMatrix(arg0,arg1,arg2,arg3,arg4,arg5); }
QVariant new_QMatrix(const QMatrix & arg0) { return QMatrix(arg0); }
void setMatrix(QMatrix* obj,qreal arg0,qreal arg1,qreal arg2,qreal arg3,qreal arg4,qreal arg5)  {obj->setMatrix(arg0,arg1,arg2,arg3,arg4,arg5); }
qreal m11(QMatrix* obj)  const  {return obj->m11(); }
qreal m12(QMatrix* obj)  const  {return obj->m12(); }
qreal m21(QMatrix* obj)  const  {return obj->m21(); }
qreal m22(QMatrix* obj)  const  {return obj->m22(); }
qreal dx(QMatrix* obj)  const  {return obj->dx(); }
qreal dy(QMatrix* obj)  const  {return obj->dy(); }
void map(QMatrix* obj,int arg0,int arg1,int * arg2,int * arg3)  const  {obj->map(arg0,arg1,arg2,arg3); }
void map(QMatrix* obj,qreal arg0,qreal arg1,qreal * arg2,qreal * arg3)  const  {obj->map(arg0,arg1,arg2,arg3); }
QRect mapRect(QMatrix* obj,const QRect & arg0)  const  {return obj->mapRect(arg0); }
QRectF mapRect(QMatrix* obj,const QRectF & arg0)  const  {return obj->mapRect(arg0); }
QPoint map(QMatrix* obj,const QPoint & arg0)  const  {return obj->map(arg0); }
QPointF map(QMatrix* obj,const QPointF & arg0)  const  {return obj->map(arg0); }
QLine map(QMatrix* obj,const QLine & arg0)  const  {return obj->map(arg0); }
QLineF map(QMatrix* obj,const QLineF & arg0)  const  {return obj->map(arg0); }
QPolygonF map(QMatrix* obj,const QPolygonF & arg0)  const  {return obj->map(arg0); }
QPolygon map(QMatrix* obj,const QPolygon & arg0)  const  {return obj->map(arg0); }
QRegion map(QMatrix* obj,const QRegion & arg0)  const  {return obj->map(arg0); }
QPainterPath map(QMatrix* obj,const QPainterPath & arg0)  const  {return obj->map(arg0); }
QPolygon mapToPolygon(QMatrix* obj,const QRect & arg0)  const  {return obj->mapToPolygon(arg0); }
void reset(QMatrix* obj)  {obj->reset(); }
bool isIdentity(QMatrix* obj)  const  {return obj->isIdentity(); }
void translate(QMatrix* obj,qreal arg0,qreal arg1)  {obj->translate(arg0,arg1); }
void scale(QMatrix* obj,qreal arg0,qreal arg1)  {obj->scale(arg0,arg1); }
void shear(QMatrix* obj,qreal arg0,qreal arg1)  {obj->shear(arg0,arg1); }
void rotate(QMatrix* obj,qreal arg0)  {obj->rotate(arg0); }
bool isInvertible(QMatrix* obj)  const  {return obj->isInvertible(); }
qreal det(QMatrix* obj)  const  {return obj->det(); }
QMatrix inverted(QMatrix* obj,bool * arg0)  const  {return obj->inverted(arg0); }
QMatrix inverted(QMatrix* obj)  const  {return obj->inverted(); }

};

