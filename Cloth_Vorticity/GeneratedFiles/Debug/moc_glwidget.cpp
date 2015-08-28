/****************************************************************************
** Meta object code from reading C++ file 'glwidget.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../glwidget.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'glwidget.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_GLWidget_t {
    QByteArrayData data[25];
    char stringdata[511];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GLWidget_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GLWidget_t qt_meta_stringdata_GLWidget = {
    {
QT_MOC_LITERAL(0, 0, 8), // "GLWidget"
QT_MOC_LITERAL(1, 9, 18), // "frameNumberUpdated"
QT_MOC_LITERAL(2, 28, 0), // ""
QT_MOC_LITERAL(3, 29, 9), // "frame_nos"
QT_MOC_LITERAL(4, 39, 22), // "animFrameNumberUpdated"
QT_MOC_LITERAL(5, 62, 23), // "setRenderingFrameNumber"
QT_MOC_LITERAL(6, 86, 12), // "frame_number"
QT_MOC_LITERAL(7, 99, 14), // "startAnimation"
QT_MOC_LITERAL(8, 114, 14), // "pauseAnimation"
QT_MOC_LITERAL(9, 129, 14), // "resetAnimation"
QT_MOC_LITERAL(10, 144, 9), // "nextFrame"
QT_MOC_LITERAL(11, 154, 18), // "bodyDisplayToggled"
QT_MOC_LITERAL(12, 173, 7), // "checked"
QT_MOC_LITERAL(13, 181, 19), // "clothDisplayToggled"
QT_MOC_LITERAL(14, 201, 19), // "renderGroundToggled"
QT_MOC_LITERAL(15, 221, 17), // "renderAxesToggled"
QT_MOC_LITERAL(16, 239, 22), // "setBodyRenderInShading"
QT_MOC_LITERAL(17, 262, 24), // "setBodyRenderInWireframe"
QT_MOC_LITERAL(18, 287, 23), // "setClothRenderInShading"
QT_MOC_LITERAL(19, 311, 25), // "setClothRenderInWireframe"
QT_MOC_LITERAL(20, 337, 31), // "setClothRenderInHeatmapVelocity"
QT_MOC_LITERAL(21, 369, 35), // "setClothRenderInHeatmapAccele..."
QT_MOC_LITERAL(22, 405, 35), // "setClothRenderInHeatmapBendin..."
QT_MOC_LITERAL(23, 441, 35), // "setClothRenderInHeatmapDampin..."
QT_MOC_LITERAL(24, 477, 33) // "setClothRenderInHeatmapShearF..."

    },
    "GLWidget\0frameNumberUpdated\0\0frame_nos\0"
    "animFrameNumberUpdated\0setRenderingFrameNumber\0"
    "frame_number\0startAnimation\0pauseAnimation\0"
    "resetAnimation\0nextFrame\0bodyDisplayToggled\0"
    "checked\0clothDisplayToggled\0"
    "renderGroundToggled\0renderAxesToggled\0"
    "setBodyRenderInShading\0setBodyRenderInWireframe\0"
    "setClothRenderInShading\0"
    "setClothRenderInWireframe\0"
    "setClothRenderInHeatmapVelocity\0"
    "setClothRenderInHeatmapAcceleration\0"
    "setClothRenderInHeatmapBendingForce\0"
    "setClothRenderInHeatmapDampingForce\0"
    "setClothRenderInHeatmapShearForce"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GLWidget[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      22,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,  124,    2, 0x06 /* Public */,
       1,    1,  127,    2, 0x06 /* Public */,
       4,    1,  130,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    1,  133,    2, 0x0a /* Public */,
       5,    1,  136,    2, 0x0a /* Public */,
       7,    0,  139,    2, 0x0a /* Public */,
       8,    0,  140,    2, 0x0a /* Public */,
       9,    0,  141,    2, 0x0a /* Public */,
      10,    0,  142,    2, 0x0a /* Public */,
      11,    1,  143,    2, 0x0a /* Public */,
      13,    1,  146,    2, 0x0a /* Public */,
      14,    1,  149,    2, 0x0a /* Public */,
      15,    1,  152,    2, 0x0a /* Public */,
      16,    1,  155,    2, 0x0a /* Public */,
      17,    1,  158,    2, 0x0a /* Public */,
      18,    1,  161,    2, 0x0a /* Public */,
      19,    1,  164,    2, 0x0a /* Public */,
      20,    1,  167,    2, 0x0a /* Public */,
      21,    1,  170,    2, 0x0a /* Public */,
      22,    1,  173,    2, 0x0a /* Public */,
      23,    1,  176,    2, 0x0a /* Public */,
      24,    1,  179,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::Double,    3,
    QMetaType::Void, QMetaType::Int,    3,
    QMetaType::Void, QMetaType::Int,    3,

 // slots: parameters
    QMetaType::Void, QMetaType::Double,    6,
    QMetaType::Void, QMetaType::Int,    6,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,
    QMetaType::Void, QMetaType::Bool,   12,

       0        // eod
};

void GLWidget::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GLWidget *_t = static_cast<GLWidget *>(_o);
        switch (_id) {
        case 0: _t->frameNumberUpdated((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: _t->frameNumberUpdated((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 2: _t->animFrameNumberUpdated((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 3: _t->setRenderingFrameNumber((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 4: _t->setRenderingFrameNumber((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 5: _t->startAnimation(); break;
        case 6: _t->pauseAnimation(); break;
        case 7: _t->resetAnimation(); break;
        case 8: _t->nextFrame(); break;
        case 9: _t->bodyDisplayToggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 10: _t->clothDisplayToggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 11: _t->renderGroundToggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 12: _t->renderAxesToggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 13: _t->setBodyRenderInShading((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 14: _t->setBodyRenderInWireframe((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 15: _t->setClothRenderInShading((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 16: _t->setClothRenderInWireframe((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 17: _t->setClothRenderInHeatmapVelocity((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 18: _t->setClothRenderInHeatmapAcceleration((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 19: _t->setClothRenderInHeatmapBendingForce((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 20: _t->setClothRenderInHeatmapDampingForce((*reinterpret_cast< bool(*)>(_a[1]))); break;
        case 21: _t->setClothRenderInHeatmapShearForce((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (GLWidget::*_t)(double );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLWidget::frameNumberUpdated)) {
                *result = 0;
            }
        }
        {
            typedef void (GLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLWidget::frameNumberUpdated)) {
                *result = 1;
            }
        }
        {
            typedef void (GLWidget::*_t)(int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&GLWidget::animFrameNumberUpdated)) {
                *result = 2;
            }
        }
    }
}

const QMetaObject GLWidget::staticMetaObject = {
    { &QGLWidget::staticMetaObject, qt_meta_stringdata_GLWidget.data,
      qt_meta_data_GLWidget,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *GLWidget::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GLWidget::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_GLWidget.stringdata))
        return static_cast<void*>(const_cast< GLWidget*>(this));
    return QGLWidget::qt_metacast(_clname);
}

int GLWidget::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 22)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 22;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 22)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 22;
    }
    return _id;
}

// SIGNAL 0
void GLWidget::frameNumberUpdated(double _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void GLWidget::frameNumberUpdated(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void GLWidget::animFrameNumberUpdated(int _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}
QT_END_MOC_NAMESPACE
