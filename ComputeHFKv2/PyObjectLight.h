#ifndef PY_OBJECT_LIGHT_H
#define PY_OBJECT_LIGHT_H

namespace py
{

class object
{
public:
    object(bool b)
      : _obj(PyBool_FromLong(b))
    { }

    object(int i)
      : _obj(PyLong_FromLong(i))
    { }

    object(size_t i)
      : _obj(PyLong_FromLong(i))
    { }

    object(const char * c)
      : _obj(PyUnicode_FromString(c))
    { }

    object(const std::string &s)
      : object(s.c_str())
    { }

    object(const object &o)
      : _obj(o._obj)
    {
        Py_INCREF(_obj);
    }
    
    template<typename T>
    object(const std::vector<T> &v)
      : _obj(PyList_New(v.size()))
    {
        for (size_t i = 0; i < v.size(); i++) {
            object t(v[i]);
            PyList_SetItem(_obj, i, t.StealObject());
        }
    }

    template<typename K, typename V>
    object(const std::map<K, V> &m)
      : _obj(PyDict_New())
    {
        for (const auto &item : m) {
            object k(item.first);
            object v(item.second);
            PyDict_SetItem(_obj, k.GetObject(), v.GetObject());
        }
    }

    template<typename T, typename U, typename ...V>
    object(const T & t, const U & u, const V & ...v)
      : _obj(_tuple_helper(object(t), object(u), object(v)...))
    {
    }

    template<typename T, typename U>
    object(const std::pair<T, U> &p)
      : object(p.first, p.second)
    {
    }

    ~object() {
        Py_DECREF(_obj);
    }

    PyObject* GetObject() const {
        return _obj;
    }

    PyObject* StealObject() const {
        Py_INCREF(_obj); return _obj;
    }

protected:
    template<typename ...T>
    static PyObject* _tuple_helper(const T &...t) {
        return PyTuple_Pack(sizeof...(T), t.GetObject()...);
    }

private:
    PyObject* _obj;
};

inline void RaiseValueError(const char * msg)
{
    Py_INCREF(PyExc_ValueError);
    PyErr_SetString(PyExc_ValueError, msg);
}

inline void RaiseValueError(const std::string &msg)
{
    RaiseValueError(msg.c_str());
}

}

#endif
