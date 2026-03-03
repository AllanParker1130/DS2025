#ifndef ALGOS_H
#define ALGOS_H

#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <stdexcept> // 添加异常处理支持

using namespace std;

// 起泡排序算法（版本1A）
void bubblesort1A(int A[], int n) {
    bool sorted = false;
    while (!sorted) {
        sorted = true;
        for (int i = 1; i < n; i++) {
            if (A[i - 1] > A[i]) {
                std::swap(A[i - 1], A[i]);
                sorted = false;
            }
        }
        n--;
    }
}

// 统计整数二进制展开中数位1的总数
int countOnes(unsigned int n) {
    int ones = 0;
    while (n) {
        ones += (1 & n);
        n >>= 1;
    }
    return ones;
}

// 幂函数2^n算法（蛮力迭代版）
__int64 power2BF_I(int n) {
    __int64 pow = 1;
    while (n-- > 0)
        pow <<= 1;
    return pow;
}

// 幂函数2^n算法（优化递归版）
__int64 sqr(__int64 a) { return a * a; }
__int64 power2(int n) {
    if (n == 0)
        return 1;
    return (n & 1) ? sqr(power2(n >> 1)) << 1 : sqr(power2(n >> 1));
}

// 数组倒置
void reverse(int* A, int lo, int hi) {
    while (lo < hi)
        std::swap(A[lo++], A[hi--]);
}

// Vector模板类定义
typedef int Rank;
#define DEFAULT_CAPACITY 3

template <typename T>
class Vector {
    protected:
        Rank _size;
        int _capacity;
        T* _elem;

        void copyFrom(T const* A, Rank lo, Rank hi);
        void expand();
        void shrink();
        bool bubble(Rank lo, Rank hi);
        void merge(Rank lo, Rank mi, Rank hi);
        Rank partition(Rank lo, Rank hi);
        void quickSort(Rank lo, Rank hi);
        void heapSort(Rank lo, Rank hi);
        void heapify(Rank lo, Rank hi, Rank root);

    public:
        // 排序算法接口
        void bubbleSort(Rank lo, Rank hi);
        void selectionSort(Rank lo, Rank hi);
        void mergeSort(Rank lo, Rank hi);
        // 构造函数、复制构造函数、赋值运算符和析构函数
        Vector(int c = DEFAULT_CAPACITY, int s = 0, T v = T());
        Vector(T const* A, Rank n);
        Vector(T const* A, Rank lo, Rank hi);
        Vector(Vector<T> const& V);
        Vector(Vector<T> const& V, Rank lo, Rank hi);
        ~Vector();

        // 向量操作
        Rank size() const { return _size; }
        bool empty() const { return !_size; }
        int disordered() const;
        Rank find(T const& e) const { return find(e, 0, _size); }
        Rank find(T const& e, Rank lo, Rank hi) const;
        Rank search(T const& e) const { return (0 >= _size) ? -1 : search(e, 0, _size); }
        Rank search(T const& e, Rank lo, Rank hi) const;

        // 运算符重载
        T& operator[](Rank r) const;
        Vector<T>& operator=(Vector<T> const&);
        
        // 元素插入与删除
        T remove(Rank r);
        int remove(Rank lo, Rank hi);
        Rank insert(Rank r, T const& e);
        Rank insert(T const& e) { return insert(_size, e); }
        
        // 向量排序与置乱
        void sort(Rank lo, Rank hi);
        void sort() { sort(0, _size); }
        void unsort(Rank lo, Rank hi);
        void unsort() { unsort(0, _size); }
        
        // 向量去重与唯一化
        int deduplicate();
        int uniquify();

        // 遍历操作
        void traverse(void (*)(T&));
        template <typename VST>
        void traverse(VST&);
};

// Vector析构函数
template <typename T>
Vector<T>::~Vector() {
    delete[] _elem;
}

// Vector模板类成员函数实现
template <typename T>
Vector<T>::Vector(int c, int s, T v) : _size(s), _capacity(c) {
    _elem = new T[_capacity];
    for (int i = 0; i < _size; i++) {
        _elem[i] = v;
    }
}

template <typename T>
void Vector<T>::copyFrom(T const* A, Rank lo, Rank hi) {
    _elem = new T[_capacity = 2 * (hi - lo)];
    _size = 0;
    while (lo < hi)
        _elem[_size++] = A[lo++];
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector<T> const& V) {
    if (_elem)
        delete[] _elem;
    copyFrom(V._elem, 0, V._size);
    return *this;
}

template <typename T>
void Vector<T>::expand() {
    if (_size < _capacity)
        return;
    if (_capacity < DEFAULT_CAPACITY)
        _capacity = DEFAULT_CAPACITY;
    T* oldElem = _elem;
    _elem = new T[(_capacity <<= 1)];
    for (int i = 0; i < _size; i++)
        _elem[i] = oldElem[i];
    delete[] oldElem;
}

template <typename T>
void Vector<T>::shrink() {
    if (_capacity < (DEFAULT_CAPACITY << 1))
        return;
    if ((_size << 2) > _capacity)
        return;
    T* oldElem = _elem;
    _elem = new T[(_capacity >>= 1)];
    for (int i = 0; i < _size; i++)
        _elem[i] = oldElem[i];
    delete[] oldElem;
}

template <typename T>
T& Vector<T>::operator[](Rank r) const {
    return _elem[r];
}

template <typename T>
void Vector<T>::unsort(Rank lo, Rank hi) {
    T* V = _elem + lo;
    for (Rank i = hi - lo; i > 0; i--)
        std::swap(V[i - 1], V[rand() % i]);
}

template <typename T>
int Vector<T>::deduplicate() {
    int oldSize = _size;
    Rank i = 1;
    while (i < _size) {
        if (find(_elem[i], 0, i) < 0)
            i++;
        else
            remove(i);
    }
    return oldSize - _size;
}

template <typename T>
void Vector<T>::traverse(void (*visit)(T&)) {
    for (int i = 0; i < _size; i++)
        visit(_elem[i]);
}

template <typename T>
template <typename VST>
void Vector<T>::traverse(VST& visit) {
    for (int i = 0; i < _size; i++)
        visit(_elem[i]);
}

template <typename T>
int Vector<T>::disordered() const {
    int n = 0;
    for (int i = 1; i < _size; i++)
        if (_elem[i - 1] > _elem[i])
            n++;
    return n;
}

template <typename T>
Rank Vector<T>::find(T const& e, Rank lo, Rank hi) const {
    while (lo < hi--) {
        if (e == _elem[hi]) {
            return hi;
        }
    }
    return -1;
}

template <typename T>
int Vector<T>::uniquify() {
    int oldSize = _size;
    Rank i = 0, j = 0;
    while (++j < _size) {
        if (_elem[i] != _elem[j]) {
            _elem[++i] = _elem[j];
        }
    }
    _size = i + 1;
    shrink();
    return oldSize - _size;
}

// 修正后的排序函数声明
template <typename T>
void Vector<T>::sort(Rank lo, Rank hi) {
    switch (rand() % 4) { 
        case 0: selectionSort(lo, hi); break;
        case 1: mergeSort(lo, hi); break;
        case 2: heapSort(lo, hi); break;
        case 3: quickSort(lo, hi); break;
    }
}

template <typename T>
void Vector<T>::selectionSort(Rank lo, Rank hi) {
    for (Rank i = lo; i < hi; ++i) {
        Rank min_idx = i;
        for (Rank j = i + 1; j < hi; ++j) {
            if (_elem[j] < _elem[min_idx]) {
                min_idx = j;
            }
        }
        if (min_idx != i) {
            std::swap(_elem[i], _elem[min_idx]);
        }
    }
}

template <typename T>
void Vector<T>::mergeSort(Rank lo, Rank hi) {
    if (hi - lo < 2)
        return;
    int mi = (lo + hi) >> 1;
    mergeSort(lo, mi);
    mergeSort(mi, hi);
    merge(lo, mi, hi);
}

template <typename T>
void Vector<T>::merge(Rank lo, Rank mi, Rank hi) {
    T* A = _elem + lo;
    int lb = mi - lo;
    T* B = new T[lb];
    for (Rank i = 0; i < lb; i++) {
        B[i] = A[i];
    }

    int lc = hi - mi;
    T* C = _elem + mi;

    for (Rank i = 0, j = 0, k = 0; (j < lb) || (k < lc);) {
        if ((j < lb) && (!(k < lc) || (B[j] <= C[k])))
            A[i++] = B[j++];
        if ((k < lc) && (!(j < lb) || (C[k] < B[j])))
            A[i++] = C[k++];
    }

    delete[] B;
}

template <typename T>
void Vector<T>::heapSort(Rank lo, Rank hi) {
    // 构建最大堆
    for (Rank i = (hi - 1) / 2; i >= lo; --i) {
        heapify(lo, hi, i);
    }
    // 交换堆顶元素到最后的位置
    for (Rank i = hi - 1; i > lo; --i) {
        std::swap(_elem[lo], _elem[i]);
        heapify(lo, i, lo);
    }
}

template <typename T>
void Vector<T>::heapify(Rank lo, Rank hi, Rank root) {
    Rank pos = root; // 假设当前根节点是最大值的位置
    Rank lchild = 2 * (root - lo) + 1 + lo; // 左孩子
    Rank rchild = 2 * (root - lo) + 2 + lo; // 右孩子

    if (lchild < hi && _elem[lchild] > _elem[pos]) {
        pos = lchild;
    }
    if (rchild < hi && _elem[rchild] > _elem[pos]) {
        pos = rchild;
    }
    if (pos != root) { // 如果最大值不是根节点，交换并递归调整
        std::swap(_elem[root], _elem[pos]);
        heapify(lo, hi, pos);
    }
}

template <typename T>
void Vector<T>::quickSort(Rank lo, Rank hi) {
    if (hi - lo < 2) return; // 基线条件：子数组长度为1或0，直接返回

    Rank pivot = partition(lo, hi); // 分区操作，返回基准值位置
    quickSort(lo, pivot);          // 递归排序左半部分
    quickSort(pivot + 1, hi);      // 递归排序右半部分
}

template <typename T>
Rank Vector<T>::partition(Rank lo, Rank hi) {
    Rank pivot = lo + rand() % (hi - lo); // 随机选择基准值
    std::swap(_elem[pivot], _elem[lo]);   // 将基准值放到起始位置
    pivot = lo;

    for (Rank i = lo + 1; i < hi; ++i) {
        if (_elem[i] < _elem[pivot]) { // 将小于基准值的元素放到左边
            std::swap(_elem[++pivot], _elem[i]);
        }
    }
    return pivot;
}

// 修正后的冒泡排序
template <typename T>
void Vector<T>::bubbleSort(Rank lo, Rank hi) {
    if (lo >= hi) {
        throw std::out_of_range("Invalid range for bubbleSort");
    }
    while (!bubble(lo, hi)) {
        hi--;
    }
}

template <typename T>
bool Vector<T>::bubble(Rank lo, Rank hi) {
    if (lo >= hi) {
        throw std::out_of_range("Invalid range for bubble");
    }
    bool sorted = true;
    for (Rank i = lo; i < hi - 1; i++) { // 修正循环条件
        if (_elem[i] > _elem[i + 1]) {
            sorted = false;
            std::swap(_elem[i], _elem[i + 1]);
        }
    }
    return sorted;
}

template <typename T>
int Vector<T>::remove(Rank lo, Rank hi) {
    for (Rank i = hi; i < _size; i++) {
        _elem[i - (hi - lo)] = _elem[i];
    }
    _size -= (hi - lo);
    shrink();
    return hi - lo;
}

template <typename T>
Rank Vector<T>::insert(Rank r, T const& e) {
    expand();
    for (int i = _size; i > r; i--) {
        _elem[i] = _elem[i - 1];
    }
    _elem[r] = e;
    _size++;
    return r;
}

template <typename T>
T Vector<T>::remove(Rank r) {
    T e = _elem[r];
    for (int i = r; i < _size - 1; i++) {
        _elem[i] = _elem[i + 1];
    }
    _size--;
    return e;
}

// 复制构造函数的实现
template <typename T>
Vector<T>::Vector(Vector<T> const& V) {
    copyFrom(V._elem, 0, V._size);
}

// 复制区间构造函数的实现
template <typename T>
Vector<T>::Vector(Vector<T> const& V, Rank lo, Rank hi) {
    copyFrom(V._elem, lo, hi);
}

// 从数组构造的实现
template <typename T>
Vector<T>::Vector(T const* A, Rank n) {
    copyFrom(A, 0, n);
}

// 从数组区间构造的实现
template <typename T>
Vector<T>::Vector(T const* A, Rank lo, Rank hi) {
    copyFrom(A, lo, hi);
}

#endif // ALGOS_H
