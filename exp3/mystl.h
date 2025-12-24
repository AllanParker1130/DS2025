#ifndef MYSTL_H
#define MYSTL_H

#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <stdexcept>
#include <utility>
#include <climits>
#include <cstring>
#include <random>
#include <stack>
#include <queue>

namespace mystl {

// ==================== 基础工具与类型定义 ====================

using Rank = int;
const int DEFAULT_CAPACITY = 3;

// 统计整数二进制展开中1的个数
inline int countOnes(unsigned int n) {
    int ones = 0;
    while (n) {
        n &= (n - 1); // 清除最低位的1
        ++ones;
    }
    return ones;
}

// 幂函数2^n（迭代版）
inline int64_t power2BF_I(int n) {
    if (n < 0) throw std::invalid_argument("Exponent must be non-negative");
    int64_t pow = 1;
    while (n-- > 0) pow <<= 1;
    return pow;
}

// 幂函数2^n（快速幂版）
inline int64_t power2(int n) {
    if (n < 0) throw std::invalid_argument("Exponent must be non-negative");
    if (n == 0) return 1;
    int64_t half = power2(n >> 1);
    return (n & 1) ? (half * half * 2) : (half * half);
}

// 数组倒置
template <typename T>
void reverse(T* A, Rank lo, Rank hi) {
    while (lo < hi) std::swap(A[lo++], A[hi--]);
}

// 起泡排序
template <typename T>
void bubblesort(T* A, Rank n) {
    bool sorted = false;
    while (!sorted) {
        sorted = true;
        for (Rank i = 1; i < n; ++i) {
            if (A[i - 1] > A[i]) {
                std::swap(A[i - 1], A[i]);
                sorted = false;
            }
        }
        --n;
    }
}

// ==================== 向量类模板 ====================

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
    void heapSort(Rank lo, Rank hi);
    void heapify(Rank lo, Rank hi, Rank root);

public:
    // 构造函数
    Vector(int c = DEFAULT_CAPACITY, Rank s = 0, T v = T());
    Vector(T const* A, Rank n);
    Vector(T const* A, Rank lo, Rank hi);
    Vector(Vector const& V);
    Vector(Vector const& V, Rank lo, Rank hi);
    ~Vector();
    
    // 重载赋值运算符
    Vector& operator=(Vector const& V);
    
    // 元素访问
    T& operator[](Rank r) { return _elem[r]; }
    T const& operator[](Rank r) const { return _elem[r]; }
    
    // 基本信息
    Rank size() const { return _size; }
    bool empty() const { return !_size; }
    int disordered() const;
    
    // 查找与搜索
    Rank find(T const& e) const { return find(e, 0, _size); }
    Rank find(T const& e, Rank lo, Rank hi) const;
    Rank search(T const& e) const { return search(e, 0, _size); }
    Rank search(T const& e, Rank lo, Rank hi) const;
    
    // 修改操作
    Rank insert(Rank r, T const& e);
    Rank insert(T const& e) { return insert(_size, e); }
    T remove(Rank r);
    int remove(Rank lo, Rank hi);
    
    // 排序算法
    void sort(Rank lo, Rank hi);
    void sort() { sort(0, _size); }
    void bubbleSort(Rank lo, Rank hi);
    void selectionSort(Rank lo, Rank hi);
    void mergeSort(Rank lo, Rank hi);
    void quickSort(Rank lo, Rank hi);
    
    // 置乱
    void unsort(Rank lo, Rank hi);
    void unsort() { unsort(0, _size); }
    
    // 去重
    int deduplicate();
    int uniquify();
    
    // 遍历
    void traverse(void (*visit)(T&));
    template <typename VST>
    void traverse(VST& visit);
};

// Vector 实现
template <typename T>
Vector<T>::~Vector() {
    delete[] _elem;
}

template <typename T>
Vector<T>::Vector(int c, Rank s, T v) : _size(s), _capacity(c) {
    _elem = new T[_capacity];
    for (Rank i = 0; i < _size; ++i) _elem[i] = v;
}

template <typename T>
void Vector<T>::copyFrom(T const* A, Rank lo, Rank hi) {
    _capacity = 2 * (hi - lo);
    _elem = new T[_capacity];
    _size = 0;
    while (lo < hi) _elem[_size++] = A[lo++];
}

template <typename T>
Vector<T>& Vector<T>::operator=(Vector const& V) {
    if (this == &V) return *this;
    delete[] _elem;
    copyFrom(V._elem, 0, V._size);
    return *this;
}

template <typename T>
void Vector<T>::expand() {
    if (_size < _capacity) return;
    if (_capacity < DEFAULT_CAPACITY) _capacity = DEFAULT_CAPACITY;
    T* oldElem = _elem;
    _elem = new T[_capacity <<= 1];
    for (Rank i = 0; i < _size; ++i) _elem[i] = oldElem[i];
    delete[] oldElem;
}

template <typename T>
void Vector<T>::shrink() {
    if (_capacity < (DEFAULT_CAPACITY << 1)) return;
    if ((_size << 2) > _capacity) return;
    T* oldElem = _elem;
    _elem = new T[_capacity >>= 1];
    for (Rank i = 0; i < _size; ++i) _elem[i] = oldElem[i];
    delete[] oldElem;
}

template <typename T>
int Vector<T>::disordered() const {
    int n = 0;
    for (Rank i = 1; i < _size; ++i)
        if (_elem[i - 1] > _elem[i]) ++n;
    return n;
}

template <typename T>
Rank Vector<T>::find(T const& e, Rank lo, Rank hi) const {
    while (lo < hi--) if (e == _elem[hi]) return hi;
    return -1;
}

template <typename T>
Rank Vector<T>::search(T const& e, Rank lo, Rank hi) const {
    // 二分查找（假设已排序）
    while (lo < hi) {
        Rank mi = (lo + hi) >> 1;
        if (e < _elem[mi]) hi = mi;
        else if (_elem[mi] < e) lo = mi + 1;
        else return mi;
    }
    return -1;
}

template <typename T>
void Vector<T>::unsort(Rank lo, Rank hi) {
    static std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    for (Rank i = hi - lo; i > 0; --i) {
        std::uniform_int_distribution<> dis(0, i - 1);
        std::swap(_elem[lo + i - 1], _elem[lo + dis(gen)]);
    }
}

template <typename T>
int Vector<T>::deduplicate() {
    int oldSize = _size;
    Rank i = 1;
    while (i < _size) {
        if (find(_elem[i], 0, i) < 0) ++i;
        else remove(i);
    }
    return oldSize - _size;
}

template <typename T>
int Vector<T>::uniquify() {
    if (_size < 2) return 0;
    int oldSize = _size;
    Rank i = 0, j = 0;
    while (++j < _size) {
        if (_elem[i] != _elem[j]) _elem[++i] = _elem[j];
    }
    _size = i + 1;
    shrink();
    return oldSize - _size;
}

template <typename T>
void Vector<T>::traverse(void (*visit)(T&)) {
    for (Rank i = 0; i < _size; ++i) visit(_elem[i]);
}

template <typename T>
template <typename VST>
void Vector<T>::traverse(VST& visit) {
    for (Rank i = 0; i < _size; ++i) visit(_elem[i]);
}

template <typename T>
Rank Vector<T>::insert(Rank r, T const& e) {
    expand();
    for (Rank i = _size; i > r; --i) _elem[i] = _elem[i - 1];
    _elem[r] = e;
    ++_size;
    return r;
}

template <typename T>
T Vector<T>::remove(Rank r) {
    if (r < 0 || r >= _size) throw std::out_of_range("Index out of range");
    T e = _elem[r];
    for (Rank i = r; i < _size - 1; ++i) _elem[i] = _elem[i + 1];
    --_size;
    shrink();
    return e;
}

template <typename T>
int Vector<T>::remove(Rank lo, Rank hi) {
    if (lo >= hi) return 0;
    if (lo < 0 || hi > _size) throw std::out_of_range("Range out of bounds");
    while (hi < _size) _elem[lo++] = _elem[hi++];
    _size = lo;
    shrink();
    return hi - lo;
}

template <typename T>
void Vector<T>::sort(Rank lo, Rank hi) {
    if (hi - lo < 2) return;
    static std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    switch (gen() % 4) {
        case 0: selectionSort(lo, hi); break;
        case 1: mergeSort(lo, hi); break;
        case 2: heapSort(lo, hi); break;
        case 3: quickSort(lo, hi); break;
    }
}

template <typename T>
void Vector<T>::bubbleSort(Rank lo, Rank hi) {
    if (hi - lo < 2) return;
    while (!bubble(lo, hi)) --hi;
}

template <typename T>
bool Vector<T>::bubble(Rank lo, Rank hi) {
    bool sorted = true;
    for (Rank i = lo; i < hi - 1; ++i) {
        if (_elem[i] > _elem[i + 1]) {
            sorted = false;
            std::swap(_elem[i], _elem[i + 1]);
        }
    }
    return sorted;
}

template <typename T>
void Vector<T>::selectionSort(Rank lo, Rank hi) {
    for (Rank i = lo; i < hi - 1; ++i) {
        Rank minIdx = i;
        for (Rank j = i + 1; j < hi; ++j) {
            if (_elem[j] < _elem[minIdx]) minIdx = j;
        }
        if (minIdx != i) std::swap(_elem[i], _elem[minIdx]);
    }
}

template <typename T>
void Vector<T>::mergeSort(Rank lo, Rank hi) {
    if (hi - lo < 2) return;
    Rank mi = (lo + hi) >> 1;
    mergeSort(lo, mi);
    mergeSort(mi, hi);
    merge(lo, mi, hi);
}

template <typename T>
void Vector<T>::merge(Rank lo, Rank mi, Rank hi) {
    T* A = _elem + lo;
    int lb = mi - lo;
    T* B = new T[lb];
    for (Rank i = 0; i < lb; ++i) B[i] = A[i];

    int lc = hi - mi;
    T* C = _elem + mi;

    for (Rank i = 0, j = 0, k = 0; (j < lb) || (k < lc);) {
        if ((j < lb) && (!(k < lc) || (B[j] <= C[k]))) A[i++] = B[j++];
        if ((k < lc) && (!(j < lb) || (C[k] < B[j]))) A[i++] = C[k++];
    }

    delete[] B;
}

template <typename T>
void Vector<T>::heapSort(Rank lo, Rank hi) {
    for (Rank i = (hi - 1 - lo) / 2 + lo; i >= lo; --i) heapify(lo, hi, i);
    for (Rank i = hi - 1; i > lo; --i) {
        std::swap(_elem[lo], _elem[i]);
        heapify(lo, i, lo);
    }
}

template <typename T>
void Vector<T>::heapify(Rank lo, Rank hi, Rank root) {
    Rank largest = root;
    Rank left = lo + 2 * (root - lo) + 1;
    Rank right = lo + 2 * (root - lo) + 2;

    if (left < hi && _elem[left] > _elem[largest]) largest = left;
    if (right < hi && _elem[right] > _elem[largest]) largest = right;
    if (largest != root) {
        std::swap(_elem[root], _elem[largest]);
        heapify(lo, hi, largest);
    }
}

template <typename T>
void Vector<T>::quickSort(Rank lo, Rank hi) {
    if (hi - lo < 2) return;
    Rank pivot = partition(lo, hi);
    quickSort(lo, pivot);
    quickSort(pivot + 1, hi);
}

template <typename T>
Rank Vector<T>::partition(Rank lo, Rank hi) {
    static std::mt19937 gen(static_cast<unsigned>(time(nullptr)));
    Rank pivot = lo + gen() % (hi - lo);
    std::swap(_elem[pivot], _elem[lo]);
    T pivotVal = _elem[lo];
    Rank i = lo, j = hi;
    
    while (true) {
        while (++i < hi && _elem[i] < pivotVal);
        while (--j > lo && pivotVal < _elem[j]);
        if (i >= j) break;
        std::swap(_elem[i], _elem[j]);
    }
    std::swap(_elem[lo], _elem[j]);
    return j;
}

// Vector 构造函数重载实现
template <typename T>
Vector<T>::Vector(T const* A, Rank n) {
    copyFrom(A, 0, n);
}

template <typename T>
Vector<T>::Vector(T const* A, Rank lo, Rank hi) {
    copyFrom(A, lo, hi);
}

template <typename T>
Vector<T>::Vector(Vector const& V) {
    copyFrom(V._elem, 0, V._size);
}

template <typename T>
Vector<T>::Vector(Vector const& V, Rank lo, Rank hi) {
    copyFrom(V._elem, lo, hi);
}

// ==================== 栈类模板 ====================

template <typename T>
class Stack : public Vector<T> {
public:
    void push(T const& e) { this->insert(this->size(), e); }
    T pop() { return this->remove(this->size() - 1); }
    T& top() { return (*this)[this->size() - 1]; }
    bool empty() const { return this->size() == 0; }
};

// ==================== 队列类模板 ====================

template <typename T>
class Queue {
protected:
    Vector<T> _elem;
    Rank _head;
    Rank _tail;
    
public:
    Queue() : _head(0), _tail(0) {}
    
    void enqueue(T const& e) { _elem.insert(_tail++, e); }
    T dequeue() {
        if (empty()) throw std::underflow_error("Queue is empty");
        return _elem[_head++];
    }
    T& front() {
        if (empty()) throw std::underflow_error("Queue is empty");
        return _elem[_head];
    }
    bool empty() const { return _head >= _tail; }
    Rank size() const { return _tail - _head; }
    void clear() { _head = _tail = 0; }
};

// ==================== 二叉树结构 ====================

template <typename T>
struct BinNode {
    T data;
    BinNode<T>* parent;
    BinNode<T>* lc;
    BinNode<T>* rc;
    int height;
    int npl;
    enum Color { RED, BLACK } color;

    BinNode(T const& e = T(), BinNode<T>* p = nullptr, 
            BinNode<T>* lc = nullptr, BinNode<T>* rc = nullptr,
            int h = 0, Color c = RED)
        : data(e), parent(p), lc(lc), rc(rc), height(h), npl(1), color(c) {}
    
    BinNode<T>* insertAsLC(T const& e) { return lc = new BinNode(e, this); }
    BinNode<T>* insertAsRC(T const& e) { return rc = new BinNode(e, this); }
    
    int size() const {
        int s = 1;
        if (lc) s += lc->size();
        if (rc) s += rc->size();
        return s;
    }
};

// 二叉树宏定义（用于便捷访问）
#define IsRoot(x) (!(x).parent)
#define IsLChild(x) (!IsRoot(x) && (&(x) == (x).parent->lc))
#define IsRChild(x) (!IsRoot(x) && (&(x) == (x).parent->rc))
#define HasLChild(x) ((x).lc)
#define HasRChild(x) ((x).rc)
#define HasChild(x) (HasLChild(x) || HasRChild(x))
#define HasBothChild(x) (HasLChild(x) && HasRChild(x))
#define IsLeaf(x) (!HasChild(x))

template <typename T>
class BinTree {
protected:
    int _size;
    BinNode<T>* _root;
    
    virtual int updateHeight(BinNode<T>* x) {
        return x->height = 1 + std::max(
            x->lc ? x->lc->height : -1,
            x->rc ? x->rc->height : -1
        );
    }
    
    void updateHeightAbove(BinNode<T>* x) {
        while (x) { updateHeight(x); x = x->parent; }
    }
    
    static int removeAt(BinNode<T>* x) {
        if (!x) return 0;
        int n = 1 + removeAt(x->lc) + removeAt(x->rc);
        delete x;
        return n;
    }

    // 遍历辅助函数
    template <typename VST>
    void travPre(BinNode<T>* x, VST& visit) {
        if (!x) return;
        visit(x->data);
        travPre(x->lc, visit);
        travPre(x->rc, visit);
    }
    
    template <typename VST>
    void travIn(BinNode<T>* x, VST& visit) {
        if (!x) return;
        travIn(x->lc, visit);
        visit(x->data);
        travIn(x->rc, visit);
    }
    
    template <typename VST>
    void travPost(BinNode<T>* x, VST& visit) {
        if (!x) return;
        travPost(x->lc, visit);
        travPost(x->rc, visit);
        visit(x->data);
    }

public:
    BinTree() : _size(0), _root(nullptr) {}
    ~BinTree() { if (_root) remove(_root); }
    
    int size() const { return _size; }
    bool empty() const { return !_root; }
    BinNode<T>* root() const { return _root; }
    
    BinNode<T>* insertAsRoot(T const& e) {
        _size = 1;
        return _root = new BinNode<T>(e);
    }
    
    BinNode<T>* insertAsLC(BinNode<T>* x, T const& e) {
        ++_size;
        x->insertAsLC(e);
        updateHeightAbove(x);
        return x->lc;
    }
    
    BinNode<T>* insertAsRC(BinNode<T>* x, T const& e) {
        ++_size;
        x->insertAsRC(e);
        updateHeightAbove(x);
        return x->rc;
    }
    
    int remove(BinNode<T>* x) {
        if (!x) return 0;
        if (IsRoot(*x)) _root = nullptr;
        else if (IsLChild(*x)) x->parent->lc = nullptr;
        else x->parent->rc = nullptr;
        
        updateHeightAbove(x->parent);
        int n = removeAt(x);
        _size -= n;
        return n;
    }
    
    template <typename VST>
    void travPre(VST& visit) { travPre(_root, visit); }
    
    template <typename VST>
    void travIn(VST& visit) { travIn(_root, visit); }
    
    template <typename VST>
    void travPost(VST& visit) { travPost(_root, visit); }
};

// ==================== 图类模板 ====================

enum class VStatus { UNDISCOVERED, DISCOVERED, VISITED };
enum class EType { UNDETERMINED, TREE, CROSS, FORWARD, BACKWARD };

template <typename Tv, typename Te>
class Graph {
protected:
    int n, e;
    
    void reset() {
        for (int i = 0; i < n; ++i) {
            status(i) = VStatus::UNDISCOVERED;
            dTime(i) = -1;
            fTime(i) = -1;
            parent(i) = -1;
            priority(i) = INT_MAX;
            for (int j = 0; j < n; ++j)
                if (exists(i, j)) type(i, j) = EType::UNDETERMINED;
        }
    }

    void BFS(int v, int& clock);
    void DFS(int v, int& clock);

public:
    // 顶点接口
    virtual int insert(Tv const&) = 0;
    virtual Tv remove(int) = 0;
    virtual Tv& vertex(int) = 0;
    virtual int inDegree(int) = 0;
    virtual int outDegree(int) = 0;
    virtual int firstNbr(int) = 0;
    virtual int nextNbr(int, int) = 0;
    virtual VStatus& status(int) = 0;
    virtual int& dTime(int) = 0;
    virtual int& fTime(int) = 0;
    virtual int& parent(int) = 0;
    virtual int& priority(int) = 0;
    
    // 边接口
    virtual bool exists(int, int) = 0;
    virtual void insert(Te const&, int, int, int) = 0;
    virtual Te remove(int, int) = 0;
    virtual EType& type(int, int) = 0;
    virtual Te& edge(int, int) = 0;
    virtual int& weight(int, int) = 0;
    
    // 遍历算法
    void bfs(int s);
    void dfs(int s);
};

template <typename Tv, typename Te>
void Graph<Tv, Te>::bfs(int s) {
    reset();
    int clock = 0;
    int v = s;
    do {
        if (VStatus::UNDISCOVERED == status(v))
            BFS(v, clock);
    } while (s != (v = (++v % n)));
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::BFS(int v, int& clock) {
    Queue<int> Q;
    status(v) = VStatus::DISCOVERED;
    Q.enqueue(v);
    while (!Q.empty()) {
        v = Q.dequeue();
        dTime(v) = ++clock;
        for (int u = firstNbr(v); -1 < u; u = nextNbr(v, u)) {
            if (VStatus::UNDISCOVERED == status(u)) {
                status(u) = VStatus::DISCOVERED;
                Q.enqueue(u);
                type(v, u) = EType::TREE;
                parent(u) = v;
            } else {
                type(v, u) = EType::CROSS;
            }
        }
        status(v) = VStatus::VISITED;
    }
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::dfs(int s) {
    reset();
    int clock = 0;
    int v = s;
    do {
        if (VStatus::UNDISCOVERED == status(v))
            DFS(v, clock);
    } while (s != (v = (++v % n)));
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::DFS(int v, int& clock) {
    dTime(v) = ++clock;
    status(v) = VStatus::DISCOVERED;
    for (int u = firstNbr(v); -1 < u; u = nextNbr(v, u)) {
        switch (status(u)) {
            case VStatus::UNDISCOVERED:
                type(v, u) = EType::TREE;
                parent(u) = v;
                DFS(u, clock);
                break;
            case VStatus::DISCOVERED:
                type(v, u) = EType::BACKWARD;
                break;
            default:
                type(v, u) = (dTime(v) < dTime(u)) ? EType::FORWARD : EType::CROSS;
                break;
        }
    }
    status(v) = VStatus::VISITED;
    fTime(v) = ++clock;
}

// ==================== 辅助算法与数据结构 ====================

// 八皇后结构
struct Queen {
    int x, y;
    Queen(int xx = 0, int yy = 0) : x(xx), y(yy) {}
    bool operator==(Queen const& q) const {
        return (x == q.x) || (y == q.y) || (x + y == q.x + q.y) || (x - y == q.x - q.y);
    }
    bool operator!=(Queen const& q) const { return !(*this == q); }
};

// 进制转换（使用栈）
inline void convert(Stack<char>& S, int64_t n, int base) {
    static const char digit[] = "0123456789ABCDEF";
    if (n < 0) throw std::invalid_argument("Negative numbers not supported");
    if (base < 2 || base > 16) throw std::invalid_argument("Base must be between 2 and 16");
    while (n > 0) {
        S.push(digit[n % base]);
        n /= base;
    }
}

// 括号匹配检查
inline bool paren(const char exp[], Rank lo, Rank hi) {
    Stack<char> S;
    for (Rank i = lo; i <= hi; ++i) {
        switch (exp[i]) {
            case '(': case '[': case '{': S.push(exp[i]); break;
            case ')': if (S.empty() || '(' != S.pop()) return false; break;
            case ']': if (S.empty() || '[' != S.pop()) return false; break;
            case '}': if (S.empty() || '{' != S.pop()) return false; break;
            default: break;
        }
    }
    return S.empty();
}

// ==================== 高级算法 ====================

// 众数验证
template <typename T>
bool majEleCheck(Vector<T> const& A, T const& maj) {
    int occurrence = 0;
    for (Rank i = 0; i < A.size(); ++i)
        if (A[i] == maj) ++occurrence;
    return 2 * occurrence > A.size();
}

// 众数候选（摩尔投票法）
template <typename T>
T majEleCandidate(Vector<T> const& A) {
    T maj;
    for (int c = 0, i = 0; i < A.size(); ++i) {
        if (c == 0) { maj = A[i]; c = 1; }
        else maj == A[i] ? ++c : --c;
    }
    return maj;
}

// 中位数（蛮力版）
template <typename T>
T trivialMedian(Vector<T>& S1, int lo1, int n1, Vector<T>& S2, int lo2, int n2) {
    int hi1 = lo1 + n1, hi2 = lo2 + n2;
    Vector<T> S;
    while ((lo1 < hi1) && (lo2 < hi2)) {
        if (S1[lo1] <= S2[lo2]) S.insert(S1[lo1++]);
        else S.insert(S2[lo2++]);
    }
    while (lo1 < hi1) S.insert(S1[lo1++]);
    while (lo2 < hi2) S.insert(S2[lo2++]);
    return S[(n1 + n2) / 2];
}

// 中位数（高效版）
template <typename T>
T median(Vector<T>& S1, int lo1, Vector<T>& S2, int lo2, int n) {
    if (n < 3) return trivialMedian(S1, lo1, n, S2, lo2, n);
    
    int mi1 = lo1 + n / 2, mi2 = lo2 + (n - 1) / 2;
    
    if (S1[mi1] < S2[mi2])
        return median(S1, mi1, S2, lo2, n + lo1 - mi1);
    else if (S1[mi1] > S2[mi2])
        return median(S1, lo1, S2, mi2, n + lo2 - mi2);
    else
        return S1[mi1];
}

} // namespace mystl

#endif // MYSTL_H
