#ifndef ALGOS_H
#define ALGOS_H

#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <stdexcept>
#include <utility>
#include <climits>
#include <cstring>

// ==================== 基础算法函数 ====================

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

// ==================== 向量类模板 ====================

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
    
    // 构造函数、析构函数
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

// Vector成员函数实现
template <typename T>
Vector<T>::~Vector() {
    delete[] _elem;
}

template <typename T>
Vector<T>::Vector(int c, int s, T v) : _size(s), _capacity(c) {
    _elem = new T[_capacity];
    for (int i = 0; i < _size; i++) {
        _elem[i] = v;
    }
}

template <typename T>
void Vector<T>::copyFrom(T const* A, Rank lo, Rank hi) {
    _capacity = 2 * (hi - lo);
    _elem = new T[_capacity];
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
    for (Rank i = (hi - 1) / 2; i >= lo; --i) {
        heapify(lo, hi, i);
    }
    for (Rank i = hi - 1; i > lo; --i) {
        std::swap(_elem[lo], _elem[i]);
        heapify(lo, i, lo);
    }
}

template <typename T>
void Vector<T>::heapify(Rank lo, Rank hi, Rank root) {
    Rank pos = root;
    Rank lchild = 2 * (root - lo) + 1 + lo;
    Rank rchild = 2 * (root - lo) + 2 + lo;

    if (lchild < hi && _elem[lchild] > _elem[pos]) {
        pos = lchild;
    }
    if (rchild < hi && _elem[rchild] > _elem[pos]) {
        pos = rchild;
    }
    if (pos != root) {
        std::swap(_elem[root], _elem[pos]);
        heapify(lo, hi, pos);
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
    Rank pivot = lo + rand() % (hi - lo);
    std::swap(_elem[pivot], _elem[lo]);
    pivot = lo;

    for (Rank i = lo + 1; i < hi; ++i) {
        if (_elem[i] < _elem[pivot]) {
            std::swap(_elem[++pivot], _elem[i]);
        }
    }
    return pivot;
}

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
    for (Rank i = lo; i < hi - 1; i++) {
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

template <typename T>
Vector<T>::Vector(Vector<T> const& V) {
    copyFrom(V._elem, 0, V._size);
}

template <typename T>
Vector<T>::Vector(Vector<T> const& V, Rank lo, Rank hi) {
    copyFrom(V._elem, lo, hi);
}

template <typename T>
Vector<T>::Vector(T const* A, Rank n) {
    copyFrom(A, 0, n);
}

template <typename T>
Vector<T>::Vector(T const* A, Rank lo, Rank hi) {
    copyFrom(A, lo, hi);
}

// ==================== 栈与队列 ====================

template <typename T>
class Stack: public Vector<T> {
public:
    void push(T const& e) { this->insert(this->size(), e); }
    T pop() { return this->remove(this->size() - 1); }
    T& top() { return (*this)[this->size() - 1]; }
};

// 队列模板类（基于向量实现）
template <typename T>
class Queue {
protected:
    Vector<T> _elem;
    Rank _head;
    Rank _tail;
    
public:
    Queue() : _head(0), _tail(0) {}
    void enqueue(T const& e) { _elem.insert(_tail++, e); }
    T dequeue() { return _elem.remove(_head++); }
    T& front() { return _elem[_head]; }
    bool empty() const { return _head >= _tail; }
    Rank size() const { return _tail - _head; }
};

// ==================== 二叉树 ====================

#define BinNodePosi(T) BinNode<T>*
#define stature(p) ((p) ? (p)->height : -1)

typedef enum { RB_RED, RB_BLACK } RBColor;

template <typename T> struct BinNode {
    T data;
    BinNodePosi(T) parent;
    BinNodePosi(T) lc;
    BinNodePosi(T) rc;
    int height;
    int npl;
    RBColor color;

    BinNode() : parent(NULL), lc(NULL), rc(NULL), height(0), npl(1), color(RB_RED) {}
    BinNode(T e, BinNodePosi(T) p = NULL, BinNodePosi(T) lc = NULL, BinNodePosi(T) rc = NULL,
            int h = 0, int l = 1, RBColor c = RB_RED) :
        data(e), parent(p), lc(lc), rc(rc), height(h), npl(l), color(c) {}

    BinNodePosi(T) insertAsLC(T const& e) { return lc = new BinNode(e, this); }
    BinNodePosi(T) insertAsRC(T const& e) { return rc = new BinNode(e, this); }
    int size() {
        int s = 1;
        if (lc) s += lc->size();
        if (rc) s += rc->size();
        return s;
    }
};

#define IsRoot(x) (!((x).parent))
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
    BinNodePosi(T) _root;
    
    virtual int updateHeight(BinNodePosi(T) x) {
        return x->height = 1 + std::max(stature(x->lc), stature(x->rc));
    }
    
    void updateHeightAbove(BinNodePosi(T) x) {
        while (x) { updateHeight(x); x = x->parent; }
    }
    
    static int removeAt(BinNodePosi(T) x) {
        if (!x) return 0;
        int n = 1 + removeAt(x->lc) + removeAt(x->rc);
        delete x;
        return n;
    }

public:
    BinTree() : _size(0), _root(NULL) {}
    ~BinTree() { if (_root) remove(_root); }
    
    int size() const { return _size; }
    bool empty() const { return !_root; }
    BinNodePosi(T) root() const { return _root; }
    
    BinNodePosi(T) insertAsRoot(T const& e) {
        _size = 1;
        return _root = new BinNode<T>(e);
    }
    
    BinNodePosi(T) insertAsLC(BinNodePosi(T) x, T const& e) {
        _size++;
        x->insertAsLC(e);
        updateHeightAbove(x);
        return x->lc;
    }
    
    BinNodePosi(T) insertAsRC(BinNodePosi(T) x, T const& e) {
        _size++;
        x->insertAsRC(e);
        updateHeightAbove(x);
        return x->rc;
    }
    
    int remove(BinNodePosi(T) x) {
        if (IsRoot(*x)) _root = NULL;
        else if (IsLChild(*x)) x->parent->lc = NULL;
        else x->parent->rc = NULL;
        
        updateHeightAbove(x->parent);
        int n = removeAt(x);
        _size -= n;
        return n;
    }
    
    template <typename VST>
    void travPre(VST& visit) {
        if (_root) {
            visit(_root->data);
            if (_root->lc) { BinTree<T> left; left._root = _root->lc; left.travPre(visit); }
            if (_root->rc) { BinTree<T> right; right._root = _root->rc; right.travPre(visit); }
        }
    }
    
    template <typename VST>
    void travIn(VST& visit) {
        if (_root) {
            if (_root->lc) { BinTree<T> left; left._root = _root->lc; left.travIn(visit); }
            visit(_root->data);
            if (_root->rc) { BinTree<T> right; right._root = _root->rc; right.travIn(visit); }
        }
    }
    
    template <typename VST>
    void travPost(VST& visit) {
        if (_root) {
            if (_root->lc) { BinTree<T> left; left._root = _root->lc; left.travPost(visit); }
            if (_root->rc) { BinTree<T> right; right._root = _root->rc; right.travPost(visit); }
            visit(_root->data);
        }
    }
};

// ==================== 图类模板 ====================

typedef enum { UNDISCOVERED, DISCOVERED, VISITED } VStatus;
typedef enum { UNDETERMINED, TREE, CROSS, FORWARD, BACKWARD } EType;

template <typename Tv, typename Te>
class Graph {
private:
    void reset() {
        for (int i = 0; i < n; i++) {
            status(i) = UNDISCOVERED;
            dTime(i) = -1;
            fTime(i) = -1;
            parent(i) = -1;
            priority(i) = INT_MAX;
            for (int j = 0; j < n; j++)
                if (exists(i, j)) type(i, j) = UNDETERMINED;
        }
    }

protected:
    void BFS(int v, int& clock);
    void DFS(int v, int& clock);

public:
    int n, e;
    
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
        if (UNDISCOVERED == status(v))
            BFS(v, clock);
    } while (s != (v = (++v % n)));
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::BFS(int v, int& clock) {
    Queue<int> Q;
    status(v) = DISCOVERED;
    Q.enqueue(v);
    while (!Q.empty()) {
        v = Q.dequeue();
        dTime(v) = ++clock;
        for (int u = firstNbr(v); -1 < u; u = nextNbr(v, u)) {
            if (UNDISCOVERED == status(u)) {
                status(u) = DISCOVERED;
                Q.enqueue(u);
                type(v, u) = TREE;
                parent(u) = v;
            } else {
                type(v, u) = CROSS;
            }
        }
        status(v) = VISITED;
    }
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::dfs(int s) {
    reset();
    int clock = 0;
    int v = s;
    do {
        if (UNDISCOVERED == status(v))
            DFS(v, clock);
    } while (s != (v = (++v % n)));
}

template <typename Tv, typename Te>
void Graph<Tv, Te>::DFS(int v, int& clock) {
    dTime(v) = ++clock;
    status(v) = DISCOVERED;
    for (int u = firstNbr(v); -1 < u; u = nextNbr(v, u)) {
        switch (status(u)) {
            case UNDISCOVERED:
                type(v, u) = TREE;
                parent(u) = v;
                DFS(u, clock);
                break;
            case DISCOVERED:
                type(v, u) = BACKWARD;
                break;
            default:
                type(v, u) = (dTime(v) < dTime(u)) ? FORWARD : CROSS;
                break;
        }
    }
    status(v) = VISITED;
    fTime(v) = ++clock;
}

// ==================== 辅助结构与算法 ====================

struct Queen {
    int x, y;
    Queen(int xx = 0, int yy = 0) : x(xx), y(yy) {}
    bool operator==(Queen const& q) const {
        return (x == q.x) || (y == q.y) || (x + y == q.x + q.y) || (x - y == q.x - q.y);
    }
    bool operator!=(Queen const& q) const { return !(*this == q); }
};

// 简化的十进制转任意进制（迭代版）
void convert(Stack<char>& S, __int64 n, int base) {
    static char digit[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
    while (n > 0) {
        S.push(digit[n % base]);
        n /= base;
    }
}

// 括号匹配检查（迭代版）
bool paren(const char exp[], int lo, int hi) {
    Stack<char> S;
    for (int i = lo; i <= hi; i++) {
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

#endif // ALGOS_H
