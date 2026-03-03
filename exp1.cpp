#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <limits>
#include "vector.h"

using namespace std;

// 定义复数类
class Complex {
public:
    double real; // 实部
    double imag; // 虚部

    Complex(double r = 0, double i = 0) : real(r), imag(i) {}

    // 计算模
    double modulus() const {
        return sqrt(real * real + imag * imag);
    }

    // 比较操作符
    bool operator==(const Complex& other) const {
        return (real == other.real) && (imag == other.imag);
    }

    bool operator!=(const Complex& other) const {
        return !(*this == other);
    }

    bool operator<(const Complex& other) const {
        if (modulus() != other.modulus()) {
            return modulus() < other.modulus();
        }
        return real < other.real;
    }

    bool operator>(const Complex& other) const {
        return modulus() > other.modulus();
    }

    bool operator<=(const Complex& other) const {
        return modulus() <= other.modulus();
    }

    bool operator>=(const Complex& other) const {
        return modulus() >= other.modulus();
    }
};

// 比较排序效率
void testSortEfficiency(Vector<Complex>& vec) {
    // 复制向量用于不同排序算法的测试
    Vector<Complex> bubbleVec = vec;
    Vector<Complex> mergeVec = vec;

    // 记录起泡排序时间
    clock_t start = clock();
    bubbleVec.bubbleSort(0, bubbleVec.size());
    clock_t end = clock();
    double bubbleTime = static_cast<double>(end - start);
    cout << "Bubble Sort Time: " << bubbleTime << " milliseconds" << endl;

    // 记录归并排序时间
    start = clock();
    mergeVec.mergeSort(0, mergeVec.size());
    end = clock();
    double mergeTime = static_cast<double>(end - start);
    cout << "Merge Sort Time: " << mergeTime << " milliseconds" << endl;
}

// 区间查找算法
Vector<Complex> findComplexInRange(const Vector<Complex>& vec, double m1, double m2) {
    Vector<Complex> result;
    for (int i = 0; i < vec.size(); i++) {
        double mod = vec[i].modulus();
        if (mod >= m1 && mod < m2) {
            result.insert(result.size(), vec[i]);
        }
    }
    return result;
}

// 字符串计算器（手动输入表达式，支持加减乘除、括号、对数、三角函数）
double calculateExpression(const string& expression) {
    Vector<Complex> stack;
    string::size_type i = 0;
    while (i < expression.length()) {
        if (isdigit(expression[i])) {
            double num = 0;
            while (i < expression.length() && isdigit(expression[i])) {
                num = num * 10 + static_cast<double>(expression[i] - '0');
                i++;
            }
            stack.insert(stack.size(), Complex(num, 0));
        } else if (expression[i] == '(') {
            stack.insert(stack.size(), Complex(nan(""), 0)); // 使用特殊值表示左括号
            i++;
        } else if (expression[i] == ')') {
            // 处理括号内的表达式
            while (stack.size() > 1 && !isnan(stack[stack.size() - 2].real)) {
                Complex b = stack.remove(stack.size() - 1);
                Complex a = stack.remove(stack.size() - 1);
                char op = static_cast<char>(stack.remove(stack.size() - 1).real);
                Complex result;
                switch (op) {
                    case '+':
                        result.real = a.real + b.real;
                        result.imag = a.imag + b.imag;
                        break;
                    case '-':
                        result.real = a.real - b.real;
                        result.imag = a.imag - b.imag;
                        break;
                    case '*':
                        result.real = a.real * b.real - a.imag * b.imag;
                        result.imag = a.real * b.imag + a.imag * b.real;
                        break;
                    case '/':
                        result.real = (a.real * b.real + a.imag * b.imag) / (b.real * b.real + b.imag * b.imag);
                        result.imag = (a.imag * b.real - a.real * b.imag) / (b.real * b.real + b.imag * b.imag);
                        break;
                }
                stack.insert(stack.size(), result);
            }
            stack.remove(stack.size() - 1); // 移除左括号
            i++;
        } else if (expression[i] == 'l' && i + 2 < expression.length() && expression.substr(i, 3) == "log") {
            i += 3;
            Complex a = stack.remove(stack.size() - 1);
            if (a.real <= 0) {
                cout << "Error: Invalid argument for log function" << endl;
                return 0.0;
            }
            Complex result;
            result.real = log10(a.real); // 使用底数为10的对数
            stack.insert(stack.size(), result);
        } else if (expression[i] == 's' && i + 2 < expression.length() && expression.substr(i, 3) == "sin") {
            i += 3;
            Complex a = stack.remove(stack.size() - 1);
            // 将角度转换为弧度
            double rad = a.real * M_PI / 180.0;
            Complex result;
            result.real = sin(rad);
            stack.insert(stack.size(), result);
        } else if (expression[i] == 'c' && i + 2 < expression.length() && expression.substr(i, 3) == "cos") {
            i += 3;
            Complex a = stack.remove(stack.size() - 1);
            // 将角度转换为弧度
            double rad = a.real * M_PI / 180.0;
            Complex result;
            result.real = cos(rad);
            stack.insert(stack.size(), result);
        } else if (expression[i] == 't' && i + 2 < expression.length() && expression.substr(i, 3) == "tan") {
            i += 3;
            Complex a = stack.remove(stack.size() - 1);
            // 将角度转换为弧度
            double rad = a.real * M_PI / 180.0;
            if (fmod(a.real, 90.0) == 0.0 && fmod(a.real / 90.0, 2) != 0.0) {
                cout << "Error: Invalid argument for tan function" << endl;
                return 0.0;
            }
            Complex result;
            result.real = tan(rad);
            stack.insert(stack.size(), result);
        } else {
            char op = expression[i];
            i++;
            if (stack.size() < 2) {
                cout << "Error: Invalid expression" << endl;
                return 0.0;
            }
            Complex b = stack.remove(stack.size() - 1);
            Complex a = stack.remove(stack.size() - 1);
            Complex result;
            switch (op) {
                case '+':
                    result.real = a.real + b.real;
                    result.imag = a.imag + b.imag;
                    break;
                case '-':
                    result.real = a.real - b.real;
                    result.imag = a.imag - b.imag;
                    break;
                case '*':
                    result.real = a.real * b.real - a.imag * b.imag;
                    result.imag = a.real * b.imag + a.imag * b.real;
                    break;
                case '/':
                    if (b.real == 0 && b.imag == 0) {
                        cout << "Error: Division by zero" << endl;
                        return 0.0;
                    }
                    result.real = (a.real * b.real + a.imag * b.imag) / (b.real * b.real + b.imag * b.imag);
                    result.imag = (a.imag * b.real - a.real * b.imag) / (b.real * b.real + b.imag * b.imag);
                    break;
            }
            stack.insert(stack.size(), result);
        }
    }
    return stack[stack.size() - 1].real;
}

// 计算矩形最大面积（手动输入柱子高度）
int largestRectangleArea(Vector<int>& heights) {
    Vector<int> stack;
    int maxArea = 0;
    for (int i = 0; i < heights.size(); i++) {
        while (!stack.empty() && heights[i] < heights[stack[stack.size() - 1]]) {
            int height = heights[stack.remove(stack.size() - 1)];
            int width = stack.empty() ? i : i - stack[stack.size() - 1] - 1;
            maxArea = max(maxArea, height * width);
        }
        stack.insert(stack.size(), i);
    }
    while (!stack.empty()) {
        int height = heights[stack.remove(stack.size() - 1)];
        int width = stack.empty() ? heights.size() : heights.size() - stack[stack.size() - 1] - 1;
        maxArea = max(maxArea, height * width);
    }
    return maxArea;
}

int main() {
    // 测试复数向量操作
    Vector<Complex> complexVector;
    
    // 扩充数据量到不少于20个
    for (int i = 0; i < 20; i++) {
        double real = rand() % 100; // 随机实部
        double imag = rand() % 100; // 随机虚部
        complexVector.insert(i, Complex(real, imag));
    }
    
    // 置乱操作
    complexVector.unsort(0, complexVector.size());
    cout << "Unsorted complex vector:" << endl;
    for (int i = 0; i < complexVector.size(); i++) {
        cout << complexVector[i].real << "+" << complexVector[i].imag << "i" << endl;
    }

    // 查找操作（实部和虚部均相同）
    int index = complexVector.find(Complex(3, 4), 0, complexVector.size());
    cout << "Index of 3+4i: " << index << endl;

    // 插入操作
    complexVector.insert(1, Complex(0, 0));
    cout << "After inserting 0+0i at position 1:" << endl;
    for (int i = 0; i < complexVector.size(); i++) {
        cout << complexVector[i].real << "+" << complexVector[i].imag << "i" << endl;
    }

    // 删除操作
    complexVector.remove(1);
    cout << "After removing position 1:" << endl;
    for (int i = 0; i < complexVector.size(); i++) {
        cout << complexVector[i].real << "+" << complexVector[i].imag << "i" << endl;
    }

    // 唯一化操作
    complexVector.deduplicate();
    cout << "After deduplication:" << endl;
    for (int i = 0; i < complexVector.size(); i++) {
        cout << complexVector[i].real << "+" << complexVector[i].imag << "i" << endl;
    }

    // 比较排序效率（计时单位为毫秒）
    testSortEfficiency(complexVector);

    // 区间查找算法
    Vector<Complex> rangeVec = findComplexInRange(complexVector, 2.0, 5.0);
    cout << "Complex numbers with modulus in [2.0, 5.0):" << endl;
    for (int i = 0; i < rangeVec.size(); i++) {
        cout << rangeVec[i].real << "+" << rangeVec[i].imag << "i" << endl;
    }

    // 基于栈数据结构实现字符串计算器（手动输入表达式）
    string expression;
    cout << "Please enter an expression (e.g., log(100)+sin(30)-tan(45)/cos(60)): ";
    cin >> expression;
    try {
        double result = calculateExpression(expression);
        cout << "Result of expression " << expression << ": " << result << endl;
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    // 计算矩形最大面积（手动输入柱子高度）
    Vector<int> heights;
    cout << "Please enter heights of the bars (press Enter when done):" << endl;
    while (true) {
        int height;
        cout << "Enter a height (or press Enter to finish): ";
        if (cin.peek() == '\n') {
            cin.ignore();
            break;
        }
        cin >> height;
        heights.insert(heights.size(), height);
    }
    cout << "Largest rectangle area: " << largestRectangleArea(heights) << endl;

    return 0;
}
