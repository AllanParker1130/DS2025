#include "mystl.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <map>

// ==================== Bitmap类（来自实验说明） ====================
class Bitmap {
private:
    unsigned char* M;
    int N;
    
public:
    Bitmap(int n = 8) {
        M = new unsigned char[N = (n + 7) / 8];
        memset(M, 0, N);
    }
    
    ~Bitmap() { delete[] M; }
    
    void set(int k) {
        if (k >= 8 * N) return;
        M[k >> 3] |= (0x80 >> (k & 0x07));
    }
    
    void clear(int k) {
        if (k >= 8 * N) return;
        M[k >> 3] &= ~(0x80 >> (k & 0x07));
    }
    
    bool test(int k) {
        return k < 8 * N && M[k >> 3] & (0x80 >> (k & 0x07));
    }
    
    char* bits2string(int n) {
        char* s = new char[n + 1];
        s[n] = '\0';
        for (int i = 0; i < n; i++) s[i] = test(i) ? '1' : '0';
        return s;
    }
};

// ==================== Huffman编码结构 ====================
struct HuffChar {
    char ch;
    int weight;
    HuffChar(char c = '^', int w = 0) : ch(c), weight(w) {}
    bool operator<(const HuffChar& other) const { return weight > other.weight; }
};

// 全局编码表
std::map<char, std::string> huffmanTable;

// 统计字母频率（不区分大小写）
int* statistics(const char* filename) {
    int* freq = new int[26]();
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }
    
    char ch;
    while (file.get(ch)) {
        if (std::isalpha(ch)) {
            freq[std::toupper(ch) - 'A']++;
        }
    }
    file.close();
    return freq;
}

// 递归生成Huffman编码
void generateCodes(BinNode<HuffChar>* node, Bitmap& code, int depth) {
    if (!node) return;
    
    // 叶子节点，生成编码
    if (!node->lc && !node->rc) {
        huffmanTable[node->data.ch] = code.bits2string(depth);
        return;
    }
    
    // 左子树编码为0
    if (node->lc) {
        code.clear(depth);
        generateCodes(node->lc, code, depth + 1);
    }
    
    // 右子树编码为1
    if (node->rc) {
        code.set(depth);
        generateCodes(node->rc, code, depth + 1);
    }
}

// ==================== 主函数 ====================
int main() {
    const char* filename = "I_have_a_dream.txt";
    
    // 检查文件
    std::ifstream test(filename);
    if (!test) {
        std::cerr << "Error: Please create file " << filename << " with the speech text" << std::endl;
        return 1;
    }
    test.close();
    
    std::cout << "========== Huffman Coding Experiment ==========\n\n";
    
    // 1. 统计字母频率
    std::cout << "Step 1: Counting letter frequencies..." << std::endl;
    int* freq = statistics(filename);
    
    // 统计有效字母总数
    int totalLetters = 0;
    for (int i = 0; i < 26; i++) {
        if (freq[i] > 0) {
            totalLetters += freq[i];
        }
    }
    
    std::cout << "Letter frequency statistics:" << std::endl;
    for (int i = 0; i < 26; i++) {
        if (freq[i] > 0) {
            std::cout << (char)('A' + i) << ": " << freq[i] << "  ";
        }
    }
    std::cout << "\nTotal letters: " << totalLetters << "\n\n";
    
    // 2. 构建Huffman森林
    std::cout << "Step 2: Building Huffman forest..." << std::endl;
    Vector<BinNode<HuffChar>*> forest;
    
    for (int i = 0; i < 26; i++) {
        if (freq[i] > 0) {
            forest.insert(forest.size(), new BinNode<HuffChar>(HuffChar('A' + i, freq[i])));
        }
    }
    std::cout << "Initial forest size: " << forest.size() << " trees\n\n";
    
    // 3. 构建Huffman树（选择权重最小的两棵树合并）
    std::cout << "Step 3: Building Huffman tree..." << std::endl;
    while (forest.size() > 1) {
        // 找到权重最小的两个节点
        Rank min1 = 0, min2 = 1;
        if (forest[1]->data.weight < forest[0]->data.weight) {
            min1 = 1; min2 = 0;
        }
        
        for (Rank i = 2; i < forest.size(); i++) {
            int w = forest[i]->data.weight;
            if (w < forest[min1]->data.weight) {
                min2 = min1; min1 = i;
            } else if (w < forest[min2]->data.weight) {
                min2 = i;
            }
        }
        
        // 合并为新树
        BinNode<HuffChar>* parent = new BinNode<HuffChar>(
            HuffChar('^', forest[min1]->data.weight + forest[min2]->data.weight)
        );
        parent->lc = forest[min1];
        parent->rc = forest[min2];
        forest[min1]->parent = parent;
        forest[min2]->parent = parent;
        
        // 替换min1为父节点，删除min2
        forest[min1] = parent;
        forest.remove(min2);
    }
    std::cout << "Huffman tree built successfully!\n\n";
    
    // 4. 生成编码表
    std::cout << "Step 4: Generating Huffman code table..." << std::endl;
    BinTree<HuffChar> huffTree;
    huffTree.insertAsRoot(forest[0]->data);
    huffTree.root()->lc = forest[0]->lc;
    huffTree.root()->rc = forest[0]->rc;
    
    Bitmap code;
    generateCodes(huffTree.root(), code, 0);
    
    std::cout << "Code table generated!\n\n";
    
    // 5. 输出编码表
    std::cout << "========== Huffman Code Table ==========" << std::endl;
    for (auto& pair : huffmanTable) {
        std::cout << pair.first << ": " << pair.second << std::endl;
    }
    
    // 6. 单词编码测试（使用自定义 Vector）
    std::cout << "\n========== Word Encoding Test ==========" << std::endl;
    Vector<std::string> testWords;
    testWords.insert(testWords.size(), std::string("dream"));
    testWords.insert(testWords.size(), std::string("freedom"));
    testWords.insert(testWords.size(), std::string("equality"));
    testWords.insert(testWords.size(), std::string("hope"));
    testWords.insert(testWords.size(), std::string("justice"));
    
    for (int i = 0; i < testWords.size(); i++) {
        const std::string& word = testWords[i];
        std::string encoded;
        for (char c : word) {
            encoded += huffmanTable[std::toupper(c)];
        }
        std::cout << word << " -> " << encoded << std::endl;
        std::cout << "  (Original ASCII: " << word.length() * 8 << " bits, Huffman: " << encoded.length() << " bits)\n";
    }
    
    // 7. 整体压缩率统计
    std::cout << "\n========== Overall Compression Statistics ==========" << std::endl;
    int originalBits = totalLetters * 8;
    int huffmanBits = 0;
    for (int i = 0; i < 26; i++) {
        if (freq[i] > 0) {
            char ch = 'A' + i;
            huffmanBits += freq[i] * huffmanTable[ch].length();
        }
    }
    double ratio = (1.0 - (double)huffmanBits / originalBits) * 100;
    std::cout << "Total letters: " << totalLetters << std::endl;
    std::cout << "Original size: " << originalBits << " bits" << std::endl;
    std::cout << "Huffman encoded: " << huffmanBits << " bits" << std::endl;
    std::cout << "Compression ratio: " << ratio << "%" << std::endl;
    
    delete[] freq;
    std::cout << "\nExperiment completed!" << std::endl;
    return 0;
}
