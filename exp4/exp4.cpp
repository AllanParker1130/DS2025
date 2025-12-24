#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>
#include "mystl.h"  // 确保此文件路径正确

// ==================== 1. 边界框数据结构 ====================

struct BoundingBox {
    float x1, y1, x2, y2;
    float score;
    
    // 关键修复1：默认构造函数
    BoundingBox() : x1(0), y1(0), x2(0), y2(0), score(0) {}
    
    BoundingBox(float x1_, float y1_, float x2_, float y2_, float score_)
        : x1(x1_), y1(y1_), x2(x2_), y2(y2_), score(score_) {}
    
    float area() const {
        return std::max(0.0f, x2 - x1) * std::max(0.0f, y2 - y1);
    }
    
    // 关键修复2：完整比较运算符
    bool operator<(const BoundingBox& other) const { return score > other.score; }
    bool operator>(const BoundingBox& other) const { return score < other.score; }
    bool operator<=(const BoundingBox& other) const { return score >= other.score; }
    bool operator>=(const BoundingBox& other) const { return score <= other.score; }
    bool operator==(const BoundingBox& other) const { return score == other.score; }
    bool operator!=(const BoundingBox& other) const { return score != other.score; }
};

// ==================== 2. IoU计算 ====================

float calculateIoU(const BoundingBox& a, const BoundingBox& b) {
    float interX1 = std::max(a.x1, b.x1);
    float interY1 = std::max(a.y1, b.y1);
    float interX2 = std::min(a.x2, b.x2);
    float interY2 = std::min(a.y2, b.y2);
    
    float interArea = std::max(0.0f, interX2 - interX1) * std::max(0.0f, interY2 - interY1);
    float unionArea = a.area() + b.area() - interArea;
    return unionArea > 0 ? interArea / unionArea : 0.0f;
}

// ==================== 3. 数据生成器 ====================

class DataGenerator {
private:
    std::mt19937 gen{std::random_device{}()};
    std::uniform_real_distribution<float> coordDist{0.0f, 1.0f};
    std::uniform_real_distribution<float> sizeDist{0.05f, 0.2f};
    std::uniform_real_distribution<float> scoreDist{0.0f, 1.0f};

public:
    // 必须返回mystl::Vector
    mystl::Vector<BoundingBox> generateRandom(int count) {
        mystl::Vector<BoundingBox> boxes;
        for (int i = 0; i < count; ++i) {
            float x = coordDist(gen), y = coordDist(gen);
            float w = sizeDist(gen), h = sizeDist(gen);
            boxes.insert(boxes.size(), BoundingBox(x, y, x+w, y+h, scoreDist(gen)));
        }
        return boxes;
    }
    
    mystl::Vector<BoundingBox> generateClustered(int count) {
        mystl::Vector<BoundingBox> boxes;
        int numClusters = std::max(1, count / 20);
        
        mystl::Vector<std::pair<float, float>> centers;
        for (int i = 0; i < numClusters; ++i) {
            centers.insert(centers.size(), {coordDist(gen), coordDist(gen)});
        }
        
        std::normal_distribution<float> clusterDist{0.0f, 0.1f};
        for (int i = 0; i < count; ++i) {
            int idx = std::uniform_int_distribution<>{0, numClusters-1}(gen);
            float x = std::clamp(centers[idx].first + clusterDist(gen), 0.0f, 1.0f);
            float y = std::clamp(centers[idx].second + clusterDist(gen), 0.0f, 1.0f);
            float w = sizeDist(gen) * 0.5f, h = sizeDist(gen) * 0.5f;
            boxes.insert(boxes.size(), BoundingBox(x, y, x+w, y+h, scoreDist(gen)));
        }
        return boxes;
    }
};

// ==================== 4. 排序枚举器 ====================

enum class SortAlgorithm { BUBBLE, SELECTION, MERGE, HEAP, QUICK };

class NMSSorter {
public:
    template <typename T>
    static void sort(mystl::Vector<T>& vec, SortAlgorithm algo) {
        switch (algo) {
            case SortAlgorithm::BUBBLE: vec.bubbleSort(0, vec.size()); break;
            case SortAlgorithm::SELECTION: vec.selectionSort(0, vec.size()); break;
            case SortAlgorithm::MERGE: vec.mergeSort(0, vec.size()); break;
            case SortAlgorithm::HEAP: vec.heapSort(0, vec.size()); break;
            case SortAlgorithm::QUICK: vec.quickSort(0, vec.size()); break;
        }
    }
    
    static const char* name(SortAlgorithm algo) {
        switch (algo) {
            case SortAlgorithm::BUBBLE: return "BubbleSort";
            case SortAlgorithm::SELECTION: return "SelectionSort";
            case SortAlgorithm::MERGE: return "MergeSort";
            case SortAlgorithm::HEAP: return "HeapSort";
            case SortAlgorithm::QUICK: return "QuickSort";
        }
        return "Unknown";
    }
};

// ==================== 5. NMS算法 ====================

class NMSAlgorithm {
private:
    float iouThreshold;
    
public:
    explicit NMSAlgorithm(float threshold = 0.5f) : iouThreshold(threshold) {}
    
    mystl::Vector<BoundingBox> apply(mystl::Vector<BoundingBox>& boxes, SortAlgorithm algo) {
        if (boxes.empty()) return boxes;
        
        // 深拷贝
        mystl::Vector<BoundingBox> sortedBoxes;
        for (int i = 0; i < boxes.size(); ++i) {
            sortedBoxes.insert(sortedBoxes.size(), boxes[i]);
        }
        
        // 排序
        NMSSorter::sort(sortedBoxes, algo);
        
        // NMS
        mystl::Vector<BoundingBox> keep;
        mystl::Vector<bool> suppressed;
        for (int i = 0; i < sortedBoxes.size(); ++i) {
            suppressed.insert(suppressed.size(), false);
        }
        
        for (int i = 0; i < sortedBoxes.size(); ++i) {
            if (suppressed[i]) continue;
            keep.insert(keep.size(), sortedBoxes[i]);
            
            for (int j = i + 1; j < sortedBoxes.size(); ++j) {
                if (!suppressed[j] && calculateIoU(sortedBoxes[i], sortedBoxes[j]) > iouThreshold) {
                    suppressed[j] = true;
                }
            }
        }
        return keep;
    }
};

// ==================== 6. 性能测试 ====================

class PerformanceTester {
private:
    DataGenerator generator;
    NMSAlgorithm nms{0.5f};
    
    template <typename Func>
    double measureTime(Func&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start).count();
    }

public:
    void runTest(const std::string& outputFile = "nms_performance.csv") {
        std::ofstream file(outputFile);
        if (!file) {
            std::cerr << "无法创建输出文件\n";
            return;
        }
        
        file << "数据规模,分布类型,排序算法,总时间(ms),排序时间(ms),NMS时间(ms),保留框数\n";
        
        mystl::Vector<int> sizes;
        sizes.insert(sizes.size(), 100); sizes.insert(sizes.size(), 500); 
        sizes.insert(sizes.size(), 1000); sizes.insert(sizes.size(), 5000); 
        sizes.insert(sizes.size(), 10000);
        
        mystl::Vector<SortAlgorithm> algorithms;
        algorithms.insert(algorithms.size(), SortAlgorithm::BUBBLE);
        algorithms.insert(algorithms.size(), SortAlgorithm::SELECTION);
        algorithms.insert(algorithms.size(), SortAlgorithm::MERGE);
        algorithms.insert(algorithms.size(), SortAlgorithm::HEAP);
        algorithms.insert(algorithms.size(), SortAlgorithm::QUICK);
        
        mystl::Vector<mystl::Vector<BoundingBox>> randomData;
        mystl::Vector<mystl::Vector<BoundingBox>> clusterData;
        
        for (int i = 0; i < sizes.size(); ++i) {
            randomData.insert(randomData.size(), generator.generateRandom(sizes[i]));
            clusterData.insert(clusterData.size(), generator.generateClustered(sizes[i]));
        }
        
        for (int i = 0; i < sizes.size(); ++i) {
            for (int j = 0; j < algorithms.size(); ++j) {
                SortAlgorithm algo = algorithms[j];
                testDistribution(file, randomData[i], "Random", algo);
                testDistribution(file, clusterData[i], "Clustered", algo);
            }
            std::cout << "完成规模 " << sizes[i] << " 的测试\n";
        }
        
        file.close();
        std::cout << "结果已保存至 " << outputFile << "\n";
    }

private:
    void testDistribution(std::ofstream& file, mystl::Vector<BoundingBox>& boxes,
                         const std::string& distType, SortAlgorithm algo) {
        double sortTime = 0, nmsTime = 0;
        
        mystl::Vector<BoundingBox> boxesCopy;
        for (int i = 0; i < boxes.size(); ++i) {
            boxesCopy.insert(boxesCopy.size(), boxes[i]);
        }
        
        sortTime = measureTime([&]() {
            NMSSorter::sort(boxesCopy, algo);
        });
        
        mystl::Vector<BoundingBox> results;
        double totalTime = measureTime([&]() {
            results = nms.apply(boxes, algo);
        });
        nmsTime = totalTime - sortTime;
        
        file << boxes.size() << ","
             << distType << ","
             << NMSSorter::name(algo) << ","
             << std::fixed << std::setprecision(2)
             << totalTime << ","
             << sortTime << ","
             << nmsTime << ","
             << results.size() << "\n";
    }
};

// ==================== 7. 主函数 ====================

int main() {
    std::cout << "=== NMS排序算法性能测试 ===\n";
    std::cout << "确保已修改mystl.h将排序方法设为public!\n\n";
    
    PerformanceTester tester;
    tester.runTest();
    
    return 0;
}
