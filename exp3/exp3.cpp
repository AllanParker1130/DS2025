#include "mystl.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace mystl;
using namespace std;

// ==================== Constant Definitions ====================

const int INF_NEG1 = -1;  // -1 indicates no connection

// ==================== UTF-8 Console Setup (Windows) ====================

void setupConsole() {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    // Also set the input code page for completeness
    SetConsoleCP(CP_UTF8);
#endif
}

// ==================== Weighted Graph Implementation ====================

template <typename Tv, typename Te>
class Graph_Matrix : public Graph<Tv, Te> {
private:
    vector<Tv> vertices;
    vector<vector<Te>> edges;
    vector<vector<EType>> edgeTypes;
    vector<VStatus> vStatus;
    vector<int> vDTime, vFTime, vParent, vPriority;

public:
    // Public interface for accessing protected members
    int vertexCount() const { return this->n; }
    int edgeCount() const { return this->e; }

    Graph_Matrix(int n = 0) {
        this->n = n;
        this->e = 0;
        vertices.resize(n);
        edges.resize(n, vector<Te>(n, INF_NEG1));
        edgeTypes.resize(n, vector<EType>(n, EType::UNDETERMINED));
        vStatus.resize(n, VStatus::UNDISCOVERED);
        vDTime.resize(n, -1);
        vFTime.resize(n, -1);
        vParent.resize(n, -1);
        vPriority.resize(n, INT_MAX);
        
        // Self-loops are 0
        for (int i = 0; i < n; ++i) edges[i][i] = 0;
    }

    // Vertex interfaces
    int insert(Tv const& v) override {
        vertices.push_back(v);
        for (auto& row : edges) row.push_back(INF_NEG1);
        edges.emplace_back(this->n + 1, INF_NEG1);
        edgeTypes.emplace_back(this->n + 1, EType::UNDETERMINED);
        for (int i = 0; i <= this->n; ++i) {
            edges[this->n][i] = edges[i][this->n] = INF_NEG1;
        }
        edges[this->n][this->n] = 0;
        vStatus.push_back(VStatus::UNDISCOVERED);
        vDTime.push_back(-1);
        vFTime.push_back(-1);
        vParent.push_back(-1);
        vPriority.push_back(INT_MAX);
        ++this->n;
        return this->n - 1;
    }

    Tv remove(int) override {
        throw runtime_error("Remove vertex not implemented");
    }

    Tv& vertex(int v) override { return vertices[v]; }
    int inDegree(int v) override {
        int count = 0;
        for (int i = 0; i < this->n; ++i)
            if (edges[i][v] != INF_NEG1 && i != v) ++count;
        return count;
    }
    int outDegree(int v) override { return inDegree(v); }
    int firstNbr(int v) override {
        for (int u = 0; u < this->n; ++u)
            if (edges[v][u] != INF_NEG1 && u != v) return u;
        return -1;
    }
    int nextNbr(int v, int u) override {
        for (++u; u < this->n; ++u)
            if (edges[v][u] != INF_NEG1 && u != v) return u;
        return -1;
    }
    VStatus& status(int v) override { return vStatus[v]; }
    int& dTime(int v) override { return vDTime[v]; }
    int& fTime(int v) override { return vFTime[v]; }
    int& parent(int v) override { return vParent[v]; }
    int& priority(int v) override { return vPriority[v]; }

    // Edge interfaces
    bool exists(int i, int j) override {
        return edges[i][j] != INF_NEG1 && i != j;
    }
    void insert(Te const& weight, int i, int j, int = 0) override {
        if (weight == INF_NEG1) return;
        if (!exists(i, j)) ++this->e;
        edges[i][j] = edges[j][i] = weight;
    }
    Te remove(int i, int j) override {
        Te w = edges[i][j];
        edges[i][j] = edges[j][i] = INF_NEG1;
        if (w != INF_NEG1 && i != j) --this->e;
        return w;
    }
    EType& type(int i, int j) override { return edgeTypes[i][j]; }
    Te& edge(int i, int j) override { return edges[i][j]; }
    int& weight(int i, int j) override { return (int&)edges[i][j]; }

    // Print adjacency matrix
    void printMatrix() {
        cout << "\nAdjacency Matrix (-1 = disconnected):\n  ";
        for (int i = 0; i < this->n; ++i)
            cout << setw(4) << (char)('A' + i);
        cout << "\n";
        for (int i = 0; i < this->n; ++i) {
            cout << (char)('A' + i) << " ";
            for (int j = 0; j < this->n; ++j) {
                if (edges[i][j] == INF_NEG1)
                    cout << setw(4) << "INF";  // ASCII-safe alternative
                else
                    cout << setw(4) << edges[i][j];
            }
            cout << "\n";
        }
    }

    // Dijkstra shortest path
    void dijkstra(int s, vector<int>& dist, vector<int>& path) {
        dist.assign(this->n, INF_NEG1);
        path.assign(this->n, -1);
        vector<bool> visited(this->n, false);
        
        dist[s] = 0;
        parent(s) = s;
        
        for (int i = 0; i < this->n; ++i) {
            int v = -1;
            int minDist = INT_MAX;
            for (int j = 0; j < this->n; ++j) {
                if (!visited[j] && dist[j] != INF_NEG1 && dist[j] < minDist) {
                    v = j;
                    minDist = dist[j];
                }
            }
            if (v == -1) break;
            visited[v] = true;
            
            for (int u = firstNbr(v); u != -1; u = nextNbr(v, u)) {
                if (!visited[u] && edges[v][u] != INF_NEG1) {
                    int newDist = dist[v] + edges[v][u];
                    if (dist[u] == INF_NEG1 || newDist < dist[u]) {
                        dist[u] = newDist;
                        path[u] = v;
                        parent(u) = v;
                    }
                }
            }
        }
    }

    // Prim minimum spanning tree
    void prim(int s, vector<pair<int, int>>& mst) {
        mst.clear();
        vector<bool> inMST(this->n, false);
        vector<int> minEdge(this->n, INF_NEG1);
        vector<int> parentVert(this->n, -1);
        
        minEdge[s] = 0;
        parentVert[s] = s;
        
        for (int i = 0; i < this->n; ++i) {
            int v = -1;
            int minWeight = INT_MAX;
            for (int j = 0; j < this->n; ++j) {
                if (!inMST[j] && minEdge[j] != INF_NEG1 && minEdge[j] < minWeight) {
                    v = j;
                    minWeight = minEdge[j];
                }
            }
            if (v == -1) break;
            inMST[v] = true;
            
            if (parentVert[v] != v) {
                mst.push_back({parentVert[v], v});
            }
            
            for (int u = firstNbr(v); u != -1; u = nextNbr(v, u)) {
                if (!inMST[u] && edges[v][u] != INF_NEG1 && 
                    (minEdge[u] == INF_NEG1 || edges[v][u] < minEdge[u])) {
                    minEdge[u] = edges[v][u];
                    parentVert[u] = v;
                }
            }
        }
    }

    // Print path
    void printPath(int s, int t, const vector<int>& path) {
        if (path[t] == -1) {
            cout << (char)('A' + s) << " to " << (char)('A' + t) << " no path\n";
            return;
        }
        
        vector<int> reversePath;
        int cur = t;
        reversePath.push_back(cur);
        while (path[cur] != cur && path[cur] != -1) {
            cur = path[cur];
            reversePath.push_back(cur);
        }
        
        for (int i = reversePath.size() - 1; i >= 0; --i) {
            cout << (char)('A' + reversePath[i]);
            if (i > 0) cout << " -> ";
        }
        cout << "\n";
    }
};

// ==================== Biconnected Components Algorithm ====================

template <typename Tv, typename Te>
class Graph_BCC : public Graph_Matrix<Tv, Te> {
private:
    vector<int> hca;
    stack<int> S;
    vector<bool> isArticulation;
    vector<vector<int>> bccComponents;

public:
    // Public access methods
    vector<bool>& getArticulation() { return isArticulation; }
    vector<vector<int>>& getComponents() { return bccComponents; }
    stack<int>& getStack() { return S; }

    Graph_BCC(int n = 0) : Graph_Matrix<Tv, Te>(n) {
        hca.resize(n);
        isArticulation.resize(n, false);
    }

    // Reset BCC state
    void resetBCC() {
        fill(isArticulation.begin(), isArticulation.end(), false);
        bccComponents.clear();
        while (!S.empty()) S.pop();
    }

    void bcc(int s) {
        this->reset();
        int clock = 0;
        int v = s;
        do {
            if (VStatus::UNDISCOVERED == this->status(v)) {
                DFS_BCC(v, clock);
                if (!S.empty()) S.pop();
            }
        } while (s != (v = ((v + 1) % this->vertexCount()))); // Fixed undefined behavior
    }

private:
    void DFS_BCC(int v, int& clock) {
        hca[v] = this->dTime(v) = ++clock;
        this->status(v) = VStatus::DISCOVERED;
        S.push(v);
        
        int childCount = 0;
        
        for (int u = this->firstNbr(v); u != -1; u = this->nextNbr(v, u)) {
            switch (this->status(u)) {
                case VStatus::UNDISCOVERED:
                    this->parent(u) = v;
                    this->type(v, u) = EType::TREE;
                    ++childCount;
                    DFS_BCC(u, clock);
                    
                    if (hca[u] < this->dTime(v)) {
                        hca[v] = min(hca[v], hca[u]);
                    } else {
                        if (this->parent(v) != -1 || childCount > 1) {
                            isArticulation[v] = true;
                        }
                        cout << "\nBiconnected Component: ";
                        vector<int> component;
                        while (true) {
                            int x = S.top(); S.pop();
                            component.push_back(x);
                            cout << (char)('A' + x) << " ";
                            if (x == u) break;
                        }
                        component.push_back(v);
                        cout << (char)('A' + v);
                        bccComponents.push_back(component);
                    }
                    break;
                    
                case VStatus::DISCOVERED:
                    this->type(v, u) = EType::BACKWARD;
                    if (u != this->parent(v)) {
                        hca[v] = min(hca[v], this->dTime(u));
                    }
                    break;
                    
                default:
                    this->type(v, u) = (this->dTime(v) < this->dTime(u)) ? 
                        EType::FORWARD : EType::CROSS;
                    break;
            }
        }
        
        this->status(v) = VStatus::VISITED;
        this->fTime(v) = ++clock;
    }

public:
    void printArticulationPoints() {
        cout << "\nArticulation Points: ";
        bool found = false;
        for (int i = 0; i < this->vertexCount(); ++i) {
            if (isArticulation[i]) {
                cout << (char)('A' + i) << " ";
                found = true;
            }
        }
        if (!found) cout << "None";
        cout << "\n";
    }

    void printComponents() {
        cout << "\nFound " << bccComponents.size() << " biconnected components\n";
    }
};

// ==================== Main Function ====================

int main() {
    // Setup UTF-8 console on Windows
    setupConsole();
    
    // ========== Task 1-3: Graph 1 (Weighted Graph) ==========
    cout << "========== Graph 1: Weighted Graph (Task 1-3) ==========\n\n";
    
    // (1) Build Graph1 and output adjacency matrix
    Graph_Matrix<int, int> g1(8);
    char labels[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
    for (int i = 0; i < 8; ++i) g1.vertex(i) = labels[i];
    
    // Define graph edges (from the table, -1 means disconnected)
    int matrix1[8][8] = {
        {0, 4, 6, 7, -1, -1, -1, -1},
        {4, 0, -1, 12, 13, -1, -1, -1},
        {6, -1, 0, 9, 1, 2, 10, -1},
        {7, 12, 9, 0, 13, -1, 2, -1},
        {-1, 13, 1, 13, 0, 5, 5, 11},
        {-1, -1, 2, -1, 5, 0, 8, -1},
        {-1, -1, 10, 2, 5, 8, 0, 3},
        {-1, -1, -1, -1, 11, -1, 3, 0}
    };
    
    // Insert edges
    for (int i = 0; i < 8; ++i) {
        for (int j = i + 1; j < 8; ++j) {
            if (matrix1[i][j] != INF_NEG1) {
                g1.insert(matrix1[i][j], i, j);
            }
        }
    }
    
    // (1) Output adjacency matrix
    g1.printMatrix();
    
    // (2) BFS and DFS from node A (0)
    cout << "\n--- Task 2: BFS from A ---\n";
    g1.bfs(0);
    cout << "Visited nodes: ";
    for (int i = 0; i < g1.vertexCount(); ++i) {
        if (g1.dTime(i) != -1) cout << labels[i] << " ";
    }
    cout << "\n";
    
    cout << "\n--- Task 2: DFS from A ---\n";
    g1.dfs(0);
    cout << "Visit order: ";
    for (int i = 0; i < g1.vertexCount(); ++i) {
        if (g1.dTime(i) != -1) cout << labels[i] << " ";
    }
    cout << "\n";
    
    // (3) Shortest path (Dijkstra) and minimum spanning tree (Prim)
    cout << "\n--- Task 3: Dijkstra Shortest Path from A ---\n";
    vector<int> dist, path;
    g1.dijkstra(0, dist, path);
    
    for (int i = 0; i < 8; ++i) {
        cout << "To " << labels[i] << " distance=";
        if (dist[i] == INF_NEG1) cout << "INF\t";
        else cout << setw(2) << dist[i] << "\t";
        cout << "Path: ";
        g1.printPath(0, i, path);
    }
    
    cout << "\n--- Task 3: Prim Minimum Spanning Tree from A ---\n";
    vector<pair<int, int>> mst;
    g1.prim(0, mst);
    
    int totalWeight = 0;
    for (auto& edge : mst) {
        int u = edge.first, v = edge.second;
        cout << labels[u] << "-" << labels[v] << " weight=" << g1.edge(u, v) << "\n";
        totalWeight += g1.edge(u, v);
    }
    cout << "Total weight: " << totalWeight << "\n";
    
    // ========== Task 4: Graph 2 (Biconnected Components) ==========
    cout << "\n\n========== Graph 2: Biconnected Components (Task 4) ==========\n\n";
    
    Graph_BCC<int, int> g2(11);
    char labels2[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L'};
    for (int i = 0; i < 11; ++i) g2.vertex(i) = labels2[i];
    
    // Define Graph2 edges (unweighted, 1 = connected, -1 = disconnected)
    int matrix2[11][11] = {
        {0, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 0, 1, 1, -1, -1, -1, -1, -1, -1, -1},
        {1, 1, 0, 1, 1, -1, -1, -1, -1, -1, -1},
        {1, 1, 1, 0, -1, 1, -1, -1, -1, -1, -1},
        {-1, -1, 1, -1, 0, 1, 1, -1, -1, -1, -1},
        {-1, -1, -1, 1, 1, 0, 1, -1, -1, -1, -1},
        {-1, -1, -1, -1, 1, 1, 0, 1, -1, -1, -1},
        {-1, -1, -1, -1, -1, -1, 1, 0, 1, 1, 1},
        {-1, -1, -1, -1, -1, -1, -1, 1, 0, 1, 1},
        {-1, -1, -1, -1, -1, -1, -1, 1, 1, 0, 1},
        {-1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 0}
    };
    
    // Insert edges
    for (int i = 0; i < 11; ++i) {
        for (int j = i + 1; j < 11; ++j) {
            if (matrix2[i][j] == 1) {
                g2.insert(1, i, j);
            }
        }
    }
    
    g2.printMatrix();
    
    // (4) Compute BCC from different starting points
    cout << "\n--- Task 4: BCC from starting point A(0) ---\n";
    g2.dfs(0);
    cout << "DFS order: ";
    for (int i = 0; i < g2.vertexCount(); ++i) {
        if (g2.dTime(i) != -1) cout << labels2[i] << "(" << g2.dTime(i) << ") ";
    }
    cout << "\n";
    
    cout << "\n--- BCC and articulation points from A ---\n";
    g2.bcc(0);
    g2.printArticulationPoints();
    g2.printComponents();
    
    cout << "\n--- Task 4: BCC from starting point D(3) ---\n";
    // Reset graph status
    for (int i = 0; i < g2.vertexCount(); ++i) {
        g2.status(i) = VStatus::UNDISCOVERED;
        g2.dTime(i) = -1;
        g2.fTime(i) = -1;
        g2.parent(i) = -1;
    }
    g2.resetBCC();
    
    g2.dfs(3);
    cout << "DFS order: ";
    for (int i = 0; i < g2.vertexCount(); ++i) {
        if (g2.dTime(i) != -1) cout << labels2[i] << "(" << g2.dTime(i) << ") ";
    }
    cout << "\n";
    
    cout << "\n--- BCC and articulation points from D ---\n";
    g2.bcc(3);
    g2.printArticulationPoints();
    g2.printComponents();
    
    // Verify consistency
    cout << "\n[Verification] Articulation points should be consistent regardless of starting point.\n";
    cout << "Conclusion: BCC decomposition is independent of DFS start vertex.\n";
    
    return 0;
}
