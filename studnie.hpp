#ifndef STUDNIE_H
#define STUDNIE_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <random>

typedef long long int ll;

struct point
{
    ll idx;
    ll x;
    ll y;
};

struct task {
    ll mult;
    std::vector<point> studnie;
    std::vector<point> domy; 
};

struct result {
    std::vector< std::pair<ll, ll> > matching;
    long double total_cost;
};

long double distance(const point& p1, const point &p2) noexcept;

std::vector<task> load_file(std::string filename);
void write_results(std::string filename, std::vector<result> &res);
result mockProgram(const task &problem);

typedef std::pair<point, size_t> Studnia;
typedef std::pair<point, size_t> Dom;
// typedef std::pair<point, size_t> Vertex;
class BipartiteMatching {
public:
    struct Edge {
        size_t studnia;
        size_t dom;
    };

    //static const bool (*edge_comp)(Edge e1, Edge e2);
    static bool edge_comp(Edge &e1, Edge &e2);
private:
    const size_t N;

    std::vector<Edge> m_matching;
    std::vector<bool> m_matched_studnie;
    std::vector<bool> m_matched_domy;

public:
    BipartiteMatching(size_t N, std::vector<Edge> edges);

    bool is_perfect(void) const;
    // BipartiteMatching symm_diff(const BipartiteMatching &m);
    void symm_diff(std::vector<Edge> &m);

    bool is_matched_studnia(size_t s) const;
    bool is_matched_dom(size_t d) const;

    const std::vector<Edge>& matching(void);

};

struct eq_vert {
    bool in_Ge;
    bool to_D;
};

class Graf
{
private:
    const std::vector<Studnia>m_studnie;
    const std::vector<Dom> m_domy;
    const std::vector<std::vector<long double>> m_weights;
    const ll m_mult;
    const size_t N;

    std::vector<long double> m_etyk_domy;
    std::vector<long double> m_etyk_studnie;

    // G_e + direcions in G_{M,e}
    std::vector<std::vector<eq_vert>> m_eq_subraph;

    explicit Graf(
        const std::vector<Studnia> s,
        const std::vector<Dom> d,
        const ll mult,
        const std::vector<std::vector<long double> > c
    );
public:
    static Graf nowy_graf(
        const task &p
    );
    
    // funkcja etykietowania dla studni s {in} S'
    long double& es(size_t s);
    // funkcja eykietowania dla domów d {in} D
    long double& ed(size_t d);

    // wielkość zbioru domów D == wielkość rozszerzonego zbioru studni S'
    const size_t size(void) const;

    // Funkcja wag zadana przez problem
    const long double& c(size_t studnia, size_t dom) const;

    eq_vert& Ge(size_t studnia, size_t dom);

    const Studnia& studnie(size_t s) const;
    const Dom& domy(size_t d) const;

    void regenerate_Ge(std::vector<BipartiteMatching::Edge> * new_edges);
    void regenerate_Gem(BipartiteMatching &M);

    std::vector<BipartiteMatching::Edge> from_S_to_D(std::vector<BipartiteMatching::Edge> edges);

    void print(void) const noexcept;

    ~Graf();
};

result studnie(const task &problem);

result studnie2(const task &problem);

std::vector<BipartiteMatching::Edge> Find_Augmneting_Path(Graf &G, BipartiteMatching &M);



#endif // STUDNIE_H