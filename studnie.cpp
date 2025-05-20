#include "studnie.hpp"
#include <iostream>
#include <fstream>
#include <queue>

Studnia (&nowa_stunia)(const point &__p, const size_t &__s) = 
    std::make_pair<const point &, const size_t &>;

Dom (&nowy_dom)(const point &__p, const size_t &__s) = 
    std::make_pair<const point &, const size_t &>;

long double distance(const point& p1, const point &p2) noexcept {
    long double dx = (long double)( p2.x - p1.x);
    long double dy = (long double)( p2.y - p1.y);
    return std::sqrt( dx * dx + dy * dy );
}

result mockProgram(const task &problem) {
    ll mult = problem.mult;
    std::vector<point> domy = problem.domy;
    std::vector<point> studnie = problem.studnie;
    
    long double total_cost;

    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(domy.begin(), domy.end(), g);

    result res;

    res.total_cost = 0.0;

    for (size_t i = 0; i < studnie.size(); i++)
    {
        auto studnia = studnie[i];
        for (size_t j = 0; j < mult; j++)
        {
            auto dom = domy[i*mult + j];

            auto dystans = distance(studnia, dom);
            res.total_cost += dystans;
            res.matching.push_back(std::make_pair(studnia.idx, dom.idx) );
        }
    }
    


    return res;
}

std::vector<task> load_file(std::string filename) {
    
    std::ifstream in(filename, std::ios::in);

    if (! in.is_open()) throw std::runtime_error("Failed to open input file");

    std::string line;

    std::getline(in, line);

    ll task_count = 0;

    std::sscanf(line.c_str(), "%lld", &task_count);

    std::cout << task_count << std::endl;

    std::vector<task> tasks(task_count);

    std::getline(in, line); //skip separator

    for (ll i = 0; i < task_count; i++) {
        ll no_wells, multiplicity;
        std::getline(in, line);
        sscanf(line.c_str(), "%lld %lld", &no_wells, &multiplicity);

        tasks[i].mult = multiplicity;

        std::getline(in, line); //skip separator
        
        // wells
        for (ll j = 0; j < no_wells; j++)
        {
            ll x, y;
            std::getline(in, line);
            
            sscanf(line.c_str(), "(%lld, %lld)", &x, &y);

            point p{
                .idx = j,
                .x = x,
                .y = y
            };

            tasks[i].studnie.push_back(p);
        }
       
        std::getline(in, line); //skip separator

        // homes
        for (ll j = 0; j < no_wells * multiplicity; j++)
        {
            ll x, y;
            std::getline(in, line);
            
            sscanf(line.c_str(), "(%lld, %lld)", &x, &y);

            point p{
                .idx = j,
                .x = x,
                .y = y
            };

            tasks[i].domy.push_back(p);
        }

        std::getline(in, line);
    }

    in.close();

    return tasks;
}

void write_results(std::string filename, std::vector<result> &res) {

    std::ofstream out(filename, std::ios::out | std::ios::trunc);

    out << res.size() << std::endl << std::endl;

    for (ll i = 0 ; i < res.size(); i++){

        out << res[i].total_cost << std::endl;
        
        for (const auto &m : res[i].matching) {
            out << "(" << m.first << ", " << m.second << ")" << std::endl;
        }
        out << std::endl;
    }
    
    out.close();
}

Graf::Graf(
        const std::vector<Studnia> s,
        const std::vector<Dom> d,
        const ll mult,
        const std::vector<std::vector<long double> > c
) : m_studnie(s), m_domy(d), m_mult(mult), N(s.size()), m_weights(c) {
    m_etyk_domy = std::vector<long double> (N, 0);
    m_etyk_studnie = std::vector<long double> (N, 0);
    m_eq_subraph = std::vector<std::vector<eq_vert>> (N, std::vector<eq_vert>(N, eq_vert{.in_Ge = false, .to_D = false}));
}

Graf Graf::nowy_graf(const task &p)
{
    std::vector<Studnia> s = std::vector<Studnia>();
    std::vector<Dom> d = std::vector<Dom>();

    for (size_t i = 0; i < p.studnie.size(); i++)
    {
        for (size_t j = 0; j < p.mult; j++)
        {
            s.push_back(nowa_stunia(p.studnie[i], i));
        }
    }

    for (size_t i = 0; i < p.domy.size(); i++){
        d.push_back(nowy_dom(p.domy[i], i));
    }
     
    const size_t N = s.size();

    std::vector<std::vector<long double>> weights 
        = std::vector<std::vector<long double>>(N, std::vector<long double>(N, 0));

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            weights[i][j] = - distance(s[i].first, d[j].first);
        }
    } 

    return Graf(
        s,
        d,
        p.mult,
        weights
    );
}

void Graf::print(void) const noexcept
{
    const auto _prev_cout_w = std::cout.width(4);
    std::cout << "Studnie: ";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << ": {(" << m_studnie[i].first.x << ", " 
            << m_studnie[i].first.x << "), " << m_studnie[i].second << "}, ";
    }
    
    std::cout << "\nDomy: ";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << ": {(" << m_domy[i].first.x << ", " 
            << m_domy[i].first.x << "), " << m_domy[i].second << "}, ";
    }
    std::cout << "\n\n c(s,d) = \n";

    
    std::cout << " d| ";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << " ";
    }
    
    std::cout << "\ns---";
    for (size_t i = 0; i < N; i++)
    {
        std::cout << "--";
    }
    std::cout << '\n';
    for (size_t i = 0; i < N; i++)
    {
        std::cout << i << " | ";
        for (size_t j = 0; j < N; j++)
        {
            std::cout << m_weights[i][j] << " ";
        }
        std::cout << '\n';
    }
    
    std::cout.width(_prev_cout_w);
    std::cout.flush();
}

Graf::~Graf() {
}

const Studnia& Graf::studnie(size_t s) const {
    return m_studnie[s];
}

const Dom& Graf::domy(size_t d) const {
    return m_domy[d];
}

const long double& Graf::c(size_t studnia, size_t dom) const {
    return m_weights[studnia][dom];
}

long double& Graf::es(size_t s) { return m_etyk_studnie[s]; }
long double& Graf::ed(size_t s) { return m_etyk_domy[s]; }
eq_vert& Graf::Ge(size_t studnia, size_t dom) {
    return m_eq_subraph[studnia][dom];
}

const size_t Graf::size(void) const {
    return N;
}

void Graf::regenerate_Ge(std::vector<BipartiteMatching::Edge> * new_edges = nullptr) {
    for (size_t s = 0; s < N; s++)
    {
        for (size_t d = 0; d < N; d++)
        {
            if( es(s) + ed(d) == c(s, d) ) {
                if(new_edges != nullptr && m_eq_subraph[s][d].in_Ge == false) {
                    new_edges->push_back(BipartiteMatching::Edge{.studnia = s, .dom = d});
                }
                m_eq_subraph[s][d].in_Ge = true;
            }
            else {
                m_eq_subraph[s][d].in_Ge = false;
            }
        }
    }
}



void Graf::regenerate_Gem(BipartiteMatching &M)
{
    for (size_t s = 0; s < N; s++)
    {
        for (size_t d = 0; d < N; d++)
        {
            m_eq_subraph[s][d].to_D = false;
        }
    }
    
    for ( const BipartiteMatching::Edge & e : M.matching()) {
        m_eq_subraph[e.studnia][e.dom].to_D = true;
    }
}

std::vector<BipartiteMatching::Edge> Graf::from_S_to_D(std::vector<BipartiteMatching::Edge> edges)
{
    std::vector<BipartiteMatching::Edge> new_dir_edges;
    for (const BipartiteMatching::Edge &e : edges) {
        if(m_eq_subraph[e.studnia][e.dom].to_D) {
            new_dir_edges.push_back(e);
        }
    }
    return new_dir_edges;
}

const std::vector<BipartiteMatching::Edge> & BipartiteMatching::matching(void) {
    return m_matching;
}

bool BipartiteMatching::is_matched_dom(size_t d) const
{
    return m_matched_domy[d];
}

bool BipartiteMatching::is_matched_studnia(size_t s) const
{
    return m_matched_studnie[s];
}

BipartiteMatching::BipartiteMatching(size_t N, std::vector<Edge> edges)
 : N(N) {
    m_matched_domy = std::vector<bool> (N, false);
    m_matched_studnie = std::vector<bool> (N, false);
    m_matching = edges;
    std::sort(m_matching.begin(), m_matching.end(), 
    edge_comp);

    for (const Edge & edge : m_matching) {
        // checks
        if (m_matched_domy[edge.dom] == true || m_matched_studnie[edge.studnia] == true) {
            // not a matching
            throw std::runtime_error("Attempted to create an invalid matching");
        }
        m_matched_domy[edge.dom] = true;
        m_matched_studnie[edge.studnia] = true;
    }

}

bool BipartiteMatching::edge_comp(Edge &e1, Edge &e2) {
    if (e1.studnia == e2.studnia) {
        return (e1.dom > e2.dom);
    }
    return e1.studnia > e2.studnia;
};

void BipartiteMatching::symm_diff(std::vector<Edge> &m) {
    
    std::sort(m.begin(), m.end(), edge_comp);

    std::vector<Edge> new_matching;
    
    std::set_symmetric_difference(
        m_matching.begin(), m_matching.end(),
        m.begin(), m.end(),
        std::back_inserter(new_matching),
        edge_comp
    );

    m_matching = new_matching;

    m_matched_domy = std::vector<bool> (N, false);
    m_matched_studnie = std::vector<bool> (N, false);
    std::sort(m_matching.begin(), m_matching.end(), 
    edge_comp);

    for (const Edge & edge : m_matching) {
        // checks
        if (m_matched_domy[edge.dom] == true || m_matched_studnie[edge.studnia] == true) {
            // not a matching
            throw std::runtime_error("Attempted to create an invalid matching");
        }
        m_matched_domy[edge.dom] = true;
        m_matched_studnie[edge.studnia] = true;
    }
}

bool BipartiteMatching::is_perfect(void) const {
    // There are 2N vertices
    // A matching is perfect if |M| == |V| / 2
    // therefore this mathcing is perfect if |M| == N
    return this->m_matching.size() == N;
}

result studnie(const task &problem) {
    Graf G = Graf::nowy_graf(problem);

    for (size_t i = 0; i < G.size(); i++)
    {
        long double max = std::numeric_limits<long double>::lowest();
        for (size_t j = 0; j < G.size(); j++)
        {
            max = std::max(G.c(i, j), max);
        }
        G.es(i) = max;
    }
    
    for (size_t i = 0; i < G.size(); i++)
    {
        G.ed(i) = 0.0;
    }
    
    G.regenerate_Ge();
    BipartiteMatching M(G.size(), {});
    G.regenerate_Gem(M);
    while (! M.is_perfect() )
    {
        std::vector<BipartiteMatching::Edge> augmenting_path = Find_Augmneting_Path(G, M);
        M.symm_diff(augmenting_path);
        G.regenerate_Ge();
        G.regenerate_Gem(M);
    }
    result Res;
    Res.total_cost = 0;
    Res.matching.clear();
    for (const BipartiteMatching::Edge & e : M.matching()) {
        Res.total_cost += G.c(e.studnia, e.dom);
        Res.matching.push_back(std::pair<ll, ll>(G.studnie(e.studnia).second, G.domy(e.dom).second ));
    }

    G.print();
    return Res;
}

struct Vertex {
    size_t idx;
    bool is_studnia;
};

std::vector<BipartiteMatching::Edge> Find_Augmneting_Path(Graf &G, BipartiteMatching &M) {

    std::vector<ll> prev_studnie(G.size(), -1);
    std::vector<ll> prev_domy(G.size(), -1);

    std::queue<Vertex> Q;
    std::vector<size_t> Fs;
    std::vector<size_t> Fd;

    size_t path_end;

    for (size_t i=0; i < G.size(); i++) { // s {in} S' not in M
        if (M.is_matched_studnia(i)) {
            continue;
        }
        Q.push(Vertex{.idx = i, .is_studnia = true});
        Fs.push_back(i);
    }

    do
    {
        if (Q.empty()) {
            long double delta = std::numeric_limits<long double>::max();
            for (const size_t &s : Fs)
            {
                size_t it = 0;
                for (size_t j = 0; j < G.size(); j++)
                {
                    if (it < Fd.size() && Fd[it] == j) // d {in} Fd
                    {
                        it++;
                        continue;
                    }
                    long double tmp1 = G.es(s);
                    long double tmp2 = G.ed(j);
                    long double tmp3 = G.c(s,j);
                    delta = std::min(delta, tmp1 + tmp2 - tmp3 );
                }
            }
            for (const size_t &s : Fs) {
                G.es(s) = G.es(s) - delta;
            }
            for (const size_t &d : Fd) {
                G.ed(d) = G.ed(d) + delta;
            }
            std::vector<BipartiteMatching::Edge> new_edges;
            G.regenerate_Ge(&new_edges);
            G.regenerate_Gem(M);
            new_edges = G.from_S_to_D(new_edges);

            for (BipartiteMatching::Edge & e: new_edges)
            {
                bool inFd = false;
                for ( const size_t & d : Fd) {
                    if(d == e.dom) {
                        inFd = true;
                        break;
                    }
                }
                if (inFd == true) continue;
                // d {not in} Fd
                prev_domy[e.dom] = e.studnia;
                if (M.is_matched_dom(e.dom)) {
                    path_end = e.dom;
                    goto found_path;
                }
                else{
                    Q.push(Vertex{.idx = e.dom, .is_studnia = false});
                    Fd.push_back(e.dom);
                }
            }
        }
        Vertex u = Q.front();
        Q.pop();

        // sąsiedzi u w Gem

        if(!u.is_studnia) {
            for (size_t i = 0; i < G.size(); i++)
            {
                const auto v = G.Ge(i, u.idx);
                if (v.in_Ge && (!v.to_D)) { // edge in Ge && target of direction
                    prev_studnie[i] = u.idx;
                    Fs.push_back(i);
                    Q.push(Vertex{.idx = i, .is_studnia = true});
                }
            }
        }
        else {
            for (size_t i = 0; i < G.size(); i++)
            {
                const auto v = G.Ge(u.idx, i);
                if (!(v.in_Ge && v.to_D)) {
                    // not in Ge || not in Gem
                    continue;
                }
                bool inFd = false;
                for ( const size_t & d : Fd) {
                    if(d == i) {
                        inFd = true;
                        break;
                    }
                }
                if (inFd == true) continue;
                prev_domy[i] = u.idx;
                if (M.is_matched_dom(i)) {
                    path_end = i;
                    goto found_path;
                }
                else {
                    Q.push(Vertex{.idx = i, .is_studnia = false});
                    Fd.push_back(i);
                }
            }
        }
    } while (true);
found_path:
    std::vector<BipartiteMatching::Edge> aug_path;
    size_t current = path_end;
    size_t prev = prev_domy[current];
    bool on_house = true;
    
    do {
        if (on_house) {
            aug_path.push_back(BipartiteMatching::Edge{.studnia = prev, .dom = current});
            current = prev;
            prev = prev_studnie[current];
            on_house = false;
        }
        else {
            aug_path.push_back(BipartiteMatching::Edge{.studnia = current, .dom = prev});
            current = prev;
            prev = prev_domy[current];
            on_house = true;
        }
    } while (prev != -1);
    return aug_path;
}

enum Mark {
    NONE  = 0b000,
    STAR  = 0b001,
    PRIME = 0b010,
    COVER = 0b100
};

Mark operator|(Mark m1, Mark m2) {
    return static_cast<Mark>(static_cast<int>(m1) | static_cast<int>(m2));
}

Mark operator&(Mark m1, Mark m2) {
    return static_cast<Mark>(static_cast<int>(m1) & static_cast<int>(m2));
}

Mark operator^(Mark m1, Mark m2) {
    return static_cast<Mark>(static_cast<int>(m1) ^ static_cast<int>(m2));
}

Mark operator~(Mark m) {
    return static_cast<Mark>(~ static_cast<int>(m));
}

Mark& operator|=(Mark &m1, Mark m2) {
    return (m1 = m1 | m2);
}

Mark& operator&=(Mark &m1, Mark m2) {
    return (m1 = m1 & m2);
}

Mark& operator^=(Mark &m1, Mark m2) {
    return (m1 = m1 ^ m2);
}

struct Zero {
    size_t stud;
    size_t dom;
};

result studnie2(const task &t) {
    const size_t N = t.domy.size();

    std::vector<std::vector<long double> > c;
    // c[studnia][dom] == waga krawędzi

    for (size_t s = 0; s < t.studnie.size(); s++)
    {
        std::vector<long double> stud_lens;
        for (size_t d = 0; d < N; d++)
        {
            stud_lens.push_back( distance(t.studnie[s], t.domy[d]));
        }
        for (size_t i = 0; i < t.mult; i++)
        {
            c.push_back(stud_lens);
        }
    }
    const std::vector<std::vector<long double> > weights(c);
    // c -> matryca wag

    std::vector<std::vector<Mark>> marks(N, std::vector<Mark>(N, NONE));
    
    std::vector<Mark> mark_rows(N, NONE);
    std::vector<Mark> mark_cols(N, NONE);

    for (size_t stud = 0; stud < N; stud++)
    // odjemujemy minimum rzędu w każdym jego elemencie
    {
        long double min = c[stud][0];
        for (size_t dom = 1; dom < N; dom++)
        {
            min = std::min<long double>(min, c[stud][dom]);
        }
        for (size_t dom = 0; dom < N; dom++)
        {
            c[stud][dom] -= min;
        }
    }
    
    for (size_t dom = 0; dom < N; dom++)
    // odejmujemy minimum kolumny w każdym jej elemencie
    {
        long double min = c[0][dom];
        for (size_t stud = 1; stud < N; stud++)
        {
            min = std::min<long double>(min, c[stud][dom]);
        }
        for (size_t stud = 0; stud < N; stud++)
        {
            c[stud][dom] -= min;
        }
    }

    std::vector<Zero> zeros;

    for (size_t stud = 0; stud < N; stud++)
    {
        for (size_t dom = 0; dom < N; dom++)
        {
            bool is_zero = false;
            if (c[stud][dom] == 0.0l) {
                is_zero = true;
                zeros.push_back(Zero{.stud = stud, .dom = dom});
            }
            if ((mark_cols[dom]  & Mark::STAR) != 0 ||
                (mark_rows[stud] & Mark::STAR) != 0){
                continue;
            }
            if(is_zero) {
                marks[stud][dom] |= Mark::STAR;
                mark_rows[stud] |= Mark::STAR;
                mark_cols[dom] |= Mark::STAR | Mark::COVER;
            }
        }
    }
    
    size_t covered_lines = 0;
    Zero uncovered_zero;
    std::vector<Zero> zero_seq;
    long double h;
    step1:

    for (const Zero & z : zeros) {
        if( (mark_cols[z.dom]  & Mark::COVER) != 0 ||
            (mark_rows[z.stud] & Mark::COVER) != 0) {
            //zero is covered
            continue;
        }
        //zero is uncovered -> prime it
        marks[z.stud][z.dom] |= Mark::PRIME;

        

        // if((mark_rows[z.stud] & Mark::STAR) == 0) {
        //     uncovered_zero = z;
        //     // there is no starred zero in this row
        //     goto step2; 
        // }
        //there is a starred zero in this row
        bool has_star = false;
        for (size_t col = 0; col < N; col++)
        {
            if ((marks[z.stud][col] & Mark::STAR) != 0)
            {
                // there is a starred zero in this row
                // cover the starred zero's row and uncover column
                mark_rows[z.stud] |= Mark::COVER;
                mark_cols[col] &= ~Mark::COVER;
                has_star = true;
                break;
            }
        }
        if (!has_star) {
            uncovered_zero = z;
            goto step2;
        }
    }
    // we have covered all zeros
    // if we use exactly N lines to cover all these zeros, we have arrived at our result
    // The starred zeros form the perfect matching
    covered_lines = 0;
    for (size_t row = 0; row < N; row++)
    {
        auto t1 = mark_rows[row];
        auto t2 = Mark::COVER;
        auto t3 = t1 & t2;
        auto t4 = (t1 & t2) != 0;
        if ((mark_rows[row] & Mark::COVER) != 0) {
            covered_lines++;
        }
    }
    
    for (size_t col = 0; col < N; col++)
    {
        if((mark_cols[col] & Mark::COVER) != 0) {
            covered_lines ++;
        }
    }
    
    if (covered_lines == N) {
        goto fin;
    }

    goto step3;

step2:

    zero_seq.clear();
    zero_seq.push_back(uncovered_zero);

    while (true)
    {
        bool found_0st = false;
        size_t row_0st;
        // find 0* in current head's column
        // if none -> sequence is complete
        for (size_t row = 0; row < N; row++)
        {
            if ((marks[row][zero_seq[zero_seq.size()-1].dom] & Mark::STAR) != 0)
            {
                found_0st = true;
                row_0st = row;
                break;
            }
        }
        if (!found_0st) break;

        zero_seq.push_back(Zero{.stud = row_0st, .dom = zero_seq[zero_seq.size()-1].dom});

        // find 0' in current head's row
        // it must exist
        
        bool found_0pr = false;
        size_t col_0pr;

        for (size_t col = 0; col < N; col++)
        {
            if ((marks[zero_seq[zero_seq.size()-1].stud][col] & Mark::PRIME) != 0)
            {
                found_0pr = true;
                col_0pr = col;
                break;
            }
        }
        if (!found_0pr) throw std::runtime_error("Algorithm error: impossible zero sequence");

        zero_seq.push_back(Zero{.stud = zero_seq[zero_seq.size()-1].stud, .dom = col_0pr });
    }
    
    for ( Zero & z : zero_seq)
    {
        marks[z.stud][z.dom] &= ~ Mark::STAR;
        if ((marks[z.stud][z.dom] & Mark::PRIME) != 0)
        {
            marks[z.stud][z.dom] |= Mark::STAR;
        }
    }
    
    for (size_t col = 0; col < N; col++)
    {
        mark_cols[col] = Mark::NONE;
    }
    for (size_t row = 0; row < N; row++)
    {
        mark_rows[row] = Mark::NONE;
    }

    covered_lines = 0;
    for ( Zero & z : zeros) {
        // unprime all zeros and cover rows w/ stared 0s
        marks[z.stud][z.dom] &= ~Mark::PRIME;
        if ((marks[z.stud][z.dom] & Mark::STAR) != 0) {
            mark_cols[z.dom] = Mark::COVER;
            covered_lines++;
        }
    }

    if (covered_lines == N) {
        goto fin;
    }

    goto step1;

step3:

    // all zeros covered, size of cover < N ==> there are uncovered fields
    h = std::numeric_limits<long double>::max();
    // find smallest uncovered element
    for (size_t row = 0; row < N; row++)
    {
        if ((mark_rows[row] & Mark::COVER) != 0) {
            continue;
        }
        for (size_t col = 0; col < N; col++)
        {
            if ((mark_cols[col] & Mark::COVER) != 0) {
                continue;
            }
            // uncovered element
            h = std::min(h, c[row][col]);
        }
    }
    for (size_t row = 0; row < N; row++)
    {
        if ((mark_rows[row] & Mark::COVER) == 0) {
            //uncovered -> skip
            continue;
        }
        for (size_t col = 0; col < N; col++)
        {
            c[row][col] += h;
        }
    }
    
    for (size_t col = 0; col < N; col++)
    {
        if ((mark_cols[col] & Mark::COVER) != 0) {
            //uncovered -> skip
            continue;
        }
        for (size_t row = 0; row < N; row++)
        {
            c[row][col] -= h;
        }
    }

    //reindex zeros

    zeros.clear();

    for (size_t row = 0; row < N; row++)
    {
        for (size_t col = 0; col < N; col++)
        {
            if (c[row][col] == 0.0l)
            {
                zeros.push_back(Zero{.stud = row, .dom = col});
            }
        }
    }
    

    goto step1;

fin:

    result Res;
    Res.total_cost = 0;
    Res.matching.clear();

    for ( const Zero &z : zeros)
    {
        if ((marks[z.stud][z.dom] & Mark::STAR) != 0) {
            Res.matching.push_back(std::pair<ll, ll>(z.stud / t.mult, z.dom));
            Res.total_cost += weights[z.stud][z.dom];
        }
    }
    
    return Res;
}