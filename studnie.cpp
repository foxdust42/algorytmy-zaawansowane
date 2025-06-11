#include "studnie.hpp"
#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include <iomanip>

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

    out.precision(10);

    out << res.size() << std::endl << std::endl;

    for (ll i = 0 ; i < res.size(); i++){

        out << std::fixed << res[i].total_cost << std::endl;
        
        for (const auto &m : res[i].matching) {
            out << "(" << m.first << ", " << m.second << ")" << std::endl;
        }
        out << std::endl;
    }
    
    out.close();
}

void write_tasks(std::string filename,std::vector<task> &tasks) {
    std::ofstream out(filename, std::ios::out | std::ios::trunc);

    out << tasks.size() << std::endl << std::endl;

    for (ll i=0; i< tasks.size(); ++i) {
        out << tasks[i].studnie.size() << " " << tasks[i].mult << std::endl << std::endl;

        for(const auto &s : tasks[i].studnie) {
            out << "(" << s.x << ", " << s.y << ")" << std::endl;
        }
        out << std::endl;

        for (const auto &d : tasks[i].domy) {
            out << "(" << d.x << ", " << d.y << ")" << std::endl;
        }
        out << std::endl;
    }
    out.close();
}

template <class IntType = int>
class _generator {
    public:
    virtual IntType get_next(void) = 0; 
};

template <class IntType = int>
class _static_generator : public _generator<IntType> {
    private:
    const IntType value;
    public: 
    _static_generator(IntType val) : value(val) {
    }
    IntType get_next(void) override {
        return value;
    }
};

template <class IntType = int>
class _uniform_generator : public _generator<IntType> {
    private:
    std::uniform_int_distribution<IntType> rand_dist;
    std::mt19937 * r_source;
    public: 
    _uniform_generator(std::mt19937 * rand_source, IntType val_low, IntType val_high) {
        r_source = rand_source;
        rand_dist = std::uniform_int_distribution<IntType>(val_low, val_high);
    }
    IntType get_next(void) override {
        return rand_dist(*r_source);
    }
};

std::vector<task> _generate_tasks(
    unsigned int  count,            unsigned int  coord_lim,
    unsigned int  stud_count,       unsigned int  mult,
    unsigned int *stud_count_high,  unsigned int *mult_high)
{
    std::vector<task> tasks;

    std::random_device rd;
    std::mt19937 mt(rd());

    std::uniform_int_distribution<ll> r_coords(- ll(coord_lim), ll(coord_lim));

    _generator<unsigned int> * stud_count_gen;
    _generator<unsigned int> * mult_gen;

    if(stud_count_high != nullptr && *stud_count_high != stud_count) {
        uint cl = stud_count;
        uint ch = *stud_count_high;
        if (cl > ch ) std::swap(cl ,ch);
        stud_count_gen = new _uniform_generator<unsigned int>(&mt, cl, ch);    
    }
    else {
        stud_count_gen = new _static_generator<unsigned int>(stud_count);
    }

    if(mult_high != nullptr && *mult_high != mult) {
        uint ml = mult;
        uint mh = *mult_high;
        if (ml > mh ) std::swap(ml ,mh);
        mult_gen = new _uniform_generator<unsigned int>(&mt, ml, mh);    
    }
    else {
        mult_gen = new _static_generator<unsigned int>(mult);
    }

    for (ll i = 0; i < count; i++) {

        ll studnie = stud_count_gen->get_next();
        ll mult = mult_gen->get_next();

        task t;

        t.mult = mult;

        for (ll s = 0; s < studnie; s++)
        {
            t.studnie.push_back(point{.idx = s, .x = r_coords(mt), .y = r_coords(mt)});
        }
        
        for (ll d = 0; d < studnie * t.mult; d++)
        {
            t.domy.push_back(point{.idx = d, .x = r_coords(mt), .y = r_coords(mt)});
        }

        tasks.push_back(t);
    }

    delete stud_count_gen;
    delete mult_gen;

    return tasks;
}

std::vector<task> generate_tasks_int(void) {

    ll to_create;
    uint coord_lim;
    uint well_count;
    uint well_count_upper;
    uint mult;
    uint mult_upper;

    std::string ok;

    while (true)
    {
        std::cout << "Specify number of tasks to create (<=0 to abort): ";
        std::cin >> to_create;

        if (std::cin.fail()) {
            std::cout << "Input is not a whole number\n";
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<unsigned int>::max(), '\n');
            continue;
        }

        if (to_create <= 0) {
            std::cout << "aborting generation...\n";
            return {};
        }

        std::cout << "specify coordinate limit (-lim <= x, y <= lim): ";
        std::cin >> coord_lim;

        std::cout << "specify well count lower bound (inclusive): ";
        std::cin >> well_count;

        std::cout << "specify well count upper bound (inclusive): ";
        std::cin >> well_count_upper;

        std::cout << "specify mult lower bound (inclusive)  (# houses = {well_count} * {mult}): ";
        std::cin >> mult;

        std::cout << "specify mult upper bound (inclusive): ";
        std::cin >> mult_upper;

        std::cout << "\nWill generate " << to_create << " tasks\n"
            << "With coordinates in [-" << coord_lim << ", " << coord_lim << "]^2 \n"
            << "With number of wells in [" << well_count <<  ", " << well_count_upper << "] \n"
            << "And number of houses per well in [" << mult << ", " << mult_upper << "]\n"
            << "Proceed with generation (y/n)?: ";

        ok.clear();
        std::cin >> ok;
        if (ok == "Y" || ok == "y") { break; }

        std::cout << "Starting over...\n\n";
    }
    
    return _generate_tasks(to_create, coord_lim,
            well_count, mult,
            &well_count_upper, &mult_upper);
}

std::vector<task> generate_tasks_str(std::string s) {
    return {};
}

void run_benchmark(void) {

    uint no_tasks_per = 1000;
    
    uint low_wells = 2;
    uint high_wells = 30;
    uint mult = 5;

    
    std::vector<std::pair<uint, std::chrono::milliseconds>> timer_results;
    for (uint i = low_wells; i <= high_wells; i++)
    {
        const auto tasks = _generate_tasks(no_tasks_per, 100, 
            i, mult,
            &i, &mult
        );

        std::vector<result> res;
        auto t0 = std::chrono::high_resolution_clock::now(); 
        for(const task &t : tasks) {
            res.push_back(studnie2(t));
        }
        auto t1 = std::chrono::high_resolution_clock::now();

        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);

        timer_results.push_back(std::make_pair(i, ms));      
        std::cout << "Done with " << std::to_string(i) << "\n"; 
    }
    
    std::string out;

    out = "\"No. wells\", \"total time (ms)\"\n";

    for (auto &&res : timer_results)
    {
        out += std::to_string(res.first) + ", " + std::to_string(res.second.count()) + "\n";
    }
    
    std::cout << out;
            
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