#include "studnie.h"
#include <iostream>
#include <fstream>

long double distance(const point& p1, const point &p2) {
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