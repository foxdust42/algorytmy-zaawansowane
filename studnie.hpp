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
void write_tasks(std::string filename, std::vector<task> &tasks);

std::vector<task> generate_tasks_str(std::string config);

std::vector<task> generate_tasks_int(void);

std::vector<task> _generate_tasks(unsigned int count, unsigned int coord_lim, 
    unsigned int   stud_count,      unsigned int   mult,
    unsigned int * stud_count_high, unsigned int * mult_high);

void run_benchmark(void);

result mockProgram(const task &problem);

result studnie2(const task &problem);



#endif // STUDNIE_H