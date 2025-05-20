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


result studnie(const task &problem);

result studnie2(const task &problem);




#endif // STUDNIE_H