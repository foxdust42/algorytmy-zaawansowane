#include <vector>
#include <string>
#include <ctime>
#include <iostream>

#include "external_libraries/include/cargs.h"

#include "studnie.hpp"

static struct cag_option options[] = {
    {
        .identifier = 'i',
        .access_letters = "i",
        .access_name = "input",
        .value_name = "INPUT_FILE",
        .description = "Input file"
    },
    {
        .identifier = 'o',
        .access_letters = "o",
        .access_name = "output",
        .value_name = "OUTPUT_FILE",
        .description = "Output file"
    },
    {
        .identifier = 'f',
        .access_letters = nullptr,
        .access_name = "help-format",
        .description = "Print file format information and exit"
    },
    {
        .identifier = 'h',
        .access_letters = "h",
        .access_name = "help",
        .description = "Show help text"
    },
    {
        .identifier = 'g',
        .access_letters = "g",
        .access_name = "generate",
        .description = "Run problem generator\n",
    },
    {
        .identifier = 'b',
        .access_letters = "b",
        .access_name = "benchmark",
        .description = "run program benchamrking"
    }
};

static const char * const format_help = 
    "----- Input File Format -----\n"
    "{Number of tasks}\n"
    "\n"
    "{Number of wells} {multiplicity}   # Task 1\n"
    "\n"
    "({well_1_x}, {well_2_y})           # repeated once for each well in task\n"
    "...\n"
    "({well_n_x}, {well_n_y})\n"
    "\n"
    "({house_1_x}, {house_1_y})         # repeated (wells * multiplicity) times\n"
    "...\n" 
    "({house_nk_x}, {house_nk_y}})"
    "\n"
    "{Number of wells} {multiplicity}   # Task 2 - repeats pattern\n"
    "...\n"
    "\n"
    "--- Example ---\n"
    "2\n"
    "\n"
    "1 3\n"
    "\n"
    "(1, 2)\n"
    "\n"
    "(3, 4)\n"
    "(5, 6)\n"
    "(7, 8)\n"
    "\n"
    "2 2\n"
    "\n"
    "(-1, -1)\n"
    "(1, 1)\n"
    "\n"
    "(-2, -2)\n"
    "(0, 0)\n"
    "(1, 2)\n"
    "(2, 1)\n"
    "\n"
    "----- Output File Format -----\n"
    "{number of tasks}                  # Tasks will print in the same order as input\n"
    "\n"
    "{minimal weight}                   # Task 1\n"
    "({well_index_1}, {house_index_1})\n"
    "({well_index_1}, {house_index_2})\n"
    "...\n"
    "({well_index_n}, {house_index_kn})\n"
    "\n"
    "{minimal weight}                   # Task 2 - repeats pattern\n"
    "\n"
    "--- Example ---\n"
    "2\n"
    "\n"
    "16.9706\n"
    "(0, 0)\n"
    "(0, 1)\n"
    "(0, 2)\n"
    "\n"
    "4.82843\n"
    "(0, 0)\n"
    "(0, 1)\n"
    "(1, 2)\n"
    "(1, 3)\n"
    ;

int main(int argc, char *argv[])
{
    const auto prog_name = argv[0];
    const char * value = NULL;
    cag_option_context context;
    cag_option_init(&context, options, CAG_ARRAY_SIZE(options), argc, argv);
    const char * tmp;

    std::string infile;
    std::string outfile;

    bool run_gen = false;
    bool run_bench = false;

    bool infile_set = false;
    bool outfile_set = false;

    while (cag_option_fetch(&context))
    {
        switch (cag_option_get_identifier(&context))
        {
        case 'i':
            if (infile_set){
                printf("Only one infile is allowed, will write to: %s", infile.c_str());
                break;
            }
            infile = std::string(cag_option_get_value(&context));
            infile_set = true;
            break;
        case 'o':
            if (outfile_set){
                printf("Only one outfile is allowed, will write to: %s", outfile.c_str());
                break;
            }
            outfile = std::string(cag_option_get_value(&context));
            outfile_set = true;
            break;
        case 'f':
            printf("%s", format_help);
            return 0;
        
        case 'g':
            run_gen = true;
            break;

        case 'b':
            run_bench = true;
            break;

        case 'h':
            printf("Usage: %s [OPTION]\n", prog_name);
            cag_option_print(options, CAG_ARRAY_SIZE(options), stdout);
            return 0;

        case '?':
        default:
            break;
        }
    }
    
    if (! outfile_set){
        char time_string[50];
        std::time_t now = std::time(nullptr);
        std::strftime(time_string, 50, "%F_%H-%M-%S", std::localtime(&now));
        outfile = "Studnie_" + std::string(time_string) + ".txt";
    }
    
    if (run_gen) {

        std::vector<task> t;

        std::cout << "Runnig generator in interactive mode" << std::endl;
        t = generate_tasks_int();

        if(t.empty()) {
            std::cout << "Nothing to output\n";
            return 0;
        }

        if (!outfile_set) {
            outfile = "Problem_" + outfile;
        }

        write_tasks(outfile, t);

        std::cout << "Wrote " << t.size() << " problems to " << outfile << "\n";

        return 0;
    }

    if(run_bench) {
        run_benchmark();
        return 0;
    }

    std::cout << "outfile:" << std::endl; 
    std::cout << outfile << std::endl;
    std::cout << "infile:" << std::endl;
    std::cout << infile << std::endl;
    std::vector<task> tasks;
    
    try
    {
        tasks = load_file(infile);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
    
    std::vector<result> results;

    for (auto &&task : tasks)
    {
        results.push_back(studnie2(task));
    }

    write_results(outfile, results);

    return 0;
}
