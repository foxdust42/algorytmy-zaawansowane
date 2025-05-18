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
    }
};

static const char * const format_help = "format help\n";



int main(int argc, char *argv[])
{
    const auto prog_name = argv[0];
    const char * value = NULL;
    cag_option_context context;
    cag_option_init(&context, options, CAG_ARRAY_SIZE(options), argc, argv);

    std::string infile;
    std::string outfile;

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
        std::strftime(time_string, 50, "%F-%T", std::localtime(&now));
        outfile = "Studnie_" + std::string(time_string) + ".txt";
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
