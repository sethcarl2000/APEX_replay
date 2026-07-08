
#include "../include/argparse.hpp"

#include <THaRun.h>

#include <memory>
#include <iostream>
#include <fstream> 
#include <string>
#include <stdexcept> 

namespace {
    std::fstream outfile;
};

int error_exit(); 

int main(int argc, char* argv[])
{
    argparse::ArgumentParser program("count_coda_events", "1.0"); 

    program.add_argument("infile_path")
        .help("Input CODA file path")
        .required();

    program.add_argument("outfile_path")
        .help("Output file path")
        .required();

    try {
        program.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cerr << err.what() << std::endl;       
        std::cerr << program << std::endl;
        return -1; 
    }

    std::string path_input = program.get<std::string>("infile_path");
    std::string path_output = program.get<std::string>("outfile_path");
    
    outfile.open(path_output.c_str(), std::ios::out | std::ios::trunc);
    if (outfile.is_open() == false) {
        std::cerr << "<count_coda_events>: Failed to open output file '" << path_output << "'\n";
        return -1;
    }
    
    auto run = std::make_unique<THaRun>(path_input.c_str(), "CODA input file");
    run->SetDataRequired(0);
    
    auto st = run->Init();
    if( st != THaRunBase::READ_OK ) {
        std::cerr << "<count_coda_events> Error initializing" << std::endl;
        return error_exit(); 
    }

    st = run->Open();
    if( st != THaRunBase::READ_OK ) {
        std::cerr << "<count_coda_events> Error opening" << std::endl;
        return error_exit(); 
    }

    ULong64_t iev = 0;
    while( (st = run->ReadEvent()) == THaRunBase::READ_OK ) { ++iev; }
    
    std::cout << "<count_coda_events> Number of events read = " << iev << std::endl;
    
    std::cout << "<count_coda_events> Finished file scan, final status = ";
    switch( st ) {
    case THaRunBase::READ_EOF:
        std::cout << "EOF"; break;
    case THaRunBase::READ_ERROR:
        std::cout << "ERROR"; return error_exit(); break;
        
    case THaRunBase::READ_FATAL:
        std::cout << "FATAL"; return error_exit(); break;
    default:
        return error_exit(); 
        std::cout << "UNKNOWN? = " << st; return error_exit(); break;
    }
    std::cout << std::endl;
    if( st != THaRunBase::READ_EOF ) {
       
        std::cerr << "<count_coda_events> Error reading ev " << iev << std::endl;
        return error_exit(); 
    }

    st = run->Close();
    if( st != THaRunBase::READ_OK ) {
        
        std::cerr << "<count_coda_events> Error closing?" << std::endl;
        return error_exit(); 
    }
    std::cout << "<count_coda_events> All successful" << std::endl;
    
    outfile << iev << "\n";
    outfile.close();

    return 0; 
}

int error_exit() {
    outfile << "bad-read\n";
    outfile.close();
    return -1;
}
