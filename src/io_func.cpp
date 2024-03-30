// io_func.cpp

#include "io_func.hpp"
#include <cstring>
#include <cstdint>
#include <fstream>
#include <algorithm>

size_t count_line(std::istream &is)
{
    
    if( is.bad() ) return 0;  
    
    is.clear();

    is.seekg(0);

    size_t line_cnt;

    size_t lf_cnt = std::count(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>(), '\n');

    line_cnt = lf_cnt;

    if( is.get() == '\n' ) { --line_cnt ; }

    
    is.clear() ; 


    return line_cnt;
}

void toLower(char* str) {
    for (; *str; ++str) {
        *str = std::tolower(static_cast<unsigned char>(*str));
    }
}

int getFileFormat(const char* filename) {
    const char* dot = std::strrchr(filename, '.');
    if (!dot || dot == filename) {
        std::cerr << "Error: Unable to determine file type." << std::endl;
        return -1;
    }

    char* ext = new char[std::strlen(dot + 1) + 1];
    std::strcpy(ext, dot + 1);
    toLower(ext);

    int linesPerRecord = -1;

    if (std::strcmp(ext, "fasta") == 0 || std::strcmp(ext, "fa") == 0) {
        linesPerRecord = 2;
    } else if (std::strcmp(ext, "fastq") == 0 || std::strcmp(ext, "fq") == 0) {
        linesPerRecord = 4;
    } else if (std::strcmp(ext, "gz") == 0) {
        const char* prevDot = dot;
        while (prevDot != filename && *prevDot != '.') {
            --prevDot;
        }
        if (prevDot == filename) {
            std::cerr << "Error: Unable to determine reads file type." << std::endl;
            delete[] ext;
            return -1;
        }

        char* prevExt = new char[std::strlen(prevDot + 1) + 1];
        std::strcpy(prevExt, prevDot + 1);
        toLower(prevExt);

        if (std::strcmp(prevExt, "fastq") == 0 || std::strcmp(prevExt, "fq") == 0) {
            linesPerRecord = 4;
        } else if (std::strcmp(prevExt, "fasta") == 0 || std::strcmp(prevExt, "fa") == 0) {
            linesPerRecord = 2;
        }

        delete[] prevExt;
    }

    delete[] ext;
    return linesPerRecord;
}