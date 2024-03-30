#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <vector>
#include <string>
#include <mutex>

extern const int file_end;
extern const char *reFnames[];
extern const int recSize[];
extern std::vector<std::vector<std::string>> clstrNames;
extern int clstrSizes[];
extern char **formatted_Ropt;
extern int kmer_l;
extern int kmer_2l;
extern int interval;
extern int window_s;

void initializeConstants(); 

#endif  // CONSTANTS_H
