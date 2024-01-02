#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <cstring>
#include <cstdint>
#include <fstream>

#include "hashK.hpp"
#include "constants.hpp"
#include "search.hpp"

off_t mapping(const int reads_cnt, int linesPerRecord, off_t reads_len, off_t offset, int reads_fd, std::ofstream &outputTsv, std::ofstream &outputReads);
