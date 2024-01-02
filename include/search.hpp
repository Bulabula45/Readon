#include "hashK.hpp"
#include "constants.hpp"


// std::mutex op_Mutex;

void searchRefs(uint64_t **miQQuery, int *query_mi_cnts, int query_cnt, const char *reFnames[], const int recSize[],
                int p_start, int p_end,
                res_record *res, res_record *res_r);