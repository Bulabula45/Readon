
#include "io_func.hpp"
#include "readon.hpp"

using std::cout;
using std::endl;


int main(int argc, char *argv[])
{

    const char *reads_file = nullptr;
    const char *output_file = nullptr;
    const char *mode = "all";
    const char *helpMessage = R"(
Usage: <program_name> -i <reads_file> -o <output_file> [-m <mode>]
-m <mode> : Specify the mode ('error' to select error mode with sequencing data errors; otherwise,
a larger k-mer size will be used by default)
)";

    int opt;
    while ((opt = getopt(argc, argv, "i:o:m:")) != -1)
    {
        switch (opt)
        {
        case 'i':
            reads_file = optarg;
            break;
        case 'o':
            output_file = optarg;
            break;
        case 'm':
            if (optarg[0] != '\0')
            {
                mode = optarg;
            }
            break;
        default:
            std::cerr << helpMessage << endl;
            return 1;
        }
    }

    if (reads_file == nullptr || output_file == nullptr)
    {
        std::cerr << helpMessage << endl;
        return 1;
    }

    if (strcmp(mode, "e") == 0)
    {
        kmer_l = 15;
        interval = 6;
        kmer_2l = kmer_l * 2;
        window_s = kmer_l - interval;
    }
    else
    {
        kmer_l = 25;
        interval = 10;
        kmer_2l = kmer_l * 2;
        window_s = kmer_l - interval;
    }

    initializeConstants();

    std::string output_tsv_str = std::string(output_file) + ".tsv";
    const char* output_tsv = output_tsv_str.c_str();
    std::string output_fa_str = std::string(output_file) + ".reads.fa";
    const char* output_fa = output_fa_str.c_str();

    cout << "kmer length: " << kmer_l << endl;
    cout << "window size: " << window_s << endl;

    int reads_fd = open(reads_file, O_RDONLY);
    off_t reads_len = lseek(reads_fd, 0, SEEK_END);

    int linesPerRecord = getFileFormat(reads_file);

    std::ifstream inFile(reads_file);
    int lines_cnt = count_line(inFile);

    if (reads_len < file_size_MAX)
    {
        
        int reads_cnt =  int((lines_cnt + linesPerRecord - 1) / linesPerRecord);
        std::ofstream outputTsv(output_tsv);
        std::ofstream outputReads(output_fa);
        mapping(reads_cnt, linesPerRecord, reads_len, 0, reads_fd, outputTsv, outputReads);


    }
    else {

        int file_splits = (reads_len + file_size_MAX - 1) / file_size_MAX;
        int reads_cnt =  (lines_cnt + linesPerRecord - 1) / linesPerRecord;
        const int reads_per = (reads_cnt + file_splits - 1) / file_splits;
        const int last_reads_per = reads_cnt - (file_splits - 1) * reads_per;
        std::ofstream outputTsv(output_tsv);
        std::ofstream outputReads(output_fa);
        off_t tmp_reads_len = file_size_MAX*1.5;
        off_t offset = mapping(reads_per, linesPerRecord, tmp_reads_len, 0, reads_fd, outputTsv, outputReads);
        for (int i = 1; i < file_splits-1; ++i) {
            offset = mapping(reads_per, linesPerRecord, tmp_reads_len, offset, reads_fd, outputTsv, outputReads);
        }
        mapping(last_reads_per, linesPerRecord, reads_len - offset, offset, reads_fd, outputTsv, outputReads);
        close(reads_fd);

    
    
    }


    

    return 0;
}
