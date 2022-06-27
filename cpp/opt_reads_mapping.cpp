#include <fstream>
#include <iostream>
#include <sstream>

#include <map>
#include <string>
#include <vector>
#include <omp.h>

size_t ham_dist(const char* s1, const char* s2, size_t len) {
    size_t dist = 0;
    for (size_t i = 0; i < len; ++i) {
        dist += (s1[i] != s2[i]);
    }

    return dist;
}

std::map<std::string, std::vector<size_t>> build_kmers_map(const std::string& genome, size_t kmer_len) {
    std::map<std::string, std::vector<size_t>> kmers;
    for (size_t i = 0; i < genome.size() - kmer_len + 1; ++i) {
        std::string kmer = genome.substr(i, kmer_len);
        if (kmers.find(kmer) == kmers.end()) {
            kmers.insert({ kmer, {i} });
        }
        else {
            kmers[kmer].push_back(i);
        }
    }

    return kmers;
}

std::vector<int> simple_reads_mapping(const std::string& genome, const std::vector<std::string>& reads, size_t kmer_len) {
    std::vector<int> mapped_reads(reads.size());
    std::map<std::string, std::vector<size_t>> kmers_map = build_kmers_map(genome, kmer_len);
    auto genome_cstr = genome.c_str();

    size_t count = 0;
    for (const auto& read : reads) {
        bool found = false;

        std::string seed = read.substr(0, kmer_len);
        auto seed_it = kmers_map.find(seed);
        std::vector<size_t> poss = {};

        if (seed_it != kmers_map.end()) {
            poss = seed_it->second;
        }

        for (size_t pos : poss) {
            if (pos + read.size() < genome.size()) {
                if (std::strncmp(read.c_str(), genome_cstr + pos, read.size()) == 0) {
                    mapped_reads[count] = int(pos);
                    found = true;
                }
            }
        }

        if (!found) {
            mapped_reads[count] = -1; // no perfect match for read in genome
        }
        ++count;
    }

    return mapped_reads;
}

int main() {
    std::ifstream genome_file;
    std::ifstream reads_file;

    //std::string genome_file_name = "test_example/genome.fasta";
    std::string genome_file_name = "../data/genome.fasta";
    genome_file.open(genome_file_name);
    if (!genome_file) {
        throw std::runtime_error("Unable to open file " + genome_file_name);
    }

    //std::string reads_file_name = "test_example/reads.txt";
    std::string reads_file_name = "../data/reads.txt";
    reads_file.open(reads_file_name);
    if (!reads_file) {
        throw std::runtime_error("Unable to open file " + reads_file_name);
    }

    std::string line;
    std::string genome = "";
    std::getline(genome_file, line); // skip header
    while (std::getline(genome_file, line)) {
        genome += line;
    }
    genome_file.close();

    std::vector<std::string> reads;
    while (std::getline(reads_file, line)) {
        reads.push_back(line);
    }
    reads_file.close();
    auto mapped_reads = simple_reads_mapping(genome, reads, 3);


    std::ofstream res("res.txt");
    for (size_t i = 0; i < mapped_reads.size(); ++i) {
        res << reads[i] << " " << mapped_reads[i] << '\n';
    }

    res.close();

    return 0;
}
