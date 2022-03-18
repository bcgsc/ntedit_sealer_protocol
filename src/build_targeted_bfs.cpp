#include "btllib/status.hpp"
#include "btllib/seq_writer.hpp"
#include "btllib/counting_bloom_filter.hpp"

#include <algorithm>
#include <random>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <memory>
#include <thread>
#include <chrono>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

static const double MAX_READS_PER_CONTIG_10KBP = 120.0;
static const std::string INPIPE = "input";
static const std::string CONFIRMPIPE = "confirm";
static const std::string SEPARATOR = "_";
static const std::string BF_EXTENSION = ".bf";
static const std::string END_SYMBOL = "x";

namespace opt {
  std::string contigs_filepath;
  std::string contigs_index_filepath;
  std::string contigs_reads_filepath;
  std::string reads_filepath;
  std::string reads_index_filepath;
  size_t cbf_bytes = 4ULL * 1024ULL * 1024ULL;
  size_t bf_bytes = 1ULL * 1024ULL * 1024ULL;
  unsigned kmer_threshold = 5;
  unsigned mx_threshold = 5;
  std::string prefix = "targeted";
  std::vector<unsigned> ks = { 32, 28, 24, 20 };
  bool ks_set = false;
  unsigned hash_num = 4;
  unsigned threads = 1;
}

using SeqId = std::string;

struct SeqCoordinates {
  size_t seq_start, seq_len;

  SeqCoordinates(const size_t seq_start, const size_t seq_len)
    : seq_start(seq_start)
    , seq_len(seq_len) {}
};

using ReadsMapping = std::unordered_map<SeqId, std::vector<SeqId>>;
using Index = std::unordered_map<SeqId, SeqCoordinates>;

void load_contigs_reads(ReadsMapping& contigs_reads, const std::string& filepath, const unsigned mx_threshold) {
  btllib::log_info(std::string("Loading contig reads from ") + filepath + "... ");
  std::ifstream ifs(filepath);
  std::string token, read_id, contig_id;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 3) {
      case 0: read_id = std::move(token); break;
      case 1: contig_id = std::move(token); break;
      case 2: {
        const auto minimizers = std::stoul(token);
        if (minimizers >= mx_threshold) {
          contigs_reads[contig_id].push_back(read_id);
        }
        break;
      }
      default: {
        btllib::log_error("Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("Done!");
}

void load_index(Index& index, const std::string& filepath) {
  btllib::log_info(std::string("Loading index from ") + filepath + "... ");
  std::ifstream ifs(filepath);
  std::string token, id;
  unsigned long id_start, id_end, seq_end;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 4) {
      case 0: id = std::move(token); break;
      case 1: id_start = std::stoul(token); break;
      case 2: id_end = std::stoul(token); break;
      case 3: {
        seq_end = std::stoul(token);
        index.emplace(std::piecewise_construct, std::make_tuple(id), std::make_tuple(id_end + 1, seq_end - id_end - 1));
        break;
      }
      default: {
        btllib::log_error("Invalid switch branch.");
        std::exit(EXIT_FAILURE);
      }
    }
    i++;
  }
  btllib::log_info("Done!");
}

std::vector<size_t> get_random_indices(const size_t total_size, const size_t count) {
  btllib::check_error(count > total_size, "get_random_indices: count cannot be larger than total_size.");

  static std::random_device dev;
  static std::mt19937 rng(dev());

  std::vector<size_t> all_indices(total_size);
  std::iota(all_indices.begin(), all_indices.end(), 0);
  std::shuffle(all_indices.begin(), all_indices.end(), rng);
  decltype(all_indices) indices(all_indices.begin(), all_indices.begin() + count);

  return indices;
}

char* get_seq_with_index(const std::string& id, std::ifstream& seqfile, const Index& index) {
  static const size_t max_seqlen = 300'000UL;
  static char* seq = new char[max_seqlen];
  
  const auto& coords = index.at(id);
  btllib::check_error(coords.seq_len >= max_seqlen, "Read size over max.");
  seqfile.seekg(coords.seq_start);
  seqfile.read(seq, coords.seq_len);
  seq[coords.seq_len] = '\0';
  return seq;
}

void serve(const std::string& contigs_filepath,
           const Index& contigs_index,
           const ReadsMapping& contigs_reads,
           const std::string& reads_filepath,
           const Index& reads_index,
           const std::string& bf_prefix,
           const size_t cbf_bytes,
           const size_t bf_bytes,
           const unsigned kmer_threshold,
           const std::string& input_pipepath,
           const std::string& confirm_pipepath,
           const unsigned hash_num,
           const std::vector<unsigned>& ks) {
  btllib::check_error(kmer_threshold == 0, "k-mer threshold must be >0.");
  btllib::check_error(mkfifo(input_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");
  btllib::check_error(mkfifo(confirm_pipepath.c_str(), S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");

  std::vector<std::string> bf_paths_base;
  for (const auto k : ks) {
    bf_paths_base.push_back(std::string(bf_prefix) + "k" + std::to_string(k) + BF_EXTENSION);
  }
  std::vector<std::string> bf_paths;

  std::ifstream contigs_file(contigs_filepath);
  std::ifstream reads_file(reads_filepath);

  std::string contig_set_prefix;
  std::string contig_id;

  btllib::log_info(std::string("Accepting contig IDs at ") + input_pipepath);
  while (true) {
    std::vector<std::unique_ptr<btllib::KmerCountingBloomFilter8>> cbfs;
    std::vector<std::unique_ptr<btllib::KmerBloomFilter>> bfs;
    for (const auto k : ks) {
      cbfs.push_back(std::make_unique<btllib::KmerCountingBloomFilter8>(cbf_bytes, hash_num, k));
      bfs.push_back(std::make_unique<btllib::KmerBloomFilter>(bf_bytes, hash_num, k));
    }
    std::ifstream inputstream(input_pipepath);
    if (!(inputstream >> contig_set_prefix)) { break; }
    if (contig_set_prefix == END_SYMBOL) { break; }
    bf_paths.clear();
    for (const auto& bf_path_suffix : bf_paths_base) {
      bf_paths.push_back(contig_set_prefix + bf_path_suffix);
    }
    while (inputstream >> contig_id) {
      if (contig_id == END_SYMBOL) { break; }

      const auto contig_seq = get_seq_with_index(contig_id, contigs_file, contigs_index);
      const auto& contig_reads_vector = contigs_reads.at(contig_id);
      const auto contig_reads_num = contig_reads_vector.size();
      const auto contig_reads_num_adjusted = std::min(contig_reads_num, decltype(contig_reads_num)(double(std::strlen(contig_seq)) / 10'000.0 * MAX_READS_PER_CONTIG_10KBP));

      const auto random_indices = get_random_indices(contig_reads_num, contig_reads_num_adjusted);

#pragma omp parallel
#pragma omp single
      for (const auto read_id_idx : random_indices) {
        const auto read_id = contig_reads_vector[read_id_idx];
        const auto seq = get_seq_with_index(read_id, reads_file, reads_index);
        const std::string seqcopy(seq);
#pragma omp task firstprivate(seqcopy)
        for (size_t i = 0; i < ks.size(); i++)
        {
          const auto& cbf = cbfs[i];
          const auto& bf = bfs[i];
          btllib::NtHash nthash(seqcopy, hash_num, ks[i]);
          while (nthash.roll()) {
            if (cbf->insert_thresh_contains(nthash.hashes(), kmer_threshold) >= kmer_threshold) {
              bf->insert(nthash.hashes());
            }
          }
        }
      }
    }
    if (contig_id == END_SYMBOL) { break; }
    for (size_t i = 0; i < bfs.size(); i++) {
      bfs[i]->save(bf_paths[i]);
    }
    std::ofstream outputstream(confirm_pipepath);
    outputstream << "1" << std::endl;
  }
  btllib::log_info("Targeted BF builder done!");

  std::remove(input_pipepath.c_str());
  std::remove(confirm_pipepath.c_str());
}

void start_watchdog() {
  (new std::thread([] () {
    if (getppid() == 1) {
      std::exit(-1);
    }
    std::this_thread::sleep_for(std::chrono::seconds(1));
  }))->detach();
}

int main(int argc, char** argv) {
  btllib::check_error(argc != 8, "Wrong args.");

  start_watchdog();

  unsigned arg = 1;
  opt::contigs_filepath = argv[arg++];
  opt::contigs_index_filepath = argv[arg++];
  opt::contigs_reads_filepath = argv[arg++];
  opt::reads_filepath = argv[arg++];
  opt::reads_index_filepath = argv[arg++];
  opt::kmer_threshold = std::stoi(argv[arg++]);
  opt::mx_threshold = std::stoi(argv[arg++]);

  omp_set_num_threads(opt::threads);

  ReadsMapping contigs_reads;
  Index reads_index, contigs_index;
  
  load_contigs_reads(contigs_reads, opt::contigs_reads_filepath, opt::mx_threshold);
  load_index(contigs_index, opt::contigs_index_filepath);
  load_index(reads_index, opt::reads_index_filepath);

  serve(opt::contigs_filepath,
        contigs_index,
        contigs_reads,
        opt::reads_filepath,
        reads_index,
        opt::prefix + SEPARATOR,
        opt::cbf_bytes,
        opt::bf_bytes,
        opt::kmer_threshold,
        opt::prefix + SEPARATOR + INPIPE,
        opt::prefix + SEPARATOR + CONFIRMPIPE,
        opt::hash_num,
        opt::ks);

  return 0;
}