#include "btllib/status.hpp"
#include "btllib/seq_writer.hpp"

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

#include <sys/types.h>
#include <sys/stat.h>

static const double MAX_READS_PER_CONTIG_10KBP = 120.0;
static const int MINIMIZER_THRESHOLD = 3;

using SeqId = std::string;

struct SeqCoordinates {
  size_t seq_start, seq_len;

  SeqCoordinates(const size_t seq_start, const size_t seq_len)
    : seq_start(seq_start)
    , seq_len(seq_len) {}
};

using ReadsMapping = std::unordered_map<SeqId, std::vector<SeqId>>;
using Index = std::unordered_map<SeqId, SeqCoordinates>;

void load_contigs_reads(ReadsMapping& contigs_reads, const std::string& filepath) {
  std::cerr << "Loading contig reads from " << filepath << "... " << std::flush;
  std::ifstream ifs(filepath);
  std::string token, read_id, contig_id;
  unsigned long i = 0;
  while (ifs >> token) {
    switch (i % 3) {
      case 0: read_id = std::move(token); break;
      case 1: contig_id = std::move(token); break;
      case 2: {
        const auto minimizers = std::stoul(token);
        if (minimizers >= MINIMIZER_THRESHOLD) {
          contigs_reads[contig_id].push_back(read_id);
        }
        break;
      }
      default: std::cerr << "Invalid switch branch.\n"; std::exit(EXIT_FAILURE);
    }
    i++;
  }
  std::cerr << "Done!" << std::endl;
}

void load_index(Index& index, const std::string& filepath) {
  std::cerr << "Loading index from " << filepath << "... " << std::flush;
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
      default: std::cerr << "Invalid switch branch.\n"; std::exit(EXIT_FAILURE);
    }
    i++;
  }
  std::cerr << "Done!" << std::endl;
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

void serve(const char* contigs_filepath,
           const Index& contigs_index,
           const ReadsMapping& contigs_reads,
           const char* reads_filepath,
           const Index& reads_index,
           const char* output_filepath,
           const char* input_pipepath,
           const char* confirm_pipepath) {
  btllib::check_error(mkfifo(input_pipepath, S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");
  btllib::check_error(mkfifo(confirm_pipepath, S_IRUSR | S_IWUSR) != 0, "mkfifo failed.");

  std::ifstream contigs_file(contigs_filepath);
  std::ifstream reads_file(reads_filepath);

  std::string contig_id;

  std::cerr << "Accepting contig IDs at " << input_pipepath << " and outputting sequences at " << output_filepath << std::endl;
  while (contig_id != "x") {
    std::cerr << "Waiting for a new contig ID..." << std::endl;
    std::ifstream inputstream(input_pipepath);
    while (inputstream >> contig_id) {
      if (contig_id == "x") { break; }

      std::cerr << "Received contig ID " << contig_id << ", outputting reads..." << std::endl;

      const auto contig_seq = get_seq_with_index(contig_id, contigs_file, contigs_index);
      const auto& contig_reads_vector = contigs_reads.at(contig_id);
      const auto contig_reads_num = contig_reads_vector.size();
      const auto contig_reads_num_adjusted = std::min(contig_reads_num, decltype(contig_reads_num)(double(std::strlen(contig_seq)) / 10'000.0 * MAX_READS_PER_CONTIG_10KBP));

      std::cerr << "Total number of contig reads = " << contig_reads_num << std::endl;
      std::cerr << "Adjusted number of contig reads = " << contig_reads_num_adjusted << std::endl;

      const auto random_indices = get_random_indices(contig_reads_num, contig_reads_num_adjusted);

      std::cerr << "Writing..." << std::endl;
      btllib::SeqWriter writer(output_filepath);
      unsigned long counter = 0;
      for (const auto read_id_idx : random_indices) {
        const auto read_id = contig_reads_vector[read_id_idx];
        const auto seq = get_seq_with_index(read_id, reads_file, reads_index);

        writer.write(std::to_string(counter), "", seq);
        counter++;
      }
      writer.close();
      std::cerr << "Done writing " << counter << " reads!" << std::endl;

      std::ofstream outputstream(confirm_pipepath);
      outputstream << "1" << std::endl;
      std::cerr << output_filepath << " is ready!" << std::endl;
    }
  }
  std::cerr << "Done!" << std::endl;

  std::remove(input_pipepath);
  std::remove(confirm_pipepath);
}

int main(int argc, char** argv) {
  if (argc != 9) {
    std::cerr << "Wrong args.\n";
    std::exit(EXIT_FAILURE);
  }

  unsigned arg = 1;
  const auto contigs_filepath = argv[arg++];
  const auto contigs_index_filepath = argv[arg++];
  const auto contigs_reads_filepath = argv[arg++];
  const auto reads_filepath = argv[arg++];
  const auto reads_index_filepath = argv[arg++];
  const auto tmpreads_filepath = argv[arg++];
  const auto input_pipepath = argv[arg++];
  const auto confirm_pipepath = argv[arg++];

  ReadsMapping contigs_reads;
  Index reads_index, contigs_index;
  
  load_contigs_reads(contigs_reads, contigs_reads_filepath);
  load_index(contigs_index, contigs_index_filepath);
  load_index(reads_index, reads_index_filepath);

  serve(contigs_filepath,
        contigs_index,
        contigs_reads,
        reads_filepath,
        reads_index,
        tmpreads_filepath,
        input_pipepath,
        confirm_pipepath);

  return 0;
}