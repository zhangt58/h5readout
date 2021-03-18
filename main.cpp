#include <cstring>

#include <iostream>

#include <CDataSource.h>
#include <CDataSourceFactory.h>
#include <CErrnoException.h>
#include <CPortManagerException.h>

#include "processor.h"
#include "misc.h"
#include "h5readout.h"

int main(int argc, char **argv) {

  // argument parser
  auto argparser = ArgumentParser();
  argparser.parse(argc, argv);

  // --version
  if (argparser.print_version_if_possible()) {
    return EXIT_SUCCESS;
  }

  // -h or missing -i
  if (argparser.print_help_if_possible()) {
    return EXIT_FAILURE;
  }

  // invalid input data source?
  if (!argparser.validate_input_data_source()) {
    fprintf(stderr, "Failed to find the real path for '%s': %s\n",
            argparser.get_input_data_source().c_str(), strerror(errno));
    return EXIT_FAILURE;
  }

  // get input data source as URI
  std::string input_uri = argparser.get_input_data_source_as_uri();

  // get output full path
  std::string output_filepath = argparser.get_output_filepath();

  // give warning if to overwrite existing file.
  char outfilepath_char[output_filepath.length() + 1];
  strcpy(outfilepath_char, output_filepath.c_str());
  if (is_file(outfilepath_char)) {
    printf("Warning: Overwrite output file '%s'? ([Y]/n)", outfilepath_char);
    char c = getc(stdin);
    if (c == 'n' || c == 'N') {
      return EXIT_SUCCESS;
    }
  }

  // data gzip compress level and compress method
  int gfactor = argparser.get_gzip_compress_level();
  std::string comp_meth = argparser.get_compress_method();

  // maximum number of events to process
  uint64_t max_evt_cnt = argparser.get_max_evt();

  // chunk dims
  hsize_t *chunk_dims = argparser.get_chunk_dims();

  // verbose?
  bool verbose = argparser.is_verbose();

  // debug
  argparser.print_all_args();


  // Ring Item types that can be sampled
  std::vector<std::uint16_t> sample;
  // Ring Item types that can be filtered out
  std::vector<std::uint16_t> exclude; // = {PHYSICS_EVENT};

  CDataSource *data_source;
  try {
    data_source = CDataSourceFactory::makeSource(input_uri, sample, exclude);
  } catch (CPortManagerException &ex) {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  } catch (CErrnoException &ex) {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }

  // container for all fragments (exclude trace)
  std::vector<FragmentData> *pfragdata = new std::vector<FragmentData>();
  // container for all trace data
  std::vector<uint16_t> *ptracedata = new std::vector<uint16_t>();

  CRingItem *pItem;
  uint64_t event_id = 0;  // event count
  uint64_t frag_cnt = 0;  // total framgnets count

  RunMetaData run_metadata;
  CRingItemProcessor processor;

  try {
    // start reading
    fprintf(stdout, "Reading data from: %s\n", input_uri);
    // reading fragment data
    while (event_id < max_evt_cnt && (pItem = data_source->getItem())) {
        processRingItem(processor, pItem, run_metadata, event_id, frag_cnt, pfragdata, ptracedata, verbose);
    }
  } catch (CErrnoException &ex) {
    fprintf(stderr, "Failed to read data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }

  // update metadata
  run_metadata.n_frags = frag_cnt;
  run_metadata.n_events = event_id;
  std::cout << "Run #: " << run_metadata.number << "\n"
            << " Title: " << run_metadata.title << "\n"
            << " Timestamp: " << run_metadata.ts << "\n"
            << " Datetime: " << run_metadata.date << "\n"
            << " RingFormat: " << run_metadata.fmt << "\n"
            << " Read events: " << run_metadata.n_events << "\n"
            << " Read fragments: " << run_metadata.n_frags
            << std::endl;


  return EXIT_SUCCESS;
}
