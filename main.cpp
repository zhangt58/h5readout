#include <cstring>

#include <iostream>

#include <CDataSource.h>
#include <CDataSourceFactory.h>
#include <CErrnoException.h>
#include <CPortManagerException.h>

#include "processor.h"
#include "misc.h"
#include "h5readout.h"
#include "modules.h"

int main(int argc, char **argv)
{

  // argument parser
  ArgumentParser argparser = ArgumentParser();
  argparser.parse(argc, argv);

  // --version
  if (argparser.print_version_if_possible())
  {
    return EXIT_SUCCESS;
  }

  // -h or missing -i
  if (argparser.print_help_if_possible())
  {
    return EXIT_FAILURE;
  }

  // invalid input data source?
  if (!argparser.validate_input_data_source())
  {
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
  if (is_file(outfilepath_char))
  {
    printf("Warning: Overwrite output file '%s'? ([Y]/n)", outfilepath_char);
    char c = getc(stdin);
    if (c == 'n' || c == 'N')
    {
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

  // verbosity level
  int verbosity = argparser.get_verbosity();

  // controller type
  std::string ctrl_type = argparser.get_ctrl_type();

  // module list for non-DDAS controller type
  std::vector<std::string> *module_list = argparser.get_module_list();

  // exclude type list
  std::vector<uint16_t> *pexclude = argparser.get_exclude_types();

  // debug
  // argparser.print_all_args();

  // Ring Item types that can be sampled
  std::vector<std::uint16_t> sample;
  // Ring Item types that can be filtered out
  // by passing --exclude argument, which accepts a string of item types, separated by ','
  // std::vector<std::uint16_t> exclude; // = {PHYSICS_EVENT};
  std::vector<std::uint16_t> exclude = {pexclude->begin(), pexclude->end()};

  CDataSource *data_source;
  try
  {
    data_source = CDataSourceFactory::makeSource(input_uri, sample, exclude);
  }
  catch (CPortManagerException &ex)
  {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }
  catch (CErrnoException &ex)
  {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }

  // (DDAS) container for all fragments (exclude trace)
  std::vector<FragmentData> *pfragdata = new std::vector<FragmentData>();

  // (DDAS) container for all trace data
  std::vector<uint16_t> *ptracedata = new std::vector<uint16_t>();

  // container for all scaler data
  std::vector<uint32_t> *pscalerdata = new std::vector<uint32_t>();

  // scaler counter
  std::vector<uint32_t> *pscalerlen = new std::vector<uint32_t>();

  // and timestamp
  std::vector<time_t> *pscalerts = new std::vector<time_t>();

  // container for fragments (VME/V785)
  std::vector<uint64_t> *pfragdata_v785 = new std::vector<uint64_t>();

  CRingItem *pItem;
  uint64_t event_id = 0; // event count
  uint64_t frag_cnt = 0; // total framgnets count

  RunMetaData run_metadata;
  CRingItemProcessor processor;

  try
  {
    // start reading
    fprintf(stdout, "Reading data from: %s\n", input_uri.c_str());
    // reading fragment data
    while (event_id < max_evt_cnt && (pItem = data_source->getItem()))
    {
      processRingItem(processor, pItem, run_metadata, event_id, frag_cnt,
                      pfragdata, ptracedata,
                      pscalerts, pscalerlen, pscalerdata,
                      verbosity,
                      ctrl_type, module_list, pfragdata_v785);
    }
  }
  catch (CErrnoException &ex)
  {
    fprintf(stderr, "Failed to read data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }

  // update metadata
  run_metadata.n_frags = frag_cnt;
  run_metadata.n_events = event_id;
  strcpy(run_metadata.ctrl_type, ctrl_type.c_str());
  std::cout << "Run #: " << run_metadata.number << "\n"
            << " Title: " << run_metadata.title << "\n"
            << " Begin: " << run_metadata.date0 << "\n"
            << " End: " << run_metadata.date1 << "\n"
            << " Duration: " << run_metadata.dt << "\n"
            << " Ring Format: " << run_metadata.fmt << "\n"
            << " Controller: " << run_metadata.ctrl_type << "\n"
            << " Read physics events: " << run_metadata.n_events << "\n"
            << " Read fragments: " << run_metadata.n_frags
            << std::endl;

  // start writing
  fprintf(stdout, "Writing data to: %s\n", outfilepath_char);

  try
  {
    //H5::Exception::dontPrint();
    // create an h5 file handle
    H5::H5File *h5file = new H5::H5File(output_filepath, H5F_ACC_TRUNC);

    // create groups for events and scalers
    H5::Group *evt_grp = new H5::Group(h5file->createGroup(PHYSICS_EVENT_GROUP_NAME));
    H5::Group *scl_grp = new H5::Group(h5file->createGroup(SCALER_GROUP_NAME));

    // meta data
    bool meta_is_written = write_metadata(run_metadata, h5file);
    if (meta_is_written)
    {
      fprintf(stdout, "Writing meta data is done.\n");
    }

    bool frag_is_written;
    if (ctrl_type == "DDAS")
    {
      // fragment data
      frag_is_written = write_fragdata(pfragdata, evt_grp);
    }
    else
    {
      // for (auto it = pfragdata_v785->begin(); it != pfragdata_v785->end(); ++it)
      // {
      //   for (int j = 0; j < 2; ++j)
      //   {
      //     for (int i = 0; i < 100; ++i)
      //     {
      //       std::cout << *it++ << " ";
      //     }
      //     std::cout << std::endl;
      //   }
      //   break;
      // }

      frag_is_written = write_fragdata_vme(pfragdata_v785, evt_grp, frag_cnt,
                                           comp_meth, gfactor);
    }
    if (frag_is_written)
    {
      fprintf(stdout, "Writing fragment data is done.\n");
    }

    // trace data dataset
    // long tracedata_cnt = ptracedata->size();
    // if (tracedata_cnt != 0)
    // { // skip PHYSICS_EVENT
    //   fprintf(stdout, "Trace data size: %i (%g KB)\n", tracedata_cnt, (float)(tracedata_cnt * sizeof(uint16_t) / 1024));
    //   fprintf(stdout, "Total fragments: %i (%g KB)\n", frag_cnt, (float)(frag_cnt * sizeof(FragmentData) / 1024));
    // }

    bool trace_is_written = write_tracedata(ptracedata, evt_grp, frag_cnt, chunk_dims, comp_meth, gfactor);
    if (trace_is_written)
    {
      fprintf(stdout, "Writing trace data is done.\n");
    }

    bool scaler_is_written = write_scalerdata(pscalerdata, scl_grp, pscalerlen, pscalerts, ctrl_type);
    if (scaler_is_written)
    {
      fprintf(stdout, "Writing scaler data is done.\n");
    }

    // clean up
    delete evt_grp;
    delete scl_grp;
    delete h5file;
  }
  catch (H5::FileIException &error)
  {
    error.printErrorStack();
    return EXIT_FAILURE;
  }
  catch (H5::DataSetIException &error)
  {
    error.printErrorStack();
    return EXIT_FAILURE;
  }
  catch (H5::DataSpaceIException &error)
  {
    error.printErrorStack();
    return EXIT_FAILURE;
  }
  catch (H5::DataTypeIException &error)
  {
    error.printErrorStack();
    return EXIT_FAILURE;
  }
  catch (...)
  {
    fprintf(stderr, "Failed to read and export event stream to HDF5.\n");
    return EXIT_FAILURE;
  }

  // fprintf(stdout, "Read and write %i fragments in %i events in total.\n", frag_cnt, event_id);

  return EXIT_SUCCESS;
}
