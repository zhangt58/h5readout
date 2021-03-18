#include <cstring>

#include <iostream>

#include <CDataSource.h>
#include <CDataSourceFactory.h>
#include <CErrnoException.h>
#include <CPortManagerException.h>

#include "processor.h"
#include "misc.h"
#include "h5readout.h"

// rows of trace data written at once.
const int T_NROWS_PER_WRITE = 1;

// rows of scaler data written at once.
const int S_NROWS_PER_WRITE = 1;

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
  auto chunk_dims = argparser.get_chunk_dims();

  // verbose?
  bool verbose = argparser.is_verbose();

  // debug
  // argparser.print_all_args();

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
  // container for all scaler data
  std::vector<uint32_t> *pscalerdata = new std::vector<uint32_t>();
  // scaler counter
  std::vector<uint32_t> *pscalerlen = new std::vector<uint32_t>();
  // and timestamp
  std::vector<time_t> *pscalerts = new std::vector<time_t>();

  CRingItem *pItem;
  uint64_t event_id = 0;  // event count
  uint64_t frag_cnt = 0;  // total framgnets count

  RunMetaData run_metadata;
  CRingItemProcessor processor;

  try {
    // start reading
    fprintf(stdout, "Reading data from: %s\n", input_uri.c_str());
    // reading fragment data
    while (event_id < max_evt_cnt && (pItem = data_source->getItem())) {
        processRingItem(processor, pItem, run_metadata, event_id, frag_cnt,
                        pfragdata, ptracedata,
                        pscalerts, pscalerlen, pscalerdata,
                        verbose);
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
            << " Begin: " << run_metadata.date0 << "\n"
            << " End: " << run_metadata.date1 << "\n"
            << " Duration: " << run_metadata.dt << "\n"
            << " Ring Format: " << run_metadata.fmt << "\n"
            << " Read physics events: " << run_metadata.n_events << "\n"
            << " Read fragments: " << run_metadata.n_frags
            << std::endl;

  // start writing
  fprintf(stdout, "Writing data to: %s\n", outfilepath_char);

  try {
    // H5::Exception::dontPrint();
    // create an h5 file handle
    auto h5file = new H5::H5File(output_filepath, H5F_ACC_TRUNC);

    // create mem data type for FragmentData
    const H5::CompType frag_dtype(sizeof(FragmentData));
    frag_dtype.insertMember(FRAGMENT_DATA_FRAGMENT_ID,         HOFFSET(FragmentData, fragment_id),        H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_EVENT_ID,            HOFFSET(FragmentData, event_id),           H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_TIMESTAMP,           HOFFSET(FragmentData, timestamp),          H5::PredType::NATIVE_DOUBLE);
    frag_dtype.insertMember(FRAGMENT_DATA_COARSE_TIME,         HOFFSET(FragmentData, coarse_time),        H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_ENERGY,              HOFFSET(FragmentData, energy),             H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TRACE_LENGTH,        HOFFSET(FragmentData, trace_length),       H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CRATE_ID,            HOFFSET(FragmentData, crate_id),           H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_SLOT_ID,             HOFFSET(FragmentData, slot_id),            H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_ID,          HOFFSET(FragmentData, channel_id),         H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_FINISH_CODE,         HOFFSET(FragmentData, finish_code),        H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_OVER_UNDER_FLOW, HOFFSET(FragmentData, adc_over_underflow), H5::PredType::NATIVE_SHORT);
    frag_dtype.insertMember(FRAGMENT_DATA_CFD_FAIL_BIT,        HOFFSET(FragmentData, cfd_fail_bit),       H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_OVERFLOW_CODE,       HOFFSET(FragmentData, overflow_code),      H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_MODMSPS,             HOFFSET(FragmentData, adc_frequency),      H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_RESOLUTION,      HOFFSET(FragmentData, adc_resolution),     H5::PredType::NATIVE_INT);

    // create a new group "PhysicsEvent"
    auto grp = new H5::Group(h5file->createGroup(PHYSICS_EVENT_GROUP_NAME));

    // meta data
    // run number
    auto attr_run_number = new H5::Attribute(grp->createAttribute(META_DATA_RUN_NUMBER,
                                             H5::PredType::NATIVE_INT,
                                             H5::DataSpace(H5S_SCALAR)));
    attr_run_number->write(H5::PredType::NATIVE_INT, &run_metadata.number);

    // title
    H5::StrType stype(H5::PredType::C_S1, 80);
    auto attr_title = new H5::Attribute(grp->createAttribute(META_DATA_TITLE,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_title->write(stype, run_metadata.title);

    // ts
    auto attr_ts0 = new H5::Attribute(grp->createAttribute(META_DATA_TIMESTAMP_0,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_ts0->write(H5::PredType::NATIVE_INT, &run_metadata.ts0);

    auto attr_ts1 = new H5::Attribute(grp->createAttribute(META_DATA_TIMESTAMP_1,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_ts1->write(H5::PredType::NATIVE_INT, &run_metadata.ts1);

    // date
    auto attr_date0 = new H5::Attribute(grp->createAttribute(META_DATA_DATETIME_0,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_date0->write(stype, run_metadata.date0);

    auto attr_date1 = new H5::Attribute(grp->createAttribute(META_DATA_DATETIME_1,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_date1->write(stype, run_metadata.date1);

    // dt
    auto attr_dt = new H5::Attribute(grp->createAttribute(META_DATA_ELAPSEDTIME,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_dt->write(H5::PredType::NATIVE_INT, &run_metadata.dt);

    // ring format
    H5::StrType stype1(H5::PredType::C_S1, 6);
    auto attr_fmt = new H5::Attribute(grp->createAttribute(META_DATA_RING_FORMAT,
                                      stype1,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_fmt->write(stype1, run_metadata.fmt);

    // total events
    auto attr_n_events = new H5::Attribute(grp->createAttribute(META_DATA_TOTAL_EVENTS,
                                           H5::PredType::NATIVE_LONG,
                                           H5::DataSpace(H5S_SCALAR)));
    attr_n_events->write(H5::PredType::NATIVE_LONG, &run_metadata.n_events);

    // total fragments
    auto attr_n_frags = new H5::Attribute(grp->createAttribute(META_DATA_TOTAL_FRAGMENTS,
                                          H5::PredType::NATIVE_LONG,
                                          H5::DataSpace(H5S_SCALAR)));
    attr_n_frags->write(H5::PredType::NATIVE_LONG, &run_metadata.n_frags);

    // create a dataspace for fragments data
    hsize_t dim[] = {pfragdata->size()};
    auto dspace = new H5::DataSpace(FRAGMENT_DATA_RANK, dim);

    // create a dataset under the new group
    auto dset = new H5::DataSet(grp->createDataSet(FRAGMENTS_DSET_NAME, frag_dtype, *dspace));

    // write dataset: fragments
    dset->write(pfragdata->data(), frag_dtype);

    // trace data dataset
    long tracedata_cnt = ptracedata->size();
    fprintf(stdout, "Trace data size: %i (%g KB)\n", tracedata_cnt, (float) (tracedata_cnt * sizeof(uint16_t) / 1024));
    fprintf(stdout, "Total fragments: %i (%g KB)\n", frag_cnt, (float) (frag_cnt * sizeof(FragmentData) / 1024));

    // create 2d array for trace data
    int dim0, dim1;
    int trace_length = tracedata_cnt / frag_cnt; // frag_cnt (nrows), trace_length (ncolumns)
    dim0 = frag_cnt;
    dim1 = trace_length;
    int nrows_sub = T_NROWS_PER_WRITE; // how many rows write per time
    int current_row_id = 0;          // starting at the first fragment
    uint16_t tracedata_subarr[nrows_sub][dim1];

    // auto chunk if 0x0
    if (chunk_dims[0] == 0) {
      chunk_dims[1] = trace_length;
      chunk_dims[0] = 1.0 * 1024 * 1024 / sizeof(uint16_t) / trace_length; // 1MB
    }

    // create a new dataset under the defined group, TraceData
    H5::IntType trace_dtype(H5::PredType::NATIVE_SHORT);
    trace_dtype.setOrder(H5T_ORDER_LE);

    // extensible dataset for Traces
    hsize_t trace_dims[2]; // initial dset shape
    trace_dims[0] = nrows_sub;
    trace_dims[1] = dim1;
    // hsize_t max_trace_dims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t max_trace_dims[2] = {(hsize_t) frag_cnt, (hsize_t) trace_length};
    H5::DataSpace trace_dspace(TRACE_DATA_RANK, trace_dims, max_trace_dims);

    // modify dataset creation properties, e.g. enable chunking
    H5::DSetCreatPropList cprops;
    cprops.setChunk(TRACE_DATA_RANK, chunk_dims);

    if (comp_meth == "szip") {
        cprops.setSzip(H5_SZIP_NN_OPTION_MASK, 16);
    } else if (comp_meth == "gzip") {
        cprops.setDeflate(gfactor);
    }

    // create Traces dataset
    H5::DataSet trace_dset = grp->createDataSet(TRACES_DSET_NAME, trace_dtype, trace_dspace, cprops);

    H5::DataSpace fspace;
    hsize_t size[2];
    size[0] = 0;
    size[1] = dim1;
    hsize_t offset[2];
    offset[1] = 0;
    hsize_t trace_dims1[2];
    trace_dims1[1] = dim1;

    while (current_row_id < dim0) {
        // prep data
        for(int i = 0; i < nrows_sub; i++) {
            for(int j = 0; j < dim1; j++) {
                tracedata_subarr[i][j] = (*ptracedata) [(i + current_row_id) * dim1 + j];
            }
        }

        // extend size along dim0
        size[0] += nrows_sub;

        // extend the dataset
        trace_dset.extend(size);

        // select a hyperslab
        fspace = trace_dset.getSpace();
        offset[0] = current_row_id;
        trace_dims1[0] = nrows_sub;
        fspace.selectHyperslab(H5S_SELECT_SET, trace_dims1, offset);

        trace_dset.write(tracedata_subarr, trace_dtype, trace_dspace, fspace);

        current_row_id += nrows_sub;
    }

    // create 2d array for scaler data
    int s_dim0 = pscalerts->size(), s_dim1 = (*pscalerlen)[0] + 1; // assume all scaler of the same length
    int s_nrows_sub = S_NROWS_PER_WRITE;                 //  |
    int s_current_row_id = 0;                            //  |
    uint32_t scalerdata_subarr[s_nrows_sub][s_dim1];     //  -> put first element with timestamp

    // apply chunk-dims of tracedata
    chunk_dims[1] = s_dim1;
    chunk_dims[0] = 1024 / sizeof(uint32_t) / s_dim1; // 1KB

    std::cout << s_dim0 << " " << s_dim1 << std::endl;
    std::cout << chunk_dims[0] << "x" << chunk_dims[1] << std::endl;

    // create a new dataset under the defined group, Scalers
    H5::IntType scaler_dtype(H5::PredType::NATIVE_INT);
    scaler_dtype.setOrder(H5T_ORDER_LE);

    // extensible dataset for scalers
    hsize_t scaler_dims[2]; // initial dset shape
    scaler_dims[0] = s_nrows_sub;
    scaler_dims[1] = s_dim1;
    // hsize_t max_trace_dims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    hsize_t max_scaler_dims[2] = {(hsize_t) s_dim0, (hsize_t) s_dim1};
    H5::DataSpace scaler_dspace(SCALER_DATA_RANK, scaler_dims, max_scaler_dims);

    // modify dataset creation properties, e.g. enable chunking
    H5::DSetCreatPropList s_cprops;
    s_cprops.setChunk(SCALER_DATA_RANK, chunk_dims);

    if (comp_meth == "szip") {
        s_cprops.setSzip(H5_SZIP_NN_OPTION_MASK, 16);
    } else if (comp_meth == "gzip") {
        s_cprops.setDeflate(gfactor);
    }

    // create scalers dataset
    H5::DataSet scaler_dset = grp->createDataSet(SCALERS_DSET_NAME, scaler_dtype, scaler_dspace, s_cprops);

    H5::DataSpace s_fspace;
    hsize_t s_size[2];
    s_size[0] = 0;
    s_size[1] = s_dim1;
    hsize_t s_offset[2];
    s_offset[1] = 0;
    hsize_t scaler_dims1[2];
    scaler_dims1[1] = s_dim1;

    while (s_current_row_id < s_dim0) {
        // prep data
        for(int i = 0; i < s_nrows_sub; i++) {
            // scalerdata_subarr[i][0] = (uint32_t) ((*pscalerts) [i]);
            scalerdata_subarr[i][0] = 0;
            for(int j = 0; j < s_dim1 - 1; j++) {
                scalerdata_subarr[i][j+1] = (*pscalerdata) [(i + s_current_row_id) * (s_dim1 - 1) + j];
            }
        }

        // extend size along dim0
        s_size[0] += s_nrows_sub;

        // extend the dataset
        scaler_dset.extend(s_size);

        // select a hyperslab
        s_fspace = scaler_dset.getSpace();
        s_offset[0] = s_current_row_id;
        scaler_dims1[0] = s_nrows_sub;
        s_fspace.selectHyperslab(H5S_SELECT_SET, scaler_dims1, s_offset);

        scaler_dset.write(scalerdata_subarr, scaler_dtype, scaler_dspace, s_fspace);

        s_current_row_id += s_nrows_sub;
    }

  } catch (H5::FileIException &error) {
    error.printErrorStack();
    return EXIT_FAILURE;
  } catch (H5::DataSetIException &error) {
    error.printErrorStack();
    return EXIT_FAILURE;
  } catch (H5::DataSpaceIException &error) {
    error.printErrorStack();
    return EXIT_FAILURE;
  } catch (H5::DataTypeIException &error) {
    error.printErrorStack();
    return EXIT_FAILURE;
  } catch (...) {
    fprintf(stderr, "Failed to read and export event stream to HDF5.\n");
    return EXIT_FAILURE;
  }

  fprintf(stdout, "Read and write %i fragments in %i events in total.\n", frag_cnt, event_id);

  return EXIT_SUCCESS;
}
