#include <cstring>
#include <ctime>

#include <iostream>
#include <algorithm>

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
  ArgumentParser argparser = ArgumentParser();
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
  hsize_t* chunk_dims = argparser.get_chunk_dims();

  // verbosity level
  int verbosity = argparser.get_verbosity();

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
                        verbosity);
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
    H5::H5File* h5file = new H5::H5File(output_filepath, H5F_ACC_TRUNC);

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
    H5::Group* grp = new H5::Group(h5file->createGroup(PHYSICS_EVENT_GROUP_NAME));

    // meta data
    // run number
    H5::Attribute* attr_run_number = new H5::Attribute(grp->createAttribute(META_DATA_RUN_NUMBER,
                                             H5::PredType::NATIVE_INT,
                                             H5::DataSpace(H5S_SCALAR)));
    attr_run_number->write(H5::PredType::NATIVE_INT, &run_metadata.number);

    // title
    H5::StrType stype(H5::PredType::C_S1, 80);
    H5::Attribute* attr_title = new H5::Attribute(grp->createAttribute(META_DATA_TITLE,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_title->write(stype, run_metadata.title);

    // ts
    H5::Attribute* attr_ts0 = new H5::Attribute(grp->createAttribute(META_DATA_TIMESTAMP_0,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_ts0->write(H5::PredType::NATIVE_INT, &run_metadata.ts0);

    H5::Attribute* attr_ts1 = new H5::Attribute(grp->createAttribute(META_DATA_TIMESTAMP_1,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_ts1->write(H5::PredType::NATIVE_INT, &run_metadata.ts1);

    // date
    H5::Attribute* attr_date0 = new H5::Attribute(grp->createAttribute(META_DATA_DATETIME_0,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_date0->write(stype, run_metadata.date0);

    H5::Attribute* attr_date1 = new H5::Attribute(grp->createAttribute(META_DATA_DATETIME_1,
                                        stype,
                                        H5::DataSpace(H5S_SCALAR)));
    attr_date1->write(stype, run_metadata.date1);

    // dt
    H5::Attribute* attr_dt = new H5::Attribute(grp->createAttribute(META_DATA_ELAPSEDTIME,
                                      H5::PredType::NATIVE_INT,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_dt->write(H5::PredType::NATIVE_INT, &run_metadata.dt);

    // ring format
    H5::StrType stype1(H5::PredType::C_S1, 6);
    H5::Attribute* attr_fmt = new H5::Attribute(grp->createAttribute(META_DATA_RING_FORMAT,
                                      stype1,
                                      H5::DataSpace(H5S_SCALAR)));
    attr_fmt->write(stype1, run_metadata.fmt);

    // total events
    H5::Attribute* attr_n_events = new H5::Attribute(grp->createAttribute(META_DATA_TOTAL_EVENTS,
                                           H5::PredType::NATIVE_LONG,
                                           H5::DataSpace(H5S_SCALAR)));
    attr_n_events->write(H5::PredType::NATIVE_LONG, &run_metadata.n_events);

    // total fragments
    H5::Attribute* attr_n_frags = new H5::Attribute(grp->createAttribute(META_DATA_TOTAL_FRAGMENTS,
                                          H5::PredType::NATIVE_LONG,
                                          H5::DataSpace(H5S_SCALAR)));
    attr_n_frags->write(H5::PredType::NATIVE_LONG, &run_metadata.n_frags);

    // create a dataspace for fragments data
    hsize_t dim[] = {pfragdata->size()};
    H5::DataSpace* dspace = new H5::DataSpace(FRAGMENT_DATA_RANK, dim);

    // create a dataset under the new group
    H5::DataSet* dset = new H5::DataSet(grp->createDataSet(FRAGMENTS_DSET_NAME, frag_dtype, *dspace));

    // write dataset: fragments
    dset->write(pfragdata->data(), frag_dtype);

    // clean up
    delete pfragdata;

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

    // clean up
    delete ptracedata;

    // create mem data type for ScalerData
    const H5::CompType ddas_scaler_dtype(sizeof(DDASScalerData));
    ddas_scaler_dtype.insertMember(DDAS_SCALER_TIMESTAMP, HOFFSET(DDASScalerData, timestamp),    H5::PredType::NATIVE_LONG);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_DATETIME,  HOFFSET(DDASScalerData, datetime),     H5::StrType(H5::PredType::C_S1, 80));
    ddas_scaler_dtype.insertMember(DDAS_SCALER_CRATE_ID,  HOFFSET(DDASScalerData, crate_id),     H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_SLOT_ID,   HOFFSET(DDASScalerData, slot_id),      H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH00,  HOFFSET(DDASScalerData, inc_raw_ch00), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH00,  HOFFSET(DDASScalerData, inc_val_ch00), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH01,  HOFFSET(DDASScalerData, inc_raw_ch01), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH01,  HOFFSET(DDASScalerData, inc_val_ch01), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH02,  HOFFSET(DDASScalerData, inc_raw_ch02), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH02,  HOFFSET(DDASScalerData, inc_val_ch02), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH03,  HOFFSET(DDASScalerData, inc_raw_ch03), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH03,  HOFFSET(DDASScalerData, inc_val_ch03), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH04,  HOFFSET(DDASScalerData, inc_raw_ch04), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH04,  HOFFSET(DDASScalerData, inc_val_ch04), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH05,  HOFFSET(DDASScalerData, inc_raw_ch05), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH05,  HOFFSET(DDASScalerData, inc_val_ch05), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH06,  HOFFSET(DDASScalerData, inc_raw_ch06), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH06,  HOFFSET(DDASScalerData, inc_val_ch06), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH07,  HOFFSET(DDASScalerData, inc_raw_ch07), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH07,  HOFFSET(DDASScalerData, inc_val_ch07), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH08,  HOFFSET(DDASScalerData, inc_raw_ch08), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH08,  HOFFSET(DDASScalerData, inc_val_ch08), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH09,  HOFFSET(DDASScalerData, inc_raw_ch09), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH09,  HOFFSET(DDASScalerData, inc_val_ch09), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH10,  HOFFSET(DDASScalerData, inc_raw_ch10), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH10,  HOFFSET(DDASScalerData, inc_val_ch10), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH11,  HOFFSET(DDASScalerData, inc_raw_ch11), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH11,  HOFFSET(DDASScalerData, inc_val_ch11), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH12,  HOFFSET(DDASScalerData, inc_raw_ch12), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH12,  HOFFSET(DDASScalerData, inc_val_ch12), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH13,  HOFFSET(DDASScalerData, inc_raw_ch13), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH13,  HOFFSET(DDASScalerData, inc_val_ch13), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH14,  HOFFSET(DDASScalerData, inc_raw_ch14), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH14,  HOFFSET(DDASScalerData, inc_val_ch14), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_RAW_CH15,  HOFFSET(DDASScalerData, inc_raw_ch15), H5::PredType::NATIVE_INT);
    ddas_scaler_dtype.insertMember(DDAS_SCALER_VAL_CH15,  HOFFSET(DDASScalerData, inc_val_ch15), H5::PredType::NATIVE_INT);


    // create 2d array for scaler data
    int s_dim0 = pscalerts->size(), s_dim1 = (*pscalerlen)[0]; // assume all scaler of the same length
    uint32_t max_slot_num = 1 + s_dim1 / 33;  // max slot number
    uint32_t slot_num = 2;                    // first slot number
    time_t ts_tmp;
    char date_tmp[80];
    int ind0;

    // process scalerdata
    std::vector<DDASScalerData> *pddasscaler = new std::vector<DDASScalerData>();
    for(int i = 0; i < s_dim0; i++) {
        ts_tmp = (*pscalerts) [i];
        std::strftime(date_tmp, sizeof(date_tmp), "%c", std::localtime(&ts_tmp));
        slot_num = 2;
        while (slot_num <= max_slot_num) {
            // create a new DDASScalerData every 33, 33, ...
            // crate_id + 16 * 2 channel readings
            DDASScalerData i_ddas_scaler_data {
                .timestamp = ts_tmp,
                .slot_id = slot_num
            };
            strcpy(i_ddas_scaler_data.datetime, date_tmp);
            ind0 = (slot_num - 2) * 33 + i * s_dim1;
            i_ddas_scaler_data.crate_id = (*pscalerdata)[ind0];
            i_ddas_scaler_data.inc_raw_ch00 = (*pscalerdata)[ind0 + 1];
            i_ddas_scaler_data.inc_val_ch00 = (*pscalerdata)[ind0 + 2];
            i_ddas_scaler_data.inc_raw_ch01 = (*pscalerdata)[ind0 + 3];
            i_ddas_scaler_data.inc_val_ch01 = (*pscalerdata)[ind0 + 4];
            i_ddas_scaler_data.inc_raw_ch02 = (*pscalerdata)[ind0 + 5];
            i_ddas_scaler_data.inc_val_ch02 = (*pscalerdata)[ind0 + 6];
            i_ddas_scaler_data.inc_raw_ch03 = (*pscalerdata)[ind0 + 7];
            i_ddas_scaler_data.inc_val_ch03 = (*pscalerdata)[ind0 + 8];
            i_ddas_scaler_data.inc_raw_ch04 = (*pscalerdata)[ind0 + 9];
            i_ddas_scaler_data.inc_val_ch04 = (*pscalerdata)[ind0 + 10];
            i_ddas_scaler_data.inc_raw_ch05 = (*pscalerdata)[ind0 + 11];
            i_ddas_scaler_data.inc_val_ch05 = (*pscalerdata)[ind0 + 12];
            i_ddas_scaler_data.inc_raw_ch06 = (*pscalerdata)[ind0 + 13];
            i_ddas_scaler_data.inc_val_ch06 = (*pscalerdata)[ind0 + 14];
            i_ddas_scaler_data.inc_raw_ch07 = (*pscalerdata)[ind0 + 15];
            i_ddas_scaler_data.inc_val_ch07 = (*pscalerdata)[ind0 + 16];
            i_ddas_scaler_data.inc_raw_ch08 = (*pscalerdata)[ind0 + 17];
            i_ddas_scaler_data.inc_val_ch08 = (*pscalerdata)[ind0 + 18];
            i_ddas_scaler_data.inc_raw_ch09 = (*pscalerdata)[ind0 + 19];
            i_ddas_scaler_data.inc_val_ch09 = (*pscalerdata)[ind0 + 20];
            i_ddas_scaler_data.inc_raw_ch10 = (*pscalerdata)[ind0 + 21];
            i_ddas_scaler_data.inc_val_ch10 = (*pscalerdata)[ind0 + 22];
            i_ddas_scaler_data.inc_raw_ch11 = (*pscalerdata)[ind0 + 23];
            i_ddas_scaler_data.inc_val_ch11 = (*pscalerdata)[ind0 + 24];
            i_ddas_scaler_data.inc_raw_ch12 = (*pscalerdata)[ind0 + 25];
            i_ddas_scaler_data.inc_val_ch12 = (*pscalerdata)[ind0 + 26];
            i_ddas_scaler_data.inc_raw_ch13 = (*pscalerdata)[ind0 + 27];
            i_ddas_scaler_data.inc_val_ch13 = (*pscalerdata)[ind0 + 28];
            i_ddas_scaler_data.inc_raw_ch14 = (*pscalerdata)[ind0 + 29];
            i_ddas_scaler_data.inc_val_ch14 = (*pscalerdata)[ind0 + 30];
            i_ddas_scaler_data.inc_raw_ch15 = (*pscalerdata)[ind0 + 31];
            i_ddas_scaler_data.inc_val_ch15 = (*pscalerdata)[ind0 + 32];
            ++slot_num;
            pddasscaler->push_back(i_ddas_scaler_data);
        }
    }

    struct sort_key {
        inline bool operator() (const DDASScalerData& left, const DDASScalerData& right) {
            if (left.slot_id < right.slot_id) return true;
            if (left.slot_id > right.slot_id) return false;
            if (left.timestamp < right.timestamp) return true;
            if (left.timestamp > right.timestamp) return false;
            return false;
        }
    };
    // sort pddasscaler with slot ID
    std::sort(pddasscaler->begin(), pddasscaler->end(), sort_key());

    // create a dataspace for scaler data
    hsize_t s_dim[] = {pddasscaler->size()};
    H5::DataSpace* s_dspace = new H5::DataSpace(SCALER_DATA_RANK, s_dim);

    // create a dataset under the new group
    H5::DataSet* s_dset = new H5::DataSet(grp->createDataSet(SCALERS_DSET_NAME, ddas_scaler_dtype, *s_dspace));

    // write dataset: scalerdata
    s_dset->write(pddasscaler->data(), ddas_scaler_dtype);

    // clean up
    delete pscalerdata;
    delete pscalerts;
    delete pscalerlen;
    delete pddasscaler;

    //
    delete attr_run_number;
    delete attr_title;
    delete attr_ts0;
    delete attr_ts1;
    delete attr_date0;
    delete attr_date1;
    delete attr_dt;
    delete attr_fmt;
    delete attr_n_events;
    delete attr_n_frags;
    delete dspace;
    delete dset;
    delete grp;
    delete s_dspace;
    delete s_dset;
    delete h5file;

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
