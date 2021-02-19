// C imports
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <linux/limits.h> // PATH_MAX
#include <unistd.h>

// C++ imports
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

// NSCLDAQ imports
#include <CDataSource.h>
#include <CDataSourceFactory.h>
#include <CErrnoException.h>
#include <CPhysicsEventItem.h>
#include <CPortManagerException.h>
#include <CRingItemFactory.h>
#include <DataFormat.h>
#include <FragmentIndex.h>

// DDAS imports
#include <DDASHitUnpacker.h>

// HDF5
#include "h5readout.h"

using namespace DAQ::DDAS;

static void print_trace(std::vector<uint16_t> &trace) {
  for (auto &i : trace) {
    std::cout << i << " ";
  }
  std::cout << std::endl;
}

uint16_t process_item(uint64_t event_id, CPhysicsEventItem &item,
                         std::vector<FragmentData> *fragdata,
                         std::vector<uint16_t> *tracedata) {
  FragmentIndex indexer{static_cast<uint16_t *>(item.getBodyPointer())};
  int total_fragmts = indexer.getNumberFragments();

  CRingItemFactory factory;
  DDASHitUnpacker unpacker;

  uint64_t frag_i = 0;
  // Parse fragments into hits
  for (auto &fragment : indexer) {
    // data container for fragment
    FragmentData i_fragdata{.event_id = event_id};
    frag_i++;
    //
    // fprintf(stdout, "Processing event %i, fragment %i.\n", event_id, fragment.s_sourceId);
    //
    CRingItem *frag_raw =
        factory.createRingItem(reinterpret_cast<uint8_t *>(fragment.s_itemhdr));
    std::unique_ptr<CPhysicsEventItem> frag_item(
        dynamic_cast<CPhysicsEventItem *>(frag_raw));

    uint32_t *begin = reinterpret_cast<uint32_t *>(frag_item->getBodyPointer());
    uint32_t *sentinel = begin + frag_item->getBodySize() / sizeof(uint32_t);

    DDASHit hit;
    unpacker.unpack(begin, sentinel, hit);

    uint32_t crate_id = static_cast<uint32_t>(hit.GetCrateID());
    i_fragdata.crate_id = crate_id;

    uint32_t slot_id = static_cast<uint32_t>(hit.GetSlotID());
    i_fragdata.slot_id = slot_id;

    uint32_t channel_id = static_cast<uint32_t>(hit.GetChannelID());
    i_fragdata.channel_id = channel_id;

    double ts = hit.GetTime();
    i_fragdata.timestamp = ts;

    uint64_t ts_coarse = static_cast<uint64_t>(hit.GetCoarseTime());
    i_fragdata.coarse_time = ts_coarse;

    uint32_t energy = static_cast<uint32_t>(hit.GetEnergy());
    i_fragdata.energy = energy;

    uint32_t time_high = static_cast<uint32_t>(hit.GetTimeHigh());
    i_fragdata.time_high = time_high;

    uint32_t time_low = static_cast<uint32_t>(hit.GetTimeLow());
    i_fragdata.time_low = time_low;

    uint32_t time_cfd = static_cast<uint32_t>(hit.GetTimeCFD());
    i_fragdata.time_cfd = time_cfd;

    uint32_t pileup = static_cast<uint32_t>(hit.GetFinishCode());
    i_fragdata.finish_code = pileup;

    uint32_t chan_len = hit.GetChannelLength();
    i_fragdata.channel_length = chan_len;

    uint32_t chan_hlen = hit.GetChannelLengthHeader();
    i_fragdata.channel_header_length = chan_hlen;

    uint32_t overflow_code = hit.GetOverflowCode();
    i_fragdata.overflow_code = overflow_code;

    uint32_t cfd_trig_source_bit = hit.GetCFDTrigSource();
    i_fragdata.cfd_trig_source_bit = cfd_trig_source_bit;

    uint32_t cfd_fail_bit = static_cast<uint32_t>(hit.GetCFDFailBit());
    i_fragdata.cfd_fail_bit = cfd_fail_bit;

    uint32_t trace_len = hit.GetTraceLength();
    i_fragdata.trace_length = trace_len;

    uint32_t frequency = hit.GetModMSPS();
    i_fragdata.adc_frequency = frequency;

    uint32_t hw_rev = static_cast<uint32_t>(hit.GetHardwareRevision());
    i_fragdata.hardware_revision = hw_rev;

    uint32_t resolution = hit.GetADCResolution();
    i_fragdata.adc_resolution = resolution;

    uint16_t adc_over_underflow =
        static_cast<uint16_t>(hit.GetADCOverflowUnderflow());
    uint16_t flags =
        (adc_over_underflow << 2) | (pileup << 1) | (cfd_fail_bit << 0);
    i_fragdata.adc_over_underflow = adc_over_underflow;

    /*
    fprintf(stdout, " Get %20s: %i\n", "CrateID"             , crate_id);
    fprintf(stdout, " Get %20s: %i\n", "SlotID"              , slot_id);
    fprintf(stdout, " Get %20s: %i\n", "ChannelID"           , channel_id);
    fprintf(stdout, " Get %20s: %f\n", "Time"                , ts);
    fprintf(stdout, " Get %20s: %d\n", "CoarseTime"          , ts_coarse);
    fprintf(stdout, " Get %20s: %i\n", "Energy"              , energy);
    fprintf(stdout, " Get %20s: %i\n", "Time High"           , time_high);
    fprintf(stdout, " Get %20s: %i\n", "Time Low"            , time_low);
    fprintf(stdout, " Get %20s: %i\n", "Time CFD"            , time_cfd);
    fprintf(stdout, " Get %20s: %i\n", "FinishCode"          , pileup);
    fprintf(stdout, " Get %20s: %i\n", "ChannelLength"       , chan_len);
    fprintf(stdout, " Get %20s: %i\n", "ChannelHeaderLength" , chan_hlen);
    fprintf(stdout, " Get %20s: %i\n", "Overflow Code"       , overflow_code);
    fprintf(stdout, " Get %20s: %i\n", "CFD Trig Source"     , cfd_trig_source_bit);
    fprintf(stdout, " Get %20s: %i\n", "CFDFailBit"          , cfd_fail_bit);
    fprintf(stdout, " Get %20s: %i\n", "TraceLength"         , trace_len);
    fprintf(stdout, " Get %20s: %i\n", "ADC Frequency"       , frequency);
    fprintf(stdout, " Get %20s: %i\n", "HardwareRevision"    , hw_rev);
    fprintf(stdout, " Get %20s: %i\n", "ADCResolution"       , resolution);
    fprintf(stdout, " Get %20s: %i\n", "ADCOverflowUnderflow", adc_over_underflow);
    */
    std::vector<uint16_t> &trace = hit.GetTrace();
    uint16_t *ptr = trace.data();
    for(int i = 0; i < trace_len; i++) {
         tracedata->push_back(*ptr++);
    }

    fragdata->push_back(i_fragdata);
  }
  return frag_i;
}

int main(int argc, char *argv[]) {

  if (argc != 2) {
    fprintf(stderr, "Usage: %s source.evt\n", argv[0]);
    fprintf(stderr, "   Output will be named source.evt.evtdata\n");
    return EXIT_FAILURE;
  }

  const char *source = argv[1];

  char fullsource[PATH_MAX + 1];
  if (realpath(source, fullsource) == nullptr) {
    fprintf(stderr, "Failed to find the real path for '%s': %s", source,
            strerror(errno));
    return EXIT_FAILURE;
  }

  char uri[PATH_MAX + 8];
  snprintf(uri, sizeof(uri), "file://%s", fullsource);

  // Ring Item types that can be sampled
  std::vector<std::uint16_t> sample;
  // Ring Item types that can be filtered out
  std::vector<std::uint16_t> exclude; // = {PHYSICS_EVENT};

  CDataSource *data_source;
  try {
    data_source = CDataSourceFactory::makeSource(uri, sample, exclude);
  } catch (CPortManagerException &ex) {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  } catch (CErrnoException &ex) {
    fprintf(stderr, "Failed to open data source: %s\n", ex.ReasonText());
    return EXIT_FAILURE;
  }

  printf("Starting event stream for %s\n", uri);

  char output[PATH_MAX + 9] = {};
  snprintf(output, sizeof(output), "%s.h5", fullsource);
  printf("Dumping events to file %s\n", output);

  // container for all fragments (exclude trace)
  std::vector<FragmentData> *pfragdata = new std::vector<FragmentData>();
  // container for all trace data
  std::vector<uint16_t> *ptracedata = new std::vector<uint16_t>();

  CRingItem *pItem;
  uint64_t event_id = 0;
  uint64_t frag_id = 0;
  uint64_t frag_i = 0;
  uint64_t max_evt_cnt = INT_MAX;
  try {
    // H5::Exception::dontPrint();
    // create an h5 file handle
    auto *h5file = new H5::H5File(output, H5F_ACC_TRUNC);

    // create mem data type for FragmentData
    const H5::CompType frag_dtype(sizeof(FragmentData));
    frag_dtype.insertMember(FRAGMENT_DATA_EVENT_ID,              HOFFSET(FragmentData, event_id),              H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_TIMESTAMP,             HOFFSET(FragmentData, timestamp),             H5::PredType::NATIVE_DOUBLE);
    frag_dtype.insertMember(FRAGMENT_DATA_COARSE_TIME,           HOFFSET(FragmentData, coarse_time),           H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_ENERGY,                HOFFSET(FragmentData, energy),                H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TRACE_LENGTH,          HOFFSET(FragmentData, trace_length),          H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CRATE_ID,              HOFFSET(FragmentData, crate_id),              H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_SLOT_ID,               HOFFSET(FragmentData, slot_id),               H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_ID,            HOFFSET(FragmentData, channel_id),            H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TIME_HIGH,             HOFFSET(FragmentData, time_high),             H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TIME_LOW,              HOFFSET(FragmentData, time_low),              H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TIME_CFD,              HOFFSET(FragmentData, time_cfd),              H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_FINISH_CODE,           HOFFSET(FragmentData, finish_code),           H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_LENGTH,        HOFFSET(FragmentData, channel_length),        H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_HEADER_LENGTH, HOFFSET(FragmentData, channel_header_length), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_OVERFLOW_CODE,         HOFFSET(FragmentData, overflow_code),         H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CFD_TRIG_SOURCE_BIT,   HOFFSET(FragmentData, cfd_trig_source_bit),   H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CFD_FAIL_BIT,          HOFFSET(FragmentData, cfd_fail_bit),          H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_MODMSPS,               HOFFSET(FragmentData, adc_frequency),         H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_HW_REVISION,           HOFFSET(FragmentData, hardware_revision),     H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_RESOLUTION,        HOFFSET(FragmentData, adc_resolution),        H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_OVER_UNDER_FLOW,   HOFFSET(FragmentData, adc_over_underflow),    H5::PredType::NATIVE_SHORT);

    //
    while (event_id < max_evt_cnt && (pItem = data_source->getItem())) {
      std::unique_ptr<CRingItem> item(pItem);
      std::unique_ptr<CRingItem> castableItem(
          CRingItemFactory::createRingItem(*item));

      // We're only interested in Physics Events
      if (castableItem->type() != PHYSICS_EVENT)
        continue;

      // each event data goes to one group named "Event<event_id>" under root
      frag_i = process_item(event_id++, dynamic_cast<CPhysicsEventItem &>(*castableItem),
                            pfragdata, ptracedata);
      frag_id += frag_i;
    }

    // create a new group "PhysicsEvent"
    auto *grp = new H5::Group(h5file->createGroup(PHYSICS_EVENT_GROUP_NAME));

    // create a dataspace for fragments data
    hsize_t dim[] = {pfragdata->size()};
    auto *dspace = new H5::DataSpace(FRAGMENT_DATA_RANK, dim);

    // create a dataset under the new group
    auto *dset = new H5::DataSet(
        grp->createDataSet(FRAGMENTS_DSET_NAME, frag_dtype, *dspace));
    // write dataset
    dset->write(pfragdata->data(), frag_dtype);

    // trace data as another dataset
    fprintf(stdout, "Trace data size: %i\n", ptracedata->size());
    fprintf(stdout, "Total fragments: %i\n", frag_id);

    // create 2d array for trace data
    int dim0, dim1;
    int trace_length = ptracedata->size() / frag_id; // frag_id (nrows), trace_length (ncolumns)
    dim0 = frag_id;
    dim1 = trace_length;
    int nrows_sub = 10; // how many rows write per time
    int current_row_id = 0;
    uint16_t tracedata_subarr[nrows_sub][dim1];

    // create a new dataset under the defined group, TraceData
    H5::IntType trace_dtype(H5::PredType::NATIVE_SHORT);
    trace_dtype.setOrder(H5T_ORDER_LE);

    // extensible dataset for Traces
    hsize_t trace_dims[2] = {nrows_sub, dim1}; // initial dset shape
    hsize_t max_trace_dims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
    H5::DataSpace trace_dspace(TRACE_DATA_RANK, trace_dims, max_trace_dims);

    // modify dataset creation properties, e.g. enable chunking
    H5::DSetCreatPropList cprops;
    hsize_t chunk_dims[2] = {10, 1000};
    cprops.setChunk(TRACE_DATA_RANK, chunk_dims);

    // create Traces dataset
    H5::DataSet trace_dset = grp->createDataSet(TRACES_DSET_NAME, trace_dtype, trace_dspace, cprops);

    H5::DataSpace fspace;
    hsize_t size[2];
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

        /*
        std::cout << "Extend size to: " << size[0] << ", " << size[1] << std::endl;
        std::cout << "Starting row: " << current_row_id << std::endl;
        std::cout << "Offset: " << offset[0] << ", " << offset[1] << std::endl;
        std::cout << "Sub dataset size: " << trace_dims1[0] << ", " << trace_dims[1] << std::endl;
        */

        trace_dset.write(tracedata_subarr, trace_dtype, trace_dspace, fspace);
        current_row_id += nrows_sub;
    }

    delete pfragdata;
    delete ptracedata;
    delete dset;
    delete dspace;
    delete grp;
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
  fprintf(stdout, "Read and process %i events in total.\n", frag_id);
  return EXIT_SUCCESS;
}
