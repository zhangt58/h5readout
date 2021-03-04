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
#include <locale>

#include "cxxopts.hpp" // args parser

// NSCLDAQ imports
#include <CDataFormatItem.h>
#include <CDataSource.h>
#include <CDataSourceFactory.h>
#include <CErrnoException.h>
#include <CRingStateChangeItem.h>
#include <CRingPhysicsEventCountItem.h>
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

// rows of trace data written at once.
const int NROWS_PER_WRITE = 1;

// version
const std::string VERSION = "1.1";

/**
 * Process ringtime to extract fragment data from physics event.
 *
 */
void process_item(uint64_t event_id, uint64_t &frag_cnt, CPhysicsEventItem &item,
                  std::vector<FragmentData> *pfragdata,
                  std::vector<uint16_t> *ptracedata, bool verbose) {

//    if (item.hasBodyHeader()) {
//        std::cout << "Timestamp: " << item.getEventTimestamp() << "\n"
//                  << "Source ID: " << item.getSourceId() << std::endl;
//    }
    FragmentIndex indexer(
            static_cast<uint16_t *>(item.getBodyPointer())
    );

    int total_fragmts = indexer.getNumberFragments();
//    std::cout << total_fragmts << std::endl;

   CRingItemFactory factory;
   DDASHitUnpacker unpacker;

  // Parse fragments into hits
  for (auto &fragment : indexer) {
    // data container for fragment
    FragmentData i_fragdata {
      .fragment_id = frag_cnt++,
      .event_id = event_id
    };

    if (verbose) {
    fprintf(stdout, "Processing evt-id: %u, frgmt: %u, ts: %u, accu frgmts: %u\n",
            event_id, fragment.s_sourceId, fragment.s_timestamp, frag_cnt);
    }

    CRingItem *frag_raw = factory.createRingItem(reinterpret_cast<uint8_t *>(fragment.s_itemhdr));
    std::unique_ptr<CPhysicsEventItem> frag_item(dynamic_cast<CPhysicsEventItem *>(frag_raw));

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

    uint32_t pileup = static_cast<uint32_t>(hit.GetFinishCode());
    i_fragdata.finish_code = pileup;

    uint32_t chan_len = static_cast<uint32_t>(hit.GetChannelLength());
    i_fragdata.channel_length = chan_len;

    uint32_t chan_hlen = static_cast<uint32_t>(hit.GetChannelLengthHeader());
    i_fragdata.channel_header_length = chan_hlen;

    uint32_t overflow_code = static_cast<uint32_t>(hit.GetOverflowCode());
    i_fragdata.overflow_code = overflow_code;

    uint32_t cfd_fail_bit = static_cast<uint32_t>(hit.GetCFDFailBit());
    i_fragdata.cfd_fail_bit = cfd_fail_bit;

    uint32_t trace_len = static_cast<uint32_t>(hit.GetTraceLength());
    i_fragdata.trace_length = trace_len;

    uint32_t frequency = static_cast<uint32_t>(hit.GetModMSPS());
    i_fragdata.adc_frequency = frequency;

    uint32_t hw_rev = static_cast<uint32_t>(hit.GetHardwareRevision());
    i_fragdata.hardware_revision = hw_rev;

    uint32_t resolution = static_cast<uint32_t>(hit.GetADCResolution());
    i_fragdata.adc_resolution = resolution;

    uint16_t adc_over_underflow = static_cast<uint16_t>(hit.GetADCOverflowUnderflow());
    uint16_t error_flag = (adc_over_underflow << 2) | (pileup << 1) | (cfd_fail_bit << 0);
    i_fragdata.adc_over_underflow = adc_over_underflow;
    i_fragdata.error_flag = error_flag;

    std::vector<uint16_t> &trace = hit.GetTrace();
    uint16_t *ptr = trace.data();
    for(int i = 0; i < trace_len; i++) {
         ptracedata->push_back(*ptr++);
    }
    pfragdata->push_back(i_fragdata);
  }
}


int main(int argc, char** argv) {

  cxxopts::Options options("h5readout", "Readout DDAS event data to H5 format.");
  options.add_options()
    ("i,input", "URI for input event data", cxxopts::value<std::string>())
    ("o,output", "File path for HDF5 data", cxxopts::value<std::string>()->default_value("<input>.h5"))
    ("n,events", "Number of events to readout", cxxopts::value<int>()->default_value(std::to_string(INT_MAX)))
    ("s,chunk-size", "Chunk size MxN for HDF5 data", cxxopts::value<std::string>()->default_value("0x0"))
    ("c,compress", "Compression method, 'szip' or 'gzip'", cxxopts::value<std::string>()->default_value("gzip"))
    ("compress-level", "GZip Compression level(0-9)", cxxopts::value<int>()->default_value("8"))
    ("v,verbose", "Show verbose message", cxxopts::value<bool>()->default_value("false"))
    ("version", "Show version info", cxxopts::value<bool>())
    ("h,help", "Print this message")
   ;
  auto result = options.parse(argc, argv);
  if (result.count("version")) {
    std::cout << "Version: " << VERSION << std::endl;
    return 0;
  }

  if (result.count("help") || ! result.count("input")) {
    std::cout << options.help() << std::endl;
    return EXIT_FAILURE;
  }

  int gfactor = result["compress-level"].as<int>();
  std::string comp_meth = result["compress"].as<std::string>();

  std::string ifname = result["input"].as<std::string>();
  std::string ofname = result["output"].as<std::string>();
  if (ofname == "<input>.h5") {
    ofname = ifname + ".h5";
  }
  uint64_t max_evt_cnt = result["events"].as<int>();
  std::string chunk_size = result["chunk-size"].as<std::string>();
  std::stringstream ss(chunk_size);
  hsize_t chunk_dims[2];
  int i;
  int idx = 0;
  while (ss >> i) {
    chunk_dims[idx++] = i;
    if (ss.peek() == 'x') ss.ignore();
  }
  bool verbose = result["verbose"].as<bool>();

  char fullsource[PATH_MAX + 1];
  if (realpath(ifname.c_str(), fullsource) == nullptr) {
    fprintf(stderr, "Failed to find the real path for '%s': %s\n", ifname.c_str(),
            strerror(errno));
    return EXIT_FAILURE;
  }
  char uri[PATH_MAX + 8];
  snprintf(uri, sizeof(uri), "file://%s", fullsource);

  char outfile[PATH_MAX + 1];
  if (realpath(ofname.c_str(), outfile) != nullptr) {
    fprintf(stdout, "Warning: Overwrite output file: '%s' ([Y]/n?) ", outfile);
    char c = getc(stdin);
    if (c == 'n' || c == 'N') {
        return 0;
    }
  }

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

  // container for all fragments (exclude trace)
  std::vector<FragmentData> *pfragdata = new std::vector<FragmentData>();
  // container for all trace data
  std::vector<uint16_t> *ptracedata = new std::vector<uint16_t>();

  CRingItem *pItem;
  uint64_t event_id = 0;  // event count
  uint64_t frag_cnt = 0;  // total framgnets count

  RunMetaData run_metadata;

try {
  // start to read
  fprintf(stdout, "Reading data from: %s\n", uri);
  // readout fragment data
  while (event_id < max_evt_cnt && (pItem = data_source->getItem())) {
    std::unique_ptr<CRingItem> item(pItem);
    std::unique_ptr<CRingItem> castableItem(CRingItemFactory::createRingItem(*item));

    switch (castableItem->type()) {
      case BEGIN_RUN:
      {
        CRingStateChangeItem &item2 = dynamic_cast<CRingStateChangeItem &>(*castableItem);
        run_metadata.number = item2.getRunNumber();
        strcpy(run_metadata.title, item2.getTitle().c_str());
        run_metadata.ts = item2.getTimestamp();
        std::strftime(run_metadata.date, sizeof(run_metadata.date), "%c", std::localtime(&run_metadata.ts));
        break;
      }
      case PHYSICS_EVENT: // Physics Event
      {
        // each event data goes to one group named "Event<event_id>" under root
        CPhysicsEventItem &phyItem = dynamic_cast<CPhysicsEventItem &>(*castableItem);
        // std::cout << " " << "Process physics event" << std::endl;
        // std::cout << " " << "BodySize "  << phyItem.getBodySize() << std::endl;
        // std::cout << " " << "toString " << phyItem.toString() << std::endl;

        process_item(event_id, frag_cnt, phyItem, pfragdata, ptracedata, verbose);
        break;
      }
      case PHYSICS_EVENT_COUNT: // Physics Event Count
      {
        // CRingPhysicsEventCountItem &item1 = dynamic_cast<CRingPhysicsEventCountItem &>(*castableItem);
        // std::cout << "Physics Event Count: " << item1.getEventCount() << std::endl;
        break;
      }
      case RING_FORMAT:
      {
        CDataFormatItem &fmtItem = dynamic_cast<CDataFormatItem &>(*castableItem);
        std::stringstream fmtss;
        fmtss << fmtItem.getMajor() << "." << fmtItem.getMinor();
        strcpy(run_metadata.fmt, fmtss.str().c_str());
        break;
      }
      default:
        break;
    }
    ++event_id;
  }

} catch (std::string &ex) {
    std::cout << ex << std::endl;
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

  // start to write
  fprintf(stdout, "Writing data to: %s\n", outfile);

  try {
    // H5::Exception::dontPrint();
    // create an h5 file handle
    auto h5file = new H5::H5File(outfile, H5F_ACC_TRUNC);

    // create mem data type for FragmentData
    const H5::CompType frag_dtype(sizeof(FragmentData));
    frag_dtype.insertMember(FRAGMENT_DATA_FRAGMENT_ID,           HOFFSET(FragmentData, fragment_id),           H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_EVENT_ID,              HOFFSET(FragmentData, event_id),              H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_TIMESTAMP,             HOFFSET(FragmentData, timestamp),             H5::PredType::NATIVE_DOUBLE);
    frag_dtype.insertMember(FRAGMENT_DATA_COARSE_TIME,           HOFFSET(FragmentData, coarse_time),           H5::PredType::NATIVE_LONG);
    frag_dtype.insertMember(FRAGMENT_DATA_ENERGY,                HOFFSET(FragmentData, energy),                H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_TRACE_LENGTH,          HOFFSET(FragmentData, trace_length),          H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CRATE_ID,              HOFFSET(FragmentData, crate_id),              H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_SLOT_ID,               HOFFSET(FragmentData, slot_id),               H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_ID,            HOFFSET(FragmentData, channel_id),            H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_LENGTH,        HOFFSET(FragmentData, channel_length),        H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_CHANNEL_HEADER_LENGTH, HOFFSET(FragmentData, channel_header_length), H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_FINISH_CODE,           HOFFSET(FragmentData, finish_code),           H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_OVER_UNDER_FLOW,   HOFFSET(FragmentData, adc_over_underflow),    H5::PredType::NATIVE_SHORT);
    frag_dtype.insertMember(FRAGMENT_DATA_CFD_FAIL_BIT,          HOFFSET(FragmentData, cfd_fail_bit),          H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_OVERFLOW_CODE,         HOFFSET(FragmentData, overflow_code),         H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ERROR_FLAG,            HOFFSET(FragmentData, error_flag),            H5::PredType::NATIVE_SHORT);
    frag_dtype.insertMember(FRAGMENT_DATA_MODMSPS,               HOFFSET(FragmentData, adc_frequency),         H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_ADC_RESOLUTION,        HOFFSET(FragmentData, adc_resolution),        H5::PredType::NATIVE_INT);
    frag_dtype.insertMember(FRAGMENT_DATA_HW_REVISION,           HOFFSET(FragmentData, hardware_revision),     H5::PredType::NATIVE_INT);

    // create a new group "PhysicsEvent"
    auto *grp = new H5::Group(h5file->createGroup(PHYSICS_EVENT_GROUP_NAME));

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
    auto attr_ts = new H5::Attribute(grp->createAttribute(META_DATA_TIMESTAMP,
                                     H5::PredType::NATIVE_INT,
                                     H5::DataSpace(H5S_SCALAR)));
    attr_ts->write(H5::PredType::NATIVE_INT, &run_metadata.ts);

    // date
    auto attr_date = new H5::Attribute(grp->createAttribute(META_DATA_DATETIME,
                                       stype,
                                       H5::DataSpace(H5S_SCALAR)));
    attr_date->write(stype, run_metadata.date);

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

    // write dataset
    dset->write(pfragdata->data(), frag_dtype);

    // trace data as another dataset
    long tracedata_cnt = ptracedata->size();
    fprintf(stdout, "Trace data size: %i (%g KB)\n", tracedata_cnt, (float) (tracedata_cnt * sizeof(uint16_t) / 1024));
    fprintf(stdout, "Total fragments: %i (%g KB)\n", frag_cnt, (float) (frag_cnt * sizeof(FragmentData) / 1024));

    // create 2d array for trace data
    int dim0, dim1;
    int trace_length = tracedata_cnt / frag_cnt; // frag_cnt (nrows), trace_length (ncolumns)
    dim0 = frag_cnt;
    dim1 = trace_length;
    int nrows_sub = NROWS_PER_WRITE; // how many rows write per time
    int current_row_id = 0;          // starting at the first fragment
    uint16_t tracedata_subarr[nrows_sub][dim1];

    // auto chunk if 0x0
    if (chunk_dims[0] == 0) {
      chunk_dims[1] = trace_length;
      chunk_dims[0] = 1000 * 1000 / sizeof(uint16_t) / trace_length;
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

        /*
        std::cout << "Extend size to: " << size[0] << ", " << size[1] << std::endl;
        std::cout << "Starting row: " << current_row_id << std::endl;
        std::cout << "Offset: " << offset[0] << ", " << offset[1] << std::endl;
        std::cout << "Sub dataset size: " << trace_dims1[0] << ", " << trace_dims[1] << std::endl;
        */

        trace_dset.write(tracedata_subarr, trace_dtype, trace_dspace, fspace);

        current_row_id += nrows_sub;
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
