#pragma once

#include <cstdint>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include "H5Cpp.h"

// version
const std::string VERSION = "1.3";

// type for run meta
typedef struct RunMetaData {
  uint32_t number;    // run number
  char title[80];     // title
  time_t ts;          // posix timestamp
  char date[80];      // datetime from ts
  char fmt[6];      // ring format
  uint64_t n_events;  // total events
  uint64_t n_frags;   // total fragments
} RunMetaData;

// runmetadata compound data type
const std::string META_DATA_RUN_NUMBER("Run_Number");           // uint32_t
const std::string META_DATA_TITLE("Title");                     // char[80]
const std::string META_DATA_TIMESTAMP("Timestamp");             // time_t
const std::string META_DATA_DATETIME("Datetime");               // char[80]
const std::string META_DATA_TOTAL_EVENTS("Total_Events");       // uint64_t
const std::string META_DATA_TOTAL_FRAGMENTS("Total_Fragments"); // uint64_t
const std::string META_DATA_RING_FORMAT("Ring_Format");         // char[6]


// type for physics fragment data
typedef struct FragmentData {
  uint64_t fragment_id;
  uint64_t event_id;
  uint32_t crate_id;
  uint32_t slot_id;
  uint32_t channel_id;
  double timestamp;
  uint64_t coarse_time;
  uint32_t energy;
  uint32_t finish_code;
  uint32_t channel_length;
  uint32_t channel_header_length;
  uint32_t overflow_code;
  uint32_t cfd_fail_bit;
  uint32_t trace_length;
  uint32_t adc_frequency;
  uint32_t hardware_revision;
  uint32_t adc_resolution;
  uint16_t adc_over_underflow;
  uint16_t error_flag;
} FragmentData;

// fragmentdata H5 compound data type
// See Also:
// http://docs.nscl.msu.edu/daq/newsite/ddas-1.1/DDASHit_8h_source.html
const std::string FRAGMENT_DATA_CRATE_ID("Crate_ID");                           // uint32_t
const std::string FRAGMENT_DATA_SLOT_ID("Slot_ID");                             // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_ID("Channel_ID");                       // uint32_t
const std::string FRAGMENT_DATA_TIMESTAMP("Time");                              // double
const std::string FRAGMENT_DATA_COARSE_TIME("Coarse_Time");                     // uint64_t
const std::string FRAGMENT_DATA_ENERGY("Energy");                               // uint32_t
const std::string FRAGMENT_DATA_FINISH_CODE("Finish_Code");                     // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_LENGTH("Channel_Length");               // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_HEADER_LENGTH("Channel_Header_Length"); // uint32_t
const std::string FRAGMENT_DATA_OVERFLOW_CODE("Overflow_Code");                 // uint32_t
const std::string FRAGMENT_DATA_CFD_FAIL_BIT("CFD_Fail_Bit");                   // uint32_t
const std::string FRAGMENT_DATA_TRACE_LENGTH("Trace_Length");                   // uint32_t
const std::string FRAGMENT_DATA_MODMSPS("ADC_Frequency");                       // uint32_t
const std::string FRAGMENT_DATA_HW_REVISION("Hardware_Revision");               // uint32_t
const std::string FRAGMENT_DATA_ADC_RESOLUTION("ADC_Resolution");               // uint32_t
const std::string FRAGMENT_DATA_ADC_OVER_UNDER_FLOW("ADC_Over/Underflow");      // uint16_t

// others, including post-processed parameters
const std::string FRAGMENT_DATA_FRAGMENT_ID("Fragment_ID");                     // uint64_t
const std::string FRAGMENT_DATA_EVENT_ID("Event_ID");                           // uint64_t
const std::string FRAGMENT_DATA_ERROR_FLAG("Error_Flag");                       // uint16_t

// group name
const std::string PHYSICS_EVENT_GROUP_NAME("DDAS_Data");

// dataset name
const std::string FRAGMENTS_DSET_NAME("Fragments");
const std::string TRACES_DSET_NAME("Traces");

// fragment data rank
const int FRAGMENT_DATA_RANK = 1;

// trace data rank
const int TRACE_DATA_RANK = 2;

/*
data structure:
 / (root group)
   Fragments (dataset)
   Traces (dataset)
*/

/*
 * Return if a path is a exsiting directory.
 */
bool is_dir(char *path) {
    struct stat sb = {};
    stat(path, &sb);
    return (sb.st_mode & S_IFMT) == S_IFDIR;
}

/*
 * Return if a path is a exsiting file.
 */
bool is_file(char *path) {
    struct stat sb = {};
    stat(path, &sb);
    return (sb.st_mode & S_IFMT) == S_IFREG;
}
