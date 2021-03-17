#pragma once

#include <cstdint>
#include <linux/limits.h> // PATH_MAX

#include <vector>
#include <string>
#include <iostream>

#include "argh.h"      // args parser
#include "H5Cpp.h"     // HDF5

// version
const std::string VERSION = "1.4";

// type for run meta
typedef struct RunMetaData {
  uint32_t number;    // run number
  char title[80];     // title
  time_t ts;          // posix timestamp
  char date[80];      // datetime from ts
  char fmt[6];        // ring format
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
  uint32_t overflow_code;
  uint32_t cfd_fail_bit;
  uint32_t trace_length;
  uint32_t adc_frequency;
  uint32_t adc_resolution;
  uint16_t adc_over_underflow;
} FragmentData;

// fragmentdata H5 compound data type (DDAS)
// See Also: http://docs.nscl.msu.edu/daq/newsite/ddas-1.1/DDASHit_8h_source.html
// directly readout parameters
const std::string FRAGMENT_DATA_CRATE_ID("Crate_ID");                           // uint32_t
const std::string FRAGMENT_DATA_SLOT_ID("Slot_ID");                             // uint32_t
const std::string FRAGMENT_DATA_CHANNEL_ID("Channel_ID");                       // uint32_t
const std::string FRAGMENT_DATA_TIMESTAMP("Time");                              // double
const std::string FRAGMENT_DATA_COARSE_TIME("Coarse_Time");                     // uint64_t
const std::string FRAGMENT_DATA_ENERGY("Energy");                               // uint32_t
const std::string FRAGMENT_DATA_FINISH_CODE("Err_Pileup");                      // uint32_t
const std::string FRAGMENT_DATA_OVERFLOW_CODE("Overflow_Code");                 // uint32_t
const std::string FRAGMENT_DATA_CFD_FAIL_BIT("Err_CFD_Fail");                   // uint32_t
const std::string FRAGMENT_DATA_TRACE_LENGTH("Trace_Length");                   // uint32_t
const std::string FRAGMENT_DATA_MODMSPS("ADC_Frequency");                       // uint32_t
const std::string FRAGMENT_DATA_ADC_RESOLUTION("ADC_Resolution");               // uint32_t
const std::string FRAGMENT_DATA_ADC_OVER_UNDER_FLOW("Err_ADC_Saturation");      // uint16_t
// others, including post-processed parameters
const std::string FRAGMENT_DATA_FRAGMENT_ID("Fragment_ID");                     // uint64_t
const std::string FRAGMENT_DATA_EVENT_ID("Event_ID");                           // uint64_t

// group name for DDAS Data
const std::string PHYSICS_EVENT_GROUP_NAME("DDAS_Data");

// dataset names under DDAS_Data group
const std::string FRAGMENTS_DSET_NAME("Fragments");
const std::string TRACES_DSET_NAME("Traces");

// fragment data rank
const int FRAGMENT_DATA_RANK = 1;

// trace data rank
const int TRACE_DATA_RANK = 2;

/*
data structure:
 / (root group)
 /DDAS_Data
  - Fragments (dataset)
  - Traces (dataset)
*/


/**
 * Argument parser
 *
 *
 */
class ArgumentParser {

    public:
        ArgumentParser();
        ~ArgumentParser();

        // default param value
        std::string get_default(std::string name);

        // parse args
        void parse(int argc, char** argv);

        // has option?
        bool has_option(std::string opt);

        // print help, if -h or missing mandatory: -i, return True if printed
        bool print_help_if_possible();

        // print version info, if -v, return True if printed
        bool print_version_if_possible();

        void print_help();
        void print_examples();

        // version, if --version
        std::string get_version();

        // gzip compress level
        int get_gzip_compress_level();

        // compress method
        std::string get_compress_method();

        // chunk size
        hsize_t* get_chunk_dims();

        // max evt number
        uint64_t get_max_evt();

        // verbose?
        bool is_verbose();

        // input data source (original arg)
        std::string get_input_data_source();

        // valida uri input data source
        std::string get_input_data_source_as_uri();

        // validate input data source
        bool validate_input_data_source();

        // full path for output data file
        std::string get_output_filepath();

        // print all args
        void print_all_args();

        // prog name and desc
        std::string get_progname();
        std::string get_progdesc();

    private:
        std::string m_prog;
        std::string m_progname;
        std::string m_progdesc;
        argh::parser m_parser;

        std::map<std::string, std::string> m_default_params;

        int m_gzip_comp_level;
        std::string m_comp_method;
        hsize_t *m_chunk_dims;
        uint64_t m_max_evt;

        bool m_verbose;

        std::string m_ifname;
        std::string m_ofname;
        std::string m_ifname_full;
        std::string m_ofname_full;

        void init_options();
};
