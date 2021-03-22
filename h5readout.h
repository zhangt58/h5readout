#pragma once

#include <cstdint>
#include <linux/limits.h> // PATH_MAX

#include <vector>
#include <string>
#include <iostream>

#include "argh.h"      // args parser
#include "H5Cpp.h"     // HDF5

// version
const std::string VERSION = "2.0";

// type for run meta
typedef struct RunMetaData {
  uint32_t number;    // run number
  char title[80];     // title
  time_t ts0;         // begin timestamp
  time_t ts1;         // end timestamp
  char date0[80];     // datetime from ts0
  char date1[80];     // datetime from ts1
  uint32_t dt;        // elapsed time in sec
  char fmt[6];        // ring format
  uint64_t n_events;  // total events
  uint64_t n_frags;   // total fragments
} RunMetaData;

// runmetadata compound data type
const std::string META_DATA_RUN_NUMBER("Run_Number");           // uint32_t
const std::string META_DATA_TITLE("Title");                     // char[80]
const std::string META_DATA_TIMESTAMP_0("Begin_Timestamp");     // time_t
const std::string META_DATA_TIMESTAMP_1("End_Timestamp");       // time_t
const std::string META_DATA_ELAPSEDTIME("Elapsed_Seconds");     // uint32_t
const std::string META_DATA_DATETIME_0("Begin_Datetime");       // char[80]
const std::string META_DATA_DATETIME_1("End_Datetime");         // char[80]
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
const std::string SCALERS_DSET_NAME("Scalers");

// fragment data rank
const int FRAGMENT_DATA_RANK = 1;

// trace data rank
const int TRACE_DATA_RANK = 2;

// scaler data rank
const int SCALER_DATA_RANK = 1;

// type for DDAS scaler data ()
typedef struct DDASScalerData {
  time_t timestamp;
  uint32_t slot_id;      // beginning at 2
  char datetime[80];
  uint32_t crate_id;
  uint32_t inc_raw_ch00; // raw trigger incremental count at slot-id.ch0
  uint32_t inc_val_ch00; // validated trigger incremental count
  uint32_t inc_raw_ch01;
  uint32_t inc_val_ch01;
  uint32_t inc_raw_ch02;
  uint32_t inc_val_ch02;
  uint32_t inc_raw_ch03;
  uint32_t inc_val_ch03;
  uint32_t inc_raw_ch04;
  uint32_t inc_val_ch04;
  uint32_t inc_raw_ch05;
  uint32_t inc_val_ch05;
  uint32_t inc_raw_ch06;
  uint32_t inc_val_ch06;
  uint32_t inc_raw_ch07;
  uint32_t inc_val_ch07;
  uint32_t inc_raw_ch08;
  uint32_t inc_val_ch08;
  uint32_t inc_raw_ch09;
  uint32_t inc_val_ch09;
  uint32_t inc_raw_ch10;
  uint32_t inc_val_ch10;
  uint32_t inc_raw_ch11;
  uint32_t inc_val_ch11;
  uint32_t inc_raw_ch12;
  uint32_t inc_val_ch12;
  uint32_t inc_raw_ch13;
  uint32_t inc_val_ch13;
  uint32_t inc_raw_ch14;
  uint32_t inc_val_ch14;
  uint32_t inc_raw_ch15;
  uint32_t inc_val_ch15;
} DDASScalerData;

// for H5 compound data type
const std::string DDAS_SCALER_TIMESTAMP("Timestamp");// time_t
const std::string DDAS_SCALER_DATETIME("Datetime");  // char[80]
const std::string DDAS_SCALER_CRATE_ID("Crate_ID");  //uint32_t
const std::string DDAS_SCALER_SLOT_ID("Slot_ID");    //uint32_t
const std::string DDAS_SCALER_RAW_CH00("Raw_Ch00");  //uint32_t
const std::string DDAS_SCALER_VAL_CH00("Val_Ch00");  //uint32_t
const std::string DDAS_SCALER_RAW_CH01("Raw_Ch01");  //uint32_t
const std::string DDAS_SCALER_VAL_CH01("Val_Ch01");  //uint32_t
const std::string DDAS_SCALER_RAW_CH02("Raw_Ch02");  //uint32_t
const std::string DDAS_SCALER_VAL_CH02("Val_Ch02");  //uint32_t
const std::string DDAS_SCALER_RAW_CH03("Raw_Ch03");  //uint32_t
const std::string DDAS_SCALER_VAL_CH03("Val_Ch03");  //uint32_t
const std::string DDAS_SCALER_RAW_CH04("Raw_Ch04");  //uint32_t
const std::string DDAS_SCALER_VAL_CH04("Val_Ch04");  //uint32_t
const std::string DDAS_SCALER_RAW_CH05("Raw_Ch05");  //uint32_t
const std::string DDAS_SCALER_VAL_CH05("Val_Ch05");  //uint32_t
const std::string DDAS_SCALER_RAW_CH06("Raw_Ch06");  //uint32_t
const std::string DDAS_SCALER_VAL_CH06("Val_Ch06");  //uint32_t
const std::string DDAS_SCALER_RAW_CH07("Raw_Ch07");  //uint32_t
const std::string DDAS_SCALER_VAL_CH07("Val_Ch07");  //uint32_t
const std::string DDAS_SCALER_RAW_CH08("Raw_Ch08");  //uint32_t
const std::string DDAS_SCALER_VAL_CH08("Val_Ch08");  //uint32_t
const std::string DDAS_SCALER_RAW_CH09("Raw_Ch09");  //uint32_t
const std::string DDAS_SCALER_VAL_CH09("Val_Ch09");  //uint32_t
const std::string DDAS_SCALER_RAW_CH10("Raw_Ch10");  //uint32_t
const std::string DDAS_SCALER_VAL_CH10("Val_Ch10");  //uint32_t
const std::string DDAS_SCALER_RAW_CH11("Raw_Ch11");  //uint32_t
const std::string DDAS_SCALER_VAL_CH11("Val_Ch11");  //uint32_t
const std::string DDAS_SCALER_RAW_CH12("Raw_Ch12");  //uint32_t
const std::string DDAS_SCALER_VAL_CH12("Val_Ch12");  //uint32_t
const std::string DDAS_SCALER_RAW_CH13("Raw_Ch13");  //uint32_t
const std::string DDAS_SCALER_VAL_CH13("Val_Ch13");  //uint32_t
const std::string DDAS_SCALER_RAW_CH14("Raw_Ch14");  //uint32_t
const std::string DDAS_SCALER_VAL_CH14("Val_Ch14");  //uint32_t
const std::string DDAS_SCALER_RAW_CH15("Raw_Ch15");  //uint32_t
const std::string DDAS_SCALER_VAL_CH15("Val_Ch15");  //uint32_t

/*
data structure:
 / (root group)
 /DDAS_Data
  - Fragments (dataset)
  - Traces (dataset)
  - Scalers (dataset)
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

        // chunk size, pointer
        hsize_t* get_chunk_dims();

        // max evt number
        uint64_t get_max_evt();

        // verbosity level, 0 (default), 1, 2
        // 0: no messages
        // 1: scaler
        // 2: scaler, events
        // 3: all
        int get_verbosity();

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

        int m_verbosity;

        std::string m_ifname;
        std::string m_ofname;
        std::string m_ifname_full;
        std::string m_ofname_full;

        void init_options();
};
