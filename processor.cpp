#include <iostream>
#include <iomanip>
#include <memory>

#include <cstring>
#include <ctime>

#include <FragmentIndex.h>
#include <CRingItemFactory.h>

//#include <fragment.h>
//#include <DataFormat.h>
//#include "FragmentIndexXYZ.h"

#include <DDASHitUnpacker.h>

#include "processor.h"

void processRingItem(CRingItemProcessor &processor, CRingItem *item,
                     RunMetaData &run_meta,
                     uint64_t &event_id, uint64_t &frag_cnt,
                     std::vector<FragmentData> *pfragdata,
                     std::vector<uint16_t> *ptracedata,
                     std::vector<time_t> *pscalerts,
                     std::vector<uint32_t> *pscalerlen,
                     std::vector<uint32_t> *pscalerdata,
                     int &verbosity,
                     std::string &ctrl_type,
                     std::vector<std::string> *pmodules,
                     std::vector<uint64_t> *pfragdata_v785)
{
    std::unique_ptr<CRingItem> item_(item);
    std::unique_ptr<CRingItem> castableItem(CRingItemFactory::createRingItem(*item_));

    uint16_t item_type = castableItem->type();
    switch (item_type)
    {
    case BEGIN_RUN:
    case END_RUN:
    case PAUSE_RUN:
    case RESUME_RUN:
    {
        CRingStateChangeItem &state_change_item = dynamic_cast<CRingStateChangeItem &>(*castableItem);
        processor.processStateChangeItem(state_change_item, run_meta, item_type);
        break;
    }
    case PERIODIC_SCALERS:
    {

        CRingScalerItem &scaler_item = dynamic_cast<CRingScalerItem &>(*castableItem);
        processor.processScalerItem(scaler_item, pscalerts, pscalerlen, pscalerdata, verbosity);
        break;
    }
    case PACKET_TYPES:
    case MONITORED_VARIABLES:
    {
        CRingTextItem &text_item = dynamic_cast<CRingTextItem &>(*castableItem);
        processor.processTextItem(text_item, verbosity);
        break;
    }
    case PHYSICS_EVENT:
    {
        // count events
        CPhysicsEventItem &phy_item = dynamic_cast<CPhysicsEventItem &>(*castableItem);
        if (ctrl_type == "DDAS")
        {
            processor.processEvent(phy_item, event_id, frag_cnt, pfragdata, ptracedata, verbosity);
            event_id++;
        }
        else if (ctrl_type == "VME")
        {
            processor.processEvent_VME(phy_item, event_id, frag_cnt, pfragdata_v785, verbosity);
            event_id++;
        }
        break;
    }
    case PHYSICS_EVENT_COUNT:
    {
        CRingPhysicsEventCountItem &phy_cnt_item = dynamic_cast<CRingPhysicsEventCountItem &>(*castableItem);
        processor.processEventCount(phy_cnt_item);
        break;
    }
    case RING_FORMAT:
    {
        CDataFormatItem &fmt_item = dynamic_cast<CDataFormatItem &>(*castableItem);
        processor.processFormat(fmt_item, run_meta);
        break;
    }
    case EVB_GLOM_INFO:
    {
        CGlomParameters &glom_item = dynamic_cast<CGlomParameters &>(*castableItem);
        processor.processGlomParams(glom_item, verbosity);
        break;
    }
    default:
    {
        processor.processUnknownItem(*item_, verbosity);
        break;
    }
    }
}

/**
 * processStateChangeItem
 *   Extract metadata.
 *
 */
void CRingItemProcessor::processStateChangeItem(CRingStateChangeItem &item,
                                                RunMetaData &run_meta, uint16_t &item_type)
{
    if (item_type == BEGIN_RUN)
    {
        run_meta.number = item.getRunNumber();
        strcpy(run_meta.title, item.getTitle().c_str());
        run_meta.ts0 = item.getTimestamp();
        std::strftime(run_meta.date0, sizeof(run_meta.date0), "%c", std::localtime(&run_meta.ts0));
        run_meta.ts1 = item.getTimestamp();
        std::strftime(run_meta.date1, sizeof(run_meta.date1), "%c", std::localtime(&run_meta.ts1));
        run_meta.dt = item.getElapsedTime();
    }
    else if (item_type == END_RUN)
    {
        run_meta.ts1 = item.getTimestamp();
        std::strftime(run_meta.date1, sizeof(run_meta.date1), "%c", std::localtime(&run_meta.ts1));
        run_meta.dt = item.getElapsedTime();
    }
}

/**
 * processScalerItem
 *   Deal with scaler data.
 *
 */
void CRingItemProcessor::processScalerItem(CRingScalerItem &item,
                                           std::vector<time_t> *pscalerts,
                                           std::vector<uint32_t> *pscalerlen,
                                           std::vector<uint32_t> *pscalerdata,
                                           int &verbosity)
{
    time_t tm = item.getTimestamp();
    uint32_t cnt = item.getScalerCount();
    pscalerts->push_back(tm);
    pscalerlen->push_back(cnt);

    char tm_date[1000];
    std::strftime(tm_date, 1000, "%c", std::localtime(&tm));

    for (auto const &v : item.getScalers())
    {
        pscalerdata->push_back(v);
    }

    if (verbosity > 0)
    {
        std::cout << "--> " << item.typeName() << ": ";
        std::cout << tm_date << " "
                  << "(" << cnt << "):";

        for (auto const &v : item.getScalers())
        {
            std::cout << " " << v;
        }
        std::cout << std::endl;
    }
}

/**
 * processTextItem
 *
 */
void CRingItemProcessor::processTextItem(CRingTextItem &item, int &verbosity)
{
    time_t tm = item.getTimestamp();
    char tm_date[1000];
    std::strftime(tm_date, 1000, "%c", std::localtime(&tm));
    if (verbosity > 2)
    {
        std::cout << "--> " << item.typeName() << ": "
                  << tm_date << " " << item.getTimeOffset() << " seconds into the run\n";
        std::cout << "---> Here are the record strings: \n";
        for (auto const &s : item.getStrings())
        {
            std::cout << "----> " << s << "\n";
        }
    }
}

/**
 * processEvent_VME
 *   process event item from VME/V785 module
 * 
 */

void CRingItemProcessor::processEvent_VME(CPhysicsEventItem &item,
                                          uint64_t &event_id, uint64_t &frag_cnt,
                                          std::vector<uint64_t> *pfragdata,
                                          int &verbosity)

{
    CRingItem *pItem = CRingItemFactory::createRingItem(item);
    uint32_t size_in_byte = pItem->getBodySize();
    uint32_t nsize = size_in_byte / sizeof(uint16_t);
    uint16_t *p = reinterpret_cast<uint16_t *>(pItem->getBodyPointer());

    // if (item.hasBodyHeader())
    // {
    //     std::cout << " Timestamp: " << item.getEventTimestamp() << " "
    //               << " Source ID: " << item.getSourceId();
    // }
    // fragment size
    // followed by nsize-1 16bits word, every two (32bits) is one record
    uint32_t n_record = *p++;

    char number[32];
    auto precords = new std::vector<uint32_t>();

    for (int i = 0; i < ((nsize - 1) / 2); ++i)
    {
        // read two words every time
        uint32_t record = 0, temp;
        record = *p++;
        temp = *p++;
        record |= (temp << 16);
        precords->push_back(record);
    }

    uint16_t last_geo;
    uint16_t last_crate_id;
    uint16_t last_n_channels;
    uint64_t last_event_id;
    uint16_t chid;

    uint16_t temp_val[32], temp_ov[32], temp_un[32];

    for (auto it = precords->begin(); it != precords->end(); ++it)
    {
        if (*it != 0xFFFFFFFF)
        {
            VME::V785Hit hit(*it);
            // // printf("%08x\n", *it);
            // printf("GEO: %4d, Type: %4d, Crate_ID: %4d, Channels: %4d, Channel_ID: %4d, Value: %4d, UN: %4d, OV: %4d, EVT_ID: %4d\n",
            //        hit.get_geo(),
            //        hit.get_type(), hit.get_crate_id(),
            //        hit.get_channels(), hit.get_channel_id(),
            //        hit.get_value(), hit.get_threshold_code(), hit.get_overflow_code(), hit.get_event_id());
            if (hit.is_header())
            {
                last_geo = hit.get_geo();
                last_crate_id = hit.get_crate_id();
                last_n_channels = hit.get_channels();
            }
            else if (hit.is_datum())
            {
                chid = hit.get_channel_id();
                temp_val[chid] = hit.get_value();
                temp_ov[chid] = hit.get_overflow_code();
                temp_un[chid] = hit.get_threshold_code();
            }
            else if (hit.is_eob())
            {
                // event_id, geo, crate_id, n_channels, (all values ch0-32), (ov), (un)
                pfragdata->push_back(hit.get_event_id());
                pfragdata->push_back(static_cast<uint64_t>(last_geo));
                pfragdata->push_back(static_cast<uint64_t>(last_crate_id));
                pfragdata->push_back(static_cast<uint64_t>(last_n_channels));
                for(auto &v : temp_val){
                    pfragdata->push_back(static_cast<uint64_t>(v));
                }
                for(auto &v : temp_ov){
                    pfragdata->push_back(static_cast<uint64_t>(v));
                }
                for(auto &v : temp_un){
                    pfragdata->push_back(static_cast<uint64_t>(v));
                }
                ++frag_cnt;
            }
        }
    }
}

/**
 * processEvent
 *  Physics event item (DDAS).
 *
 */
void CRingItemProcessor::processEvent(CPhysicsEventItem &item,
                                      uint64_t &event_id, uint64_t &frag_cnt,
                                      std::vector<FragmentData> *pfragdata,
                                      std::vector<uint16_t> *ptracedata,
                                      int &verbosity)
{
    FragmentIndex indexer(static_cast<uint16_t *>(item.getBodyPointer()));

    int total_fragmts = indexer.getNumberFragments();
    //    std::cout << total_fragmts << std::endl;

    CRingItemFactory factory;
    DAQ::DDAS::DDASHitUnpacker unpacker;

    // Parse fragments into hits
    for (auto &fragment : indexer)
    {
        // data container for fragment
        FragmentData i_fragdata{
            .fragment_id = frag_cnt++,
            .event_id = event_id};

        if (verbosity == 2)
        {
            fprintf(stdout, "Processing evt-id: %u, frgmt: %u, ts: %u, accu frgmts: %u\n",
                    event_id, fragment.s_sourceId, fragment.s_timestamp, frag_cnt);
        }

        CRingItem *frag_raw = factory.createRingItem(reinterpret_cast<uint8_t *>(fragment.s_itemhdr));
        std::unique_ptr<CPhysicsEventItem> frag_item(dynamic_cast<CPhysicsEventItem *>(frag_raw));

        uint32_t *begin = reinterpret_cast<uint32_t *>(frag_item->getBodyPointer());
        uint32_t *sentinel = begin + frag_item->getBodySize() / sizeof(uint32_t);

        DAQ::DDAS::DDASHit hit;
        unpacker.unpack(begin, sentinel, hit);

        double ts = hit.GetTime();
        uint32_t crate_id = static_cast<uint32_t>(hit.GetCrateID());
        uint32_t slot_id = static_cast<uint32_t>(hit.GetSlotID());
        uint32_t channel_id = static_cast<uint32_t>(hit.GetChannelID());
        uint64_t ts_coarse = static_cast<uint64_t>(hit.GetCoarseTime());
        uint32_t energy = static_cast<uint32_t>(hit.GetEnergy());
        uint32_t pileup = static_cast<uint32_t>(hit.GetFinishCode());
        uint32_t overflow_code = static_cast<uint32_t>(hit.GetOverflowCode());
        uint32_t cfd_fail_bit = static_cast<uint32_t>(hit.GetCFDFailBit());
        uint32_t trace_len = static_cast<uint32_t>(hit.GetTraceLength());
        uint32_t frequency = static_cast<uint32_t>(hit.GetModMSPS());
        uint32_t resolution = static_cast<uint32_t>(hit.GetADCResolution());
        uint16_t adc_over_underflow = static_cast<uint16_t>(hit.GetADCOverflowUnderflow());
        i_fragdata.timestamp = ts;
        i_fragdata.crate_id = crate_id;
        i_fragdata.slot_id = slot_id;
        i_fragdata.channel_id = channel_id;
        i_fragdata.coarse_time = ts_coarse;
        i_fragdata.energy = energy;
        i_fragdata.finish_code = pileup;
        i_fragdata.overflow_code = overflow_code;
        i_fragdata.cfd_fail_bit = cfd_fail_bit;
        i_fragdata.trace_length = trace_len;
        i_fragdata.adc_frequency = frequency;
        i_fragdata.adc_resolution = resolution;
        i_fragdata.adc_over_underflow = adc_over_underflow;

        std::vector<uint16_t> &trace = hit.GetTrace();
        uint16_t *ptr = trace.data();
        for (int i = 0; i < trace_len; i++)
        {
            ptracedata->push_back(*ptr++);
        }
        pfragdata->push_back(i_fragdata);
    }
}

/**
 * processEventCount
 *
 *
 */
void CRingItemProcessor::processEventCount(CRingPhysicsEventCountItem &item)
{
}

/**
 * processFormat
 *
 */
void CRingItemProcessor::processFormat(CDataFormatItem &item, RunMetaData &run_meta)
{
    std::stringstream fmtss;
    fmtss << item.getMajor() << "." << item.getMinor();
    strcpy(run_meta.fmt, fmtss.str().c_str());
}

/**
 * processGlomParams
 *
 */
void CRingItemProcessor::processGlomParams(CGlomParameters &item, int &verbosity)
{
    if (verbosity > 2)
    {
        std::cout << "--> " << item.typeName() << ": ";
        if (item.isBuilding())
        {
            std::cout << "Coincidence interval: "
                      << item.coincidenceTicks() << ", ";
            std::cout << "Timestamp policy: " << glomPolicyMap[item.timestampPolicy()]
                      << std::endl;
        }
        else
        {
            std::cout << "operating in passthrough (non-building) mode" << std::endl;
        }
    }
}

/**
 * processUnknownItem
 *
 */
void CRingItemProcessor::processUnknownItem(CRingItem &item, int &verbosity)
{
    if (verbosity > 2)
    {
        std::cout << "Unknown Type: " << item.toString() << std::endl;
    }
}
