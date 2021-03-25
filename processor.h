#pragma once

#include <vector>
#include <map>
#include <unordered_map>
#include <set>

#include <CDataFormatItem.h>
#include <CRingItem.h>
#include <CRingScalerItem.h>
#include <CRingTextItem.h>
#include <CRingStateChangeItem.h>
#include <CRingPhysicsEventCountItem.h>
#include <CPhysicsEventItem.h>
#include <CGlomParameters.h>

#include "h5readout.h"


class CRingItemProcessor {
public:
    virtual void processScalerItem(CRingScalerItem& item,
                                   std::vector<time_t> *pscalerts,      // scaler timestamp
                                   std::vector<uint32_t> *pscalerlen,   // scaler length
                                   std::vector<uint32_t> *pscalerdata,  // scaler array
                                   int& verbosity);
    virtual void processStateChangeItem(CRingStateChangeItem& item,
                                        RunMetaData& run_meta,          // meta data
                                        uint16_t& item_type);
    virtual void processTextItem(CRingTextItem& item, int& verbosity);
    virtual void processEvent(CPhysicsEventItem& item,
                              uint64_t& event_id,                       // physics event id
                              uint64_t& frag_cnt,                       // accumulated fragments
                              std::vector<FragmentData> *pfragdata,     // all FragmentData
                              std::vector<uint16_t> *ptracedata,        // all tracedata
                              int& verbosity);
    virtual void processEventCount(CRingPhysicsEventCountItem& item);
    virtual void processFormat(CDataFormatItem& item, RunMetaData& run_meta);
    virtual void processGlomParams(CGlomParameters& item, int& verbosity);
    virtual void processUnknownItem(CRingItem& item, int& verbosity);
};

static std::map<CGlomParameters::TimestampPolicy, std::string> glomPolicyMap = {
    {CGlomParameters::first, "first"},
    {CGlomParameters::last, "last"},
    {CGlomParameters::average, "average"}
};


void processRingItem(CRingItemProcessor& processor, CRingItem* item,
                     RunMetaData& run_meta,
                     uint64_t& event_id, uint64_t& frag_cnt,
                     std::vector<FragmentData> *pfragdata,
                     std::vector<uint16_t> *ptracedata,
                     std::vector<time_t> *pscalerts,
                     std::vector<uint32_t> *pscalerlen,
                     std::vector<uint32_t> *pscalerdata,
                     int& verbosity);


static std::unordered_map<std::string, uint16_t> const ITEM_TYPE = {
    {"PHYSICS_EVENT", PHYSICS_EVENT},
    {"PERIODIC_SCALERS", PERIODIC_SCALERS},
    {"EVB_GLOM_INFO", EVB_GLOM_INFO}
};

// supported VME modules
static std::set<std::string> const VME_MODULES = {
    "V785"
};