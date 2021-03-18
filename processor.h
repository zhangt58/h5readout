#pragma once

#include <vector>
#include <map>

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
    virtual void processScalerItem(CRingScalerItem& item);
    virtual void processStateChangeItem(CRingStateChangeItem& item,
                                        RunMetaData& run_meta);
    virtual void processTextItem(CRingTextItem& item);
    virtual void processEvent(CPhysicsEventItem& item,
                              uint64_t& event_id, uint64_t& frag_cnt,
                              std::vector<FragmentData> *pfragdata,
                              std::vector<uint16_t> *ptracedata,
                              bool& verbose);
    virtual void processEventCount(CRingPhysicsEventCountItem& item);
    virtual void processFormat(CDataFormatItem& item, RunMetaData& run_meta);
    virtual void processGlomParams(CGlomParameters& item);
    virtual void processUnknownItem(CRingItem& item);
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
                     bool& verbose);
