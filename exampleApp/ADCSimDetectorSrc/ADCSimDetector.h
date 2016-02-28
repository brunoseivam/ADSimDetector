/* ADCS simulation driver
 * Mark Rivers
 * Feb. 28, 2016
 */

#include <epicsEvent.h>
#include "asynNDArrayDriver.h"

#define SimAcquireString        "SIM_ACQUIRE"
#define SimNumTimePointsString  "SIM_NUM_TIME_POINTS"
#define SimAmplitudeString      "SIM_AMPLITUDE"
#define SimFrequencyString      "SIM_FREQUENCY"
#define SimPhaseString          "SIM_PHASE"
#define SimNoiseString          "SIM_NOISE"

#define MAX_SIGNALS 8

/** ADC simulation driver; does 1-D waveforms on 8 channels. 
  * Inherits from asynNDArrayDriver */
class epicsShareClass ADCSimDetector : public asynNDArrayDriver {
public:
    ADCSimDetector(const char *portName, int numTimePoints, NDDataType_t dataType,
                   int maxBuffers, size_t maxMemory,
                   int priority, int stackSize);

    /* These are the methods that we override from asynNDArrayDriver */
    virtual asynStatus writeInt32(asynUser *pasynUser, epicsInt32 value);
    virtual void report(FILE *fp, int details);
    void simTask(); /**< Should be private, but gets called from C, so must be public */

protected:
    int P_Acquire;
    #define FIRST_SIM_DETECTOR_PARAM P_Acquire
    int P_NumTimePoints;
    int P_Amplitude;
    int P_Frequency;
    int P_Phase;
    int P_Noise;
    #define LAST_SIM_DETECTOR_PARAM P_Noise


private:
    /* These are the methods that are new to this class */
    template <typename epicsType> void computeArraysT();
    void computeArrays();

    /* Our data */
    epicsEventId startEventId_;
    epicsEventId stopEventId_;
    int arrayCounter_;
    int currentPoint_;
};


#define NUM_SIM_DETECTOR_PARAMS ((int)(&LAST_SIM_DETECTOR_PARAM - &FIRST_SIM_DETECTOR_PARAM + 1))

