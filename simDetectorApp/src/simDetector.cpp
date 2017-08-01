/* simDetector.cpp
 *
 * This is a driver for a simulated area detector.
 *
 * Author: Mark Rivers
 *         University of Chicago
 *
 * Created:  March 20, 2008
 *
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <epicsTime.h>
#include <epicsThread.h>
#include <epicsEvent.h>
#include <epicsMutex.h>
#include <epicsString.h>
#include <epicsStdio.h>
#include <epicsMutex.h>
#include <cantProceed.h>
#include <iocsh.h>

#include "ADDriver.h"
#include <epicsExport.h>
#include "simDetector.h"

#include <pv/sharedVector.h>

// TODO: remove this
#include <pv/ntndarrayAttribute.h>

using namespace std;
using namespace epics::pvData;

static const char *driverName = "simDetector";

/* Some systems don't define M_PI in math.h */
#ifndef M_PI
  #define M_PI 3.14159265358979323846
#endif

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeArray(int sizeX, int sizeY)
{
    int simMode=0;
    int status = asynSuccess;

    status = getIntegerParam (SimMode, &simMode);
    switch(simMode) {
        case SimModeLinearRamp:
            status = computeLinearRampArray<epicsType>(sizeX, sizeY);
            break;
        case SimModePeaks:
            status = computePeaksArray<epicsType>(sizeX, sizeY);
            break;
        case SimModeSine:
            status = computeSineArray<epicsType>(sizeX, sizeY);
            break;
    }
    return status;
}

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeLinearRampArray(int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();
    epicsType *pMono=NULL, *pRed=NULL, *pGreen=NULL, *pBlue=NULL;
    int columnStep=0, rowStep=0, colorMode;
    epicsType incMono, incRed, incGreen, incBlue;
    int status = asynSuccess;
    double exposureTime, gain, gainX, gainY, gainRed, gainGreen, gainBlue;
    int resetImage;
    int i, j;

    status = getDoubleParam (ADGain,        &gain);
    status = getDoubleParam (SimGainX,      &gainX);
    status = getDoubleParam (SimGainY,      &gainY);
    status = getDoubleParam (SimGainRed,    &gainRed);
    status = getDoubleParam (SimGainGreen,  &gainGreen);
    status = getDoubleParam (SimGainBlue,   &gainBlue);
    status = getIntegerParam(SimResetImage, &resetImage);
    status = getIntegerParam(NDColorMode,   &colorMode);
    status = getDoubleParam (ADAcquireTime, &exposureTime);

    /* The intensity at each pixel[i,j] is:
     * (i * gainX + j* gainY) + imageCounter * gain * exposureTime * 1000. */
    incMono  = (epicsType) (gain      * exposureTime * 1000.);
    incRed   = (epicsType) gainRed   * incMono;
    incGreen = (epicsType) gainGreen * incMono;
    incBlue  = (epicsType) gainBlue  * incMono;

    switch (colorMode) {
        case NDColorModeMono:
            pMono = (epicsType *)rawData;
            break;
        case NDColorModeRGB1:
            columnStep = 3;
            rowStep = 0;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData+1;
            pBlue  = (epicsType *)rawData+2;
            break;
        case NDColorModeRGB2:
            columnStep = 1;
            rowStep = 2 * sizeX;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData + sizeX;
            pBlue  = (epicsType *)rawData + 2*sizeX;
            break;
        case NDColorModeRGB3:
            columnStep = 1;
            rowStep = 0;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData + sizeX*sizeY;
            pBlue  = (epicsType *)rawData + 2*sizeX*sizeY;
            break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);

    if (resetImage) {
        for (i=0; i<sizeY; i++) {
            switch (colorMode) {
                case NDColorModeMono:
                    for (j=0; j<sizeX; j++) {
                        (*pMono++) = (epicsType) (incMono * (gainX*j + gainY*i));
                    }
                    break;
                case NDColorModeRGB1:
                case NDColorModeRGB2:
                case NDColorModeRGB3:
                    for (j=0; j<sizeX; j++) {
                        *pRed   = (epicsType) (incRed   * (gainX*j + gainY*i));
                        *pGreen = (epicsType) (incGreen * (gainX*j + gainY*i));
                        *pBlue  = (epicsType) (incBlue  * (gainX*j + gainY*i));
                        pRed   += columnStep;
                        pGreen += columnStep;
                        pBlue  += columnStep;
                    }
                    pRed   += rowStep;
                    pGreen += rowStep;
                    pBlue  += rowStep;
                    break;
            }
        }
    } else {
        for (i=0; i<sizeY; i++) {
            switch (colorMode) {
                case NDColorModeMono:
                    for (j=0; j<sizeX; j++) {
                            *pMono++ += incMono;
                    }
                    break;
                case NDColorModeRGB1:
                case NDColorModeRGB2:
                case NDColorModeRGB3:
                    for (j=0; j<sizeX; j++) {
                        *pRed   += incRed;
                        *pGreen += incGreen;
                        *pBlue  += incBlue;
                        pRed   += columnStep;
                        pGreen += columnStep;
                        pBlue  += columnStep;
                    }
                    pRed   += rowStep;
                    pGreen += rowStep;
                    pBlue  += rowStep;
                    break;
            }
        }
    }

    pRaw_->setData<epicsType>(data);
    return(status);
}

/** Compute array for array of peaks */
template <typename epicsType> int simDetector::computePeaksArray(int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();
    epicsType *pMono=NULL, *pRed=NULL;
    epicsType *pMono2=NULL, *pRed2=NULL, *pGreen2=NULL, *pBlue2=NULL;
    int columnStep=0, colorMode;
    int peaksStartX, peaksStartY, peaksStepX, peaksStepY;
    int peaksNumX, peaksNumY, peaksWidthX, peaksWidthY;
    int status = asynSuccess;
    int i,j,k,l;
    int minX, maxX, minY,maxY;
    int offsetX, offsetY;
    int peakVariation, noisePct;
    double gainVariation, noise;
    double gain, gainX, gainY, gainRed, gainGreen, gainBlue;
    double gaussX, gaussY;
    double tmpValue;

    status = getIntegerParam(NDColorMode,   &colorMode);
    status = getDoubleParam (ADGain,        &gain);
    status = getDoubleParam (SimGainX,      &gainX);
    status = getDoubleParam (SimGainY,      &gainY);
    status = getDoubleParam (SimGainRed,    &gainRed);
    status = getDoubleParam (SimGainGreen,  &gainGreen);
    status = getDoubleParam (SimGainBlue,   &gainBlue);
    status = getIntegerParam (SimPeakStartX,  &peaksStartX);
    status = getIntegerParam (SimPeakStartY,  &peaksStartY);
    status = getIntegerParam (SimPeakStepX,  &peaksStepX);
    status = getIntegerParam (SimPeakStepY,  &peaksStepY);
    status = getIntegerParam (SimPeakNumX,  &peaksNumX);
    status = getIntegerParam (SimPeakNumY,  &peaksNumY);
    status = getIntegerParam (SimPeakWidthX,  &peaksWidthX);
    status = getIntegerParam (SimPeakWidthY,  &peaksWidthY);
    status = getIntegerParam (SimPeakHeightVariation,  &peakVariation);
    status = getIntegerParam (SimNoise,  &noisePct);

       switch (colorMode) {
        case NDColorModeMono:
            pMono = (epicsType *)rawData;
            break;
        case NDColorModeRGB1:
            columnStep = 3;
            pRed   = (epicsType *)rawData;
            break;
        case NDColorModeRGB2:
            columnStep = 1;
            pRed   = (epicsType *)rawData;
            break;
        case NDColorModeRGB3:
            columnStep = 1;
            pRed   = (epicsType *)rawData;
            break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);
    switch (colorMode) {
        case NDColorModeMono:
            // Clear the Image
            pMono2 = pMono;
            for (i = 0; i<sizeY; i++) {
                for (j = 0; j<sizeX; j++) {
                    (*pMono2++) = (epicsType)0;
                }
            }
            for (i = 0; i<peaksNumY; i++) {
                for (j = 0; j<peaksNumX; j++) {
                    gaussX = 0;
                    gaussY = 0;
                    if (peakVariation !=0) {
                        gainVariation = 1.0 + (rand()%peakVariation+1)/100.0;
                    }
                    else{
                        gainVariation = 1.0;
                    }
                    offsetY = i * peaksStepY + peaksStartY;
                    offsetX = j * peaksStepX + peaksStartX;
                    minX = (offsetX>4*peaksWidthX) ?(offsetX -4*peaksWidthX):0;
                    maxX = (offsetX+4*peaksWidthX<sizeX) ?(offsetX + 4*peaksWidthX):sizeX;
                    minY = (offsetY>4*peaksWidthY) ?(offsetY -4*peaksWidthY):0;
                    maxY = (offsetY+4*peaksWidthY<sizeY) ?(offsetY + 4*peaksWidthY):sizeY;
                    for (k =minY; k<maxY; k++) {
                        pMono2 = pMono + (minX + k*sizeX);
                        for (l=minX; l<maxX; l++) {
                            if (noisePct !=0) {
                                noise = 1.0 + (rand()%noisePct+1)/100.0;
                            }
                            else {
                                noise = 1.0;
                            }
                            gaussY = gainY * exp( -pow((double)(k-offsetY)/(double)peaksWidthY,2.0)/2.0 );
                            gaussX = gainX * exp( -pow((double)(l-offsetX)/(double)peaksWidthX,2.0)/2.0 );
                            tmpValue =  gainVariation*gain * gaussX * gaussY*noise;
                            (*pMono2) += (epicsType)tmpValue;
                            pMono2++;
                        }
                    }
                }
            }
            break;
        case NDColorModeRGB1:
        case NDColorModeRGB2:
        case NDColorModeRGB3:
            // Clear the Image
            pRed2 = pRed;
            for (i = 0; i<sizeY; i++) {
                for (j = 0; j<sizeX; j++) {
                    (*pRed2++) = (epicsType)0;  //Since we are just clearing the field we will do this with one pointer
                    (*pRed2++) = (epicsType)0;
                    (*pRed2++) = (epicsType)0;
                }
            }
            for (i = 0; i<peaksNumY; i++) {
                for (j = 0; j<peaksNumX; j++) {
                    if (peakVariation !=0) {
                        gainVariation = 1.0 + (rand()%peakVariation+1)/100.0;
                    }
                    else{
                        gainVariation = 1.0;
                    }
                    offsetY = i * peaksStepY + peaksStartY;
                    offsetX = j * peaksStepX + peaksStartX;
                    minX = (offsetX>4*peaksWidthX) ?(offsetX -4*peaksWidthX):0;
                    maxX = (offsetX+4*peaksWidthX<sizeX) ?(offsetX + 4*peaksWidthX):sizeX;
                    minY = (offsetY>4*peaksWidthY) ?(offsetY -4*peaksWidthY):0;
                    maxY = (offsetY+4*peaksWidthY<sizeY) ?(offsetY + 4*peaksWidthY):sizeY;
                    for (k =minY; k<maxY; k++) {
                        //Move to the starting point for this peak
                        switch (colorMode) {
                            case NDColorModeRGB1:
                                pRed2 = pRed + (minX*columnStep + k*sizeX*columnStep);
                                pGreen2 = pRed2 + 1;
                                pBlue2 = pRed2 + 2;
                                break;
                            case NDColorModeRGB2:
                                pRed2 = pRed + (minX*columnStep + k*3*sizeX*columnStep);
                                pGreen2 = pRed2 + sizeX;
                                pBlue2 = pRed2 + 2*sizeX;
                                break;
                            case NDColorModeRGB3:
                                pRed2 = pRed + (minX*columnStep + k*sizeX*columnStep);
                                pGreen2 = pRed2 + sizeX*sizeY;
                                pBlue2 = pRed2 + 2*sizeX*sizeY;
                                break;
                        }
                        //Fill in a row for this peak
                        for (l=minX; l<maxX; l++) {
                            if (noisePct !=0) {
                                noise = 1.0 + (rand()%noisePct+1)/100.0;
                            }
                            else {
                                noise = 1.0;
                            }
                            gaussY = gainY * exp( -pow((double)(k-offsetY)/(double)peaksWidthY,2.0)/2.0 );
                            gaussX = gainX * exp( -pow((double)(l-offsetX)/(double)peaksWidthX,2.0)/2.0 );
                            tmpValue =  gainVariation*gain * gaussX * gaussY*noise;
                            (*pRed2) += (epicsType)(gainRed*tmpValue);
                            (*pGreen2) += (epicsType)(gainGreen*tmpValue);
                            (*pBlue2) += (epicsType)(gainBlue*tmpValue);

                            pRed2 += columnStep;
                            pGreen2 += columnStep;
                            pBlue2 += columnStep;
                        }
                    }
                }
            }



            break;

    }
    pRaw_->setData<epicsType>(data);
    return status;
}

/** Template function to compute the simulated detector data for any data type */
template <typename epicsType> int simDetector::computeSineArray(int sizeX, int sizeY)
{
    shared_vector<epicsType> data(pRaw_->getData<epicsType>());
    epicsType *rawData = data.data();
    epicsType *pMono=NULL, *pRed=NULL, *pGreen=NULL, *pBlue=NULL;
    int columnStep=0, rowStep=0, colorMode;
    int status = asynSuccess;
    int xSineOperation, ySineOperation;   
    double exposureTime, gain, gainX, gainY, gainRed, gainGreen, gainBlue;
    double sineOffset, sineNoise;
    double xSine1Amplitude, xSine1Frequency, xSine1Phase;
    double xSine2Amplitude, xSine2Frequency, xSine2Phase;
    double ySine1Amplitude, ySine1Frequency, ySine1Phase;
    double ySine2Amplitude, ySine2Frequency, ySine2Phase;
    double rndm;
    double xTime, yTime;
    int resetImage;
    int i, j;

    status = getDoubleParam (ADGain,            &gain);
    status = getDoubleParam (SimGainX,          &gainX);
    status = getDoubleParam (SimGainY,          &gainY);
    status = getDoubleParam (SimGainRed,        &gainRed);
    status = getDoubleParam (SimGainGreen,      &gainGreen);
    status = getDoubleParam (SimGainBlue,       &gainBlue);
    status = getIntegerParam(SimResetImage,     &resetImage);
    status = getIntegerParam(NDColorMode,       &colorMode);
    status = getDoubleParam (ADAcquireTime,     &exposureTime);
    status = getDoubleParam(SimSineOffset,      &sineOffset);
    status = getDoubleParam(SimSineNoise,       &sineNoise);
    status = getIntegerParam(SimXSineOperation, &xSineOperation);
    status = getDoubleParam(SimXSine1Amplitude, &xSine1Amplitude);
    status = getDoubleParam(SimXSine1Frequency, &xSine1Frequency);
    status = getDoubleParam(SimXSine1Phase,     &xSine1Phase);
    status = getDoubleParam(SimXSine2Amplitude, &xSine2Amplitude);
    status = getDoubleParam(SimXSine2Frequency, &xSine2Frequency);
    status = getDoubleParam(SimXSine2Phase,     &xSine2Phase);
    status = getIntegerParam(SimYSineOperation, &ySineOperation);
    status = getDoubleParam(SimYSine1Amplitude, &ySine1Amplitude);
    status = getDoubleParam(SimYSine1Frequency, &ySine1Frequency);
    status = getDoubleParam(SimYSine1Phase,     &ySine1Phase);
    status = getDoubleParam(SimYSine2Amplitude, &ySine2Amplitude);
    status = getDoubleParam(SimYSine2Frequency, &ySine2Frequency);
    status = getDoubleParam(SimYSine2Phase,     &ySine2Phase);

    switch (colorMode) {
        case NDColorModeMono:
            pMono = (epicsType *)rawData;
            break;
        case NDColorModeRGB1:
            columnStep = 3;
            rowStep = 0;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData+1;
            pBlue  = (epicsType *)rawData+2;
            break;
        case NDColorModeRGB2:
            columnStep = 1;
            rowStep = 2 * sizeX;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData + sizeX;
            pBlue  = (epicsType *)rawData + 2*sizeX;
            break;
        case NDColorModeRGB3:
            columnStep = 1;
            rowStep = 0;
            pRed   = (epicsType *)rawData;
            pGreen = (epicsType *)rawData + sizeX*sizeY;
            pBlue  = (epicsType *)rawData + 2*sizeX*sizeY;
            break;
    }
    //pRaw_->pAttributeList->add("ColorMode", "Color mode", NDAttrInt32, &colorMode);

    if (resetImage) {
      if (xSine1_) free(xSine1_);
      if (xSine2_) free(xSine2_);
      if (ySine1_) free(ySine1_);
      if (ySine2_) free(ySine2_);
      xSine1_ = (double *)calloc(sizeX, sizeof(double));
      xSine2_ = (double *)calloc(sizeX, sizeof(double));
      ySine1_ = (double *)calloc(sizeY, sizeof(double));
      ySine2_ = (double *)calloc(sizeY, sizeof(double));
      xSineCounter_ = 0;
      ySineCounter_ = 0;
    } 
    
    for (i=0; i<sizeX; i++) {
        xTime = xSineCounter_++ * gainX / sizeX;
        xSine1_[i] = xSine1Amplitude * sin((xTime  * xSine1Frequency + xSine1Phase/360.) * 2. * M_PI);
        xSine2_[i] = xSine2Amplitude * sin((xTime  * xSine2Frequency + xSine2Phase/360.) * 2. * M_PI);
    }
    for (i=0; i<sizeY; i++) {
        yTime = ySineCounter_++ * gainY / sizeY;
        ySine1_[i] = ySine1Amplitude * sin((yTime  * ySine1Frequency + ySine1Phase/360.) * 2. * M_PI);
        ySine2_[i] = ySine2Amplitude * sin((yTime  * ySine2Frequency + ySine2Phase/360.) * 2. * M_PI);
    }                             
    
    if (colorMode == NDColorModeMono) {
        if (xSineOperation == SimSineOperationAdd) {
            for (i=0; i<sizeX; i++) {
                xSine1_[i] = xSine1_[i] + xSine2_[i];
            }
        }
        else {
            for (i=0; i<sizeX; i++) {
                xSine1_[i] = xSine1_[i] * xSine2_[i];
            }
        }
        if (ySineOperation == SimSineOperationAdd) {
            for (i=0; i<sizeY; i++) {
                ySine1_[i] = ySine1_[i] + ySine2_[i];
            }
        }
        else {
            for (i=0; i<sizeY; i++) {
                ySine1_[i] = ySine1_[i] * ySine2_[i];
            }
        }
    }
    for (i=0; i<sizeY; i++) {
        switch (colorMode) {
            case NDColorModeMono:
                for (j=0; j<sizeX; j++) {
                    rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                    *pMono++ = (epicsType) (gain * (sineOffset + sineNoise*rndm + ySine1_[i] + xSine1_[j]));
                }
                break;
            case NDColorModeRGB1:
            case NDColorModeRGB2:
            case NDColorModeRGB3:
                for (j=0; j<sizeX; j++) {
                    rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                    *pRed   = (epicsType)(gain * gainRed   * (sineOffset + sineNoise*rndm + xSine1_[j]));
                    rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                    *pGreen = (epicsType)(gain * gainGreen * (sineOffset + sineNoise*rndm + ySine1_[i]));
                    rndm = 2.*(rand()/(double)RAND_MAX - 0.5);
                    *pBlue  = (epicsType)(gain * gainBlue  * (sineOffset + sineNoise*rndm + (xSine2_[j] + ySine2_[i])/2.));
                    pRed   += columnStep;
                    pGreen += columnStep;
                    pBlue  += columnStep;
                }
                pRed   += rowStep;
                pGreen += rowStep;
                pBlue  += rowStep;
                break;
        }
    }
    pRaw_->setData<epicsType>(data);
    return(status);
}

/** Controls the shutter */
void simDetector::setShutter(int open)
{
    int shutterMode;

    getIntegerParam(ADShutterMode, &shutterMode);
    if (shutterMode == ADShutterModeDetector) {
        /* Simulate a shutter by just changing the status readback */
        setIntegerParam(ADShutterStatus, open);
    } else {
        /* For no shutter or EPICS shutter call the base class method */
        ADDriver::setShutter(open);
    }
}

/** Computes the new image data */
NDArrayPtr simDetector::computeImage()
{
    ScalarType dataType;
    int itemp;
    int binX, binY, minX, minY, sizeX, sizeY, reverseX, reverseY;
    int xDim=0, yDim=1, colorDim=-1;
    int resetImage;
    int maxSizeX, maxSizeY;
    int colorMode;
    const char* functionName = "computeImage";

    /* NOTE: The caller of this function must have taken the mutex */

    int status;
    status  = getIntegerParam(ADBinX,         &binX);
    status |= getIntegerParam(ADBinY,         &binY);
    status |= getIntegerParam(ADMinX,         &minX);
    status |= getIntegerParam(ADMinY,         &minY);
    status |= getIntegerParam(ADSizeX,        &sizeX);
    status |= getIntegerParam(ADSizeY,        &sizeY);
    status |= getIntegerParam(ADReverseX,     &reverseX);
    status |= getIntegerParam(ADReverseY,     &reverseY);
    status |= getIntegerParam(ADMaxSizeX,     &maxSizeX);
    status |= getIntegerParam(ADMaxSizeY,     &maxSizeY);
    status |= getIntegerParam(NDColorMode,    &colorMode);
    status |= getIntegerParam(NDDataType,     &itemp); dataType = (ScalarType)itemp;
    status |= getIntegerParam(SimResetImage,  &resetImage);
    if (status) asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
                    "%s:%s: error getting parameters\n",
                    driverName, functionName);

    /* Make sure parameters are consistent, fix them if they are not */
    if (binX < 1) {
        binX = 1;
        status |= setIntegerParam(ADBinX, binX);
    }
    if (binY < 1) {
        binY = 1;
        status |= setIntegerParam(ADBinY, binY);
    }
    if (minX < 0) {
        minX = 0;
        status |= setIntegerParam(ADMinX, minX);
    }
    if (minY < 0) {
        minY = 0;
        status |= setIntegerParam(ADMinY, minY);
    }
    if (minX > maxSizeX-1) {
        minX = maxSizeX-1;
        status |= setIntegerParam(ADMinX, minX);
    }
    if (minY > maxSizeY-1) {
        minY = maxSizeY-1;
        status |= setIntegerParam(ADMinY, minY);
    }
    if (minX+sizeX > maxSizeX) {
        sizeX = maxSizeX-minX;
        status |= setIntegerParam(ADSizeX, sizeX);
    }
    if (minY+sizeY > maxSizeY) {
        sizeY = maxSizeY-minY;
        status |= setIntegerParam(ADSizeY, sizeY);
    }

    size_t ndims = 0;
    switch (colorMode) {
        case NDColorModeMono:
            ndims = 2;
            xDim = 0;
            yDim = 1;
            break;
        case NDColorModeRGB1:
            ndims = 3;
            colorDim = 0;
            xDim     = 1;
            yDim     = 2;
            break;
        case NDColorModeRGB2:
            ndims = 3;
            colorDim = 1;
            xDim     = 0;
            yDim     = 2;
            break;
        case NDColorModeRGB3:
            ndims = 3;
            colorDim = 2;
            xDim     = 0;
            yDim     = 1;
            break;
    }

    size_t dims[3];
    if (resetImage) {
        /* Allocate the raw buffer we use to compute images. */
        dims[xDim] = maxSizeX;
        dims[yDim] = maxSizeY;
        if (ndims > 2) dims[colorDim] = 3;
        pRaw_ = NDArrayPool_.alloc(ndims, dims, dataType);

        if (!pRaw_) {
            asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
                      "%s:%s: error allocating raw buffer\n",
                      driverName, functionName);
            return pRaw_;
        }
    }

    switch (dataType) {
        case pvByte:
            status |= computeArray<int8>(maxSizeX, maxSizeY);
            break;
        case pvUByte:
            status |= computeArray<uint8>(maxSizeX, maxSizeY);
            break;
        case pvShort:
            status |= computeArray<int16>(maxSizeX, maxSizeY);
            break;
        case pvUShort:
            status |= computeArray<uint16>(maxSizeX, maxSizeY);
            break;
        case pvInt:
            status |= computeArray<int32>(maxSizeX, maxSizeY);
            break;
        case pvUInt:
            status |= computeArray<uint32>(maxSizeX, maxSizeY);
            break;
        case pvFloat:
            status |= computeArray<float>(maxSizeX, maxSizeY);
            break;
        case pvDouble:
            status |= computeArray<double>(maxSizeX, maxSizeY);
            break;
        default:
            break;
    }

    /* Extract the region of interest with binning.
     * If the entire image is being used (no ROI or binning) that's OK because
     * convertImage detects that case and is very efficient */
    NDDimension_t dimsOut[3];
    NDArray::initDimension(dimsOut[xDim], sizeX, minX, binX, reverseX);
    NDArray::initDimension(dimsOut[yDim], sizeY, minY, binY, reverseY);
    if(ndims > 2) NDArray::initDimension(dimsOut[colorDim], 3);

    vector<NDDimension_t> dimsOutVector;
    for(size_t i = 0; i < ndims; ++i)
        dimsOutVector.push_back(dimsOut[i]);

    NDArrayPtr result(NDArrayPool_.convert(pRaw_, dataType, dimsOutVector));
    if(!result) {
        asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
                    "%s:%s: error allocating buffer in convert()\n",
                    driverName, functionName);
        return result;
    }

    NDArrayInfo_t arrayInfo = result->getInfo();
    result->setCompressedSize(arrayInfo.totalBytes);
    result->setUncompressedSize(arrayInfo.totalBytes);

    status |= setIntegerParam(NDArraySize,  (int)arrayInfo.totalBytes);
    status |= setIntegerParam(NDArraySizeX, (int)arrayInfo.x.size);
    status |= setIntegerParam(NDArraySizeY, (int)arrayInfo.y.size);
    status |= setIntegerParam(SimResetImage, 0);
    if (status) asynPrint(this->pasynUserSelf, ASYN_TRACE_ERROR,
                    "%s:%s: error setting parameters\n",
                    driverName, functionName);
    return result;
}

static void simTaskC(void *drvPvt)
{
    simDetector *pPvt = (simDetector *)drvPvt;

    pPvt->simTask();
}

/** This thread calls computeImage to compute new image data and does the callbacks to send it to higher layers.
  * It implements the logic for single, multiple or continuous acquisition. */
void simDetector::simTask()
{
    int status = asynSuccess;
    int imageCounter;
    int numImages, numImagesCounter;
    int imageMode;
    int arrayCallbacks;
    int acquire=0;
    double acquireTime, acquirePeriod, delay;
    epicsTimeStamp startTime, endTime;
    double elapsedTime;
    const char *functionName = "simTask";

    this->lock();
    /* Loop forever */
    while (1) {
        /* If we are not acquiring then wait for a semaphore that is given when acquisition is started */
        if (!acquire) {
            /* Release the lock while we wait for an event that says acquire has started, then lock again */
            asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                "%s:%s: waiting for acquire to start\n", driverName, functionName);
            this->unlock();
            status = epicsEventWait(startEventId_);
            this->lock();
            acquire = 1;
            setStringParam(ADStatusMessage, "Acquiring data");
            setIntegerParam(ADNumImagesCounter, 0);
        }

        /* We are acquiring. */
        /* Get the current time */
        epicsTimeGetCurrent(&startTime);
        getIntegerParam(ADImageMode, &imageMode);

        /* Get the exposure parameters */
        getDoubleParam(ADAcquireTime, &acquireTime);
        getDoubleParam(ADAcquirePeriod, &acquirePeriod);

        setIntegerParam(ADStatus, ADStatusAcquire);

        /* Open the shutter */
        setShutter(ADShutterOpen);

        /* Call the callbacks to update any changes */
        callParamCallbacks();

        /* Simulate being busy during the exposure time.  Use epicsEventWaitWithTimeout so that
         * manually stopping the acquisition will work */

        if (acquireTime > 0.0) {
            this->unlock();
            status = epicsEventWaitWithTimeout(stopEventId_, acquireTime);
            this->lock();
        } else {
            status = epicsEventTryWait(stopEventId_);
        }        
        if (status == epicsEventWaitOK) {
            acquire = 0;
            if (imageMode == ADImageContinuous) {
              setIntegerParam(ADStatus, ADStatusIdle);
            } else {
              setIntegerParam(ADStatus, ADStatusAborted);
            }
            callParamCallbacks();
        }

        /* Update the image */
        NDArrayPtr pImage(computeImage());
        if (!pImage) continue;

        /* Close the shutter */
        setShutter(ADShutterClosed);
        
        if (!acquire) continue;

        setIntegerParam(ADStatus, ADStatusReadout);
        /* Call the callbacks to update any changes */
        callParamCallbacks();

        /* Get the current parameters */
        getIntegerParam(NDArrayCounter, &imageCounter);
        getIntegerParam(ADNumImages, &numImages);
        getIntegerParam(ADNumImagesCounter, &numImagesCounter);
        getIntegerParam(NDArrayCallbacks, &arrayCallbacks);
        imageCounter++;
        numImagesCounter++;
        setIntegerParam(NDArrayCounter, imageCounter);
        setIntegerParam(ADNumImagesCounter, numImagesCounter);

        /* Put the frame number and time stamp into the buffer */
        epicsTimeStamp epicsTS;
        updateTimeStamp(&epicsTS);

        pImage->setUniqueId(imageCounter);
        pImage->setTimeStamp(TimeStamp(startTime.secPastEpoch, startTime.nsec));
        pImage->setEpicsTimeStamp(TimeStamp(epicsTS.secPastEpoch, epicsTS.nsec));

        /* Get any attributes that have been defined for this driver */
        NDAttributeListPtr pImageAttrs(pImage->getAttributeList());
        this->getAttributes(pImageAttrs);

        if (arrayCallbacks) {
            /* Publish the NDArray */
            this->unlock();
            NDArrayData->put(pImage);
            this->lock();
            asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: putting to NDArrayData\n", driverName, functionName);
        }

        /* See if acquisition is done */
        if ((imageMode == ADImageSingle) ||
            ((imageMode == ADImageMultiple) &&
             (numImagesCounter >= numImages))) {

          /* First do callback on ADStatus. */
          setStringParam(ADStatusMessage, "Waiting for acquisition");
          setIntegerParam(ADStatus, ADStatusIdle);
          callParamCallbacks();

          acquire = 0;
          setIntegerParam(ADAcquire, acquire);
          asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: acquisition completed\n", driverName, functionName);
        }

        /* Call the callbacks to update any changes */
        callParamCallbacks();

        /* If we are acquiring then sleep for the acquire period minus elapsed time. */
        if (acquire) {
          epicsTimeGetCurrent(&endTime);
          elapsedTime = epicsTimeDiffInSeconds(&endTime, &startTime);
          delay = acquirePeriod - elapsedTime;
          asynPrint(this->pasynUserSelf, ASYN_TRACE_FLOW,
                    "%s:%s: delay=%f\n",
                    driverName, functionName, delay);
          if (delay >= 0.0) {
            /* We set the status to waiting to indicate we are in the period delay */
            setIntegerParam(ADStatus, ADStatusWaiting);
            callParamCallbacks();
            this->unlock();
            status = epicsEventWaitWithTimeout(stopEventId_, delay);
            this->lock();
            if (status == epicsEventWaitOK) {
              acquire = 0;
              if (imageMode == ADImageContinuous) {
                setIntegerParam(ADStatus, ADStatusIdle);
              } else {
                setIntegerParam(ADStatus, ADStatusAborted);
              }
              callParamCallbacks();
            }
          }
        }
    }
}

/** Called when asyn clients call pasynInt32->write().
  * This function performs actions for some parameters, including ADAcquire, ADColorMode, etc.
  * For all parameters it sets the value in the parameter library and calls any registered callbacks..
  * \param[in] pasynUser pasynUser structure that encodes the reason and address.
  * \param[in] value Value to write. */
asynStatus simDetector::writeInt32(asynUser *pasynUser, epicsInt32 value)
{
    int function = pasynUser->reason;
    int adstatus;
    int acquiring;
    int imageMode;
    asynStatus status = asynSuccess;

    /* Ensure that ADStatus is set correctly before we set ADAcquire.*/
    getIntegerParam(ADStatus, &adstatus);
    getIntegerParam(ADAcquire, &acquiring);
    getIntegerParam(ADImageMode, &imageMode);
    if (function == ADAcquire) {
      if (value && !acquiring) {
        setStringParam(ADStatusMessage, "Acquiring data");
      }
      if (!value && acquiring) {
        setStringParam(ADStatusMessage, "Acquisition stopped");
        if (imageMode == ADImageContinuous) {
          setIntegerParam(ADStatus, ADStatusIdle);
        } else {
          setIntegerParam(ADStatus, ADStatusAborted);
        }
        setIntegerParam(ADStatus, ADStatusAcquire); 
      }
    }
    callParamCallbacks();
 
    /* Set the parameter and readback in the parameter library.  This may be overwritten when we read back the
     * status at the end, but that's OK */
    status = setIntegerParam(function, value);

    /* For a real detector this is where the parameter is sent to the hardware */
    if (function == ADAcquire) {
        if (value && !acquiring) {
            /* Send an event to wake up the simulation task.
             * It won't actually start generating new images until we release the lock below */
            epicsEventSignal(startEventId_); 
        }
        if (!value && acquiring) {
            /* This was a command to stop acquisition */
            /* Send the stop event */
            epicsEventSignal(stopEventId_); 
        }
    } else if ((function == NDDataType) || 
               (function == NDColorMode) ||
               (function == SimMode)) {
        status = setIntegerParam(SimResetImage, 1);
    } else {
        /* If this parameter belongs to a base class call its method */
        if (function < FIRST_SIM_DETECTOR_PARAM) status = ADDriver::writeInt32(pasynUser, value);
    }

    /* Do callbacks so higher layers see any changes */
    callParamCallbacks();

    if (status)
        asynPrint(pasynUser, ASYN_TRACE_ERROR,
              "%s:writeInt32 error, status=%d function=%d, value=%d\n",
              driverName, status, function, value);
    else
        asynPrint(pasynUser, ASYN_TRACEIO_DRIVER,
              "%s:writeInt32: function=%d, value=%d\n",
              driverName, function, value);
    return status;
}


/** Called when asyn clients call pasynFloat64->write().
  * This function performs actions for some parameters, including ADAcquireTime, ADGain, etc.
  * For all parameters it sets the value in the parameter library and calls any registered callbacks..
  * \param[in] pasynUser pasynUser structure that encodes the reason and address.
  * \param[in] value Value to write. */
asynStatus simDetector::writeFloat64(asynUser *pasynUser, epicsFloat64 value)
{
    int function = pasynUser->reason;
    asynStatus status = asynSuccess;

    /* Set the parameter and readback in the parameter library.  This may be overwritten when we read back the
     * status at the end, but that's OK */
    status = setDoubleParam(function, value);

    /* Changing any of the following parameters requires recomputing the base image */
    if ((function == ADAcquireTime) ||
        (function == ADGain) ||
        (function == SimGainX) ||
        (function == SimGainY) ||
        (function == SimGainRed) ||
        (function == SimGainGreen) ||
        (function == SimGainBlue)) {
            status = setIntegerParam(SimResetImage, 1);
    } else {
        /* If this parameter belongs to a base class call its method */
        if (function < FIRST_SIM_DETECTOR_PARAM) status = ADDriver::writeFloat64(pasynUser, value);
    }

    /* Do callbacks so higher layers see any changes */
    callParamCallbacks();
    if (status)
        asynPrint(pasynUser, ASYN_TRACE_ERROR,
              "%s:writeFloat64 error, status=%d function=%d, value=%f\n",
              driverName, status, function, value);
    else
        asynPrint(pasynUser, ASYN_TRACEIO_DRIVER,
              "%s:writeFloat64: function=%d, value=%f\n",
              driverName, function, value);
    return status;
}

/** Report status of the driver.
  * Prints details about the driver if details>0.
  * It then calls the ADDriver::report() method.
  * \param[in] fp File pointed passed by caller where the output is written to.
  * \param[in] details If >0 then driver details are printed.
  */
void simDetector::report(FILE *fp, int details)
{

    fprintf(fp, "Simulation detector %s\n", this->portName);
    if (details > 0) {
        int nx, ny, dataType;
        getIntegerParam(ADSizeX, &nx);
        getIntegerParam(ADSizeY, &ny);
        getIntegerParam(NDDataType, &dataType);
        fprintf(fp, "  NX, NY:            %d  %d\n", nx, ny);
        fprintf(fp, "  Data type:         %d\n", dataType);
    }
    /* Invoke the base class method */
    ADDriver::report(fp, details);
}

/** Constructor for simDetector; most parameters are simply passed to ADDriver::ADDriver.
  * After calling the base class constructor this method creates a thread to compute the simulated detector data,
  * and sets reasonable default values for parameters defined in this class, asynNDArrayDriver and ADDriver.
  * \param[in] portName The name of the asyn port driver to be created.
  * \param[in] maxSizeX The maximum X dimension of the images that this driver can create.
  * \param[in] maxSizeY The maximum Y dimension of the images that this driver can create.
  * \param[in] dataType The initial data type (NDDataType_t) of the images that this driver will create.
  * \param[in] maxBuffers The maximum number of NDArray buffers that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited number of buffers.
  * \param[in] maxMemory The maximum amount of memory that the NDArrayPool for this driver is
  *            allowed to allocate. Set this to -1 to allow an unlimited amount of memory.
  * \param[in] priority The thread priority for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  * \param[in] stackSize The stack size for the asyn port driver thread if ASYN_CANBLOCK is set in asynFlags.
  */
simDetector::simDetector(const char *portName, string const & arrayName, int maxSizeX,
                         int maxSizeY, epics::pvData::ScalarType dataType,
                         int maxBuffers, size_t maxMemory, int priority, int stackSize)

    : ADDriver(portName, arrayName, 1, 0, maxBuffers, maxMemory,
               0, 0, /* No interfaces beyond those set in ADDriver.cpp */
               0, 1, /* ASYN_CANBLOCK=0, ASYN_MULTIDEVICE=0, autoConnect=1 */
               priority, stackSize),
      startEventId_(epicsEventCreate(epicsEventEmpty)),
      stopEventId_(epicsEventCreate(epicsEventEmpty)),
      pRaw_(), xSine1_(0), xSine2_(0), ySine1_(0), ySine2_(0)

{
    int status = asynSuccess;
    char versionString[20];
    const char *functionName = "simDetector";

    createParam(SimGainXString,               asynParamFloat64, &SimGainX);
    createParam(SimGainYString,               asynParamFloat64, &SimGainY);
    createParam(SimGainRedString,             asynParamFloat64, &SimGainRed);
    createParam(SimGainGreenString,           asynParamFloat64, &SimGainGreen);
    createParam(SimGainBlueString,            asynParamFloat64, &SimGainBlue);
    createParam(SimNoiseString,               asynParamInt32,   &SimNoise);
    createParam(SimResetImageString,          asynParamInt32,   &SimResetImage);
    createParam(SimModeString,                asynParamInt32,   &SimMode);
    createParam(SimPeakNumXString,            asynParamInt32,   &SimPeakNumX);
    createParam(SimPeakNumYString,            asynParamInt32,   &SimPeakNumY);
    createParam(SimPeakStepXString,           asynParamInt32,   &SimPeakStepX);
    createParam(SimPeakStepYString,           asynParamInt32,   &SimPeakStepY);
    createParam(SimPeakStartXString,          asynParamInt32,   &SimPeakStartX);
    createParam(SimPeakStartYString,          asynParamInt32,   &SimPeakStartY);
    createParam(SimPeakWidthXString,          asynParamInt32,   &SimPeakWidthX);
    createParam(SimPeakWidthYString,          asynParamInt32,   &SimPeakWidthY);
    createParam(SimPeakHeightVariationString, asynParamInt32,   &SimPeakHeightVariation);
    createParam(SimSineOffsetString,          asynParamFloat64, &SimSineOffset);
    createParam(SimSineNoiseString,           asynParamFloat64, &SimSineNoise);
    createParam(SimXSineOperationString,      asynParamInt32,   &SimXSineOperation);
    createParam(SimYSineOperationString,      asynParamInt32,   &SimYSineOperation);
    createParam(SimXSine1AmplitudeString,     asynParamFloat64, &SimXSine1Amplitude);
    createParam(SimXSine1FrequencyString,     asynParamFloat64, &SimXSine1Frequency);
    createParam(SimXSine1PhaseString,         asynParamFloat64, &SimXSine1Phase);
    createParam(SimXSine2AmplitudeString,     asynParamFloat64, &SimXSine2Amplitude);
    createParam(SimXSine2FrequencyString,     asynParamFloat64, &SimXSine2Frequency);
    createParam(SimXSine2PhaseString,         asynParamFloat64, &SimXSine2Phase);
    createParam(SimYSine1AmplitudeString,     asynParamFloat64, &SimYSine1Amplitude);
    createParam(SimYSine1FrequencyString,     asynParamFloat64, &SimYSine1Frequency);
    createParam(SimYSine1PhaseString,         asynParamFloat64, &SimYSine1Phase);
    createParam(SimYSine2AmplitudeString,     asynParamFloat64, &SimYSine2Amplitude);
    createParam(SimYSine2FrequencyString,     asynParamFloat64, &SimYSine2Frequency);
    createParam(SimYSine2PhaseString,         asynParamFloat64, &SimYSine2Phase);

    /* Set some default values for parameters */
    status =  setStringParam (ADManufacturer, "Simulated detector");
    status |= setStringParam (ADModel, "Basic simulator");
    epicsSnprintf(versionString, sizeof(versionString), "%d.%d.%d", 
                  DRIVER_VERSION, DRIVER_REVISION, DRIVER_MODIFICATION);
    setStringParam(NDDriverVersion, versionString);
    setStringParam(ADSDKVersion, versionString);
    setStringParam(ADSerialNumber, "No serial number");
    setStringParam(ADFirmwareVersion, "No firmware");

    status |= setIntegerParam(ADMaxSizeX, maxSizeX);
    status |= setIntegerParam(ADMaxSizeY, maxSizeY);
    status |= setIntegerParam(ADMinX, 0);
    status |= setIntegerParam(ADMinY, 0);
    status |= setIntegerParam(ADBinX, 1);
    status |= setIntegerParam(ADBinY, 1);
    status |= setIntegerParam(ADReverseX, 0);
    status |= setIntegerParam(ADReverseY, 0);
    status |= setIntegerParam(ADSizeX, maxSizeX);
    status |= setIntegerParam(ADSizeY, maxSizeY);
    status |= setIntegerParam(NDArraySizeX, maxSizeX);
    status |= setIntegerParam(NDArraySizeY, maxSizeY);
    status |= setIntegerParam(NDArraySize, 0);
    status |= setIntegerParam(NDDataType, dataType);
    status |= setIntegerParam(ADImageMode, ADImageContinuous);
    status |= setDoubleParam (ADAcquireTime, .001);
    status |= setDoubleParam (ADAcquirePeriod, .005);
    status |= setIntegerParam(ADNumImages, 100);
    status |= setIntegerParam(SimNoise, 3);
    status |= setIntegerParam(SimResetImage, 1);
    status |= setDoubleParam (SimGainX, 1);
    status |= setDoubleParam (SimGainY, 1);
    status |= setDoubleParam (SimGainRed, 1);
    status |= setDoubleParam (SimGainGreen, 1);
    status |= setDoubleParam (SimGainBlue, 1);
    status |= setIntegerParam(SimMode, 0);
    status |= setIntegerParam(SimPeakStartX, 1);
    status |= setIntegerParam(SimPeakStartY, 1);
    status |= setIntegerParam(SimPeakWidthX, 10);
    status |= setIntegerParam(SimPeakWidthY, 20);
    status |= setIntegerParam(SimPeakNumX, 1);
    status |= setIntegerParam(SimPeakNumY, 1);
    status |= setIntegerParam(SimPeakStepX, 1);
    status |= setIntegerParam(SimPeakStepY, 1);
    status |= setIntegerParam(SimPeakHeightVariation, 3);

    if (status) {
        printf("%s: unable to set camera parameters\n", functionName);
        return;
    }

    /* Create the thread that updates the images */
    status = (epicsThreadCreate("SimDetTask",
                                epicsThreadPriorityMedium,
                                epicsThreadGetStackSize(epicsThreadStackMedium),
                                (EPICSTHREADFUNC)simTaskC,
                                this) == NULL);
    if (status) {
        printf("%s:%s epicsThreadCreate failure for image task\n",
            driverName, functionName);
        return;
    }
}

/** Configuration command, called directly or from iocsh */
extern "C" int simDetectorConfig(const char *portName, const char *arrayName, int maxSizeX, int maxSizeY, int dataType,
                                 int maxBuffers, int maxMemory, int priority, int stackSize)
{
    new simDetector(portName, arrayName, maxSizeX, maxSizeY, (ScalarType)dataType,
                    (maxBuffers < 0) ? 0 : maxBuffers,
                    (maxMemory < 0) ? 0 : maxMemory, 
                    priority, stackSize);
    return(asynSuccess);
}

/** Code for iocsh registration */
static const iocshArg simDetectorConfigArg0 = {"Port name", iocshArgString};
static const iocshArg simDetectorConfigArg1 = {"Array name", iocshArgString};
static const iocshArg simDetectorConfigArg2 = {"Max X size", iocshArgInt};
static const iocshArg simDetectorConfigArg3 = {"Max Y size", iocshArgInt};
static const iocshArg simDetectorConfigArg4 = {"Data type", iocshArgInt};
static const iocshArg simDetectorConfigArg5 = {"maxBuffers", iocshArgInt};
static const iocshArg simDetectorConfigArg6 = {"maxMemory", iocshArgInt};
static const iocshArg simDetectorConfigArg7 = {"priority", iocshArgInt};
static const iocshArg simDetectorConfigArg8 = {"stackSize", iocshArgInt};
static const iocshArg * const simDetectorConfigArgs[] =  {&simDetectorConfigArg0,
                                                          &simDetectorConfigArg1,
                                                          &simDetectorConfigArg2,
                                                          &simDetectorConfigArg3,
                                                          &simDetectorConfigArg4,
                                                          &simDetectorConfigArg5,
                                                          &simDetectorConfigArg6,
                                                          &simDetectorConfigArg7,
                                                          &simDetectorConfigArg8};
static const iocshFuncDef configsimDetector = {"simDetectorConfig", 9, simDetectorConfigArgs};
static void configsimDetectorCallFunc(const iocshArgBuf *args)
{
    simDetectorConfig(args[0].sval, args[1].sval, args[2].ival, args[3].ival,
                      args[4].ival, args[5].ival, args[6].ival, args[7].ival,
                      args[7].ival);
}


static void simDetectorRegister(void)
{

    iocshRegister(&configsimDetector, configsimDetectorCallFunc);
}

extern "C" {
epicsExportRegistrar(simDetectorRegister);
}
