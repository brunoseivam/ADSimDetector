# Must have loaded envPaths via st.cmd.linux or st.cmd.win32

errlogInit(20000)

dbLoadDatabase("$(TOP)/dbd/simDetectorApp.dbd")
simDetectorApp_registerRecordDeviceDriver(pdbbase) 

epicsEnvSet("PREFIX", "13SIM1:")
epicsEnvSet("PORT",   "SIM1")
epicsEnvSet("QSIZE",  "20")
epicsEnvSet("XSIZE",  "1024")
epicsEnvSet("YSIZE",  "1024")
epicsEnvSet("NCHANS", "2048")
epicsEnvSet("CBUFFS", "500")
epicsEnvSet("MAX_THREADS", "8")
epicsEnvSet("EPICS_DB_INCLUDE_PATH", "$(ADCORE)/db")

asynSetMinTimerPeriod(0.001)

# Create a simDetector driver
# simDetectorConfig(const char *portName, const char *pvName, int maxSizeX, int maxSizeY, int dataType,
#                   int maxBuffers, int maxMemory, int priority, int stackSize)
simDetectorConfig("$(PORT)", $(PREFIX)cam1:Image, $(XSIZE), $(YSIZE), 1, 0, 0)
# To have the rate calculation use a non-zero smoothing factor use the following line
#dbLoadRecords("simDetector.template",     "P=$(PREFIX),R=cam1:,PORT=$(PORT),ADDR=0,TIMEOUT=1,RATE_SMOOTH=0.2")
dbLoadRecords("$(ADSIMDETECTOR)/db/simDetector.template","P=$(PREFIX),R=cam1:,PORT=$(PORT),ADDR=0,TIMEOUT=1")

NDStatsConfigure("STATS1", "$(PREFIX)Stats1:Image", $(QSIZE), 0, $(PREFIX)cam1:Image, 0, 0, 0, 0, $(MAX_THREADS=5))
dbLoadRecords("NDStats.template", "P=$(PREFIX),R=Stats1:, PORT=STATS1,ADDR=0,TIMEOUT=1,HIST_SIZE=256,XSIZE=$(XSIZE),YSIZE=$(YSIZE),NCHANS=$(NCHANS),NTNDARRAY_PV=$(PREFIX)cam1:Image")

# Create a second simDetector driver
#simDetectorConfig("SIM2", 300, 200, 1, 50, 50000000)
#dbLoadRecords("$(ADSIMDETECTOR)/db/simDetector.template","P=$(PREFIX),R=cam2:,PORT=SIM2,ADDR=0,TIMEOUT=1")

# Load all other plugins using commonPlugins.cmd
#< $(ADCORE)/iocBoot/commonPlugins.cmd

startPVAServer()

iocInit()

dbgrep *Block*
dbpf $(PREFIX)cam1:AcquirePeriod 1
dbpf $(PREFIX)Stats1:BlockingCallbacks 0
dbpf $(PREFIX)cam1:Acquire 1
dbpf $(PREFIX)Stats1:EnableCallbacks 1
#dbpf $(PREFIX)cam1:NDAttributesFile "/home/bmartins/workspace-adv4/areaDetector/ADSimDetector/iocs/simDetectorIOC/iocBoot/iocSimDetector/simDetectorAttributes.xml"
dbpf $(PREFIX)cam1:NDAttributesFile "simDetectorAttributes.xml"
