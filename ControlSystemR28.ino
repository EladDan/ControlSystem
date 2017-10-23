#include <MegunoLink.h>
#include <CommandHandler.h>
#include <CommandProcessor.h>
#include <OneWire.h>
#include <EEPROMStore.h>
#include <Filter.h>
#include <math.h>
#include <SoftwareSerial.h>

//Pins definition
# define Lout_Pin 2  // OUTPUT of TSL235R out-connected to Digital pin 2 (interrupt 0)
# define Lin_Pin 3   // OUTPUT of TSL235R in -connected to Digital pin 3 (interrupt 1)
// # define Lre_Pin 21  // OUTPUT of TSL235R reflection -connected to Digital pin 21 (interrupt 2)
# define RedPin 5    //   PWM pin 5 for DIM Red LEDS 664 nm
# define BluePin 6   //   PWM pin 6 for DIM Blue LEDS 451 nm
// PINS 0 (RX) AND 1 (TX) are not allowed: reserved for Arduino's serial monitor !!
# define rx 16           // pin rx: RX pH5 to TX2 Arduino Mega 16 (Serial2) green
# define tx 17           // pin Tx: TX pH5 to RX2 Arduino Mega 17 (Serial2) yellow
# define DS18B20_Pin 43  // Temperature sensor pin
# define Relay_CO2 45    // for Arduino Mega; for Arduino Uno, use Digital I/O pin 11
# define Relay_Off 0     // CO2 valve Closed, when reaching pHLow from above
# define Relay_On 1      // CO2 valve Open, when reaching pHHigh from below

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Declaration of parameters and variables XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// Debug parameters for reduced running time
// int debug_mode = 1; // flag
unsigned int red_light_millis;
unsigned int blue_light_millis;

// Dark Background [Hz]
double BgLi = 0.099;  // CA: Changed from 0.1
double BgLo = 0.167;  // CA: Changed from 0.157
// double BgLr = 0.1XX;

// Transmittance through filters in, out and refl for red and blue - INCLUDED IN col_NPARa1
// Usage: filter_corr_sensor = sensor/xfiltery
// const double rFilteri = 0.013402;   // 2 layers filter (Sensor in)
// const double rFiltero = 0.087755;   // 1 layer filter (Sensors out and refl)
// const double bFilteri = 0.037965;   // 2 layers filter (Sensor in)
// const double bFiltero = 0.159598;   // 1 layer filter (Sensors out and refl)

// sensor sensitivity, relative to WL 625 nm - INCLUDED IN col_NPARa1
// Usage: corr_sensor = sensor/X_sens
// const double rSens = 1.05;
// const double bSens = 0.65;

// Predict Normalized sensor using PWM; from file Calib161010.xlsx/
// Estimates for Box-Cox: Light = ((a0 + PWM * a1) * lambda + 1) ^ (1/lambda)
// Back : PWM =(((Light^lambda-1)/lambda)-a0)/a1

// Red
double rLambda = 0.6649892;
double ra1 = 0.01278822;
double ra0 = -1.823489;
// Blue
double bLambda = 0.5178513;
double ba1 = 0.01688053;
double ba0 = -2.3530692;

// Linear transform sensor [Hz] to PAR [uE] normalized values. global parameter from Calibs_fit.xlsx/Red-Blue norm 140
// Usage: PAR_norm = xNPARa1 * Hz_corr_norm
double rNPARa1 = 1.041;
double bNPARa1 = 0.9192;

// EpsR (652-674) and EpsB (440-464) for Abs  in L/cm/mg Chl from file PWMNormPWM.xlsx/packSpect (Bidigare)
// Abs = -ln(Trans) * 2.303 = -log10(Trans)

// Red
double rEpsa = 0.07523 ;  // Red epsilon for Chla specific to Chla          CA: Updated Spectra.xlsx/packSpec
double rEpsb = 0.07815 ;  // Red epsilon for Chlb specific to Chla          CA: Updated Spectra.xlsx/packSpec
double rEpsc = 0.0 ;       // Red epsilon for Car specific to Chla
double rKd = 0.002476495;     // Red coefficient for DW_; Last value from LAB
// Blue
double bEpsa = 0.1093 ;  // Blue epsilon for Chla specific to Chla
double bEpsb = 0.1705 ;  // Blue epsilon for Chlb specific to Chla
double bEpsc = 0.5189 ;  // Blue epsilon for Car specific to Chla
double bKd = 0.000896663;  // Blue coefficient for DW_; Last value from LAB
// Measured pigments and DW concentration (mg/L) MANUAL INPUT  through !SetComment
double ChlaLab;      //  MANUAL INPUT
double ChlbLab;      //  MANUAL INPUT
double CarLab;       //  MANUAL INPUT
double DWLab;        //  MANUAL INPUT

// Initial Concentrations, relative to Chla (Averages of pooled Lab Data) to be updated further
// increased light <==> decreased Chlb/Chla (Chlb_)
double Chlb_ = 0.502;    // = ChlbLab / ChlaLab
// increased stress <==> increased Car/Chla (Car_)
double Car_ = 0.721;     // = CarLab / ChlaLab    from last LAB
// increased stress <==> increased Car/Chla (Car_)
double DW_ = 43.3;       // = DWLab / ChlaLab    from last LAB

// Current concentrations for pigments and DW (mg/L)
double Chla;
double Chlb;
double Car;
double DW;

// Combined epsilons (weighted by concentration relative to Chla) for each color, independent of Chla
// The factors XKd*DW_*Chla represent the increase in optical path due to scattering, independent of Chla
double rSumeps_;           // = (rEpsa + rEpsb * Chlb_ + rEpsc * Car_)
double bSumeps_;           // = (bEpsa + bEpsb * Chlb_ + bEpsc * Car_)

const double width = 4.0;   // Light path width of the culture vessel in cm
double rOp_;                // Red multiplyer optical path  = rKd * DW_ , independent of Chla and width
double bOp_;                // Blue multiplyer optical path = bKd * DW_ , independent of Chla and width

// Blank (no added algae) values for sensors [Hz] at PWM=140, from Calib: Calibs_fit.xlsx/Red-Blue norm 140
double rLi_0 = 14576 - BgLi;
double rLo_0 = 3028 - BgLo;
//double rLr_0 = XXXX - BgLr;
double bLi_0 = 7462 - BgLi;
double bLo_0 = 683 - BgLo;
//double bLr_0 = XXXX - BgLr;

// Light measured for Con2 calculation
// sensor [Hz] at PWM=140. Measured at setup() and updated periodically
double rLi_ref;
double rLo_ref;
//double rLr_ref;
double bLi_ref;
double bLo_ref;
//double bLr_ref;

// Area under the corrected light profile in the PBR (as PAR) per chlorophyll
double rL_avail;  // = <L>r/Chla
double bL_avail;  // = <L>b/Chla
double L_avail;   // = rL_avail + bL_avail

// Calculated (split L into rL and bL, in Hz); a solution of system of 2 equations with 2 unknowns
// (PWMNormPWM.xlsx/ FullFormula)
double rLi;
double rLo;
//double rLr;
double bLi;
double bLo;
//double bLr;

// and transformed to PAR:
double rPARi;
double rPARo;
//double rPARr;
double bPARi;
double bPARo;
//double bPARr;
double PARi;   // rPARi + bPARi

// PAR values at 140 PWM with algae
double rPARi_ref;
double rPARo_ref;
//double rPARr_ref;
double bPARi_ref;
double bPARo_ref;
//double bPARr_ref;

// PAR [uE/m2/s] at PWM=140 at zero algae, from Calib: Calibs_fit.xlsx/Red-Blue norm 140
double rPARi_const = 753.0;
double bPARi_const = 228.8;
//double rPARr_const = XXX.X;
double rPARo_const = 43.30;       // (used for Transformation Sensor<=>PAR only).
double bPARo_const = 6.94;        // (used for Transformation Sensor<=>PAR only).
//double bPARr_const = XXX.X;

// factors for conversion Hz to PAR. calculate as follows: (see PWMNormPWM.xls/Hz2PAR)
// Transform Lightin [Hz] to Lightin [PAR]
// norm_Hz = (Hz-Bg) / L_0;            // normalized Hz
// norm_PAR = (Hz-Bg) * NPARa1;        // transform to normalized PAR
// PAR = norm_PAR * L_PAR_const;       // PAR

// Conversion factors from Hz to PAR
double rHz2PARi = 0.05311;      // = rNPARa1*rPARi_const/rLi_0;
double rHz2PARo = 0.01488;      // = rNPARa1*rPARo_const/rLo_0;
//double rHz2PARr = 0.XXXXX;      // = rNPARa1*rPARr_const/rLr_0;
double bHz2PARi = 0.02818;      // = bNPARa1*bPARi_const/bLi_0;
double bHz2PARo = 0.009338;      // = bNPARa1*bPARo_const/bLo_0;
//double bHz2PARr = 0.XXXXX;      // = bNPARa1*bPARr_const/bLr_0;

// PBR Absorbance (-log(PARout/PARin), in uE): fixed and global at zero Con,
// from file AbsCalc.xlsx/Param

const double rPBR = -log10(rLo_0*rHz2PARo / rLi_0 / rHz2PARi);      // 1.235  Changed from log to log10
const double bPBR = -log10(bLo_0*bHz2PARo / bLi_0 / bHz2PARi);      // 1.518  Changed from log to log10

// dividers (= Transmittance through 1 face of PBR) for col_PARi to supply.
// usage: {col_PARi needed} = {col_PARi to supply} / col_Xi
const double rTransPBR = pow(10, -0.333 * rPBR);                          // 0.3879 = red transmitted by one PBR face
const double bTransPBR = pow(10, -0.333 * bPBR);                          // 0.3122 = blue transmitted by one PBR face

// Light and Con stuff
double ConMax = 10.0;      // concentration at and beyond which PWM (at noon) is PWMmax; user adjustable
double rPARMax = 1978;     // Maximal red PAR that can be supplied by the LEDs (~ 1 Sun)
double bPARMax = 719;      // Maximal blue PAR that can be supplied by the LEDs (~ 0.3 Sun)
// double PARMax = rPARMax + bPARMax;      // Maximal PAR (R+B) that can be supplied by the LEDs (~ 1.3 Sun)

double  rSAC = (rEpsa + rEpsb * Chlb_ + rEpsc * Car_) * rKd * DW_;
double  bSAC = (bEpsa + bEpsb * Chlb_ + bEpsc * Car_) * bKd * DW_;

// Constant Specific Light supply = ratio between PARi [uE/m2/s] (including correction by (1-e)) to Con^3 [mg^3/L^3]. Taken from PWMNormPWM.xlsx/PARimax(FrR)
// Constant Specific Light supply to algae (NOT to PBR by LEDS)
double rSpecLSupply = rPARMax * rTransPBR * ( 1 - exp(-2.30285 * rSAC * pow(ConMax, 2) * width)) / pow(ConMax, 3); // 0.09591
double bSpecLSupply = bPARMax * bTransPBR * ( 1 - exp(-2.30285 * bSAC * pow(ConMax, 2) * width)) / pow(ConMax, 3); // 0.02806

//Time constants for checks (in milliseconds)
unsigned long TCon;                  // Period for Con calculation (in msec), dependent on the illum mode and time of day
const unsigned long TLmax = 20000;   // Maximal period for Light measure - 20 sec
const unsigned long TpH = 12500;     // Period for pH measure - 12.5 sec
const unsigned long TTemp = 300000;  // Period for Temp measure - 5 min

// ensure entrance to con2 at the beginning (so we don't need to wait up to 4 hours until con2 measurement)
double TCon2 = 0.0;
// initialize with zero so its related "if" condition will be true at start of run,
// until first time entering Con2 measurement,
// then, it get the value 3.9 - time between Con2 measurements.
double Con2Period = 4.0;             // Period for Con2 calculation at night (!IsDay) [hours]

unsigned long TLi = 2000;
unsigned long TLo = 2000;
//unsigned long TLr = 2000;
unsigned long TLastLi = 0;
unsigned long TLastLo = 0;
//unsigned long TLastLr = 0;
unsigned long TCurrentLi = 0;       // Recorded time just before Li check
unsigned long TCurrentLo = 0;       // Recorded time just before Lo check
//unsigned long TCurrentLr = 0;       // Recorded time just before Lr check
unsigned long TCurrentCon = 0;

unsigned long TLastCon = 0;
unsigned long TLastpH = 0;
unsigned long TLastTemp = 0;  // Time for last Temperature measure
unsigned long TLastCon2 = 0;  // Time for last Con2 measure in millis

unsigned long PerCon = 0;         // seconds spent from the begining of current TCon

volatile unsigned long TOpenCO2 = 0.0;     // *** CA: changed from 'LastCO2' for consistent nomenclature
volatile unsigned long TCurrentCO2 = 0.0;  // *** CA: changed from 'CurrentCO2' to 'TCurrentCO2'
volatile unsigned long TCloseCO2 = 0.0;    // *** CA: changed from 'CO2CycleBegin' to 'TCloseCO2'; used as start/end of CO2 cycle

// Counters Light
unsigned int TimesLo = 0;
unsigned int TimesLi = 0;
//unsigned int TimesLr = 0;

volatile unsigned long pulses_out = 0 ;
volatile unsigned long pulses_in = 0 ;
// volatile unsigned long pulses_re = 0 ;
unsigned long TempPulseI;
unsigned long TempPulseO;
// unsigned long TempPulseR;

// Fraction Red available
double FrR_target = 0.95;   // $$$ Still need to be determined adjustable between 0.75 and 1.0
double FrR;                 //

// Illumination status (switch mode)
unsigned int Illum = 3;
// 0 - continous dark (default)
// 1 - continous light (non stop)
// 2 - constant light during the day
// 3 - sinusoid light function

// Cummulative and input vars [Hz]
double Li ;
double Lo ;
// double Lr ;

// Since we are checking different times, we are going to need counters for how many tests we did to each sensor
const unsigned int ko = 5 ;    // Number of times Lo is measured for small average
const unsigned int ki = 10;    // Number of times Li is measured for small average
// const unsigned int kr = 5;  // Number of times Lr is measured for small average
double LoSu = 0.0 ;            // temporary sum in TLo before the LoSum for LoAvg in TCon
double LiSu = 0.0 ;            // temporary sum in TLi before the LiSum for LiAvg in TCon
// double LrSu = 0.0 ;         // temporary sum in TLr before the LrSum for LrAvg in TCon
double LoSum = 0.0 ;
double LiSum = 0.0 ;
// double LrSum = 0.0 ;
double LoAvg = 0.0 ;     // Time average (during TCon) for Lo
double LiAvg = 0.0 ;     // Time average (during TCon) for Li
// double LrAvg = 0.0 ;  // Time average (during TCon) for Lr

double rPARimax;         // Red in to be supplied [uE] at noon in any Illum 3
double bPARimax;         // Blue in to be supplied [uE] at noon in any Illum 3

double rAbs2;
double bAbs2;

//double rRef2;
//double bRef2;

double rAbs;
double bAbs;

//double rRef;
//double bRef;

double b_rAbs;          // ratio bAbs/rAbs
//double b_rRef;        // ratio bRef/rRef

const unsigned int Con2PWM = 140;    // Fixed PWM (light intensity) for Con2 measurement and calculation (AKA _ref)
double rCon2;
double bCon2;
double Con2Old;      // Concentration of Chla in mg Chl/L at the begining of current Con2Period

double rCon;
double bCon;

// rate calculation related variables
double Mu;
unsigned int iStartMu;
double ConMuFirst;      // value of first con calc (renew every change of parameter) previous mix
double ConRate;         // calculated rate of growth in mgChl/L/day
double SpConRate;       // specific growth rate  = 2.0*ConRate/(Chla+Con2Old) in 1/d
// specific growth rate = 24.0*log(Chla/Con2Old)/(tor - TimeMuOld) in 1/d

// unsigned long TCurrentParameter = 0;
// unsigned long TLastParameter = 0;
// unsigned long TParameter = 0;

// PWM vars
// Constants
const unsigned int rPWMmin = 27;
const unsigned int bPWMmin = 26;
const double PWMmax = 245;      // Absolute maximum for constrained PWMs

unsigned int rPWM;    // Power for RED (at actual time of day)
unsigned int bPWM;    // Power for Blue (at actual time of day)

// Time state boolean
boolean IsDay;           // if daytime (IsDay = TRUE), else (IsDay = FALSE)
boolean IsSetup;         // used in Cmd_UpdateIndicators() & ConPWMCalc()

// Time vars
double tod0;              // the hour of the day at program start in hours (24 h)
double daylength = 18.0;  // length of the day in hours out of 24 h, adjustable with DATE
double daystart = (24.0 - daylength) / 2.0;      // start of the day; Max Light at Noon
double tod;               // time of the day in hour (0.000-23.999)
double todOld;            // time measured previously
double todIP;             // time displayed previously (every 0.01 h = 36 sec)
double tor;               // time of run since the begining (h)
double TMuFirst;            // time of first con calculation (renew every change of parameter)
double TMu;               // time elapsed since TMuFirst

unsigned int dayofRun = 0; // from start
double radr;     // time of the day in radian for red leds PWM adjustment in diel variable light
double radb;     // time of the day in radian for blue leds PWM adjustment in diel variable light
double cosred;   // cosine of radr
double cosblue;  // cosine of radb

// Temperature
OneWire  ds(DS18B20_Pin); // on pin 10 (a 4.7K resistor is necessary)
// double Tmax = 27.0 ;
byte TempPresent = 0;
byte TempData[12];
byte addr[8];
double Temp;
char TempAscii[8];   // temperature in ascii format for pH serial command: XX.XX\r

// pH
byte received_from_pH = 0;
char pHData[20];      // array for ph_circuit incoming data
double pHAvg = 8.5;   //  *** CA: made adaptive to Day/Nite; effective pH = pHAvg +/- 0.05
double pHAvgOld = 8.5;         // for pHAvg update on control panel
double pHHigh = pHAvg + 0.05;  // Maximal pH allowed
double pHLow = pHAvg - 0.05;   // Minimal pH allowed -- difference = 0.1
double pH;             // used to hold a doubleing point number that is the pH
double pHMin = 6.5;
double pHRange = 2.0;  // = 8.5 - 6.5
// parameters for logistic function pHAvg vs. Li
double steepness = -0.00008; // slope dpH / dLi at LiMiddle point (inflexion point)
unsigned int LiMiddle = 2000;       // Li value at the inflexion point (Hz)

// CO2
unsigned int CO2 = 7;     // Used for solenoid status: Closed = 7; Open = 8.
unsigned int CO2Old = 7;  // used for data acquisition of CO2 status
double CO2Fraction;       // fraction time of CO2 open in last CO2 cycle
double Jc;                // CO2 Flux in mM/h
boolean IsValveOn = false;

// for parsing lab results received from MLP
String comment;
String comment_input;
int comma_index;
String Chla_Lab;
String Chlb_Lab;
String Car_Lab;
String DW_Lab;

char ToString[20];     // common temporary array for converting double/int to string
char timeInput[30];    // temporary char[] for parsing Date and Time
String Date;

// XXXXXXXXXXXXXXXXXXXXXX END of Declaration of parameters and variables XXXXXXXXXXXXXXXXXXXXXXXXXXX

/* XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX List SUBs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // List here Annotated subroutines used (include input /output if applicable)
  SetupMLP()
  TimeCalc()
  ClearIndicatorsFunction();  // cleaning the Control Panel fields
  UpdateIndicatorsFunction(); // update indicators
  SetupPins();
  SetupRelay();
  SetupTemp();
  SetupPH();
  InitLights(rPWM, bPWM);          // Send temporary PWMs to LEDs, Measure and report SD and AVG to CSV
  ConCalc(Abs, color)              // Abs(PAR), corrected for PBR; color = R or B
  PWMCalcSend();
  PrintCSV("Setup");                    // using IsSetup
  ReadTemp();
  ReadPH();
  LightIntervals();  // adaptive TLi and TLo calculation; returns proper TLi and TLo
  MeasureAbs2(false);
  PrintCSV("continuous"); // Print datum to CSV File
  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX End of  List SUBs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
*/

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX MEGUNOLINK Section XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//Declare TIMEPLOTs
TimePlot RPWMTP;
TimePlot BPWMTP;
TimePlot LoTP;
TimePlot LiTP;
TimePlot LoAvTP;
TimePlot LiAvTP;
TimePlot FrRTP;
TimePlot rAvaTP;
TimePlot bAvaTP;
TimePlot Car_TP;
TimePlot DW_TP;
TimePlot ConRTP;
TimePlot ConBTP;
TimePlot b_rAbsTP;
TimePlot rPARiTP;
TimePlot bPARiTP;
TimePlot rPARoTP;
TimePlot bPARoTP;
TimePlot PARiTP;
TimePlot MuTP;
TimePlot SpRaTP;
TimePlot TempTP;
TimePlot PHTP;
TimePlot COTP;      // CO2 status: 7 = Valve on;  or 8 = Valve off
TimePlot JcTP;      // CO2 Flow in mM/h

// Message Monitor
// "LogFile" = "CSV" = the target message channel (remember to select this in megunolink)
Message MyCSVMessage("CSV");
// "LogFile" = "Params" = the target message channel (remember to select this in megunolink)
Message MyParams("Params");

// *************************** MEGUNOLINK INTERFACE *****************************
// Communication with the Interface Panel (Panel to Arduino)
CommandHandler <10, 45> SerialCommandHandler; // memory reserved for up to 10 Commands, each of length 40

InterfacePanel Panel;                     // used for Monitor & Control
InterfacePanel LPanel("LibraryPanel");    // invoking commands from the message library,
// it has to be done in a different panel, meaning a different channel

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX End of MEGUNOLINK Section XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX SETUP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void setup()
{
  Serial.begin(9600) ;
  delay(300);
  Serial2.begin(38400);  // for pH communication
  delay(300);

  PowerLEDs(0, 0); // allow sensors to decay
  delay(5000);     // Give the user time to open megunolink

  IsSetup = true;

  ClearIndicatorsFunction();  // cleaning the Control Panel fields

  // Initialize TimePlots, CommandHandler(Interface Panel Functions) & Graphic Indicator for New Run
  MessageMonitor("Screen", "## Initialize", "Megunolink Pro functions");

  SetupMLP();

  TimeCalc();

  MessageMonitor("Screen", "## Graphic Indicator for", "Full Reset");
  for (unsigned int ii = 0  ;  ii < 31; ii = ii + 2)
  {
    BPWMTP.SendData("BPWM", ii);
  }

  MessageMonitor("Screen", "## Setup", "Light, Temp & pH Procedures");
  SetupPins();
  SetupRelay();
  SetupTemp();
  SetupPH();
  delay(100);             // waiting because of MegunoLink opening

  // print header for CSV title
  MyCSVMessage.Begin();
  Serial.println("New run.");
  Serial.print(Date);
  Serial.print(",");
  Serial.print(tod);
  MyCSVMessage.End();

  MessageMonitor("Screen", "## Initialize", "Lights, Con & PWM max calculation");

  // Initial Red light
  MessageMonitor("Screen", "## Red Light Absorbance", "10*2 sec");
  InitLights(Con2PWM, bPWMmin);          // Measure red Abs after bg subtraction
  rLi_ref = LiAvg;
  rLo_ref = LoAvg;
  // rLr_ref = LrAvg;
  UpdateIndicatorsFunction("rLi_ref");
  UpdateIndicatorsFunction("rLo_ref");
  //  UpdateIndicatorsFunction("rLr_ref");

  // Convert Units from Hz to uE/s, Calculate PAR_in and PAR_out
  rPARi_ref = rHz2PARi * rLi_ref;
  rPARo_ref = rHz2PARo * rLo_ref;
  //  rPARr_ref = rHz2PARr * rLr_ref;

  // calculate and report raw rAbs
  rAbs2 = -log10(rPARo_ref / rPARi_ref) - rPBR;           //  Changed from log to log10
  UpdateIndicatorsFunction("rAbs2");

  // Initial Blue light
  MessageMonitor("Screen", "## Blue Light Absorbance", "10*20 sec");

  // Send temporary PWMs to LEDs, Measure and report SD and AVG to CSV
  InitLights(rPWMmin, Con2PWM);
  bLi_ref = LiAvg;
  bLo_ref = LoAvg;
  //  bLr_ref = LrAvg;
  UpdateIndicatorsFunction("bLi_ref");
  UpdateIndicatorsFunction("bLo_ref");
  //  UpdateIndicatorsFunction("bLr_ref");

  // Convert Units from Hz to uE/s, Calculate PAR_in and PAR_out
  bPARi_ref = bHz2PARi * bLi_ref;
  bPARo_ref = bHz2PARo * bLo_ref;
  //  bPARr_ref = bHz2PARr * bLr_ref;

  bAbs2 = -log10(bPARo_ref / bPARi_ref) - bPBR;         //  Changed from log to log10
  UpdateIndicatorsFunction("bAbs2");

  // Calculate ratio bAbs/rAbs
  b_rAbs = bAbs2 / rAbs2;
  UpdateIndicatorsFunction("b_rAbs");

  UpdateSumepsOp();

  CalcParameters();

  //  UpdateIndicatorsFunction("rSpecLSupply");
  //  UpdateIndicatorsFunction("bSpecLSupply");

  CalcPARimax();   // Calculate total incident light [uE/m2/s] to supply to the algae (NET) at noon

  /*    Obselete
    if (rPARimax < 150)   // Li at noon < 2800 Hz. Light supply is too low
    {
      ConMax = ConMax / 2.0;              // Reduce adjustable ConMax to increase Light supply
      UpdateIndicatorsFunction("ConMax");
      MessageMonitor("Screen", "## ConMax reduced to ", String(ConMax));

      rSpecLSupply = rPARMax * rTransPBR / pow(ConMax, 3);  //   767.2/ConMax^3
    //    bSpecLSupply = bPARMax * bTransPBR / pow(ConMax, 3);  //   224.5/ConMax^3
      CalcPARimax();   // Recalculate total incident light [uE/m2/s] to supply to the algae (NET) at noon
    }

    if (rPARimax > 770)   // Li at noon > 45000 Hz. Light supply is too high
    {
      ConMax = ConMax * 1.5;              // Increase adjustable ConMax to decrease Light supply
      UpdateIndicatorsFunction("ConMax");
      MessageMonitor("Screen", "## ConMax reduced to ", String(ConMax));

      rSpecLSupply = rPARMax * rTransPBR / pow(ConMax, 3);  //   767.2/ConMax^3
    //    bSpecLSupply = bPARMax * bTransPBR / pow(ConMax, 3);  //   224.5/ConMax^3
      CalcPARimax();   // Recalculate total incident light [uE/m2/s] to supply to the algae (NET) at noon
    }
  */

  TimeCalc();

  // allow sensors to decay
  PowerLEDs(0, 0);
  delay(2000);

  // Calculate, update and send PWM
  PWMCalcSend();

  iStartMu = 0;                      // initiate Mu cycles counter
  PrintCSV("Setup");                 // using IsSetup

  ReadTemp();
  ReadPH();

  UpdateIndicatorsFunction("ALL");

  MessageMonitor("Screen", "## Waiting for sensors to decay", "20 sec");

  IsSetup = false;

  MessageMonitor("Screen", "## Setup() is done, Starting Loop()", "##################");
  Con2Old = Chla;
  TLastCon2 = millis();         // CA: Restart TCon2 timer here
  TLastCon = millis();         // CA: Restart TCon timer here

}
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX End of Setup() XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX LOOP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void loop()
{
  // Read from SerialCommandHandler, receive commands and write to Interface Panel
  SerialCommandHandler.Process();

  TimeCalc();       // Update time (hour) and IsDay status and date if its a new day
  PWMCalcSend();    // Calculate current PWMs from col_PARimax, Illum and time; sends PWM to LEDs and reports the change

  LightIntervals();  // adaptive TLi and TLo calculation; returns proper TLi and TLo

  if (tod > 23.95 || tod < 0.05)
  {
    MessageMonitor("Screen", "------------------------------------------", "");
    MessageMonitor("Screen", "close to midnight; tod ", String(tod));
    MessageMonitor("Screen", "close to midnight; first condition ( <= 0.01)", String(fmod(tod, Con2Period)));
    MessageMonitor("Screen", "close to midnight; TCon2 ", String(TCon2));
    MessageMonitor("Screen", "close to midnight; second condition ( > TCon2)", String((millis() - TLastCon2) / (60.0 * 60.0 * 1000.0)));
    MessageMonitor("Screen", "------------------------------------------", "");
  }

  // Check if (its the Con2Period th hour e.g. 0:00, 4:00, 8:00 ..., means time to measure con2)
  // && (it has been over (Con2Period-0.1) (say 0.1 hours less than Con2Period) hours since last measurement, except the first time)
  if ( (fmod(tod, Con2Period) <= 0.01) && ((millis() - TLastCon2) / (60.0 * 60.0 * 1000.0) > TCon2) )
  {
    if (TCon2 != (Con2Period - 0.1))
    {
      TCon2 = Con2Period - 0.1;     // First time calculating Con2 in sketch
    }

    // periodic
    MeasureAbs2();

    CalcParameters();

    // Alternative rates calculation and reporting once in 4 h
    // volumetric rate of growth in mgChl/L/d

    ConRate = 24.0 * 60.0 * 60.0 * 1000.0 * (Chla - Con2Old) / (millis() - TLastCon2);
    // specific non-exponential growth rate in 1/d
    SpConRate = 2.0 * ConRate / (Chla + Con2Old);
    UpdateIndicatorsFunction("ConRate");
    UpdateIndicatorsFunction("SpConRate");

    CalcPARimax();   // Calculate total incident light [uE/m2/s] to supply to the algae (NET) at noon

    TimeCalc();
    PWMCalcSend();

    PrintCSV("periodic"); // Print datum to CSV File

    Con2Old = Chla;                // Save for ConRate calculations after 4 hours
    TLastCon2 = millis();
  }    // End of periodic measure

  // Seconds in TCon
  if (PerCon != (unsigned long)(millis() - TLastCon) / 1000)
  {
    UpdateIndicatorsFunction("PerCon");                        // ????? is this if statement necessary?
  }

  PerCon = (unsigned long)(millis() - TLastCon) / 1000;

  TCurrentCon = millis();


  // -------------------- Measure Li -------------------------------------------
  TCurrentLi = millis();

  if (TCurrentLi - TLastLi >= TLi)
  {
    TempPulseI = pulses_in;  // Seize the number of pulses in counted untill now.

    if (TempPulseI < 1)   // Error handling
    {
      TempPulseI = 1;
      MessageMonitor("Screen", "## We had 0 Pulses In, Constrain to", "1");
    }

    UpdateIndicatorsFunction("TLi");
    UpdateIndicatorsFunction("TempPulseI");

    Li = TempPulseI / ((TCurrentLi - TLastLi) / 1000.0);   // Calculating Light_in in Hz.
    Li = Li - BgLi;                                        // !!! CA: Background in subtraction

    LiSu += Li;          // Sums Light to calculate their small average later.
    LiSum += Li;         // Sums Light to calculate the big average at TCon.
    pulses_in = 0;       // Reset the pulses in counting
    TLastLi = millis();  // Time at the end of each single Li measurement
    TimesLi += 1;        // Update the number of Li measurements for calculating LiAvg

    // Calculate small average for Li
    if ((TimesLi % ki) == 0)
    {
      Li = LiSu / ki;

      CleanArray(); // After small average
      String(TimesLi).toCharArray(ToString , 10);
      MessageMonitor("Screen", "                TimesLi", ToString);

      CleanArray(); // After small average
      dtostrf(Li, 4, 2, ToString);
      MessageMonitor("Screen", "                   Li", ToString);

      UpdateIndicatorsFunction("Li");

      LiSu = 0.0;
    }
  }
  // ------------------------------- End Li measure --------------------------------

  // -------------------------------  Measure Lo -----------------------------------
  TCurrentLo = millis();

  if (TCurrentLo - TLastLo >= TLo)
  {
    TempPulseO = pulses_out; // Seize the number of pulses out counted untill now.
    if (TempPulseO < 1)      // Error handling
    {
      TempPulseO = 1;
      MessageMonitor("Screen", "## We had 0 Pulses Out, Constrain to", "1");
    }

    UpdateIndicatorsFunction("TLo");
    UpdateIndicatorsFunction("TempPulseO");

    Lo = TempPulseO / ((TCurrentLo - TLastLo) / 1000.0); // Calculating Light_out in Hz.
    Lo = Lo - BgLo;                                      // !!! CA: Background out subtraction

    LoSu += Lo;         // Sums Light frequencies to calculate their (small) average later.
    LoSum += Lo;        // Sums Light frequencies to calculate the big average at TCon.
    pulses_out = 0;     // Reset the pulses out counting
    TLastLo = millis(); // Time at the end of each single Lo measurement
    TimesLo += 1;       // Update the number of Li measurements for calculating LoAvg

    // Calculate small average for Lo
    if ((TimesLo % ko) == 0)  // Calculate the small average
    {
      Lo = LoSu / ko;

      CleanArray();        // After small average
      String(TimesLo).toCharArray(ToString , 10);
      MessageMonitor("Screen", "                TimesLo", ToString);

      CleanArray();
      dtostrf(Lo, 1, 2, ToString);
      MessageMonitor("Screen", "                    Lo", ToString);

      UpdateIndicatorsFunction("Lo");

      LoSu = 0.0;
    }
  }
  // --------------------------------- End Lo measure ------------------------------------------

  /*
    // -------------------------------  Measure Lr -----------------------------------
    TCurrentLr = millis();

    if (TCurrentLr - TLastLr >= TLr)
    {
      TempPulseR = pulses_re; // Seize the number of pulses reflected counted untill now.
      if (TempPulseR < 1)      // Error handling
      {
        TempPulseR = 1;
        MessageMonitor("Screen", "## We had 0 Pulses Reflected, Constrain to", "1");
      }

      UpdateIndicatorsFunction("TLr");
      UpdateIndicatorsFunction("TempPulseR");

      Lr = TempPulseR / ((TCurrentLr - TLastLr) / 1000.0); // Calculating Light_re in Hz.
      Lr = Lr - BgLr;                                      // !!! CA: Background reflected subtraction

      LrSu += Lr;         // Sums Light frequencies to calculate their (small) average later.
      LrSum += Lr;        // Sums Light frequencies to calculate the big average at TCon.
      pulses_re = 0;     // Reset the pulses reflected counting
      TLastLr = millis(); // Time at the end of each single Lr measurement
      TimesLr += 1;       // Update the number of Lr measurements for calculating LrAvg

      UpdateIndicatorsFunction("Lr");


      // Calculate small average for Lr
      if ((TimesLr % kr) == 0)  // Calculate the small average
      {
        Lr = LrSu / kr;

        CleanArray();        // After small average
        String(TimesLr).toCharArray(ToString , 10);
        MessageMonitor("Screen", "                TimesLr", ToString);

        CleanArray();
        dtostrf(Lr, 1, 2, ToString);
        MessageMonitor("Screen", "                    Lr", ToString);

        UpdateIndicatorsFunction("Lr");

        LrSu = 0.0;
      }
    }
    // --------------------------------- End Lr measure ------------------------------------------
  */

  // ---------------------------- Con, PWMmax calculation --------------------------------------
  if ((TCurrentCon - TLastCon) >= TCon)
  {
    LiAvg = LiSum / TimesLi;
    LoAvg = LoSum / TimesLo;
    //    LrAvg = LrSum / TimesLr;

    UpdateIndicatorsFunction("LiAvg");
    UpdateIndicatorsFunction("LoAvg");

    // reset counters
    LiSu = 0.0;
    LoSu = 0.0;
    LiSum = 0.0;
    LoSum = 0.0;
    TimesLi = 0;
    TimesLo = 0;

/*   Apparently wrong calculation
    if (rPWM > 60 && bPWM > 60)
    {
      // Convert Hz, Calculate Abs, Calculate concentrations (absolute Chla and relative Chlb_, Car_ and DW_)
      // We use the transmittance ratio of reference (140)
      // In any case, the ratio b_rAbs doesn't change since no new information on that is gathered

      // Calculating rLi, rLo, bLi & bLo as solution of system of 2 equations with 2 unknowns (PWMNormPWM.xlsx/ FullFormula)
      // "Split"(Li,Lo) into red and blue components

      bLi = (LoAvg - LiAvg * rLo_ref / rLi_ref) / ((bLo_ref / bLi_ref) - (rLo_ref / rLi_ref));
      bLo = bLi * bLo_ref / bLi_ref;
      rLi = LiAvg - bLi;
      rLo = rLi * rLo_ref / rLi_ref;

      // calculate rPARi, rPARo, bPARi, bPARo by converting Units from Hz to uE/s
      rPARi = rHz2PARi * rLi;   // YIELDS NEGATIVE VALUES FOR rPARi when PWMs are lower than ~60
      rPARo = rHz2PARo * rLo;   // YIELDS NEGATIVE VALUES FOR rPARo
      bPARi = bHz2PARi * bLi;
      bPARo = bHz2PARo * bLo;
      UpdateIndicatorsFunction("rPARi");
      UpdateIndicatorsFunction("bPARi");
      UpdateIndicatorsFunction("rPARo");
      UpdateIndicatorsFunction("bPARo");

      // Calculate rAbs
      rAbs = -log10(rPARo / rPARi) - rPBR;                   //  Changed from log to log10
      UpdateIndicatorsFunction("rAbs");

      // Calculate bAbs
      bAbs = -log10(bPARo / bPARi) - bPBR;                   //  Changed from log to log10
      UpdateIndicatorsFunction("bAbs");

      // Calculate Con red
      rCon = ConCalc(rAbs , "R");
      UpdateIndicatorsFunction("rCon");
      MessageMonitor("Screen", "            new", "rCon");

      // Calculate and update Pigments and DW concentrations    ???? CA: Is that OK? ???????????
      Chla = rCon;
      UpdateIndicatorsFunction("Chla");
      UpdatePigments(Chla);

      // Calculate Con blue
      bCon = ConCalc(bAbs , "B");
      UpdateIndicatorsFunction("bCon");
      MessageMonitor("Screen", "            new", "bCon");

    } // end of if (PWMs>60)

    else // (PWMs <= 60)
    {
      MessageMonitor("Screen", "            PWMs", "<= 60");
      Chla = rCon2;
    }
*/

    CalcPARimax();   // Calculate total incident light [uE/m2/s] to supply to the algae (NET) at noon
    // update PWMs based on PWMmax, time and illumination
    TimeCalc();
    PWMCalcSend(); // calc rPARi, bPARi & PWMs

    // Calculate and update L_avails once per period
    L_availCalc(rPARi, "R");
    L_availCalc(bPARi, "B");

    L_avail = rL_avail + bL_avail;
    UpdateIndicatorsFunction("L_avail");

    if (IsDay)
    {
      // Calculate current FrR
      FrR = rL_avail / L_avail;
    }
    else
    {
      FrR = -0.15;
    }
    UpdateIndicatorsFunction("FrR");

    // Mu section
    iStartMu++;
    if (iStartMu == 2)
    {
      ConMuFirst = Chla;
      TMuFirst = tor;
    }
    else if (iStartMu > 2)
    {
      TMu = tor - TMuFirst;
      Mu = 24.0 * log(Chla / ConMuFirst) / TMu;      // specific exponential growth rate in 1/d

      UpdateIndicatorsFunction("Mu");
      MessageMonitor("Screen", "####              Calculated Mu (1/d) ", String(Mu));
      UpdateIndicatorsFunction("TMu");
      MessageMonitor("Screen", "####                 For the last (h) ", String(TMu));
    }

    PrintCSV("continuous"); // Print datum to CSV File
    TLastCon = millis();
  }
  // ------------------------------ End of Con, PWMmax calculation -------------------------------

  // ------------------------------ Temp, pH, CO2 procedure --------------------------------------

  if ((millis() - TLastTemp) >= TTemp)    // read the temperature every TTemp millis
  {
    ReadTemp();
    TLastTemp = millis();
  }

  if ((millis() - TLastpH) >= TpH)
  {
    // prevent pHAvg update during Night illumination
    if (IsDay)
    {
      pHAvg = pHMin + pHRange / (1 + exp(-steepness * (Li - LiMiddle))); // Logistic function
    }
    else // night
    {
      pHAvg = 8.5;
    }

    pHHigh = pHAvg + 0.05;
    pHLow = pHAvg - 0.05;

    UpdateIndicatorsFunction("pHAvg");

    ReadPH();
    if (!IsValveOn)    // if the valve is OFF
    {
      if (pH >= pHHigh) // reached the max pH boundary. Open the CO2 Valve to reduce the pH.
      {
        delay(200);  // CA: added to slowdown communication (doesn't work very well)

        MessageMonitor("Screen", "###  Valve is OFF", "pH >= pH_High");
        digitalWrite(Relay_CO2, Relay_On);
        delay(200);
        MessageMonitor("Screen", "###  CO2", "valve is now ON"); // XXXXX MegunoLink freezes here XXXXXX

        COTP.SendData("CO", CO2);
        CO2 = 8;
        COTP.SendData("CO", CO2);
        IsValveOn = true;
        Panel.SetCheck("CO2", IsValveOn);
        TOpenCO2 = millis();             // record time of valve opening
      }
    }
    else             // if the valve is ON
    {
      if (pH <= pHLow)     // CloseValve, end CO2 cycleF.
      {
        delay(200);  // CA: added to slowdown communication (doesn't work very well)

        MessageMonitor("Screen", "## Valve is ON", "pH <= pH_Low");
        digitalWrite(Relay_CO2, Relay_Off);
        delay(200);
        MessageMonitor("Screen", "## CO2", "valve is now OFF");

        COTP.SendData("CO", CO2);
        CO2 = 7;
        COTP.SendData("CO", CO2);
        IsValveOn = false;
        Panel.SetCheck("CO2", IsValveOn); // Valve closing => the begining of a new cycle

        TCloseCO2 = TCurrentCO2;     // Record time of cycle start (Previous cycle end)
        TCurrentCO2 = millis();      // Record time of cycle end (Time still saved, see above line)

        CO2Fraction = ((double)(TCurrentCO2 - TOpenCO2) / (double)(TCurrentCO2 - TCloseCO2));
        // Fraction = ( Cycle End - valve opening time )  / ( Cycle end - Cycle start )
        // Percent = 100*(Time of valve open) / (Whole cycle time)

        CleanArray();
        dtostrf((TCurrentCO2 - TCloseCO2) / (60.0 * 1000.0), 1, 3, ToString);
        MessageMonitor("Screen", "CO2 Cycle ended. Cycle time [min]", ToString);

        CleanArray();
        dtostrf((TCurrentCO2 - TOpenCO2) / (60.0 * 1000.0), 1, 3, ToString);
        MessageMonitor("Screen", "Valve open time [min]", ToString);

        CleanArray();
        dtostrf((CO2Fraction * 100.0) , 1, 3, ToString);
        MessageMonitor("Screen", "% time with valve open", ToString);

        /*
          Mixing using bubbling air at 6000 ml/min (~0.35 vol air:vol culture/min)
          with 100 ml CO2/min (~1.67% in air) when relay is on)
         *** Calculation of molar flux of CO2 in the reactor when relay is on
          100 ml/min (*60min/h) = 6.0 L/h (/22.4 L/mol) = 0.268 mol/h (/20L) =
          0.0134 M/h (*1000) = 13.4 mM/h when relay is on ***
        */

        Jc = (100.0 * 60 / 22.4 / 20.0) * CO2Fraction; // CO2 flow rate (mM/h) in last cycle
        JcTP.SendFloatData("Jc", Jc, 3);

      }  // End if (pH <= pHLow)
    }  // End else Valve-is-on

    TLastpH = millis();
  }
  // ------------------------ END Temp, pH and CO2 Procedure --------------------
}
// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  END of loop() XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

// ************************** SUB: SetupPins() ************************************
void SetupPins() // Setup pins 2 & 3 to serve as input and interrupts.
{
  // Light In setup (Pin 3)
  pinMode(Lin_Pin, INPUT) ;
  digitalWrite(Lin_Pin, HIGH) ;
  attachInterrupt(1, LiFreq, RISING) ;

  // Light Out setup (Pin 2)
  pinMode(Lout_Pin, INPUT) ;
  digitalWrite(Lout_Pin, HIGH) ;
  attachInterrupt(0, LoFreq, RISING) ;

  //  Light Refl setup (Pin 21)
  //  pinMode(Lre_Pin, INPUT) ;
  //  digitalWrite(Lre_Pin, HIGH) ;
  //  attachInterrupt(2, LrFreq, RISING) ;


  pinMode(RedPin, OUTPUT);  // Declare the Red DIM pin as output
  pinMode(BluePin, OUTPUT); // Declare the Blue DIM pin as output
}
// *********************************************************************************

// ************************** SUB: SetupRelay() ************************************
void SetupRelay()
{
  pinMode(Relay_CO2, OUTPUT);         // Declare Relay pin as output
  digitalWrite(Relay_CO2, Relay_Off); // relay is  inactive at reset
}
// *********************************************************************************

// *************************** SUB: SetupTemp() ************************************
void SetupTemp()
{
  // Searching for the temp sensor
  if ( !ds.search(addr))
  {
    Serial.println("Didn't find the temp sensor");
    //ds.reset_search();
    delay(50);
    return;
  }
  else
  {
    Serial.println("found the temp sensor");
  }

  if (OneWire::crc8(addr, 7) != addr[7])
  {
    Serial.println("Error with transfering the sensor address");
    return;
  }
  ds.reset();
  ds.select(addr);
  //ds.write(0x7F); // 12-bit sampling resolution (default)

  ds.reset();
  ds.select(addr);
  ds.write(0x4E);          // write scratchpad (starts at byte 2)
  // note:  set high/low temp alarms by changing the next two values
  ds.write(0x4B);    // default value of TH reg (user byte 1)
  ds.write(0x46);    // default value of TL reg (user byte 2)
  // modify scratchpad register to set temp sampling RESOLUTION (http://www.netfluvia.org/code/ds18B20_test.pde)
  // uncomment one of the following
  // ds.write(0x7F);    // 12-bit sampling resolution (default)
  ds.write(0x5F);    // 11-bit
  // ds.write(0x3F);    // 10-bit
  // ds.write(0x1F);    // 9-bit

}
// **************************** End SetupTemp() **************************************

// **************************** SUB: SetupPH() ***************************************
void SetupPH()
{
  Serial2.print("E\r"); // take the pH Circuit out of continuous mode (stand by mode).
  delay(200);           // on start up sometimes the first command is missed.
  Serial2.print("E\r"); // so, let's send it twice.
  delay(200);           // a short delay after the pH Circuit was taken out of continuous mode
}
// **************************** End SetupPH() ****************************************

// ***************************** SUB: SetupMLP() *************************************
void SetupMLP()
{
  MessageMonitor("Screen", "## Initialize", "TimePlots & CommandHandler");
  // RPWM Time Plot
  RPWMTP.SetXlabel("Time");
  RPWMTP.SetYlabel("Value");
  RPWMTP.SetSeriesProperties(F("RPWM"), F("r#0*6")); // SetSeriesProperties("RPWM",Plot::Red,Plot::NoLine,0,Plot::Star);
  // BPWM Time Plot
  BPWMTP.SetXlabel("Time");
  BPWMTP.SetYlabel("Value");
  BPWMTP.SetSeriesProperties(F("BPWM"), F("b#0*6")); // SetSeriesProperties("BPWM",Plot::Blue,Plot::NoLine,0,Plot::Star);
  // FrR Time Plot
  FrRTP.SetXlabel("Time");
  FrRTP.SetYlabel("Value");
  FrRTP.SetSeriesProperties(F("FrR"), F("r#0*3")); // SetSeriesProperties("FrR",Plot::red,Plot::NoLine,0,Plot::Star);
  // Red Light available Time Plot
  rAvaTP.SetXlabel("Time");
  rAvaTP.SetYlabel("Value");
  rAvaTP.SetSeriesProperties(F("rAva"), F("r-1o8")); // SetSeriesProperties("rAva",Plot::Red,Plot::Dashed,1,Plot::Circle);
  // Blue Light available Time Plot
  bAvaTP.SetXlabel("Time");
  bAvaTP.SetYlabel("Value");
  bAvaTP.SetSeriesProperties(F("bAva"), F("b-1o8")); // SetSeriesProperties("Lo",Plot::Blue,Plot::Dashed,1,Plot::Circle);
  // Light In Time Plot
  LiTP.SetXlabel("Time");
  LiTP.SetYlabel("Value");
  LiTP.SetSeriesProperties(F("Li"), F("m#0o6")); // SetSeriesProperties("Li",Plot::Magenta,Plot::NoLine,1,Plot::Circle);
  // Light Out Time Plot
  LoTP.SetXlabel("Time");
  LoTP.SetYlabel("Value");
  LoTP.SetSeriesProperties(F("Lo"), F("c#0o6")); // SetSeriesProperties("Lo",Plot::Cyan,Plot::NoLine,1,Plot::Circle);
  // Light In Average Time Plot
  LiAvTP.SetXlabel("Time");
  LiAvTP.SetYlabel("Value");
  LiAvTP.SetSeriesProperties(F("LiAv"), F("m_1o10")); // SetSeriesProperties("LiAv",Plot::Magenta,Plot::Solid,1,Plot::Circle);
  // Light Out Average Time Plot
  LoAvTP.SetXlabel("Time");
  LoAvTP.SetYlabel("Value");
  LoAvTP.SetSeriesProperties(F("LoAv"), F("c_1o10")); // SetSeriesProperties("LoAv",Plot::Cyan,Plot::Solid,1,Plot::Circle);
  // Concentration from Red light
  ConRTP.SetXlabel("Time");
  ConRTP.SetYlabel("Value");
  ConRTP.SetSeriesProperties(F("ConR"), F("r_2s10")); // SetSeriesProperties("ConR",Plot::Red,Plot::Solid,2,Plot::Square);
  // Concentration from Blue light
  ConBTP.SetXlabel("Time");
  ConBTP.SetYlabel("Value");
  ConBTP.SetSeriesProperties(F("ConB"), F("b_1s8")); // SetSeriesProperties("ConB",Plot::Blue,Plot::Solid,1,Plot::Square);
  // bAbs/rAbs
  b_rAbsTP.SetXlabel("Time");
  b_rAbsTP.SetYlabel("Value");
  b_rAbsTP.SetSeriesProperties(F("b_rAbs"), F("m-1v8")); // SetSeriesProperties("b_rAbs",Plot::Magenta,Plot::Dashed,1,Plot::inverted triangle);
  // Car/Chla Time Plot
  Car_TP.SetXlabel("Time");
  Car_TP.SetYlabel("Value");
  Car_TP.SetSeriesProperties(F("Car_"), F("b-2x9")); // SetSeriesProperties("Car_",Plot::Blue,Plot::Dashed,2,Plot::Croix);
  // DW/Chla Time Plot
  DW_TP.SetXlabel("Time");
  DW_TP.SetYlabel("Value");
  DW_TP.SetSeriesProperties(F("DW_"), F("k-2+9")); // SetSeriesProperties("DW_",Plot::Black,Plot::Dashed,2,Plot::NoMarker);
  // PARi = Total (R+B) light supplied [uE/m2/s] Time Plot
  PARiTP.SetXlabel("Time");
  PARiTP.SetYlabel("Value");
  PARiTP.SetSeriesProperties(F("PARi"), F("m_1*5")); // SetSeriesProperties("PARi",Plot::Magenta,Plot::Solid,1,Plot::Star);
  // rPARi = Red light supplied [uE/m2/s] Time Plot   xxxxxxxxxxxxxxxxxx
  rPARiTP.SetXlabel("Time");
  rPARiTP.SetYlabel("Value");
  rPARiTP.SetSeriesProperties(F("rPARi"), F("r_1o5")); // SetSeriesProperties("PARi",Plot::Red,Plot::Solid,1,Plot::Circle);
  // bPARi = Blue light supplied [uE/m2/s] Time Plot   xxxxxxxxxxxxxxxxx
  bPARiTP.SetXlabel("Time");
  bPARiTP.SetYlabel("Value");
  bPARiTP.SetSeriesProperties(F("bPARi"), F("b_1o5")); // SetSeriesProperties("PARi",Plot::Blue,Plot::Solid,1,Plot::Circle);
  // rPARo = Red light outcoming [uE/m2/s] Time Plot   xxxxxxxxxxxxxxxxxx
  rPARoTP.SetXlabel("Time");
  rPARoTP.SetYlabel("Value");
  rPARoTP.SetSeriesProperties(F("rPARo"), F("r_1^7")); // SetSeriesProperties("PARi",Plot::Red,Plot::Solid,1,Plot::Triangle);
  // bPARo = Blue light outcoming [uE/m2/s] Time Plot   xxxxxxxxxxxxxxxxx
  bPARoTP.SetXlabel("Time");
  bPARoTP.SetYlabel("Value");
  bPARoTP.SetSeriesProperties(F("bPARo"), F("b_1^7")); // SetSeriesProperties("PARi",Plot::Blue,Plot::Solid,1,Plot::Triangle);
  // Rate Time Plot (Mu)
  MuTP.SetXlabel("Time");
  MuTP.SetYlabel("Value");
  MuTP.SetSeriesProperties(F("Mu"), F("g_1^8")); // SetSeriesProperties("Mu",Plot::Green,Plot::Solid,1,Plot::Triangle);
  // Rate Time Plot (SpConRate)
  SpRaTP.SetXlabel("Time");
  SpRaTP.SetYlabel("Value");
  SpRaTP.SetSeriesProperties(F("Sp"), F("g_1v8")); // SetSeriesProperties("Sp",Plot::Green,Plot::Solid,1,Plot::inverted triangle);
  // Temp Time Plot
  TempTP.SetXlabel("Time");
  TempTP.SetYlabel("Value");
  TempTP.SetSeriesProperties(F("Temp"), F("r#0^8"));  // SetSeriesProperties("Temp",Plot::Red,Plot::NoLine,0,Plot::Triangle);
  // pH Time Plot
  PHTP.SetXlabel("Time");
  PHTP.SetYlabel("Value");
  PHTP.SetSeriesProperties(F("pH"), F("g#0s3")); // SetSeriesProperties("pH",Plot::Green,Plot::NoLine,0,Plot::Square);
  // CO Time Plot
  COTP.SetXlabel("Time");
  COTP.SetYlabel("Value");
  COTP.SetSeriesProperties(F("CO"), F("m_1n")); // SetSeriesProperties("CO",Plot::Magenta,Plot::Solid,1,Plot::NoMarker);
  // Flux CO2 Time Plot (Jc)
  JcTP.SetXlabel("Time");
  JcTP.SetYlabel("Value");
  JcTP.SetSeriesProperties(F("Jc"), F("k_1^8")); // SetSeriesProperties("Jc",Plot::Magenta,Plot::Solid,1,Plot::Triangle);

  // Interface Panel Functions
  // invoked from MLP when sending the right command
  SerialCommandHandler.AddCommand(F("ResetMu"), Cmd_ResetMu);
  SerialCommandHandler.AddCommand(F("SetTimeDate"), Cmd_SetTimeDate);
  SerialCommandHandler.AddCommand(F("SetDate"), Cmd_SetDate);
  SerialCommandHandler.AddCommand(F("SetIllumination"), Cmd_SetIllumination);
  SerialCommandHandler.AddCommand(F("SetComment"), Cmd_SetComment);
  SerialCommandHandler.AddCommand(F("UpdateIndicators"), Cmd_UpdateIndicators);
  SerialCommandHandler.AddCommand(F("ClearIndicators"), Cmd_ClearIndicators);

  delay(300);
}
// *************************** End of SetupMLP() ***************************************

// *************************** SUB: PrintParams() *****************************************
void PrintParams()
{
  CleanArray();
  dtostrf(tod, 1, 3, ToString);
  MessageMonitor("Screen", "## Printing to Params file (in csv format). tod", ToString);
  // format: date, tod, Chlb_, Car_, DW_, rKd, bKd
  // body
  delay(5);
  MyParams.Begin();
  delay(5);

  Serial.print(Date);
  Serial.print(",");
  Serial.print(tod, 3);
  Serial.print(",");
  Serial.print(Chlb_);
  Serial.print(",");
  Serial.print(Car_);
  Serial.print(",");
  Serial.print(DW_);
  Serial.print(",");
  Serial.print(rKd);
  Serial.print(",");
  Serial.print(bKd);

  MyParams.End();
}
// *************************** End of PrintParams() ***************************************


// *************************** SUB: PrintCSV() *****************************************
void PrintCSV(String state) // CSV File setup
{
  CleanArray();
  dtostrf(tod, 1, 3, ToString);
  MessageMonitor("Screen", "## Printing to CSV Logfile. tod", ToString);
  delay(5);
  MyCSVMessage.Begin();
  delay(5);

  if (IsSetup)
  {
    Serial.println("Global parameters (Constants)");
    // Headers and values for Initial Parameters
    Serial.print("Date, tod, TCon, rPARMax, bPARMax, ConMax, rSpecLSupply, bSpecLSupply,");
    Serial.print("FrR_target, rPBR, bPBR, rEpsa, rEpsb, rEpsc, bEpsa, bEpsb, bEpsc, rLambda, ra1, ra0, ");
    Serial.print("bLambda, ba1, ba0, rNPARa1, bNPARa1, pHMin, width, BgLo, BgLi, ");
    Serial.println("rPWMmin, bPWMmin, daylength");
    Serial.print(Date);
    Serial.print(",");
    Serial.print(tod, 3);
    Serial.print(",");
    Serial.print(TCon);
    Serial.print(",");
    Serial.print(rPARMax, 2);
    Serial.print(",");
    Serial.print(bPARMax, 2);
    Serial.print(",");
    Serial.print(ConMax, 2);
    Serial.print(",");
    Serial.print(rSpecLSupply, 4);
    Serial.print(",");
    Serial.print(bSpecLSupply, 4);
    Serial.print(",");
    Serial.print(FrR_target, 3);
    Serial.print(",");
    Serial.print(rPBR, 4);
    Serial.print(",");
    Serial.print(bPBR, 4);
    Serial.print(",");
    Serial.print(rEpsa, 5);
    Serial.print(",");
    Serial.print(rEpsb, 5);
    Serial.print(",");
    Serial.print(rEpsc, 5);
    Serial.print(",");
    Serial.print(bEpsa, 5);
    Serial.print(",");
    Serial.print(bEpsb, 5);
    Serial.print(",");
    Serial.print(bEpsc, 5);
    Serial.print(",");
    Serial.print(rLambda, 4);
    Serial.print(",");
    Serial.print(ra1, 4);
    Serial.print(",");
    Serial.print(ra0, 4);
    Serial.print(",");
    Serial.print(bLambda, 4);
    Serial.print(",");
    Serial.print(ba1, 4);
    Serial.print(",");
    Serial.print(ba0, 4);
    Serial.print(",");
    Serial.print(rNPARa1, 4);
    Serial.print(",");
    Serial.print(bNPARa1, 4);
    Serial.print(",");
    Serial.print(pHMin);
    Serial.print(",");
    Serial.print(width, 1);
    Serial.print(",");
    Serial.print(BgLo, 3);
    Serial.print(",");
    Serial.print(BgLi, 3);
    Serial.print(",");
    Serial.print(rPWMmin);
    Serial.print(",");
    Serial.print(bPWMmin);
    Serial.print(",");
    Serial.print(daylength, 1);
    Serial.println();

    Serial.println("Continuous");
    // CSV file headers
    Serial.print("Date, tod, tor, Illum, state, rLi_ref, rLo_ref, bLi_ref, bLo_ref, ");
    Serial.print("rKd, bKd, b_rAbs, rCon2, bCon2, rCon, bCon, ");
    Serial.print("rPARimax, bPARimax, FrR, rL_avail, bL_avail, L_avail, pH, pHAvg, Temp, ");
    Serial.println("rPWM, bPWM, LiAvg, LoAvg, rLi, rLo, bLi, bLo, Mu, TMu, ConRate, SpConRate, comment");
  }
  else
  {
    Serial.print(Date);
    Serial.print(",");
    Serial.print(tod, 3);
    Serial.print(",");
    Serial.print(tor, 3);
    Serial.print(",");
    Serial.print(Illum);
    Serial.print(",");
    Serial.print(state); // continuous or periodic report (measurement)
    Serial.print(",");
    Serial.print(rLi_ref);
    Serial.print(",");
    Serial.print(rLo_ref, 3);
    Serial.print(",");
    Serial.print(bLi_ref);
    Serial.print(",");
    Serial.print(bLo_ref, 3);
    Serial.print(",");
    Serial.print(rKd, 6);
    Serial.print(",");
    Serial.print(bKd, 6);
    Serial.print(",");
    Serial.print(b_rAbs, 3);
    Serial.print(",");
    Serial.print(rCon2, 3);
    Serial.print(",");
    Serial.print(bCon2, 3);
    Serial.print(",");
    Serial.print(rCon, 3);
    Serial.print(",");
    Serial.print(bCon, 3);
    Serial.print(",");
    Serial.print(rPARimax);
    Serial.print(",");
    Serial.print(bPARimax);
    Serial.print(",");
    Serial.print(FrR, 3);
    Serial.print(",");
    Serial.print(rL_avail, 2);
    Serial.print(",");
    Serial.print(bL_avail, 2);
    Serial.print(",");
    Serial.print(L_avail, 2);
    Serial.print(",");
    Serial.print(pHData);
    Serial.print(",");
    Serial.print(pHAvg);
    Serial.print(",");
    Serial.print(Temp);
    Serial.print(",");
    Serial.print(rPWM);
    Serial.print(",");
    Serial.print(bPWM);
    Serial.print(",");
    Serial.print(LiAvg, 3);
    Serial.print(",");
    Serial.print(LoAvg, 3);
    Serial.print(",");
    Serial.print(rLi, 3);
    Serial.print(",");
    Serial.print(rLo, 3);
    Serial.print(",");
    Serial.print(bLi, 3);
    Serial.print(",");
    Serial.print(bLo, 3);
    Serial.print(",");
    Serial.print(Mu, 6);
    Serial.print(",");
    Serial.print(tor - TMuFirst, 3); // TMu
    Serial.print(",");
    Serial.print(ConRate, 6);
    Serial.print(",");
    Serial.print(SpConRate, 6);
    Serial.print(",");
    Serial.println(comment);
  }
  MyCSVMessage.End();
  comment = "";     // Reset comment
}
// ***************************** End PrintCSV() *********************************

// *********************** SUB: Read Temperature ********************************
double ReadTemp()
{
  ds.reset();
  ds.select(addr);
  ds.write(0x44, 1);        // start conversion, with parasite power on at the end

  delay(450);     // 375ms according to manufacturer

  TempPresent = ds.reset();
  ds.select(addr);
  ds.write(0xBE);         // Read Scratchpad

  for ( byte i = 0; i < 9; i++) TempData[i] = ds.read();         // we need 9 bytes, so we read 9 bytes!

  int16_t raw = (TempData[1] << 8) | TempData[0];

  byte cfg = (TempData[4] & 0x60);
  // at lower res, the low bits are undefined, so let's zero them
  if (cfg == 0x00)
  {
    raw = raw & ~7; // 9 bit resolution, 93.75 ms
  }
  else if (cfg == 0x20)
  {
    raw = raw & ~3; // 10 bit res, 187.5 ms
  }
  else if (cfg == 0x40)
  {
    raw = raw & ~1; // 11 bit res, 375 ms
  }
  // default is 12 bit resolution, 750 ms conversion time

  Temp = (float)raw / 16.0; // Temp is the temperature [celsius]

  TempTP.SendFloatData("Temp", Temp, 2);

  UpdateIndicatorsFunction("Temp");
  MessageMonitor("Screen", "  Temp", String(Temp));

  strncpy(TempAscii, ToString, 5); //copies string stamp to a new variable (char[])
  TempAscii[5] = '\r';              //adding <carriage return> at the end of string temp

  return Temp;
}
// *********************************** End ReadTemp() *****************************

// ********************************* SUB: Read pH *********************************
double ReadPH()
{
  for ( unsigned int i = 0; i < sizeof(pHData);  ++i )    // Cleaning the Temporary Array
    pHData[i] = (char)0;

  Serial2.print(TempAscii); // Send "tt.tt\r" it will set the temperature, but won't send the pH reading string.

  Serial2.print("R\r");     // Send the pH reading string.
  delay(2500);              // time needeed by the ph circuit to reply

  if (Serial2.available())  // if the pH Circuit has sent a character - read it.
  {
    received_from_pH = Serial2.readBytesUntil(13, pHData, 20); // we read the data sent from pH Circuit untill we see a <CR>.
    // We also count how many character have been received.
    pHData[received_from_pH] = 0;  // we add a 0 to the spot in the array just after the last character we received.
    // This will prevent us from transmiting incorrect data that may have been left in the buffer.
  }

  pH = (double)atof(pHData);
  UpdateIndicatorsFunction("pH");
  MessageMonitor("Screen", "   pH ", pHData);

  return pH;
}
// *************************** End ReadpH() ********************************

// ******************* SUB: InitLights(rPWM,bPWM) *********************
void InitLights(unsigned int redPWM, unsigned int bluePWM)   // Initial measurement of Li & Lo for first Con2 estimates
{
  double LiSumInit = 0.0;
  double LoSumInit = 0.0;

  double LightIn[10];
  double LightOut[10];

  PowerLEDs(redPWM, bluePWM);

  RPWMTP.SendData("RPWM", redPWM);
  BPWMTP.SendData("BPWM", bluePWM);

  CleanArray();
  String(redPWM).toCharArray(ToString, 10);
  MessageMonitor("Screen", "## Initial rPWM", ToString);

  CleanArray();
  String(bluePWM).toCharArray(ToString, 10);
  MessageMonitor("Screen", "## Initial bPWM", ToString);

  red_light_millis = 2000;
  blue_light_millis = 20000;

  for (unsigned int ii = 0; ii < 10 ; ii++)
  {
    unsigned int MillisCurrent = millis();
    pulses_in = 0;
    pulses_out = 0;
    // pulses_ref = 0;

    if ((redPWM == rPWMmin) || (redPWM == 0) )   // Blue light only
    {
      delay(blue_light_millis); // time to accumulate pulses
    }
    else
    {
      delay(red_light_millis); // time to accumulate pulses
    }

    unsigned int MillisDelayed = millis() - MillisCurrent;

    TempPulseI = pulses_in;  // Seize the number of pulses in counted untill now.
    TempPulseO = pulses_out; // Seize the number of pulses out counted untill now.
    //    TempPulseR = pulses_ref;  // Seize the number of pulses reflected counted untill now.

    UpdateIndicatorsFunction("TempPulseI");
    UpdateIndicatorsFunction("TempPulseO");
    // UpdateIndicatorsFunction("TempPulseR");

    CleanArray();
    String(MillisDelayed).toCharArray(ToString, 10);
    Panel.SetText("TLo", ToString);
    Panel.SetText("TLi", ToString);

    LightIn[ii] = TempPulseI / (MillisDelayed / 1000.0) - BgLi;          // !!! CA: bg subtraction
    LightOut[ii] = TempPulseO / (MillisDelayed / 1000.0) - BgLo;         // !!! CA: bg subtraction

    LiTP.SendFloatData("Li", LightIn[ii] / 1000.0 , 6);
    LoTP.SendFloatData("Lo", LightOut[ii], 5);

    LiSumInit = LiSumInit + LightIn[ii];
    LoSumInit = LoSumInit + LightOut[ii];

    // print results
    MessageMonitor("Screen", String(ii), "");

    CleanArray();
    dtostrf(LightIn[ii], 4, 2, ToString);
    Panel.SetText("Li", ToString);
    MessageMonitor("Screen", "                  Initial Li", ToString);

    CleanArray();
    dtostrf(LightOut[ii], 4, 2, ToString);
    Panel.SetText("Lo", ToString);
    MessageMonitor("Screen", "                  Initial Lo", ToString);
  }

  LiAvg = LiSumInit / 10; // Calculating Light in, small average in Hz.
  LoAvg = LoSumInit / 10; // Calculating Light out, small average in Hz.

  UpdateIndicatorsFunction("LiAvg");
  UpdateIndicatorsFunction("LoAvg");

  MessageMonitor("Screen", "## LiAvg = ", String(LiAvg));
  MessageMonitor("Screen", "## LoAvg = ", String(LoAvg));

  // calculate s.d of pulses
  double varSumI = 0;
  double varSumO = 0;

  for (unsigned int ii = 0; ii < 10 ; ii++)
  {
    varSumI += pow((LightIn[ii] - LiAvg), 2);
    varSumO += pow((LightOut[ii] - LoAvg), 2);
  }

  double varianceI = varSumI / (9) ;
  double varianceO = varSumO / (9) ;

  MyCSVMessage.Begin();
  delay(5);

  Serial.println("rPWM,bPWM,Li_avg,Lo_avg,Li_SD,Lo_SD");
  Serial.print(redPWM);
  Serial.print(",");
  Serial.print(bluePWM);
  Serial.print(",");
  Serial.print(LiAvg, 3);
  Serial.print(",");
  Serial.print(LoAvg, 3);
  Serial.print(",");
  Serial.print(pow(varianceI, 0.5), 2);
  Serial.print(",");
  Serial.print(pow(varianceO, 0.5), 6);

  MyCSVMessage.End();

  CleanArray();
  dtostrf(pow(varianceI, 0.5), 4, 2, ToString);
  MessageMonitor("Screen", "## S.D. for initial Light In", ToString);

  CleanArray();
  dtostrf(pow(varianceO, 0.5), 4, 2, ToString);
  MessageMonitor("Screen", "## S.D. for initial Light Out", ToString);
  Serial.println();
}
// **************************** End: InitLights() **********************************

// ********************************* SUB: TimeCalc() *******************************
// input:
// output -> tod, IsDay, Date (if changed)
void TimeCalc()
{
  if (IsSetup)
  {
    LPanel.CallCommand("GetTime"); // ask megunolink to send us the time (it will activate Cmd_SetTimeDate())
    delay(300);
    SerialCommandHandler.Process();    // Determines DATE and tod0 (in hours)
    MessageMonitor("Screen", "Setup, debugging for tod0", String(tod0));
    // determine whether it is day or night
    if ((tod0 >= daystart) && (tod0 <= (24.0 - daystart)))  // Nominal daytime
    {
      IsDay = true;
      MessageMonitor("Screen", "Setup, IsDay", "true");
    }
    else // Nominal nighttime
    {
      IsDay = false;
      MessageMonitor("Screen", "Setup, IsDay", "false");
    }

    UpdateIndicatorsFunction("IsDay");
  }

  // time-of-day = start + (millis/millisInHours)
  // fmod(double x, double y) => (x MOD y) in double numbers
  tod = fmod((double)(tod0 + double (millis() / 1000.0 / 60.0 / 60.0)), 24.0);
  tor = dayofRun * 24.0 + tod - tod0;

  if (tod < todOld)                 // new day, IsDay status.
  {
    delay(50);
    dayofRun++;
    DateCalc(); // activate the function in Megunolink that sends us the date

    todIP = 0.0;     // added to enable proper update in IP after midnight
    Panel.SetText("tod", "0.000");
  }

  // detect change in day or night condition
  if (((tod >= daystart) && (todOld <= daystart)) || ((tod >= (24.0 - daystart)) && (todOld <= 24.0 - daystart)))
  { // the sun has just rise or set
    if (!IsDay)
    {
      IsDay = true;
      MessageMonitor("Screen", "## Night ended", "New day starts");
      delay(100);                 // added to prevent another entrance to the "else" below at the next loop
    }
    else
    {
      IsDay = false;
      MessageMonitor("Screen", "## Day ended", "Night starts");
      delay(100);       // E.D: added to prevent another entrance to the "else" below 1 loop ahead.
    }

    UpdateIndicatorsFunction("IsDay");

    LiSu = 0.0;
    LoSu = 0.0;
    LiSum = 0.0;
    LoSum = 0.0;
    TimesLi = 0;
    TimesLo = 0;
    pulses_in = 0;
    pulses_out = 0;
    // pulses_ref = 0;
    TLastCon = millis();
  }

  //  error handling routine - 0.1 hour after day/night end, check if IsDay status is correct
  if ((!IsDay) && ((tod >= (daystart + 0.1)) && (tod <= (daystart + 0.12))))
    // if the sun didnt rise and its 0.1 hour past daytime
  {
    IsDay = true;
    MessageMonitor("Screen", "###Error handling routine###", "day didn't begin as expected");
    MessageMonitor("Screen", "                            ", "IsDay reset to 'true'");
    UpdateIndicatorsFunction("IsDay");
  }

  if ((IsDay) && ((tod >= (24 - daystart + 0.1)) && (tod <= (24 - daystart + 0.12))))
    // if the sun didnt set and its 0.1 hour past "night time"
  {
    IsDay = false;
    MessageMonitor("Screen", "###Error handling routine###", "night didn't begin as expected");
    MessageMonitor("Screen", "                            ", "IsDay reset to 'false'");
    UpdateIndicatorsFunction("IsDay");
  }

  if ((tod - todIP) >= 0.01)    // for display in IP only every 36 sec
  {
    UpdateIndicatorsFunction("tod");
    MessageMonitor("Screen", "##    TimeCalc(),    tod", ToString);

    CleanArray();
    dtostrf(tor, 4, 3, ToString);
    MessageMonitor("Screen", "##    TimeCalc(),    tor", ToString);

    todIP = tod;
  }

  if (IsDay)
  {
    Con2Period = 2.0;
  }
  else
  {
    Con2Period = 4.0;
  }
  // update TCon2 if needed
  if ((TCon2 != 0.0) && (TCon2 != (Con2Period - 0.1)))
  {
    TCon2 = Con2Period - 0.1;     // First time calculating Con2 in sketch
  }

  todOld = tod;
}
// ******************************* End of TimeCalc() ***************************

// ****************************** SUB: Update Date *****************************
void DateCalc()
{
  LPanel.CallCommand("GetDate"); // Ask megunolink to send us the date
  delay(500);
  SerialCommandHandler.Process();
}
// ****************************** End of DateCalc() ****************************

// **************************** SUB: PWMCalcSend() *****************************
// A function that calculates, sets the proper power to the LEDs, as well as L_avail and current FrR
// input: rPARimax, bPARimax (NET to supply to algae), Illum, rTransPBR, bTransPBR
// intermediate output1: rPARi, bPARi
// intermediate output2: rPWM, bPWM and send them to LEDs

void PWMCalcSend()
{
  // To detect changes in PWM's
  unsigned int rPWMTemp = rPWM;
  unsigned int bPWMTemp = bPWM;

  double rPARi_norm, rLi_norm;
  double bPARi_norm, bLi_norm;
  unsigned long TCon_old = TCon;             // for detecting changes in TCon

  switch (Illum)
  {
    case 1:              // Continuous light; corrected formula
      rPARi = rPARimax * 0.5 * daylength / 24.0;
      bPARi = bPARimax * 0.5 * daylength / 24.0;
      IsDay = true;
      Con2Period = 2.0;
      // update TCon2
      if (TCon2 != 0.0)
      {
        TCon2 = Con2Period - 0.1; 
      }
      TCon = 900000;  // 15 min
      break;

    case 2:              // Constant light / Dark; corrected formula
      if ( IsDay )       // if the sun is shining by tod and daylength
      {
        rPARi = rPARimax * 0.5;
        bPARi = bPARimax * 0.5;
        TCon = 900000;  // 15 min
      }
      else              // The sun has set.
      {
        rPARi = 0.0;
        bPARi = 0.0;
        TCon = 1800000;  // 30 min
      }
      break;

    case 3:             // Sinusoid daylight function
      if ( IsDay )      // daytime
      {
        radr = (2 * PI * (tod - daystart) / daylength);  // convert time of day to radian (needed for cos function)
        radb = (2 * PI * (tod - daystart) / daylength);  // convert time of day to radian (needed for cos function)
        TCon = 300000;  // 5 min
      }
      else             // nighttime
      {
        radr = 0.0;
        radb = 0.0;
        TCon = 1800000;  // 30 min
      }

      cosred = 0.5 * (1 - cos(radr));           // starts 0, 0.5, 1, 0.5, 0, 0.5, 1.0, ...
      cosblue = 0.5 * (1 - cos(radb));          // starts 0, 0.5, 1, 0.5, 0, 0.5, 1.0, ...

      rPARi = rPARimax * cosred;
      bPARi = bPARimax * cosblue;
      break;

    default:             // Dark
      rPARi = 0.0;
      bPARi = 0.0;
      IsDay = false;
      break;
  } // end switch(Illum)

  PARi = rPARi + bPARi; // for reporting only

  // Report PAR (to supply to algae) at given time and given Illumination mode
  if ((rPWM != rPWMTemp) || (bPWM != bPWMTemp) || (PerCon >= (TCon / 1000)) || IsSetup)
  {
    UpdateIndicatorsFunction("rPARi");
    UpdateIndicatorsFunction("bPARi");
    UpdateIndicatorsFunction("PARi");
  }

  // We have col_PARi - the net light to be supplied to the algae.
  // since the PBR reduces some PAR, we have to increase col_PARi to be supplied to the PBR

  // calculate PWM from rPARi and bPARi
  // take into account light supplied BRUTO to the PBR (=col_PARi/col_Xi)
  rPARi_norm = (rPARi / rTransPBR) / rPARi_const;
  bPARi_norm = (bPARi / bTransPBR) / bPARi_const;

  rLi_norm = rPARi_norm / rNPARa1;
  bLi_norm = bPARi_norm / bNPARa1;

  // NEED TO BE TESTED ... AGAIN!!!
  rPWM = (unsigned int)(constrain(rPWMmin + ((((pow(rLi_norm, rLambda) - 1) / rLambda) - ra0) / ra1) - 25, rPWMmin, PWMmax));
  bPWM = (unsigned int)(constrain(bPWMmin + ((((pow(bLi_norm, bLambda) - 1) / bLambda) - ba0) / ba1) - 25, bPWMmin, PWMmax));

  // plot and update only after change or if already at maximum once per TCon
  if ((rPWM != rPWMTemp) || ((rPWM == PWMmax) && (PerCon >= (TCon / 1000))))
  {
    UpdateIndicatorsFunction("rPWM");
    MessageMonitor("Screen", "                                   rPWM", String(rPWM));
  }

  if ((bPWM != bPWMTemp) || ((bPWM == PWMmax) && (PerCon >= (TCon / 1000))))
  {
    UpdateIndicatorsFunction("bPWM");
    MessageMonitor("Screen", "                                   bPWM", String(bPWM));
  }

  if (TCon != TCon_old)
  {
    UpdateIndicatorsFunction("TCon");
  }

  PowerLEDs(rPWM, bPWM);  // Send power to LEDs.
}
// ********************************* End PWMCalcSend() **********************************

// ******************************** SUB: PowerLEDs() ************************************
void PowerLEDs(unsigned int redPWM, unsigned int bluePWM)
{
  if ((redPWM != rPWM) || (bluePWM != bPWM))
    MessageMonitor("Screen", "PowerLEDs()", "");
  if (redPWM != rPWM)
  {
    rPWM = redPWM;
    UpdateIndicatorsFunction("rPWM");
    MessageMonitor("Screen", "                                   rPWM", ToString);
  }
  if (bluePWM != bPWM)
  {
    bPWM = bluePWM;
    UpdateIndicatorsFunction("bPWM");
    MessageMonitor("Screen", "                                   bPWM", ToString);
  }

  analogWrite(RedPin, redPWM);
  delay(50);
  analogWrite(BluePin, bluePWM);
  delay(50);
}
// ******************************* End of PowerLEDs() ***********************************

// **************************** SUB: LightIntervals() ***********************************
// adaptive TLi and TLo calculation; returns proper TLi and TLo
void  LightIntervals()
{
  if (IsDay)
  {
    TLi = constrain(0.1 * TLmax + int(5 * TLmax / Li), 0.2 * TLmax, TLmax);    // 2 - 10 [sec]
    TLo = constrain(0.5 * TLmax + int(4.2 * TLmax / Lo), 0.5 * TLmax, TLmax); // 10 - 20 [sec]
  }
  else
  {
    TLi = TLmax;
    TLo = TLmax;
  }
}
// ********************** End of LightIntervals() *************************************

void LoFreq() //Counts the pulses arrived from the sensor out.
{
  pulses_out++ ;
}

void LiFreq() //Counts the pulses arrived from the sensor in.
{
  pulses_in++ ;
}

//void LrFreq() //Counts the pulses arrived from the sensor refl.
//{
//  pulses_ref++ ;
//}

void Cmd_SetIllumination(CommandParameter & Parameter)
{
  Illum = Parameter.NextParameterAsInteger();


  switch (Illum)
  {
    case 1:
      IsDay = true;
      Con2Period = 2.0;
      // update TCon2
      if ((TCon2 != 0.0) && (TCon2 != (Con2Period - 0.1)))
      {
        TCon2 = Con2Period - 0.1;     // First time calculating Con2 in sketch
      }
      TCon = 900000;  // 15 min
      break;

    case 2:
      TCon = 900000;  // 15 min
      break;

    case 3:
      TCon = 300000;  // 5 min
      break;

    default:
      TCon = 1800000;  // 30 min
      break;
  }

  UpdateIndicatorsFunction("Illum");
  UpdateIndicatorsFunction("TCon");

  MeasureAbs2();                      // Illum changed
  TLastCon = millis();

  // Graphic indicator for Illum change
  MessageMonitor("Screen", "## Graphic Indicator for", "Illum Change");
  for (unsigned int ii = 0  ;  ii < 21; ii = ii + 1)
  {
    BPWMTP.SendData("BPWM", ii);
  }

  CleanArray();
  String(Illum).toCharArray(ToString, 10);
  MessageMonitor("Screen", "## Illumination mode changed to", ToString);
  comment = "Illumination mode changed";
}

// ******************** Sub: Cmd_SetComment - MEGUNOLINK Stuff *****************
void Cmd_SetComment(CommandParameter & Parameter)
{
  comment = Parameter.NextParameter();
  MessageMonitor("Screen", "## comment set to", comment);

  if (comment == "Mix")
  {
    MessageMonitor("Screen", "## Graphic Indicator for", "Mix");

    for (unsigned int ii = 0  ;  ii < 21; ii = ii + 2)
    {
      RPWMTP.SendData("RPWM", ii);
    }
    //  "reset PerCon", "Reset Mu", "Reset ref"
    MeasureAbs2();                // Reset everything (mix)

    CalcParameters();

    MessageMonitor("Screen", "'MIX'       Chla", String(Chla));

    ResetMuFunction();
    TLastCon = millis();         // Restart Con timer
  }

  if (comment.substring(0, 3) == "LAB") // The first 3 characters of the comment (format "LAB,Chla,Chlb,Car,DW")
  {
    UpdateIndicatorsFunction("rCon2");      // Draw point before
    UpdateIndicatorsFunction("bCon2");      // Draw point before
    UpdateIndicatorsFunction("Car_");       // Draw point before
    UpdateIndicatorsFunction("DW_");        // Draw point before

    MessageMonitor("Screen", "## Graphic Indicator for", "LAB");
    for (unsigned int ii = 0  ;  ii < 36; ii = ii + 1)
    {
      RPWMTP.SendData("RPWM", ii);
    }

    MessageMonitor("Screen", "Chla calc", "");
    comment = comment.substring(4);                     // xx.xx,yy.yy,zz.zz,kkk.k
    MessageMonitor("Screen", "comment_input", comment);
    comma_index = comment.indexOf(',');
    Chla_Lab = comment.substring(0, comma_index);
    MessageMonitor("Screen", "Chla received", Chla_Lab);

    MessageMonitor("Screen", "Chlb calc", "");
    comment = comment.substring(comma_index + 1); // yy.yy,zz.zz,kkk.k
    comma_index = comment.indexOf(',');
    Chlb_Lab = comment.substring(0, comma_index);
    MessageMonitor("Screen", "Chlb received", Chlb_Lab);

    MessageMonitor("Screen", "Car calc", "");
    comment = comment.substring(comma_index + 1); // zz.zz,kkk.k
    comma_index = comment.indexOf(',');
    Car_Lab = comment.substring(0, comma_index);
    MessageMonitor("Screen", "Car received", Car_Lab);

    DW_Lab = comment.substring(comma_index + 1); // kkk.k
    MessageMonitor("Screen", "DW received", DW_Lab);


    Chla = Chla_Lab.toFloat();
    UpdateIndicatorsFunction("Chla");
    MessageMonitor("Screen", "Measured Chla", String(Chla));

    rCon2 = Chla;
    UpdateIndicatorsFunction("rCon2");

    Chlb_ = Chlb_Lab.toFloat() / Chla;
    Car_ = Car_Lab.toFloat() / Chla;
    DW_ = DW_Lab.toFloat() / Chla;
    UpdateIndicatorsFunction("Car_");
    UpdateIndicatorsFunction("DW_");

    UpdatePigments(Chla);

    // "Reset ref", "Reset Mu", "reset PerCon"
    MeasureAbs2();     // Get rAbs2, bAbs2 and b_rAbs2

    /*
      // Relevant to the case of "LAB" - rKd is the only unknown (but ajustable) parameter.
      // Calculate and report rKd from new rAbs2 and Chla, and initial Car_ and DW_
    */

    // Calculate and report rKd from red parameters [rAbs2] and new [Chla, DW_]
    rKd = rAbs2 / (pow(Chla, 2) * (rEpsa + rEpsb * Chlb_) * width * DW_);
    UpdateIndicatorsFunction("rKd");
    MessageMonitor("Screen", "Calculated rKd", String(rKd, 9));

    // In this case (LAB), there is no need to update DW_ since we get its measured value

    UpdateSumepsOp(); // for updating Sumeps_; the update of bOp_ is irrelevant.

    // Process blue data
    // Recalculate bKd  using the condition for which bCon2 = rCon2: "bAbs2 / rAbs2 = bSumeps_*bOp_ /(rSumeps_*rOp_)"
    bKd = rKd * b_rAbs * rSumeps_ / bSumeps_;
    UpdateIndicatorsFunction("bKd");
    MessageMonitor("Screen", "Calculated bKd", String(bKd, 9));

    // write parameters to file as csv
    PrintParams();

    UpdateSumepsOp(); // Update also bOp_

    CalcPARimax();

    // save for further ConRate calculation
    Con2Old = Chla;
    TLastCon2 = millis();

    ResetMuFunction();
    TLastCon = millis();   // Restart Con timer


  }   // End comment = "LAB"
}
// ******************* End of Cmd_SetComment - MEGUNOLINK Stuff ******************

void Cmd_UpdateIndicators(CommandParameter & Parameter)
{
  UpdateIndicatorsFunction("ALL");
}

// ********************* SUB:UpdateIndicatorsFunction()****************************
void UpdateIndicatorsFunction(String var)
{
  CleanArray();
  if (var == "ALL")
  {
    MessageMonitor("Screen", "## Update all Indicators", "ControlPanel");
  }

  if ((var == "Illum") || (var == "ALL"))
  {
    switch (Illum)
    {
      case 1:
        Panel.SetText("Illumination", "Continuous");
        break;

      case 2:
        Panel.SetText("Illumination", "Const D/N");
        break;

      case 3:
        Panel.SetText("Illumination", "Sin D/N");
        break;

      default:
        Panel.SetText("Illumination", "Off");
        break;
    }
  }

  if ((var == "IsDay") || (var == "ALL"))
  {
    if (IsDay)
    {
      Panel.SetText("IsDay", "True");
    }
    else
    {
      Panel.SetText("IsDay", "False");
    }
  }

  if ((var == "CO2") || (var == "ALL"))   // ???????? What if CO2 valve should be closed ???
  {
    Panel.SetCheck("CO2", IsValveOn);
    COTP.SendData("CO", CO2);
  }

  if ((var == "ConMax") || (var == "ALL"))
  {
    dtostrf(ConMax, 4, 2, ToString);
    Panel.SetText("ConMax", ToString);
  }

  if ((var == "FrR_target") || (var == "ALL"))
  {
    dtostrf(FrR_target, 2, 3, ToString);
    Panel.SetText("FrR_target", ToString);
  }

  if ((var == "rSpecLSupply") || (var == "ALL"))
  {
    dtostrf(rSpecLSupply, 2, 4, ToString);
    Panel.SetText("rSpecLSupply", ToString);
  }

  if ((var == "bSpecLSupply") || (var == "ALL"))
  {
    dtostrf(bSpecLSupply, 2, 4, ToString);
    Panel.SetText("bSpecLSupply", ToString);
  }

  if ((var == "pHMin") || (var == "ALL"))
  {
    dtostrf(pHMin, 2, 3, ToString);
    Panel.SetText("pHMin", ToString);
  }

  if ((var == "BgLi") || (var == "ALL"))
  {
    dtostrf(BgLi, 4, 3, ToString);
    Panel.SetText("BgLi", ToString);
  }

  if ((var == "BgLo") || (var == "ALL"))
  {
    dtostrf(BgLo, 4, 3, ToString);
    Panel.SetText("BgLo", ToString);
  }

  //  if ((var == "BgLr") || (var == "ALL"))
  //  {
  //    dtostrf(BgLr, 4, 3, ToString);
  //    Panel.SetText("BgLr", ToString);
  //  }

  if ((var == "Temp") || (var == "ALL"))
  {
    dtostrf(Temp, 2, 3, ToString);
    Panel.SetText("Temp", ToString);
    TempTP.SendFloatData("Temp", Temp, 2);
  }

  if ((var == "pH") || (var == "ALL"))
  {
    dtostrf(pH, 2, 3, ToString);
    Panel.SetText("pH", ToString);
    PHTP.SendFloatData("pH", pH, 3);
  }

  if ((var == "pHAvg") || (var == "ALL"))
  {
    dtostrf(pHAvg, 2, 3, ToString);
    Panel.SetText("pHAvg", ToString);
  }

  if ((var == "tod") || (var == "ALL"))
  {
    dtostrf(tod, 4, 3, ToString);
    Panel.SetText("tod", ToString);
  }

  if ((var == "Date") || (var == "ALL"))
  {
    String(Date).toCharArray(ToString, 10);
    Panel.SetText("Date", ToString);
  }

  if ((var == "TCon") || (var == "ALL"))
  {
    String(TCon / 1000).toCharArray(ToString, 15);
    Panel.SetText("TCon", ToString);
  }

  if ((var == "rPARimax") || (var == "ALL"))
  {
    dtostrf(rPARimax, 4, 3, ToString);
    Panel.SetText("rPARimax", ToString);
  }

  if ((var == "bPARimax") || (var == "ALL"))
  {
    dtostrf(bPARimax, 4, 3, ToString);
    Panel.SetText("bPARimax", ToString);
  }

  if ((var == "PARi") || (var == "ALL"))
  {
    dtostrf(PARi, 4, 2, ToString);
    Panel.SetText("PARi", ToString);
    PARiTP.SendFloatData("PARi", PARi, 2);
  }

  if ((var == "rPARi") || (var == "ALL"))
  {
    dtostrf(rPARi, 4, 2, ToString);
    Panel.SetText("rPARi", ToString);
    rPARiTP.SendFloatData("rPARi", rPARi, 2);
  }

  if ((var == "bPARi") || (var == "ALL"))
  {
    dtostrf(bPARi, 4, 2, ToString);
    Panel.SetText("bPARi", ToString);
    bPARiTP.SendFloatData("bPARi", bPARi, 2);
  }

  if ((var == "rPARo") || (var == "ALL"))
  {
    rPARoTP.SendFloatData("rPARo", rPARo, 3);
  }

  if ((var == "bPARo") || (var == "ALL"))
  {
    bPARoTP.SendFloatData("bPARo", bPARo, 3);
  }

  if ((var == "rPWM") || (var == "ALL"))
  {
    String(rPWM).toCharArray(ToString, 10);
    Panel.SetText("rPWM", ToString);
    RPWMTP.SendData("RPWM", rPWM);
  }

  if ((var == "bPWM") || (var == "ALL"))
  {
    String(bPWM).toCharArray(ToString, 10);
    Panel.SetText("bPWM", ToString);
    BPWMTP.SendData("BPWM", bPWM);
  }

  if ((var == "rLi_ref") || (var == "ALL"))
  {
    dtostrf(rLi_ref, 4, 2, ToString);
    Panel.SetText("rLi_ref", ToString);
    LiAvTP.SendFloatData("LiAv", rLi_ref / 1000.0 , 6);
  }

  if ((var == "rLo_ref") || (var == "ALL"))
  {
    dtostrf(rLo_ref, 4, 3, ToString);
    Panel.SetText("rLo_ref", ToString);
    LoAvTP.SendFloatData("LoAv", rLo_ref, 3);
  }

  //  if ((var == "rLr_ref") || (var == "ALL"))
  //  {
  //    dtostrf(rLr_ref, 4, 3, ToString);
  //    Panel.SetText("rLr_ref", ToString);
  //    LoAvTP.SendFloatData("LrAv", rLr_ref, 3);
  //  }

  if ((var == "bLi_ref") || (var == "ALL"))
  {
    dtostrf(bLi_ref, 4, 2, ToString);
    Panel.SetText("bLi_ref", ToString);
    LiAvTP.SendFloatData("LiAv", bLi_ref / 1000.0 , 6);
  }

  if ((var == "bLo_ref") || (var == "ALL"))
  {
    dtostrf(bLo_ref, 4, 3, ToString);
    Panel.SetText("bLo_ref", ToString);
    LoAvTP.SendFloatData("LoAv", bLo_ref, 3);
  }

  //  if ((var == "bLr_ref") || (var == "ALL"))
  //  {
  //    dtostrf(bLr_ref, 4, 3, ToString);
  //    Panel.SetText("bLr_ref", ToString);
  //    LoAvTP.SendFloatData("LrAv", bLr_ref, 3);
  //  }

  if ((var == "rAbs2") || (var == "ALL"))
  {
    dtostrf(rAbs2, 2, 4, ToString);
    Panel.SetText("rAbs2", ToString);
  }

  if ((var == "bAbs2") || (var == "ALL"))
  {
    dtostrf(bAbs2, 2, 4, ToString);
    Panel.SetText("bAbs2", ToString);
  }

  if ((var == "b_rAbs") || (var == "ALL"))
  {
    dtostrf(b_rAbs, 2, 4, ToString);
    Panel.SetText("b_rAbs", ToString);
    b_rAbsTP.SendFloatData("b_rAbs", b_rAbs * 10.0 , 3);
  }

  if ((var == "FrR") || (var == "ALL"))
  {
    dtostrf(FrR, 2, 3, ToString);
    Panel.SetText("FrR", ToString);
    FrRTP.SendFloatData("FrR", FrR * 10.0 , 3);
  }

  if ((var == "L_avail") || (var == "ALL")) // !!! note that updating this variable require special string var -> "L_avail" !!!
  {
    dtostrf((rL_avail + bL_avail), 4, 3, ToString);
    Panel.SetText("L_avail", ToString);
  }

  if ((var == "rL_avail") || (var == "ALL"))
  {
    dtostrf(rL_avail, 4, 3, ToString);
    Panel.SetText("rL_avail", ToString);
    if (!IsSetup)
    {
      rAvaTP.SendFloatData("rAva", rL_avail , 3);
    }
  }

  if ((var == "bL_avail") || (var == "ALL"))
  {
    dtostrf(bL_avail, 4, 3, ToString);
    Panel.SetText("bL_avail", ToString);
    if (!IsSetup)
    {
      bAvaTP.SendFloatData("bAva", bL_avail , 3);
    }
  }

  if ((var == "Lo") || (var == "ALL"))
  {
    dtostrf(Lo, 4, 3, ToString);
    Panel.SetText("Lo", ToString);
    LoTP.SendFloatData("Lo", Lo, 6);
  }

  if ((var == "Li") || (var == "ALL"))
  {
    dtostrf(Li, 4, 2, ToString);
    Panel.SetText("Li", ToString);
    LiTP.SendFloatData("Li", Li / 1000.0 , 6);
  }

  if ((var == "LoAvg") || (var == "ALL"))
  {
    dtostrf(LoAvg, 4, 4, ToString);
    Panel.SetText("LoAvg" , ToString);
    LoAvTP.SendFloatData("LoAv", LoAvg, 3);
  }

  if ((var == "LiAvg") || (var == "ALL"))
  {
    dtostrf(LiAvg, 4, 2, ToString);
    Panel.SetText("LiAvg" , ToString);
    LiAvTP.SendFloatData("LiAv", LiAvg / 1000.0 , 6);
  }

  //  if ((var == "LrAvg") || (var == "ALL"))
  //  {
  //    dtostrf(LrAvg, 4, 3, ToString);
  //    Panel.SetText("LrAvg" , ToString);
  //    LiAvTP.SendFloatData("LrAv", LrAvg, 3);
  //  }

  if ((var == "TempPulseO") || (var == "ALL"))
  {
    String(TempPulseO).toCharArray(ToString, 10);
    Panel.SetText("TempPulseO", ToString);
  }

  if ((var == "TempPulseI") || (var == "ALL"))
  {
    String(TempPulseI).toCharArray(ToString, 10);
    Panel.SetText("TempPulseI", ToString);
  }

  //  if ((var == "TempPulseR") || (var == "ALL"))
  //  {
  //    String(TempPulseR).toCharArray(ToString, 10);
  //    Panel.SetText("TempPulseR", ToString);
  //  }

  if ((var == "TLo") || (var == "ALL"))
  {
    String(TCurrentLo - TLastLo).toCharArray(ToString, 10);
    Panel.SetText("TLo" , ToString);
  }

  if ((var == "TLi") || (var == "ALL"))
  {
    String(TCurrentLi - TLastLi).toCharArray(ToString, 10);
    Panel.SetText("TLi" , ToString);
  }
  //
  //  if ((var == "TLr") || (var == "ALL"))
  //  {
  //    String(TCurrentLr - TLastLr).toCharArray(ToString, 10);
  //    Panel.SetText("TLr" , ToString);
  //  }

  if ((var == "PerCon") || (var == "ALL")) // !!! note that updating this variable require special string var -> "PerCon" !!!
  {
    String((millis() - TLastCon) / 1000).toCharArray(ToString, 15); // $$$ consider distinction between setup and loop $$$ ///
    Panel.SetText("PerCon", ToString);
  }

  if ((var == "rCon") || (var == "ALL"))
  {
    dtostrf(rCon, 4, 3, ToString);
    Panel.SetText("rCon", ToString);
    //    ConRTP.SendFloatData("ConR", rCon, 3);
  }

  if ((var == "bCon") || (var == "ALL"))
  {
    dtostrf(bCon, 4, 3, ToString);
    Panel.SetText("bCon", ToString);
    //    ConBTP.SendFloatData("ConB", bCon, 3);
  }

  if ((var == "Chla") || (var == "ALL"))
  {
    dtostrf(Chla, 4, 3, ToString);
    Panel.SetText("Chla", ToString);
  }

  if ((var == "Chlb") || (var == "ALL"))
  {
    dtostrf(Chlb, 4, 3, ToString);
    Panel.SetText("Chlb", ToString);
  }

  if ((var == "Car") || (var == "ALL"))
  {
    dtostrf(Car, 4, 3, ToString);
    Panel.SetText("Car", ToString);
  }

  if ((var == "DW") || (var == "ALL"))
  {
    dtostrf(DW, 4, 2, ToString);
    Panel.SetText("DW", ToString);
  }

  if ((var == "Car_") || (var == "ALL"))
  {
    dtostrf(Car_, 4, 3, ToString);
    Panel.SetText("Car_", ToString);
    Car_TP.SendFloatData("Car_", 100.0 * Car_ , 2);
  }

  if ((var == "DW_") || (var == "ALL"))
  {
    dtostrf(DW_, 4, 2, ToString);
    Panel.SetText("DW_", ToString);
    DW_TP.SendFloatData("DW_", DW_ , 2);
  }

  if ((var == "rCon2") || (var == "ALL"))
  {
    dtostrf(rCon2, 3, 3, ToString);
    Panel.SetText("rCon2", ToString);
    ConRTP.SendFloatData("ConR", rCon2, 3);
  }

  if ((var == "bCon2") || (var == "ALL"))
  {
    dtostrf(bCon2, 3, 3, ToString);
    Panel.SetText("bCon2", ToString);
    ConBTP.SendFloatData("ConB", bCon2, 3);
  }

  if ((var == "Mu") || (var == "ALL"))
  {
    dtostrf(Mu, 2, 4, ToString);
    Panel.SetText("Mu", ToString);
    MuTP.SendFloatData("Mu", Mu, 4);
  }

  if ((var == "TMu") || (var == "ALL")) // !!! note that updating this variable require special string var -> "TMu" !!!
  {
    dtostrf(tor - TMuFirst, 4, 4, ToString);
    Panel.SetText("TMu", ToString);
  }

  if ((var == "ConRate") || (var == "ALL"))
  {
    dtostrf(ConRate, 4, 6, ToString);
    Panel.SetText("ConRate", ToString);
  }

  if ((var == "SpConRate") || (var == "ALL"))
  {
    dtostrf(SpConRate, 4, 6, ToString);
    Panel.SetText("SpConRate", ToString);
    SpRaTP.SendFloatData("Sp", SpConRate, 6);
  }

  if ((var == "rKd") || (var == "ALL"))
  {
    dtostrf(rKd, 1, 9, ToString);
    Panel.SetText("rKd", ToString);
  }

  if ((var == "bKd") || (var == "ALL"))
  {
    dtostrf(bKd, 1, 9, ToString);
    Panel.SetText("bKd", ToString);
  }
}
// ********************* END of UpdateIndicatorsFunction()*************************

void Cmd_ClearIndicators(CommandParameter & Parameter)
{
  ClearIndicatorsFunction();
}

// ********************* SUB: ClearIndicatorsFunction() ****************************
void ClearIndicatorsFunction()
{
  MessageMonitor("Screen", "## Clearing Indicators", "Control Panel");
  Panel.SetText("Illumination", "");
  Panel.SetText("IsDay", "");
  Panel.ClearCheck("CO2");
  Panel.SetText("ConMax", "");
  Panel.SetText("FrR_target", "");
  Panel.SetText("rSpecLSupply", "");
  Panel.SetText("bSpecLSupply", "");
  Panel.SetText("pHMin", "");
  Panel.SetText("BgLi", "");
  Panel.SetText("BgLo", "");
  Panel.SetText("Temp", "");
  Panel.SetText("pH", "");
  Panel.SetText("pHAvg", "");
  Panel.SetText("tod", "");
  Panel.SetText("Date", "");
  Panel.SetText("TCon", "");
  Panel.SetText("bPARimax", "");
  Panel.SetText("rPARimax", "");
  Panel.SetText("PARi", "");
  Panel.SetText("bPARi", "");
  Panel.SetText("rPARi", "");
  Panel.SetText("bPWM", "");
  Panel.SetText("rPWM", "");
  Panel.SetText("rLi_ref", "");
  Panel.SetText("rLo_ref", "");
  Panel.SetText("bLi_ref", "");
  Panel.SetText("bLo_ref", "");
  Panel.SetText("rAbs2", "");
  Panel.SetText("bAbs2", "");
  Panel.SetText("b_rAbs", "");
  Panel.SetText("FrR", "");
  Panel.SetText("L_avail", "");
  Panel.SetText("rL_avail", "");
  Panel.SetText("bL_avail", "");
  Panel.SetText("Lo", "");
  Panel.SetText("Li", "");
  Panel.SetText("LoAvg", "" );
  Panel.SetText("LiAvg", "" );
  Panel.SetText("TempPulseO", "");
  Panel.SetText("TempPulseI", "");
  Panel.SetText("TLo", "" );
  Panel.SetText("TLi", "" );
  Panel.SetText("PerCon", "");
  Panel.SetText("rCon", "");
  Panel.SetText("bCon", "");
  Panel.SetText("Chla", "");
  Panel.SetText("Chlb", "");
  Panel.SetText("Car", "");
  Panel.SetText("DW", "");
  Panel.SetText("Car_", "");
  Panel.SetText("DW_", "");
  Panel.SetText("rCon2", "");
  Panel.SetText("bCon2", "");
  Panel.SetText("Mu", "");
  Panel.SetText("TMu", "");
  Panel.SetText("ConRate", "");
  Panel.SetText("SpConRate", "");
  Panel.SetText("rKd", "");
  Panel.SetText("bKd", "");
}
// ********************* END of ClearIndicatorsFunction()*************************

// *************************** Setup MessageMonitor ******************************
void MessageMonitor(String channelName, String VarName , String message)
{
  Serial.print("{MESSAGE:");
  Serial.print(channelName);
  Serial.print("|data|");
  Serial.print(VarName);
  Serial.print(" : ");
  Serial.print(message);
  Serial.println("}");
}
// *********************** End Setup MessageMonitor ***************************

// *************************** Sub: CleanArray() ******************************
void CleanArray()
{
  for ( unsigned int i = 0; i < sizeof(ToString);  ++i )    // Cleaning the Temporary Array
    ToString[i] = (char)0;
}
// *************************** End: CleanArray() ******************************

// ************************ Sub: Cmd_SetTimeDate() ****************************
void Cmd_SetTimeDate(CommandParameter & Parameter)
{
  String TimeDate = Parameter.NextParameter(); // get the time in format: "HHmmss_ddMMyyyy"
  TimeDate.toCharArray(timeInput, 22);
  MessageMonitor("Screen", "Time received", timeInput);

  char *pch; // holds splited strings
  char *i;   // buffer

  //read untill the "_" , meaning read time in format "HHmmss" into pch
  pch = strtok_r(timeInput, "_", &i);

  //convert the hour to int
  unsigned int H1 = pch[0] - '0';
  unsigned int H2 = pch[1] - '0';
  unsigned int H = H1 * 10 + H2;

  //convert the minutes to int
  unsigned int m1 = pch[2] - '0';
  unsigned int m2 = pch[3] - '0';
  unsigned int m = m1 * 10 + m2;

  //convert the sec to int
  unsigned int s1 = pch[4] - '0';
  unsigned int s2 = pch[5] - '0';
  unsigned int s = s1 * 10 + s2;

  // combines them all together into a double who represents the time in hours
  tod0 = H + double(m / 60.0) + double(s / 3600.0) - double(millis() / 3600000.0);
  todOld = tod0;

  CleanArray();
  dtostrf(tod0, 4, 3, ToString);
  MessageMonitor("Screen", "tod0 Updated from MLP", ToString);

  //second iteration. read date
  pch = strtok_r(NULL, "_", &i);
  Date = String(pch);

  UpdateIndicatorsFunction("Date");
  MessageMonitor("Screen", "Date Updated from MLP", Date);

  delay(100);
}
// ************************* End: Cmd_SetTimeDate() ****************************

// *************************** Sub: Cmd_SetDate() ******************************
void Cmd_SetDate(CommandParameter & Parameter)
{
  String TimeDate = Parameter.NextParameter(); // get the time in format: "HHmmss_ddMMyyyy"
  TimeDate.toCharArray(timeInput, 22);
  MessageMonitor("Screen", "Time received", timeInput);

  char *pch; // holds splited strings
  char *i;   // buffer

  // Read untill the "_" , meaning read time in format "HHmmss" into pch
  pch = strtok_r(timeInput, "_", &i);

  // Read date
  pch = strtok_r(NULL, "_", &i);
  Date = String(pch);

  UpdateIndicatorsFunction("Date");
  MessageMonitor("Screen", "## New Day!, Date updated from MLP: ", ToString);

  delay(100);
}
// *************************** End: Cmd_SetDate() ******************************


// *************************** Sub: Cmd_ResetMu() ******************************
void Cmd_ResetMu(CommandParameter & Parameter)
{
  ResetMuFunction();
}
// *************************** End: Cmd_ResetMu() ******************************


// ************************** Sub: ResetMuFunction() ***************************
void ResetMuFunction()
{
  iStartMu = 0;
  Mu = 0.0;
  ConMuFirst = Chla;
  TMuFirst = tor;
  comment = "Reset_Mu";

  MessageMonitor("Screen", "## comment set to", comment);
  UpdateIndicatorsFunction("Mu");
  UpdateIndicatorsFunction("TMu");
}
// **************************** End: ResetMuFunction() ***********************

// ****************************** Sub: MeasureAbs2() *************************
// Short reference (140) measurement for Red and Blue;

void MeasureAbs2()
{
  MessageMonitor("Screen", "          MeasureAbs2()", "");

  UpdateIndicatorsFunction("rCon2");      // Draw point before
  UpdateIndicatorsFunction("bCon2");      // Draw point before
  UpdateIndicatorsFunction("Car_");       // Draw point before
  UpdateIndicatorsFunction("DW_");        // Draw point before

  // Measure Red lights
  MessageMonitor("Screen", "******** Light measurement", "Red");
  PowerLEDs(Con2PWM, 0);

  // reset variables
  unsigned long TCurrentCon2 = millis();
  pulses_in = 0;
  pulses_out = 0;
  //   pulses_ref = 0;

  // measure time
  delay(10000);

  // seize values
  unsigned long TDelayedCon2 = millis() - TCurrentCon2;
  TempPulseI = pulses_in;
  TempPulseO = pulses_out;
  //   TempPulseR = pulses_ref;

  // Close light and allow decay
  PowerLEDs(0, 0);
  pulses_in = 0;
  pulses_out = 0;
  delay(2000);

  // Calculate and update sensor readings at 140 PWM
  rLi_ref = TempPulseI / (TDelayedCon2 / 1000.0) - BgLi;   // bg subtraction
  rLo_ref = TempPulseO / (TDelayedCon2 / 1000.0) - BgLo;   // bg subtraction
  //   rLr_ref = TempPulseR / (TDelayedCon2 / 1000.0) - BgLr;   // bg subtraction

  UpdateIndicatorsFunction("rLi_ref");
  UpdateIndicatorsFunction("rLo_ref");
  //   UpdateIndicatorsFunction("rLr_ref");

  MessageMonitor("Screen", "Red Light Measure in Hz", "in = " + String(rLi_ref) + ", out = " + String(rLo_ref));

  // Convert Units from Hz to uE/s
  // Calculate PAR_in, PAR_out and PAR_refl
  rPARi_ref = rHz2PARi * rLi_ref;
  rPARo_ref = rHz2PARo * rLo_ref;
  //   rPARr_ref = rHz2PARr * rLr_ref;

  // calculate and report rAbs
  rAbs2 = -log10(rPARo_ref / rPARi_ref) - rPBR;
  UpdateIndicatorsFunction("rAbs2");

  // Measure Blue light
  MessageMonitor("Screen", "******** Light measurement", "Blue");

  PowerLEDs(0, Con2PWM);

  // reset variables
  TCurrentCon2 = millis();
  pulses_in = 0;
  pulses_out = 0;
  //   pulses_ref = 0;

  // time to measure con2blue (longer than for red because of very low light out)
  delay(20000);
  // seize values
  TDelayedCon2 = millis() - TCurrentCon2;
  TempPulseI = pulses_in;
  TempPulseO = pulses_out;
  //  TempPulseO = pulses_re;

  // close light and allow decay
  PowerLEDs(0, 0);
  pulses_in = 0;
  pulses_out = 0;
  delay(2000);

  // Calculate and update sensor readings [Hz] at 140 PWM
  bLi_ref = TempPulseI / (TDelayedCon2 / 1000.0) - BgLi;   // bg subtraction
  bLo_ref = TempPulseO / (TDelayedCon2 / 1000.0) - BgLo;   // bg subtraction
  //   bLr_ref = TempPulseR / (TDelayedCon2 / 1000.0) - BgLr;   // bg subtraction

  // Report measurement
  MessageMonitor("Screen", "Blue Light Measure in Hz", "in = " + String(bLi_ref) + ", out = " + String(bLo_ref));

  UpdateIndicatorsFunction("bLi_ref");
  UpdateIndicatorsFunction("bLo_ref");

  // Convert Units from Hz to uE/s
  bPARi_ref = bHz2PARi * bLi_ref;
  bPARo_ref = bHz2PARo * bLo_ref;
  //   bPARr_ref = bHz2PARr * bLr_ref;

  // calculate and report bAbs from PAR ratio
  bAbs2 = -log10(bPARo_ref / bPARi_ref) - bPBR;
  UpdateIndicatorsFunction("bAbs2");

  b_rAbs = bAbs2 / rAbs2;
  UpdateIndicatorsFunction("b_rAbs");

  // reset counters
  LiSu = 0.0;
  LoSu = 0.0;
  LiSum = 0.0;
  LoSum = 0.0;
  TimesLi = 0;
  TimesLo = 0;
  TLastLo = millis(); // Time at the end of each single Lo measurement
  TLastLi = millis(); // Time at the end of each single Li measurement
  TLastCon = millis();
}
// ****************************** End of MeasureAbs2() ********************************

// ********************************** SUB: CalcPARimax() ******************************
// input: Chla
// output -> rPARimax, bPARimax
void CalcPARimax()
{
  // Calculate total incident light [uE/m2/s] to supply to the algae (NET)
  rPARimax = rSpecLSupply * pow(Chla, 3) / ( 1 - exp(-2.30285 * rSAC * pow(Chla, 2) * width));
  // bPARimax calculated from rPARimax with appropriate weighting factors (epsilons, Kd, FrR_target)
  bPARimax = rPARimax * (bSAC / rSAC) * ((1 - FrR_target) / FrR_target) * ((1 - exp(-2.30285 * rSAC * pow(Chla, 2) * width)) / (1 - exp(-2.30285 * bSAC * pow(Chla, 2) * width)));



  UpdateIndicatorsFunction("rPARimax");
  MessageMonitor("Screen", "                                  rPARimax", ToString);

  UpdateIndicatorsFunction("bPARimax");
  MessageMonitor("Screen", "                                  bPARimax", ToString);
}
// ********************************** END: CalcPARimax() ***************************

// *********************** SUB: ConCalc(Abs,color) *********************************
// input: corr Abs and color
// output: Concentration after filter and sensitivity correction
double ConCalc(double Abs, String color)
{
  double epsa, epsb, epsc, Kd, sumeps_, Op_, concalc;
  MessageMonitor("Screen", "          ConCalc()", color);  // ED: added color report

  UpdateSumepsOp();

  if (color == "R")
  {
    sumeps_ = rSumeps_;
    Op_ = rOp_;
  }
  else if (color == "B")
  {
    sumeps_ = bSumeps_;
    Op_ = bOp_;
  }

  concalc = pow((Abs / (sumeps_ * Op_ * width)), 0.5);
  return concalc;
}
// ***************************** End of ConCalc(Abs,color) *****************************

// ******************************* Sub: CalcParameters() *******************************
// inputs: rAbs2, bAbs2, rKd, bKd, Chlb_
// outputs: rCon2, Chla, DW_, bCon2, Car_, Chlb_
// Calculation at SETUP(), "MIX" and PERIODIC;  NOT at "LAB"
// After MeasureAbs2, or InitLights(R and B)

void CalcParameters()
{
  // Calc rCon2
  rCon2 = ConCalc(rAbs2, "R");
  UpdateIndicatorsFunction("rCon2");

  Chla = rCon2;
  UpdateIndicatorsFunction("Chla");

  if (IsSetup == true)
  {
    MessageMonitor("Screen", "Initial Chla", String(Chla));
  }
  else
  {
    MessageMonitor("Screen", "Periodic Chla", String(Chla));
  }

  /*        !!!!!!!!!!!!! CA: rKd considered constant, except for "LAB"
    // Calculating rKd from RED parameters ONLY; independent of Car_
    rKd = rAbs2 / (pow(Chla, 2) * (rEpsa + rEpsb * Chlb_) * width * DW_);
    UpdateIndicatorsFunction("rKd");
  */

  // Calculate DW_ from RED parameters ONLY; independent of Car_
  DW_ = rAbs2 / (rKd * pow(Chla, 2) * width * (rEpsa + rEpsb * Chlb_));
  UpdateIndicatorsFunction("DW_");
  MessageMonitor("Screen", "DW_ based on red abs", String(DW_));

  UpdatePigments(Chla);

  // Calc bCon2
  bCon2 = ConCalc(bAbs2, "B");
  UpdateIndicatorsFunction("bCon2");

  /*        !!!!!!!!!!!!! CA: bKd considered constant, except for "LAB"
    // Recalculate bKd  using the condition for which bCon2 = rCon2: "bAbs2 / rAbs2 = bSumeps_*bOp_ /(rSumeps_*rOp_)"
    bKd = rKd * b_rAbs * rSumeps_ / bSumeps_;
    UpdateIndicatorsFunction("bKd");
  */

  // Calculate Car_ from ratio b/rAbs2 and r/bKd
  Car_ = (b_rAbs * (rKd / bKd) * (rEpsa + rEpsb * Chlb_) - (bEpsa + bEpsb * Chlb_)) / (bEpsc - b_rAbs * (rKd / bKd) * rEpsc);
  UpdateIndicatorsFunction("Car_");

  // Calculate DW_ from blue abs, bKd and new Car_; for checking purposes
  DW_ = bAbs2 / (bKd * pow(Chla, 2) * width * (bEpsa + bEpsb * Chlb_ + bEpsc * Car_));
  UpdateIndicatorsFunction("DW_");
  MessageMonitor("Screen", "DW_ based on blue Abs", String(DW_));

  UpdatePigments(Chla);

  Chlb_ = Chlb / Chla; // Methodic update of Chlb_
  UpdateIndicatorsFunction("Chlb_");
}
// ******************************* End of CalcParameters() *******************************

// ******************************* Sub: UpdatePigments(Chla) *******************************
// Update and report pigments
void UpdatePigments(double chla)
{
  Chlb = chla * Chlb_;
  Car = chla * Car_;
  DW = chla * DW_;

  UpdateIndicatorsFunction("Chlb");
  UpdateIndicatorsFunction("Car");
  UpdateIndicatorsFunction("DW");
}
// ******************************* End of UpdatePigments(Chla) *******************************

// ******************************* Sub: UpdateSumepsOp() *******************************
// Update Sumeps_ and Op_, based on current Car_, DW_ and Kd
void UpdateSumepsOp()
{
  rSumeps_ = rEpsa + rEpsb * Chlb_;                // independent of Chla and Car_
  rOp_ = rKd * DW_;                                // independent of Chla and width
  bSumeps_ = bEpsa + bEpsb * Chlb_ + bEpsc * Car_; // independent of Chla
  bOp_ = bKd * DW_;                                // independent of Chla and width
}
// ******************************* End of UpdateSumepsOp() *******************************


// *********************** SUB: L_availCalc(Li, color) **********************************
// input: PARin [uE] for color at time , color
// this represents the area under the curve {Lx vs x} from zero to width (= <L>)

void L_availCalc(double PARi, String color)
{
  double epsa, epsb, epsc, Kd, temp_pow_, l_avail;
  MessageMonitor("Screen", "          L_availCalc()", "");

  if (color == "R")
  {
    epsa = rEpsa;
    epsb = rEpsb;
    epsc = rEpsc;
    Kd = rKd;
  }
  else if (color == "B")
  {
    epsa = bEpsa;
    epsb = bEpsb;
    epsc = bEpsc;
    Kd = bKd;
  }

  temp_pow_ = 2.30285 * (epsa + epsb * Chlb_ + epsc * Car_) * Kd * DW_;
  // 2.30285 = ln(10) needed for integration (transf 'log10' into 'ln')

  // Light under the curve at PBR profile [uE] (= <PAR>)
  l_avail = PARi * ( 1 - exp(-temp_pow_  * pow(Chla, 2) * width) ) / (temp_pow_ * pow(Chla, 2)) ;

  // Save and plot Light variables
  if (color == "R")
  {
    rL_avail = l_avail;
    UpdateIndicatorsFunction("rL_avail");
    MessageMonitor("Screen", "          Red L_availCalc()", String(rL_avail));
  }
  else if (color == "B")
  {
    bL_avail = l_avail;
    UpdateIndicatorsFunction("bL_avail");
    MessageMonitor("Screen", "          Blue L_availCalc()", String(bL_avail));
  }
}
// *************************** End of L_availCalc() ******************************


/*
  // *********************** SUB: AdjustParameters() *********** #Under work# **********
  void AdjustParameters()
  {

  }
  // ************************** End of AdjustParameters() ******************************


  // testing purpose
  long LoopCounter = 0;
   LoopCounter++;
   if (LoopCounter % 20 == 0)
   {
     MessageMonitor("Screen", "Loop Counter", String(LoopCounter));
     if (LoopCounter % 100 == 0)
     {
       MeasureAbs2(false);
     }
   }

  /////////////////////////////////////////////////////
   TCurrentParameter  = millis();
   if (TCurrentParameter - TLastParameter >= TParameter) want to change parameter, lets say 24 hours.
   {
     rnd = random(1, 3);
     switch (rnd)
     {
       case 1:
         change light intensity per cell
         RelLiPow = random(17.75,27.75);
         break;
       case 2:
         change pH average
         break;
       case 3:
         change red / blue proportion
          FrR_target = random(0,1);
         break;
     }
     iStartMu = 0;
     TLastParameter = millis();
   }
  //////////////////////////////////////////////////////
*/

/*
  To call external C function from a sketch in arduino
  http://arduino.stackexchange.com/questions/946/how-to-call-c-functions-from-arduino-sketch
  in addition, one may use an arduino library to simulate thread computing:
  https://github.com/ivanseidel/ArduinoThread
*/


