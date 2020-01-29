#include <Keyboard.h>
#include <lmic.h>
#include <hal/hal.h>
#include <SPI.h>

#define DISABLE_INVERT_IQ_ON_RX
#define CFG_eu868 1

#if !defined(DISABLE_INVERT_IQ_ON_RX)
#error This example requires DISABLE_INVERT_IQ_ON_RX to be set. Update \
       config.h in the lmic library to set it.
#endif

// How often to send a packet. Note that this sketch bypasses the normal
// LMIC duty cycle limiting, so when you change anything in this sketch
// (payload length, frequency, spreading factor), be sure to check if
// this interval should not also be increased.
// See this spreadsheet for an easy airtime and duty cycle calculator:
// https://docs.google.com/spreadsheets/d/1voGAtQAjC1qBmaVuP1ApNKs1ekgUjavHuVQIXyYSvNc
#define TX_INTERVAL 4000

// Pin mapping
const lmic_pinmap lmic_pins = {
    .nss = 10,
    .rxtx = LMIC_UNUSED_PIN,
    .rst = 9,
    .dio = {2, 6, 7},
};

// These callbacks are only used in over-the-air activation, so they are
// left empty here (we cannot leave them out completely unless
// DISABLE_JOIN is set in config.h, otherwise the linker will complain).
void os_getArtEui(u1_t *buf) {}
void os_getDevEui(u1_t *buf) {}
void os_getDevKey(u1_t *buf) {}

void onEvent(ev_t ev)
{
}

osjob_t txjob;
osjob_t timeoutjob;
// static void tx_func (osjob_t* job);
static void my_tx_func(osjob_t *job);

int currentPacketIndex = 0;
const int numOfPackets = 3;
char *myPackets[numOfPackets] = {
    "This is string packet 1ACK", "This is packet 2ACK", "This is packet 3",
    //  "This is packet 4", "This is packet 5ACK", "This is packet 6ACK"
};

// Transmit the given string and call the given function afterwards
void tx(const char *str, osjobcb_t func)
{
  LMIC.datarate = DR_SF9;
   LMIC.rps = updr2rps(LMIC.datarate);
  os_radio(RADIO_RST); // Stop RX first
  delay(1);            // Wait a bit, without this os_radio below asserts, apparently because the state hasn't changed yet
  LMIC.dataLen = 0;
  while (*str)
    LMIC.frame[LMIC.dataLen++] = *str++;
  LMIC.osjob.func = func;
  os_radio(RADIO_TX);
  Serial.println("TX");
}

// Enable rx mode and call func when a packet is received
void 
rx(osjobcb_t func)
{
  LMIC.datarate = DR_SF7;
  LMIC.rps = updr2rps(LMIC.datarate);
  delay(1);
  LMIC.osjob.func = func;
  LMIC.rxtime = os_getTime(); // RX _now_
  // Enable "continuous" RX (e.g. without a timeout, still stops after
  // receiving a packet)
  os_radio(RADIO_RXON);
  Serial.println("RX");
}

// static void rxtimeout_func(osjob_t *job) {
//   digitalWrite(LED_BUILTIN, LOW); // off
// }

static void my_rx_func(osjob_t *job)
{
  Serial.print("Got ");
  Serial.print(LMIC.dataLen);
  Serial.println(" bytes");
  Serial.write(LMIC.frame, LMIC.dataLen);
  Serial.println();
  Serial.println("In HEX ");
  for (int i = 0; i < LMIC.dataLen; i++)
  Serial.print(LMIC.frame[i], HEX);
  Serial.println();

  if (LMIC.dataLen == 3)
  {
    Serial.println("Got ACK");
    // if we get our ACK, start with next transmission, reschedules transmission at half TX_INTERVAL
    currentPacketIndex++;
    os_setTimedCallback(&txjob, os_getTime() + ms2osticks(TX_INTERVAL / 2), my_tx_func);
  }
  else
  {
    Serial.println("NOT AN ACK");
    // resend packet if no ACK received within 3*TX_INTERVAL, reschedules transmission in 3* TX_INTERVAL
    os_setTimedCallback(&txjob, os_getTime() + ms2osticks(3 * TX_INTERVAL), my_tx_func);

    // listen again
    rx(my_rx_func);
  }
}

// static void rx_func (osjob_t* job) {
//   // Blink once to confirm reception and then keep the led on
//   digitalWrite(LED_BUILTIN, LOW); // off
//   delay(10);
//   digitalWrite(LED_BUILTIN, HIGH); // on

//   // Timeout RX (i.e. update led status) after 3 periods without RX
//   os_setTimedCallback(&timeoutjob, os_getTime() + ms2osticks(3*TX_INTERVAL), rxtimeout_func);

//   // Reschedule TX so that it should not collide with the other side's
//   // next TX
//   os_setTimedCallback(&txjob, os_getTime() + ms2osticks(TX_INTERVAL/2), tx_func);

//   Serial.print("Got ");
//   Serial.print(LMIC.dataLen);
//   Serial.println(" bytes");
//   Serial.write(LMIC.frame, LMIC.dataLen);
//   Serial.println();

//   // Restart RX
//   rx(rx_func);
// }

// static void txdone_func (osjob_t* job) {
//   rx(rx_func);
// }

static void my_txdone_func(osjob_t *job)
{
  rx(my_rx_func);
}

static void my_txdone_no_ack_func(osjob_t *job)
{
    Serial.println("no ack func");
  currentPacketIndex++;
}

static void my_tx_func(osjob_t *job)
{
  
  if (currentPacketIndex < numOfPackets)
  {
    if (strlen(myPackets[currentPacketIndex]) > 3)
    {
      int length = strlen(myPackets[currentPacketIndex]);
      char lastThree[3];
      memcpy(lastThree, &myPackets[currentPacketIndex][length - 3], 3);
      const char ack[] = {'A', 'C', 'K'};
      if (!memcmp(lastThree, ack, 3))
      {
        //std::cout << "match" << std::endl;
        // send and start rx for receiving ACK
        Serial.print("transmitting packet with ACK, packet: ");
        Serial.println(currentPacketIndex + 1);
        tx(myPackets[currentPacketIndex], my_txdone_func);
      }
      else
      {
        //std::cout << "no match" << std::endl;
        // send and schedule next packet
        Serial.print("transmitting packet without ACK, packet: ");
        Serial.println(currentPacketIndex + 1);
        tx(myPackets[currentPacketIndex], my_txdone_no_ack_func);
      }
      Serial.print("here");
      os_setTimedCallback(&txjob, os_getTime() + ms2osticks(TX_INTERVAL), my_tx_func);
    }
    else
    {
        /* TODO 
        packet has less than 3 characters, unsafe to compare last three character for ACK
        send packet without ACK check  */
     
    }
  }
  else
  {
     Serial.println("No more packets to send, done");
  }
}

// // log text to USART and toggle LED
// static void tx_func (osjob_t* job) {
//   // say hello
//   tx("Hello World!", txdone_func);
//   // reschedule job every TX_INTERVAL (plus a bit of random to prevent
//   // systematic collisions), unless packets are received, then rx_func
//   // will reschedule at half this time.
//   os_setTimedCallback(job, os_getTime() + ms2osticks(TX_INTERVAL + random(500)), tx_func);
// }

// application entry point
void setup()
{
  Serial.begin(9600);
  Serial.println("Starting");

  LMIC_setClockError(MAX_CLOCK_ERROR * 50 / 100);

#ifdef VCC_ENABLE
  // For Pinoccio Scout boards
  pinMode(VCC_ENABLE, OUTPUT);
  digitalWrite(VCC_ENABLE, HIGH);
  delay(1000);
#endif

  pinMode(LED_BUILTIN, OUTPUT);

  // initialize runtime env
  os_init();

  // Set up these settings once, and use them for both TX and RX

#if defined(CFG_eu868)
  // Use a frequency in the g3 which allows 10% duty cycling.
  LMIC.freq = 868500000;
#elif defined(CFG_us915)
  LMIC.freq = 902300000;
#endif

  // Maximum TX power
  LMIC.txpow = 27;
  // Use a medium spread factor. This can be increased up to SF12 for
  // better range, but then the interval should be (significantly)
  // lowered to comply with duty cycle limits as well.
  LMIC.datarate = DR_SF9; 
  // This sets CR 4/5, BW125 (except for DR_SF7B, which uses BW250)
  LMIC.rps = updr2rps(LMIC.datarate);

  Serial.println("Started");
  Serial.flush();

  // setup initial job
  os_setCallback(&txjob, my_tx_func);
}

void loop()
{
  // execute scheduled jobs and events
  os_runloop_once();
}

