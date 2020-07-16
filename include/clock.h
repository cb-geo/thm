#ifndef THM_CLOCK_H_
#define THM_CLOCK_H_

#include <chrono>
#include <iostream>

namespace cbgeo {
class Clock;
}

//! \brief A Tick-Tock clock
class cbgeo::Clock {
 public:
  //! Constructor
  Clock() {
    tick_ = std::chrono::system_clock::now();
    tock_ = std::chrono::system_clock::now();
  }

  //! Tick
  void tick() { tick_ = std::chrono::system_clock::now(); }

  //! Tock
  void tock(const std::string& message) {
    tock_ = std::chrono::system_clock::now();
    // Duration
    auto duration =
        std::chrono::duration_cast<std::chrono::milliseconds>(tock_ - tick_)
            .count();
    std::cout << "\n" << message << ": " << duration << " ms";
  }

 private:
  //! Tick and Tock
  std::chrono::time_point<std::chrono::system_clock> tick_, tock_;
};

#endif
