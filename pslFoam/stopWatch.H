#include <ctime>

class StopWatch {
public:
  StopWatch() : totalTime{0}, count{0} {};
  void start() { startTime = std::clock(); }
  void stop() {
    count++;
    totalTime += (std::clock() - startTime) / (double)CLOCKS_PER_SEC;
  }
  double getTotalTime() { return totalTime; }
  int getCount() { return count; }

private:
  double totalTime;
  std::clock_t startTime;
  int count;
};
