#pragma once

#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iostream>

/**
 * Luminosity mask (constructed from a list of good lumi block ranges)
 */
class LumiMask {
public:
  using Run = unsigned int;
  using LumiBlock = unsigned int;

  /**
   * Luminosity block range - helper struct
   */
  class LumiBlockRange {
  public:
    LumiBlockRange(Run run, LumiBlock firstLumi, LumiBlock lastLumi)
      : m_run(run), m_firstLumi(firstLumi),
        m_lastLumi(lastLumi ? lastLumi : std::numeric_limits<LumiBlock>::max())
    {}

    Run run() const { return m_run; }
    LumiBlock firstLumi() const { return m_firstLumi; }
    LumiBlock lastLumi() const { return m_lastLumi; }

    bool operator<(const LumiBlockRange& other) const {
      return (m_run == other.m_run)
              ? (m_lastLumi < other.m_firstLumi)
              : m_run < other.m_run;
    }

  private:
    Run m_run;
    LumiBlock m_firstLumi;
    LumiBlock m_lastLumi;
  };

  explicit LumiMask(const std::vector<LumiBlockRange>& accept)
    : m_accept(accept)
  {
    std::sort(m_accept.begin(), m_accept.end());
  }

  bool accept(Run run, LumiBlock lumi) const {
    return std::binary_search(m_accept.begin(), m_accept.end(), LumiBlockRange(run, lumi, lumi));
  }

  // Example input string: "12345:1-10,20-30;67890:5-15"
  static LumiMask fromCustomString(const std::string& input, LumiMask::Run firstRun=0, LumiMask::Run lastRun=0) {
    const bool noRunFilter = (firstRun == 0) && (lastRun == 0);
    std::vector<LumiBlockRange> accept;

    std::istringstream runStream(input);
    std::string runSegment;

    // Split input by semicolon to get each run's data
    while (std::getline(runStream, runSegment, ';')) {
      std::istringstream segmentStream(runSegment);
      std::string runPart;
      std::string lumiPart;

      // Split each segment into run number and lumi ranges
      if (std::getline(segmentStream, runPart, ':') && std::getline(segmentStream, lumiPart)) {
        Run run = std::stoul(runPart);

        // Check run filter
        if (noRunFilter || ((firstRun <= run) && (run <= lastRun))) {
          std::istringstream lumiStream(lumiPart);
          std::string lumiRange;

          // Split lumi ranges by comma
          while (std::getline(lumiStream, lumiRange, ',')) {
            size_t dashPos = lumiRange.find('-');
            if (dashPos != std::string::npos) {
              LumiBlock firstLumi = std::stoul(lumiRange.substr(0, dashPos));
              LumiBlock lastLumi = std::stoul(lumiRange.substr(dashPos + 1));
              accept.emplace_back(run, firstLumi, lastLumi);
            } else {
              std::cerr << "ERROR: Invalid LumiBlock range format: " << lumiRange << std::endl;
            }
          }
        }
      } else {
        std::cerr << "ERROR: Invalid run segment format: " << runSegment << std::endl;
      }
    }

    return LumiMask(accept);
  }



private:
  std::vector<LumiBlockRange> m_accept;
};
