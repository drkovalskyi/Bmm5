#include "LumiMask.h"
#include <iostream>

// g++ -o test_LumiMask test_LumiMask.cpp
// ./test_LumiMask
int main() {
    std::string input = "12345:1-10,20-30;67890:5-15";
    LumiMask mask = LumiMask::fromCustomString(input);

    // Test the accept function
    if (mask.accept(12345, 5)) {
        std::cout << "LumiBlock 5 in Run 12345 is accepted." << std::endl;
    } else {
        std::cout << "LumiBlock 5 in Run 12345 is not accepted." << std::endl;
    }

    return 0;
}
