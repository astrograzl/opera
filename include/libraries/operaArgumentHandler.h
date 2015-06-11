#ifndef OPERAARGUMENTHANDLER_H
#define OPERAARGUMENTHANDLER_H

#include "libraries/ArgumentHandler.h"

class operaArgumentHandler : public ArgumentHandler {
public:
	operaArgumentHandler();
	void AddPlotFileArguments(std::string& plotfilename, std::string& datafilename, std::string& scriptfilename, bool& interactive);
	void AddOrderLimitArguments(int& ordernumber, int& minorder, int& maxorder, const int default_value);
	
	bool verbose;
	bool debug;
	bool trace;
	bool plot;
};

#endif
