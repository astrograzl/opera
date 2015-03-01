#ifndef OPERAARGUMENTHANDLER_H
#define OPERAARGUMENTHANDLER_H

#include "libraries/ArgumentHandler.h"

class operaArgumentHandler : public ArgumentHandler {
public:
	operaArgumentHandler();
	
	bool verbose;
	bool debug;
	bool trace;
	bool plot;
};

#endif
