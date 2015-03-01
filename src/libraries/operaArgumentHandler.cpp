#include "libraries/operaArgumentHandler.h"

operaArgumentHandler::operaArgumentHandler() {
	AddSwitch("verbose", verbose, "output informational messages");
	AddSwitch("debug", debug, "output debug messages");
	AddSwitch("trace", trace, "output trace messages");
	AddSwitch("plot", plot, "produce plots");
}
