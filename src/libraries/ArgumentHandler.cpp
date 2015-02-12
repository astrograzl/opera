#include "libraries/ArgumentHandler.h"
#include <iostream>

template <>
void ArgumentList<std::string>::SetToPassedString(unsigned int index, std::string value) {
	SetToPassedValue(index, value);
}

std::string ArgumentMap::ArgumentDescriptionBuilder(ArgumentInterface& argList, unsigned int index, std::string name, char typechar) {
	ArgumentFlag flag = argList.GetFlag(index);
	std::string temp;
	if (flag != REQUIRED_ARGUMENT) temp += "[";
	temp += "--";
	temp += name;
	if (flag != SWITCH) {
		if (flag == FLEXIBLE_SWITCH) temp += "[";
		if (typechar == 'd') temp += "=double";
		else if (typechar == 'i') temp += "=int";
		else if (typechar == 's') temp += "=string";
		else if (typechar == 'b') temp += "=bool";
		else if (typechar == 'u') temp += "=unsigned";
		else throw std::logic_error("invalid character in argument map");
		if (flag == FLEXIBLE_SWITCH) temp += "]";
	}
	if (flag != REQUIRED_ARGUMENT) temp += "]";
	temp += "\t";
	temp += argList.GetDescription(index);
	return temp;
}

void ArgumentMap::PassArgumentValue(std::string arg, std::string value) {
	Selector(argMap[arg].mapchar).SetToPassedString(argMap[arg].index, value);
}

void ArgumentMap::PassSwitch(std::string arg) {
	Selector(argMap[arg].mapchar).SetToSwitchValue(argMap[arg].index);
}

void ArgumentMap::PassNothing(std::string arg) {
	Selector(argMap[arg].mapchar).SetToDefaultValue(argMap[arg].index);
}

std::string ArgumentMap::GetDescription(std::string arg) {
	return ArgumentDescriptionBuilder(Selector(argMap[arg].mapchar), argMap[arg].index, arg, argMap[arg].mapchar);
}

ArgumentInterface& ArgumentMap::Selector(char typechar) {
	switch (typechar) {
	case 'd':
		return dynamic_cast<ArgumentInterface&>(argListDouble);
	case 'i':
		return dynamic_cast<ArgumentInterface&>(argListInt);
	case 's':
		return dynamic_cast<ArgumentInterface&>(argListString);
	case 'b':
		return dynamic_cast<ArgumentInterface&>(argListBool);
	case 'u':
		return dynamic_cast<ArgumentInterface&>(argListUnsigned);
	default:
		throw std::logic_error("invalid character in argument map");
	}
}

void ArgumentMap::AddMapping(std::string name, char mapchar, int index) {
	if (argMap.count(name) != 0) throw std::logic_error("multiple arguments mapped to the name " + name);
	argMap[name] = MapIndex(mapchar, index);
}

void ArgumentMap::Add(std::string name, ArgumentWrapper<double> variable) {
	AddMapping(name, 'd', argListDouble.Size());
	argListDouble.Add(variable);
}

void ArgumentMap::Add(std::string name, ArgumentWrapper<int> variable) {
	AddMapping(name, 'i', argListInt.Size());
	argListInt.Add(variable);
}

void ArgumentMap::Add(std::string name, ArgumentWrapper<std::string> variable) {
	AddMapping(name, 's', argListString.Size());
	argListString.Add(variable);
}

void ArgumentMap::Add(std::string name, ArgumentWrapper<bool> variable) {
	AddMapping(name, 'b', argListBool.Size());
	argListBool.Add(variable);
}

void ArgumentMap::Add(std::string name, ArgumentWrapper<unsigned> variable) {
	AddMapping(name, 'u', argListUnsigned.Size());
	argListUnsigned.Add(variable);
}

void ArgumentHandler::AddSwitch(std::string name, bool& variable, std::string description) {
	ArgumentWrapper<bool> temp(variable);
	temp.flag = SWITCH;
	temp.defaultValue = false;
	temp.switchValue = true;
	temp.description = description;
	expectedArgs.Add(name, temp);
	argNames.push_back(name);
}

void ArgumentHandler::ParsingReader(int argc, char* argv[]) {
	programName = argv[0];
	for (int i = 1; i < argc; i++) {
		std::string temp = argv[i];
		if (temp.length() < 3 || temp[0] != '-' || temp[1] != '-') throw std::runtime_error("expected -- followed by argument name, got " + temp);
		std::size_t splitpos = temp.find('=');
		std::string readName = temp.substr(2, splitpos - 2);
		BoolString readVal = temp.substr(splitpos + 1, std::string::npos);
		if (splitpos == std::string::npos) readVal = false;
		if (readArgs.count(readName) > 0) throw std::runtime_error("same argument name passed multiple times: " + readName);
		readArgs[readName] = readVal;
	}
}

void ArgumentHandler::Parse(int argc, char* argv[]) {
	class HelpMessage : std::exception { };
	try {
		ParsingReader(argc, argv);
		if (readArgs.count("help") > 0) {
			throw HelpMessage();
		}
		for (unsigned int i = 0; i < argNames.size(); i++) {
			try {
				if (readArgs.count(argNames[i]) == 0) expectedArgs.PassNothing(argNames[i]);
				else if (!readArgs[argNames[i]]) expectedArgs.PassSwitch(argNames[i]);
				else expectedArgs.PassArgumentValue(argNames[i], readArgs[argNames[i]]);
			}
			catch (std::runtime_error& error) {
				throw std::runtime_error(error.what() + (": " + argNames[i]));
			}
			readArgs.erase(argNames[i]);
		}
		if (readArgs.size() > 0) throw std::runtime_error("unknown argument read in: " + readArgs.begin()->first);
	}
	catch (HelpMessage) {
		PrintUsageSyntax();
		throw std::exception();
	}
	catch (std::exception& error) {
		std::cout << "Error: " << error.what() << std::endl;
		std::cout << "See --help for available options." << std::endl;
		throw error;
	}
}

void ArgumentHandler::PrintUsageSyntax() {
	std::cout << "Usage options" << (programName.empty() ? "" : " for ") << programName << ":" << std::endl;
	for (unsigned int i = 0; i < argNames.size(); i++) {
		std::cout << expectedArgs.GetDescription(argNames[i]) << std::endl;
	}
}
